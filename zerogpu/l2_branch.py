"""
l2_branch.py --- L2 演化层的全分支动力学（GPU），以及它到 L0 项链类的推前测度 ω。

**更新律（逐字取自 simulations/zero_sum_closure_exit.py，无任何随机性）**

    每个站点持有一个"活动词"集合 E（词 = ±1 步序列，词长 = 年龄）。
    每一轮：
        for 每个活动词 w:
            if len(w) >= LIFETIME:      w 进 dead（终端），退出 E
            else:
                for step in (+1, -1):           # 全分支：Z0③ 不设概率
                    child = w + (step,)
                    if sum(child) == 0:         # 闭合
                        histories += child
                        zero_layer[canonical_cycle(child)] += 1   # <-- 推到 L0 项链类
                        退出 E
                    else:
                        child 进 E

**由此得到的结构性事实（本文件的核心观察，可独立核验）**

  1. 每个站点的 E 是一棵**二叉树**：长度 n 的词只有一个长度 n-1 的前驱，故两个不同父节点的
     子节点必不相同 => 原代码里的 `set` **从不发生去重**。整套动力学是**确定性枚举**，
     没有任何随机数。
  2. 因此 E(n) = {长度 n 的词 : 在初始词之后没有任何真前缀闭合} = **首次通过词**。
  3. 闭合事件的 `canonical_cycle(child)` 是长度 n、恰含 n/2 个 +1 的平衡词的循环类
     => 它就是 **L0 的项链**（见 l0_closure.py）。故 zero_layer 的计数测度正是
     L1′ 的**推前权重 ω**，而 q_n = Σ_C ω_C² 就是 R87 的对象。
  4. 闭合事件数在站点 '+' 上等于 Catalan 数：n=2,4,6,... -> 1,1,2,5,14,42,132,429,1430,4862。

原文档只跑到 LIFETIME=4（最多 16 个词）。本文件用 GPU 把它推到 n=32。

**一处引用更正（附带发现）**：Z0 §2.5 的对照表把 `closure_exit` 一行写成 8,16,16,16,0,0,0；
逐字复跑 simulations/zero_sum_closure_exit.py 得
    continue_after_closure -> 8,16,16,16,0,0,0
    exit_after_closure     -> 8,12,10,6,0,0,0
即表里引的是**对照臂**的数字，真正的闭合退出臂是 8,12,10,6,0,0,0。两臂都熄灭，
故 §2.5 的结论（熄灭 => 逼出重播种）不受影响；Z0_check.py F4 断言的也是"末项为 0"，通过。
"""
from __future__ import annotations

import json
import os
import sys
import time
from collections import Counter

import numpy as np

from zcl import Engine

LAYER_LIFETIME_DOC = 4
EPOCHS_DOC = 6
INITIAL_WORDS = ("+", "-", "++", "--", "+++", "---", "+-+", "-+-")


# ------------------------------------------------------------ 参考实现（逐字）
def word_from_string(text):
    return tuple(1 if c == "+" else -1 for c in text)


def is_closed(word):
    return bool(word) and sum(word) == 0


def canonical_cycle(word):
    n = len(word)
    bits = 0
    for i, v in enumerate(word):
        if v == 1:
            bits |= 1 << i
    best = bits
    for _ in range(n - 1):
        bits = ((bits << 1) | (bits >> (n - 1))) & ((1 << n) - 1)
        best = min(best, bits)
    return best


def reference_simulation(lifetime=LAYER_LIFETIME_DOC, epochs=EPOCHS_DOC,
                         keep_closed_active=False):
    active = [{word_from_string(t)} for t in INITIAL_WORDS]
    histories = [set() for _ in INITIAL_WORDS]
    dead = [set() for _ in INITIAL_WORDS]
    zero_layer = Counter()
    traj = []
    for epoch in range(epochs + 1):
        allp = [w for s in active for w in s]
        traj.append({
            "epoch": epoch,
            "active_paths": len(allp),
            "active_closed": sum(1 for w in allp if is_closed(w)),
            "history_records": sum(len(h) for h in histories),
            "zero_events": sum(zero_layer.values()),
            "zero_modes": len(zero_layer),
            "dead_paths": sum(len(d) for d in dead),
        })
        if epoch == epochs:
            break
        nxt = [set() for _ in INITIAL_WORDS]
        for site, paths in enumerate(active):
            for w in paths:
                if len(w) >= lifetime:
                    dead[site].add(w)
                    continue
                for step in (1, -1):
                    c = w + (step,)
                    if is_closed(c):
                        histories[site].add(c)
                        zero_layer[canonical_cycle(c)] += 1
                        if keep_closed_active:
                            nxt[site].add(c)
                    else:
                        nxt[site].add(c)
        active = nxt
    return traj, zero_layer


# ------------------------------------------------------------------ GPU 扫描
KERNEL = r"""
// 枚举总长度 n 的"首次通过闭合词"，且前 m 位固定为初始词 init。
//   低 m 位 = 初始词（bit=1 表示 +1）；位 m..n-1 = 扩展，由 g 的 popcount 过滤后给出。
// 命中条件：对 k = m..n-2，前缀和 != 0；且全词和 == 0。
// 命中则写出该词 n 位循环类的最小旋转。
__kernel void fp_closure2(__global ulong* out, __global uint* cnt,
                          const ulong init, const int m, const int n,
                          const int need_ones, const ulong gmask, const ulong mask) {
    size_t g = get_global_id(0);
    if ((ulong)g > gmask) return;
    if ((int)popcount((ulong)g) != need_ones) return;
    ulong w = init | (((ulong)g) << (ulong)m);
    int s = 0;
    for (int k = 0; k < m; ++k) s += ((w >> k) & 1UL) ? 1 : -1;
    for (int k = m; k < n; ++k) {
        s += ((w >> k) & 1UL) ? 1 : -1;
        if (k < n - 1 && s == 0) return;      // 提前闭合 => 该分支已被剪掉
    }
    if (s != 0) return;
    ulong best = w & mask, r = w & mask;
    const int sh = n - 1;
    for (int q = 1; q < n; ++q) {
        r = ((r << 1) | (r >> sh)) & mask;
        best = min(best, r);
    }
    size_t p = atomic_inc(cnt);
    out[p] = best;
}
"""


def site_closure_events(eng: Engine, initial: str, n: int, cap: int = 1 << 27):
    """返回该站点在总长度 n 上全部首次通过闭合词的规范项链类（uint64 数组）。"""
    m = len(initial)
    if n <= m:
        return np.zeros(0, np.uint64)
    init = 0
    ones = 0
    for i, c in enumerate(initial):
        if c == "+":
            init |= 1 << i
            ones += 1
    d = n - m
    need = n // 2 - ones           # 扩展里必须有的 +1 个数
    if need < 0 or need > d:
        return np.zeros(0, np.uint64)

    import pyopencl as cl
    k = eng.kernel("l2_fp2", KERNEL, "fp_closure2")
    cnt_buf = eng.buffer(4, readonly=False)
    cl.enqueue_copy(eng.queue, cnt_buf, np.zeros(1, np.uint32)).wait()
    outbuf = eng.empty(cap, np.uint64)
    eng.run(k, (1 << d,), outbuf, cnt_buf, np.uint64(init), np.int32(m), np.int32(n),
            np.int32(need), np.uint64((1 << d) - 1), np.uint64((1 << n) - 1))
    eng.finish()
    cnt = int(eng.from_device(cnt_buf, 1, np.uint32)[0])
    if cnt == 0:
        return np.zeros(0, np.uint64)
    return eng.from_device(outbuf, cnt, np.uint64)[:cnt]


def main():
    eng = Engine()
    print("device:", eng.info(), flush=True)
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    nmax += nmax % 2               # 只取偶数

    # ---- 1. 校准：逐字复现原仿真 ----
    print("\n[1] 参考实现（逐字复现 zero_sum_closure_exit.py, LIFETIME=4, EPOCHS=6）")
    out = {"reference": {}}
    for keep, name in ((False, "exit_after_closure"), (True, "continue_after_closure")):
        traj, zl = reference_simulation(keep_closed_active=keep)
        series = [t["active_paths"] for t in traj]
        out["reference"][name] = {"active_series": series,
                                  "zero_events": int(sum(zl.values())),
                                  "zero_modes": len(zl),
                                  "dead_total": traj[-1]["dead_paths"]}
        print(f"   {name:>24}: active={series}  zero_events={sum(zl.values())} "
              f"zero_modes={len(zl)}")

    # ---- 2. GPU 大规模扫描 ----
    print(f"\n[2] GPU 扫描：每站点首次通过闭合词，总长度 n 到 {nmax}")
    rows = []
    store = {}
    for site in INITIAL_WORDS:
        for n in range(len(site) + 1, nmax + 1):
            if n % 2:
                continue
            t0 = time.time()
            hits = site_closure_events(eng, site, n)
            rows.append({"site": site, "n": n, "events": int(len(hits)),
                         "distinct_classes": int(len(np.unique(hits))) if len(hits) else 0,
                         "seconds": round(time.time() - t0, 3)})
            if len(hits):
                store[f"{site}|{n}"] = hits.tolist()
            if site == "+":
                print(f"   n={n:3d}  events={len(hits):11,d}  "
                      f"distinct={rows[-1]['distinct_classes']:9,d}  [{time.time()-t0:5.2f}s]",
                      flush=True)
    res = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, f"l2_branch_n{nmax}.json"), "w") as f:
        json.dump({"device": eng.info(), "nmax": nmax, "rows": rows,
                   "reference": out["reference"]}, f)
    print(f"\nsaved results; {len(rows)} (site,n) cells")
    # 站点 '+' 的事件数序列 = Catalan?
    seq = [r["events"] for r in rows if r["site"] == "+"]
    from math import comb
    cat = [comb(2 * j, j) // (j + 1) for j in range(len(seq) + 2)]
    print(f"site '+' events    = {seq}")
    print(f"Catalan C_0..C_?   = {cat[:len(seq)]}")
    print(f"match = {seq == cat[:len(seq)]}")


if __name__ == "__main__":
    main()
