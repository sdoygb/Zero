"""
z0_genesis.py --- **只用 Z0 的最小生成器（宇宙初开）**：从零态出发，零参数。

**为什么需要它（对前几轮工作的诚实更正）**

我前几轮的脚本（`l2_branch.py` / `l2_capacity.py` / `l0_closure.py`）**都不是最小理论**，
它们照抄了语料里已经加了料的模型。逐项列出被塞进去的东西：

| 塞进去的结构 | 出处 | 在 Z0 里的地位 |
|:--|:--|:--|
| 8 个手工初始词 `+ - ++ -- +++ --- +-+ -+-` | `zero_sum_closure_exit.py:63-72` | **无**。Z0 没有种子 |
| `LAYER_LIFETIME = 4` | 同上 `:61` | **自由参数**（`Z0:166` 自述） |
| 摧毁周期 `T = 3` | `zero_sum_periodic_destruction.py:60` | 语料自登记 `R-Z-PERIODIC-WASHOUT-RESEED`，「周期未导出」 |
| **闭合即吸收** | `zero_sum_closure_exit.py:181-188` | 语料自登记 **`R-Z-CLOSURE-EXIT-RULE`**（`G28:120`）—— **且与 Z0① 相反** |
| 容量 K | 我自己扫的参数 | 语料自述「饱和机制是 L2 的具名输入」 |

**关键**：`Z0` §1 逐字写着「① 零不停留：零（**平衡／闭合态**）不是静止态，它持续运动」。
「闭合态**不**静止」=> **回到零并不终止运动**。所以「闭合即吸收」不是 Z0 的条款。
本脚本把它当**开关**，两种读法都跑。

---

**最小生成器（零参数）**

    状态 = ±1 词（Z1 导出：运动有先后 => 步，有所在 => 位置）
    起点 = 零态 ∅              <- 不是 8 个种子，是 Z0 本身
    一步 = 全分支（Z0③）：alive <- { w+ , w- : w ∈ alive }
    记录 = 词和为 0 者（Z0① 的「闭合态」）
    终止 = 从不（Z0②）

    absorb=True   闭合词退出（语料读法）-> 记录流 = 首次通过词，计数 = **2·Catalan(n/2-1)**
    absorb=False  闭合词继续（Z0① 字面）-> 记录流 = 全部平衡词，计数 = **C(n, n/2)**

---

**表示与速度**

从单个起点出发、无重播种 => 每个词恰有一条路径到达 => **重数恒为 1** => 状态是**集合**而非计数向量。

    => 用 bitset：比 int64 计数向量**密 64 倍**
    => 一步 = 「两份拷贝 + 与关闭掩码」= 纯内存带宽操作（GPU 赢的那类负载）

    深度 n 的状态 = **恰好 2^n 个比特**（n>=6 时 2^n 是 64 的整数倍）
    关闭掩码：位 i 关闭 ⟺ popcount(i) == n/2  （词和 = 0 ⟺ +1 的个数 = 长度的一半）
    word 技巧：位下标 i = 64w+b ⟹ popcount(i) = popcount(w) + popcount(b)
              ⟹ mask(w) = PAT[target - popcount(w)]
"""
from __future__ import annotations

import json
import os
import sys
import time
from math import comb

import numpy as np

from zcl import Engine


def build_pat() -> np.ndarray:
    """PAT[t] = 那些 popcount 恰为 t 的位位置的并集（t = 0..64）。"""
    pat = np.zeros(65, dtype=np.uint64)
    for t in range(65):
        v = 0
        for b in range(64):
            if bin(b).count("1") == t:
                v |= 1 << b
        pat[t] = np.uint64(v)
    return pat


PAT = build_pat()


def cat(k: int) -> int:
    return comb(2 * k, k) // (k + 1)


def predicted(n: int, absorb: bool) -> int:
    if n < 2 or n % 2:
        return 0
    return 2 * cat(n // 2 - 1) if absorb else comb(n, n // 2)


# --------------------------------------------------- 深度 6 的显式初值
def seed_depth6(absorb: bool = True):
    """直接枚举 2^6 = 64 个长度 6 的词，返回 (alive 掩码, 本层关闭数)。"""
    alive, closed = 0, 0
    for i in range(64):
        s, ret = 0, False
        for k in range(6):
            s += 1 if (i >> k) & 1 else -1
            if s == 0 and k < 5:
                ret = True
                break
        if ret:
            continue
        if s == 0:
            closed += 1
            if absorb:
                continue
        alive |= 1 << i
    return np.array([np.uint64(alive)], dtype=np.uint64), closed


# --------------------------------------------------------------- CPU 参考
def cpu_run(N: int, absorb: bool = True):
    """uint64 bitset 参考实现。深度 n 的状态恰好 2^n 比特（n>=6）。"""
    A, cl6 = seed_depth6(absorb)
    rows = [{"n": 6, "closed": cl6, "alive": int(np.bitwise_count(A).sum())}]
    for n in range(6, N):
        B = np.concatenate([A, A])                 # 位数组恰好翻倍
        closed = 0
        if (n + 1) % 2 == 0:
            target = (n + 1) // 2
            widx = np.arange(B.size, dtype=np.uint64)
            need = target - np.bitwise_count(widx).astype(np.int64)
            valid = (need >= 0) & (need <= 64)
            mc = np.zeros(B.size, dtype=np.uint64)
            mc[valid] = PAT[need[valid]]
            closed = int(np.bitwise_count(B & mc).sum())
            if absorb:
                B = B & ~mc
        rows.append({"n": n + 1, "closed": closed,
                     "alive": int(np.bitwise_count(B).sum())})
        A = B
    return rows


# ------------------------------------------------------------------ GPU
KERNEL = r"""
// B[w] = A[w % La] & ~mc ；mc 由 PAT 查表；同时原子累加被移除的位数。
// target < 0 表示本层不关闭（奇数层）。
__kernel void step_bits(const __global ulong* A, __global ulong* B,
                        __constant ulong* PAT, __global uint* cnt,
                        const uint La, const int target) {
    size_t w = get_global_id(0);
    ulong a = A[w % La];
    ulong mc = 0UL;
    if (target >= 0) {
        int need = target - (int)popcount((ulong)w);
        if (need >= 0 && need <= 64) mc = PAT[need];
    }
    B[w] = a & ~mc;
    ulong rm = a & mc;
    if (rm) atomic_add(cnt, (uint)popcount(rm));
}
"""


def gpu_run(N: int, eng: Engine, verbose: bool = True):
    import pyopencl as cl
    ctx, q = eng.ctx, eng.queue
    mf = cl.mem_flags

    pat_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=PAT)
    cnt_buf = cl.Buffer(ctx, mf.READ_WRITE, 4)
    k = eng.kernel("z0_step", KERNEL, "step_bits")

    A, cl6 = seed_depth6(True)
    rows = [{"n": 6, "closed": cl6}]
    if N <= 6:
        return rows, 0.0

    Ab = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=A)
    La = int(A.size)
    t_step = 0.0
    for n in range(6, N):
        Lb = La * 2
        Bb = cl.Buffer(ctx, mf.READ_WRITE, Lb * 8)
        cl.enqueue_copy(q, cnt_buf, np.zeros(1, np.uint32)).wait()
        target = ((n + 1) // 2) if (n + 1) % 2 == 0 else -1
        t0 = time.time()
        k(q, (Lb,), None, Ab, Bb, pat_buf, cnt_buf,
          np.uint32(La), np.int32(target)).wait()
        t_step += time.time() - t0
        cnt = int(eng.from_device(cnt_buf, 1, np.uint32)[0])
        rows.append({"n": n + 1, "closed": cnt})
        if verbose:
            print(f"    n={n+1:3d}  closed={cnt:>13,}  "
                  f"state={Lb*8/2**20:8.1f} MB  cum {t_step:6.2f}s", flush=True)
        Ab = Bb
        La = Lb
    return rows, t_step


# ------------------------------------------------------------------ 主
def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 26
    mode = sys.argv[2] if len(sys.argv) > 2 else "cpu"
    N_cpu = min(N, 26)

    print("=" * 96)
    print("Z0 最小生成器：从零态出发（**零参数**）")
    print("=" * 96)
    print(f"\n{'n':>3} {'吸收：2·Cat(n/2-1)':>19} {'不吸收：C(n,n/2)':>18} {'倍数':>8}")
    for n in range(2, N + 1, 2):
        a, b = predicted(n, True), predicted(n, False)
        print(f"{n:>3} {a:>19,} {b:>18,} {b/a:>8.2f}")

    t0 = time.time()
    rows_cpu = cpu_run(N_cpu, absorb=True)
    tcpu = time.time() - t0
    bad = [(r["n"], r["closed"], predicted(r["n"], True))
           for r in rows_cpu if r["closed"] != predicted(r["n"], True)]
    print(f"\n[CPU bitset] 深度 6..{N_cpu}，{tcpu:.2f}s")
    print(f"  关闭计数 = 2·Catalan(n/2-1) ? "
          f"{'全部符合' if not bad else '不符: %s' % bad[:3]}")
    print(f"  {'n':>3} {'关闭':>14} {'活动':>16}")
    for r in rows_cpu:
        if r["n"] % 2 == 0:
            print(f"  {r['n']:>3} {r['closed']:>14,} {r['alive']:>16,}")

    if mode != "gpu":
        return
    eng = Engine()
    print(f"\n[GPU bitset] {eng.info()['name']}  -> 深度 {N}")
    rows_gpu, t_step = gpu_run(N, eng)
    ncmp = min(len(rows_gpu), len(rows_cpu))
    diffs = [(rows_gpu[i]["n"], rows_gpu[i]["closed"], rows_cpu[i]["closed"])
             for i in range(ncmp) if rows_gpu[i]["closed"] != rows_cpu[i]["closed"]]
    print(f"\n  GPU vs CPU（重叠 {ncmp} 层）: "
          f"{'完全一致' if not diffs else '不一致 %s' % diffs[:3]}")
    badg = [(r["n"], r["closed"], predicted(r["n"], True)) for r in rows_gpu
            if r["closed"] != predicted(r["n"], True)]
    print(f"  GPU 关闭计数 = 2·Catalan ? {'全部符合' if not badg else '不符 %s' % badg[:3]}")
    print(f"  GPU 纯计算 {t_step:.2f}s（深度 6..{N}）")
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "z0_genesis.json"), "w") as f:
        json.dump({"N": N, "rows_cpu": rows_cpu, "rows_gpu": rows_gpu,
                   "cpu_seconds": tcpu, "gpu_step_seconds": t_step}, f, indent=2)
    print("  saved results/z0_genesis.json")


if __name__ == "__main__":
    main()
