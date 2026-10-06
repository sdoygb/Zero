"""
z0_record.py --- 最小生成器的**记录流**落到 L0 项链上的分布。

上一轮建好了零参数生成器 `z0_genesis.py`（从零态出发，全分支，闭合即记录）。
本脚本把它产出的**记录流**（每层闭合词）规范化成 L0 项链类，回答：

  §1 覆盖率：记录流打到了全部 K(n) 条项链中的多少条？
  §2 单射性：不同的闭合词会不会互为循环旋转？
  §3 手征结构：闭合类在"反射"（反序+变号）下成对还是自镜像？
  §4 时间线：从零态出发，各层"沉积"了多少对象（宇宙初开的时间序列）。

**对象与约定**

  闭合词：长度 n 的 ±1 词，词和为 0（Z0① 的"闭合态"），且**在此前从未闭合过**
          （吸收读法；语料口径 R-Z-CLOSURE-EXIT-RULE）。
  位编码：位 j = 1 表示第 j 步为 +1。词和 = 0 ⟺ popcount = n/2。
  项链类：canon(w) = n 个循环旋转中的最小位模式（与 l0_closure.py 同一约定）。
  反射  ：reflect(w) = 反序 + 变号（G35 的 \bar w），位表示 = 反转 n 位并全部取反。
"""
from __future__ import annotations

import json
import os
import sys
import time
from math import comb

import numpy as np

from zcl import Engine
import l0_closure as L0

KERNEL = r"""
// 与 z0_genesis 相同的推进，但把"本层被移除的位"（= 闭合词）写到 rm 里。
__kernel void step_rec(const __global ulong* A, __global ulong* B, __global ulong* rm,
                       __constant ulong* PAT, __global uint* cnt,
                       const uint La, const int target) {
    size_t w = get_global_id(0);
    ulong a = A[w % La];
    ulong mc = 0UL;
    if (target >= 0) {
        int need = target - (int)popcount((ulong)w);
        if (need >= 0 && need <= 64) mc = PAT[need];
    }
    ulong r = a & mc;
    rm[w] = r;
    B[w] = a & ~mc;
    if (r) atomic_add(cnt, (uint)popcount(r));
}
"""


def seed_depth6():
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
            continue
        alive |= 1 << i
    return np.array([np.uint64(alive)], dtype=np.uint64), closed


def bits_to_words(rm: np.ndarray, base_bit: int) -> np.ndarray:
    """把位数组 rm（每元素 64 位）展开成置位的绝对位下标，分块以避免爆内存。"""
    out = []
    chunk = 1 << 20                     # 1M 个 uint64 = 8 MB -> 展开后 64 MB
    for s in range(0, rm.size, chunk):
        e = min(s + chunk, rm.size)
        sub = rm[s:e]
        nz = np.nonzero(sub)[0]
        if nz.size == 0:
            continue
        bits = np.unpackbits(sub[nz].view(np.uint8).reshape(-1, 8),
                             axis=1, bitorder="little")
        bi, bb = np.nonzero(bits)
        out.append(((s + nz[bi]).astype(np.uint64) * np.uint64(64)
                    + bb.astype(np.uint64)) + np.uint64(base_bit))
    return np.concatenate(out) if out else np.zeros(0, dtype=np.uint64)


def canon_rot(w: np.ndarray, n: int) -> np.ndarray:
    """所有 n 位循环旋转取最小（向量化）。"""
    mask = np.uint64((1 << n) - 1)
    best = (w & mask).copy()
    r = best.copy()
    for _ in range(n - 1):
        r = ((r << np.uint64(1)) | (r >> np.uint64(n - 1))) & mask
        np.minimum(best, r, out=best)
    return best


def reflect_word(w: np.ndarray, n: int) -> np.ndarray:
    """reflect = 反转 n 位 + 全部取反。"""
    mask = np.uint64((1 << n) - 1)
    out = np.zeros_like(w)
    r = w & mask
    for k in range(n):
        out |= ((r >> np.uint64(k)) & np.uint64(1)) << np.uint64(n - 1 - k)
    return (~out) & mask


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 28
    eng = Engine()
    print("device:", eng.info()["name"], flush=True)
    import pyopencl as cl
    ctx, q = eng.ctx, eng.queue
    mf = cl.mem_flags

    PAT = np.zeros(65, dtype=np.uint64)
    for t in range(65):
        v = 0
        for b in range(64):
            if bin(b).count("1") == t:
                v |= 1 << b
        PAT[t] = np.uint64(v)
    pat_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=PAT)
    cnt_buf = cl.Buffer(ctx, mf.READ_WRITE, 4)
    k = eng.kernel("z0_rec", KERNEL, "step_rec")

    A, _ = seed_depth6()
    Ab = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=A)
    La = 1
    # 注意：深度 n 的状态位下标 i **就是**该词的 n 位模式，所以偏移恒为 0。
    rows = []
    t_all = time.time()
    for n in range(6, N):
        Lb = La * 2
        Bb = cl.Buffer(ctx, mf.READ_WRITE, Lb * 8)
        Rb = cl.Buffer(ctx, mf.READ_WRITE, Lb * 8)
        cl.enqueue_copy(q, cnt_buf, np.zeros(1, np.uint32)).wait()
        target = ((n + 1) // 2) if (n + 1) % 2 == 0 else -1
        k(q, (Lb,), None, Ab, Bb, Rb, pat_buf, cnt_buf,
          np.uint32(La), np.int32(target)).wait()
        cnt = int(eng.from_device(cnt_buf, 1, np.uint32)[0])

        if cnt:
            t0 = time.time()
            rm = eng.from_device(Rb, Lb, np.uint64)
            words = bits_to_words(rm, 0)
            t_dl = time.time() - t0
            nn = n + 1
            canon = canon_rot(words, nn)
            K = L0.necklace_count(nn, nn // 2)
            uniq = np.unique(canon)
            rows.append({
                "n": nn, "events": int(words.size),
                "distinct_classes": int(uniq.size), "K_n": int(K),
                "coverage": round(uniq.size / K, 6),
                "injective": bool(uniq.size == words.size),
                "collisions": int(words.size - uniq.size),
                "download_seconds": round(t_dl, 2),
            })
            print(f"  n={nn:3d} events={words.size:>11,} distinct={uniq.size:>11,} "
                  f"K(n)={K:>11,} coverage={uniq.size/K:7.4f} "
                  f"{'单射' if uniq.size==words.size else '有碰撞 %d' % (words.size-uniq.size)}"
                  f"  [{t_dl:5.1f}s]", flush=True)

            # ---- 手征结构（只做到 n<=24，reflect 是 O(n) 的 Python 循环）----
            if nn <= 24:
                refl = reflect_word(uniq, nn)
                c2 = canon_rot(refl, nn)
                self_mirror = int(np.sum(c2 == uniq))
                # 手征对：c2 != uniq；对的大小 2（对合）
                chiral_pairs = int(np.sum(c2 != uniq)) // 2
                rows[-1]["self_mirror_classes"] = self_mirror
                rows[-1]["chiral_pairs"] = chiral_pairs
                print(f"          手征：自镜像类 {self_mirror}，手征对 {chiral_pairs}", flush=True)
        else:
            rows.append({"n": n + 1, "events": 0})
        Ab = Bb
        La = Lb

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "z0_record.json"), "w") as f:
        json.dump({"N": N, "rows": rows,
                   "total_seconds": round(time.time() - t_all, 1)}, f, indent=2)
    print(f"\n总耗时 {time.time()-t_all:.1f}s，saved results/z0_record.json")
    print(f"\n{'n':>3} {'闭合词':>13} {'相异类':>13} {'K(n)':>13} {'覆盖率':>9} {'自镜像':>9} {'手征对':>9}")
    for r in rows:
        if r.get("events"):
            print(f"{r['n']:>3} {r['events']:>13,} {r['distinct_classes']:>13,} "
                  f"{r['K_n']:>13,} {r['coverage']:>9.4f} "
                  f"{r.get('self_mirror_classes','-'):>9} {r.get('chiral_pairs','-'):>9}")


if __name__ == "__main__":
    main()
