"""
z0_record_gpu.py --- 把 `z0_record.py` 的两个 CPU 瓶颈搬到 GPU。

**为什么**（实测时间分解，`z0_record.py` 到 n=32，合计 12.45 s）：

    kernel(生成)   0.120 s   1.0%     <- 已经在 GPU 上
    download       0.335 s   2.7%
    expand         3.529 s  28.4%     <- bits_to_words，numpy 逐块 unpackbits
    canon          7.746 s  62.2%     <- canon_rot，numpy 做 32 次旋转迭代
    unique         0.707 s   5.7%
    reflect        0.010 s   0.1%

=> **GPU 内核只占 1%**，所以"再加一块 GPU"（如 Intel HD630）最多优化 1%。
   真正的瓶颈是算法：规范化（62%）与位展开（28%）都还在 CPU 上。
   本脚本把这两个搬到 OpenCL：压缩用原子追加，规范化用一格一词的旋转扫描。

   预期的量级：19.4M 个词 × 32 次旋转 = 620M 次旋转步，在 RX 570 上是亚秒级。
"""
from __future__ import annotations

import json
import os
import sys
import time
from math import comb

import numpy as np

from zcl import Engine

STEP = r"""
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

COMPACT = r"""
// 把 rm 中置位的位下标（= 闭合词）收集成紧凑数组。
// 最低置位下标 b 用 popcount((v & -v) - 1) 求，避免依赖 ctz 扩展。
__kernel void compact_bits(const __global ulong* rm, __global ulong* out,
                           __global uint* cnt) {
    size_t w = get_global_id(0);
    ulong v = rm[w];
    if (!v) return;
    ulong base = (ulong)w * 64UL;
    while (v) {
        ulong low = v & (~v + 1UL);
        int b = (int)popcount(low - 1UL);
        size_t p = atomic_inc(cnt);
        out[p] = base + (ulong)b;
        v ^= low;
    }
}
"""

CANON = r"""
// n 位循环旋转取最小（一格一词，每格做 n 次旋转）
__kernel void canon_all(const __global ulong* in, __global ulong* out,
                        const int n, const ulong mask) {
    size_t g = get_global_id(0);
    ulong w = in[g] & mask;
    ulong best = w, r = w;
    const int sh = n - 1;
    for (int q = 1; q < n; ++q) {
        r = ((r << 1) | (r >> sh)) & mask;
        best = min(best, r);
    }
    out[g] = best;
}
"""


def build_pat():
    pat = np.zeros(65, dtype=np.uint64)
    for t in range(65):
        v = 0
        for b in range(64):
            if bin(b).count("1") == t:
                v |= 1 << b
        pat[t] = np.uint64(v)
    return pat


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


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    verbose = "-q" not in sys.argv
    import pyopencl as cl
    eng = Engine()
    ctx, q = eng.ctx, eng.queue
    mf = cl.mem_flags
    print("device:", eng.info()["name"], flush=True)

    PAT = build_pat()
    pat_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=PAT)
    cnt_buf = cl.Buffer(ctx, mf.READ_WRITE, 4)
    k_step = eng.kernel("z0r_step", STEP, "step_rec")
    k_comp = eng.kernel("z0r_comp", COMPACT, "compact_bits")
    k_canon = eng.kernel("z0r_canon", CANON, "canon_all")

    A, _ = seed_depth6()
    Ab = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=A)
    La = 1
    rows = []
    T = {"kernel": 0.0, "compact": 0.0, "canon": 0.0, "download": 0.0, "unique": 0.0}
    t_all = time.time()
    for n in range(6, N):
        nn = n + 1
        Lb = La * 2
        Bb = cl.Buffer(ctx, mf.READ_WRITE, Lb * 8)
        Rb = cl.Buffer(ctx, mf.READ_WRITE, Lb * 8)
        cl.enqueue_copy(q, cnt_buf, np.zeros(1, np.uint32)).wait()
        target = (nn // 2) if nn % 2 == 0 else -1
        t0 = time.time()
        k_step(q, (Lb,), None, Ab, Bb, Rb, pat_buf, cnt_buf, np.uint32(La), np.int32(target)).wait()
        T["kernel"] += time.time() - t0
        cnt = int(eng.from_device(cnt_buf, 1, np.uint32)[0])

        if cnt:
            # 压缩：位数组 -> 词下标数组
            cap = cnt
            wb = cl.Buffer(ctx, mf.READ_WRITE, max(cap, 1) * 8)
            cl.enqueue_copy(q, cnt_buf, np.zeros(1, np.uint32)).wait()
            t0 = time.time()
            k_comp(q, (Lb,), None, Rb, wb, cnt_buf).wait()
            T["compact"] += time.time() - t0
            got = int(eng.from_device(cnt_buf, 1, np.uint32)[0])

            # 规范化：全部留在 GPU 上
            cb = cl.Buffer(ctx, mf.READ_WRITE, max(cap, 1) * 8)
            t0 = time.time()
            k_canon(q, (got,), None, wb, cb, np.int32(nn), np.uint64((1 << nn) - 1)).wait()
            T["canon"] += time.time() - t0

            t0 = time.time()
            can = eng.from_device(cb, got, np.uint64)
            T["download"] += time.time() - t0
            t0 = time.time()
            uniq = np.unique(can)
            T["unique"] += time.time() - t0

            import l0_closure as L0
            K = L0.necklace_count(nn, nn // 2)
            rows.append({"n": nn, "events": int(got), "distinct_classes": int(uniq.size),
                         "K_n": int(K), "coverage": round(uniq.size / K, 6),
                         "injective": bool(uniq.size == got),
                         "collisions": int(got - uniq.size)})
            if verbose:
                print(f"  n={nn:3d} events={got:>11,} distinct={uniq.size:>11,} "
                      f"coverage={uniq.size/K:7.4f}", flush=True)
        Ab = Bb
        La = Lb

    tot = sum(T.values())
    print(f"\n时间分解（GPU 版，到 n={N}，合计 {tot:.2f}s，总墙钟 {time.time()-t_all:.2f}s）")
    for kk, v in T.items():
        print(f"  {kk:10s} {v:8.3f}s  {100*v/tot:5.1f}%")
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "z0_record_gpu.json"), "w") as f:
        json.dump({"N": N, "rows": rows, "timing": T, "total": tot}, f, indent=2)
    print("saved results/z0_record_gpu.json")


if __name__ == "__main__":
    main()
