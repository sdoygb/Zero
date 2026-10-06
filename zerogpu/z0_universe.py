"""
z0_universe.py --- **零参数宇宙生成器**（全 GPU 路径）

**为什么用最小系统**：因为这样**没有任何预设**。
语料的 L2 模型塞进了 8 个人工种子、寿命 L、摧毁周期 T、闭合吸收、容量 K——
只要塞了东西，跑出来的"结构"就分不清是算出来的还是放进去的，
**也就无法用来检验推导**。最小系统只有 Z0 三款：

    Z0① 零不停留  -> 闭合态不静止 => 回到零后继续走（"闭合即吸收"是额外规则，不用）
    Z0② 从不停歇  -> 永不终止
    Z0③ 不设概率  -> 全分支、无权重

    起点 = 零态 ∅（不是 8 个种子）
    一步 = alives <- {w+, w-}
    记录 = 词和为 0 者
    终止 = 从不

**产出与可观测量**（每个都对应一个物理量）

    活动层 |E(n)|          -> 宇宙"体积"（三维体积的离散类比）
    记录数 R(n)            -> "物质"沉积量
    区域（远足）普查        -> 成核数、尺寸谱
    闭合类普查             -> 可达态比例、手征量子数
    谱密度 rho(mu^2)        -> 是否存在极点（粒子）

**实现**：状态是 bitset（1 比特/词），一步 = 拷贝两份 + 与掩码；
规范化与手征判定都在 GPU 上。见 z0_record_gpu.py 的加速结论（canon 单项 310x）。
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from math import comb

import numpy as np

from zcl import Engine
import l0_closure as L0

# ------------------------------------------------------------------ 内核
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

# 手征伙伴：反序 + 变号（G35 的 \bar w）
REFLECT = r"""
__kernel void reflect_all(const __global ulong* in, __global ulong* out,
                          const int n, const ulong mask) {
    size_t g = get_global_id(0);
    ulong w = in[g] & mask, r = 0UL;
    for (int k = 0; k < n; ++k)
        r |= ((w >> (ulong)k) & 1UL) << (ulong)(n - 1 - k);
    out[g] = (~r) & mask;
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


def cat(k):
    return comb(2 * k, k) // (k + 1)


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    import pyopencl as cl
    eng = Engine()
    ctx, q = eng.ctx, eng.queue
    mf = cl.mem_flags
    print("device:", eng.info()["name"], flush=True)

    PAT = build_pat()
    pat_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=PAT)
    cnt_buf = cl.Buffer(ctx, mf.READ_WRITE, 4)
    k_step = eng.kernel("zu_step", STEP, "step_rec")
    k_comp = eng.kernel("zu_comp", COMPACT, "compact_bits")
    k_canon = eng.kernel("zu_canon", CANON, "canon_all")
    k_refl = eng.kernel("zu_refl", REFLECT, "reflect_all")

    A, _ = seed_depth6()
    Ab = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=A)
    La = 1
    rows = []
    t0 = time.time()
    for n in range(6, N):
        nn = n + 1
        Lb = La * 2
        Bb = cl.Buffer(ctx, mf.READ_WRITE, Lb * 8)
        Rb = cl.Buffer(ctx, mf.READ_WRITE, Lb * 8)
        cl.enqueue_copy(q, cnt_buf, np.zeros(1, np.uint32)).wait()
        target = (nn // 2) if nn % 2 == 0 else -1
        k_step(q, (Lb,), None, Ab, Bb, Rb, pat_buf, cnt_buf,
               np.uint32(La), np.int32(target)).wait()
        cnt = int(eng.from_device(cnt_buf, 1, np.uint32)[0])

        row = {"n": nn, "volume": 2 ** nn, "records": cnt}
        if cnt:
            wb = cl.Buffer(ctx, mf.READ_WRITE, cnt * 8)
            cl.enqueue_copy(q, cnt_buf, np.zeros(1, np.uint32)).wait()
            k_comp(q, (Lb,), None, Rb, wb, cnt_buf).wait()
            got = int(eng.from_device(cnt_buf, 1, np.uint32)[0])
            cb = cl.Buffer(ctx, mf.READ_WRITE, max(got, 1) * 8)
            rmask = np.uint64((1 << nn) - 1)
            k_canon(q, (got,), None, wb, cb, np.int32(nn), rmask).wait()
            can = eng.from_device(cb, got, np.uint64)
            uniq, counts = np.unique(can, return_counts=True)
            row.update({"classes": int(uniq.size),
                        "K_n": int(L0.necklace_count(nn, nn // 2)),
                        "coverage": round(uniq.size / L0.necklace_count(nn, nn // 2), 6),
                        "q_class": round(float(((counts / got) ** 2).sum()), 10)})
            # 手征：对全部事件算 reflect 后 canonical，比较
            if got <= (1 << 25):
                rb = cl.Buffer(ctx, mf.READ_WRITE, got * 8)
                k_refl(q, (got,), None, cb, rb, np.int32(nn), rmask).wait()
                cb2 = cl.Buffer(ctx, mf.READ_WRITE, got * 8)
                k_canon(q, (got,), None, rb, cb2, np.int32(nn), rmask).wait()
                refl = eng.from_device(cb2, got, np.uint64)
                row["self_mirror_events"] = int(np.sum(refl == can))
                row["chiral_events"] = int(np.sum(refl != can))
        rows.append(row)
        if cnt:
            print(f"  n={nn:3d} 体积={2**nn:>13,} 记录={cnt:>12,} "
                  f"类={row.get('classes',0):>11,} 覆盖={row.get('coverage',0):7.4f}"
                  + (f" 手征事件={row.get('chiral_events',0)/max(cnt,1):.3f}"
                     if 'chiral_events' in row else ""), flush=True)
        Ab = Bb
        La = Lb

    print(f"\n总耗时 {time.time()-t0:.2f}s")
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "z0_universe.json"), "w") as f:
        json.dump({"N": N, "rows": rows}, f, indent=2)
    print("saved results/z0_universe.json")

    # ---------------- 宇宙账本 ----------------
    print("\n" + "=" * 100)
    print("宇宙账本（零参数，从 ∅ 出发）")
    print("=" * 100)
    print(f"{'年龄 n':>7} {'体积 2^n':>16} {'记录 R(n)':>14} {'R/体积':>10} "
          f"{'区域数(累计)':>14} {'可达类':>12} {'覆盖率':>9}")
    cum_reg = 0
    for r in rows:
        n = r["n"]
        if n % 2:
            continue
        cum_reg += r["records"]
        print(f"{n:>7} {r['volume']:>16,} {r['records']:>14,} "
              f"{r['records']/r['volume']:>10.3e} {cum_reg:>14,} "
              f"{r.get('classes',0):>12,} {r.get('coverage',0):>9.4f}")
    print()
    print("对照闭式：")
    print(f"  记录数（吸收读法） 2·Cat(n/2-1) ; （不吸收读法）C(n,n/2)")
    print(f"  区域数（远足）     平均 ~ sqrt(pi*m), m=n/2")
    print(f"  临界性             E(1/4)=1  =>  平均后代数 = 1")


if __name__ == "__main__":
    main()
