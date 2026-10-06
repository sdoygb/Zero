"""
l0_closure.py --- L0 闭链图（closure graph）的生成与测量。

**对象（本文件钉死的定义，可被核验）**

  节点  = 长度 T、恰含 T/2 个 `+1` 的 **二元项链**（binary necklace）= 零和 ±1 词的循环旋转类。
          计数由 Burnside 给出：  N(T) = (1/T) * Σ_{s=0}^{T-1} C(gcd(T,s), gcd(T,s)/2)
          实测精确复现文档值 10, 26, 80, 246, 810, 2704（T=8..18）。

  边    = **循环相邻对换**：位置 i 与 (i+1) mod T 的两位相异时交换之。
          这就是 Z1 定理 1 的补偿移动 x -> x + e_j - e_i（把 +1 右移一格、-1 左移一格）
          在零和词空间上的逐字实现，也是 `zero_sum_geometry_probe.py` 的 `+-` <-> `-+`。

  度    = **相异邻居的项链类个数**（同一类只计一次）。

**已知基准（文档 G28 §4 / Z0 §2.6）**：T=8..18 的平均度 3.60 4.92 6.15 7.39 8.43 9.50；
"无截断 => 平均度无界 => 非流形离散化 => 无稳定谱维数"。

本脚本用 GPU 把 T 推到文档 3.6 倍远处（T=32 有 1880 万节点），检验该 no-go 是否继续成立。
"""
from __future__ import annotations

import sys
import time
import numpy as np

from zcl import Engine, COMMON

# --------------------------------------------------------------- 组合计数
def necklace_count(T: int, k: int | None = None) -> int:
    """Burnside：长度 T、恰含 k 个 1 的二元项链数。"""
    from math import comb, gcd
    if k is None:
        k = T // 2
    tot = 0
    for s in range(T):
        g = gcd(T, s)
        if (g * k) % T:
            continue
        tot += comb(g, (g * k) // T)
    return tot // T


# ----------------------------------------------------- Phase A: 枚举项链集
def enumerate_necklaces(T: int, eng: Engine, chunk_bits: int = 24, verbose: bool = True):
    """
    枚举全部 N(T) 个项链，返回 (reps, stats)。
    reps: 排序后的 uint64 数组，每个元素是该项链的最小旋转（规范形）。

    做法：扫描所有 bit(T-1)=1 的平衡掩码（每个项链至少有一个旋转以 1 开头，
    故这一半空间必然覆盖全部项链），用 GPU 求最小旋转，分块去重。
    """
    lo = 1 << (T - 1)
    hi = 1 << T
    chunk = 1 << chunk_bits
    half = T // 2

    src = COMMON + r"""
    __kernel void unrank_scan(__global ulong* out, __global uint* cnt,
                              const ulong lo, const ulong hi, const int target) {
        size_t g = get_global_id(0);
        ulong m = lo + (ulong)g;
        if (m >= hi) return;
        if ((int)popcount(m) != target) return;
        size_t p = atomic_inc(cnt);
        out[p] = m;
    }
    """
    prog_key = "l0_enumerate"

    uniq_parts: list[np.ndarray] = []
    n_scanned = 0
    n_balanced = 0
    t0 = time.time()

    import pyopencl as cl

    # 缓冲区在循环外分配一次（大缓冲反复分配是纯浪费）
    cw = chunk_bits + 2  # 期望平衡词密度 ~ C(T-1,T/2-1)/2^(T-1) <= ~0.23
    MAXW = min(1 << cw, 1 << 27)
    cnt_buf = eng.buffer(4, readonly=False)
    outbuf = eng.empty(MAXW, np.uint64)
    zero = np.zeros(1, np.uint32)

    start = lo
    while start < hi:
        step = min(chunk, hi - start)
        cl.enqueue_copy(eng.queue, cnt_buf, zero).wait()

        k = eng.kernel(prog_key, src, "unrank_scan")
        eng.run(k, (step,), outbuf, cnt_buf, np.uint64(start), np.uint64(start + step),
                np.int32(half))
        eng.finish()
        cnt = int(eng.from_device(cnt_buf, 1, np.uint32)[0])
        if cnt == 0:
            start += step
            continue
        if cnt > MAXW:
            raise RuntimeError(f"缓冲区过小: cnt={cnt} > MAXW={MAXW}")

        words = eng.from_device(outbuf, cnt, np.uint64)[:cnt]
        n_scanned += step
        n_balanced += cnt

        # GPU 规范形（最小旋转）
        wb = eng.to_device(words)
        cb = eng.empty(cnt, np.uint64)
        kc = eng.kernel("l0_canon", COMMON, "canon_rot")
        mask = np.uint64((1 << T) - 1)
        eng.run(kc, (cnt,), wb, cb, np.int32(T), mask)
        canon = eng.from_device(cb, cnt, np.uint64)
        uniq_parts.append(np.unique(canon))
        start += step

        if verbose:
            el = time.time() - t0
            print(f"    scanned {n_scanned/1e6:8.1f}M masks | balanced {n_balanced/1e6:7.2f}M "
                  f"| uniq-so-far {sum(len(u) for u in uniq_parts)/1e6:6.3f}M | {el:6.1f}s")

    reps = np.unique(np.concatenate(uniq_parts)) if uniq_parts else np.zeros(0, np.uint64)
    stats = {
        "T": T,
        "masks_scanned": int(n_scanned),
        "balanced_words_found": int(n_balanced),
        "necklaces": int(len(reps)),
        "necklaces_formula": int(necklace_count(T)),
        "seconds": round(time.time() - t0, 2),
    }
    return reps, stats


# ------------------------------------------------------- Phase B: 建图 + 度
FUSED = COMMON + r"""
// 融合：对每个项链，算 T 个循环相邻对换，并各自取最小旋转（规范形）。
// 这样 T*T 次旋转全部留在 GPU 上，避免 T 次内核启动与主机往返。
__kernel void neighbors_canon(const __global ulong* in, __global ulong* out,
                              const int T, const ulong mask) {
    size_t g = get_global_id(0);
    ulong w = in[g] & mask;
    const int s = T - 1;
    __global ulong* o = out + g * (size_t)T;
    for (int i = 0; i < T; ++i) {
        int j = (i + 1) % T;
        ulong bi = (w >> i) & 1UL;
        ulong bj = (w >> j) & 1UL;
        ulong v = w;
        if (bi != bj) v ^= (1UL << i) | (1UL << j);
        ulong best = v, r = v;
        for (int q = 1; q < T; ++q) {
            r = ((r << 1) | (r >> s)) & mask;
            best = min(best, r);
        }
        o[i] = best;
    }
}
"""


def build_graph(T: int, reps: np.ndarray, eng: Engine, chunk: int = 1 << 19,
                verbose: bool = False):
    """
    返回 (degrees, n_edges)。
    度 = 相异邻居的项链类个数（排除自环）。n_edges = sum(deg)/2。
    """
    N = len(reps)
    mask = np.uint64((1 << T) - 1)
    kf = eng.kernel("l0_fused", FUSED, "neighbors_canon")

    degrees = np.zeros(N, np.int32)
    n_edges = 0
    t0 = time.time()
    for s in range(0, N, chunk):
        e = min(s + chunk, N)
        sub = np.ascontiguousarray(reps[s:e])
        m = e - s
        rb = eng.to_device(sub)
        cb = eng.empty(m * T, np.uint64)
        eng.run(kf, (m,), rb, cb, np.int32(T), mask)
        can = eng.from_device(cb, m * T, np.uint64).reshape(m, T)

        idx = np.searchsorted(reps, can)
        np.clip(idx, 0, N - 1, out=idx)
        valid = reps[idx] == can
        rows = np.arange(s, e, dtype=np.int64)[:, None]
        valid &= idx != rows                      # 排除自环
        big = np.where(valid, idx, np.int64(-1))
        big.sort(axis=1)                          # -1 排在最前
        deg = (big[:, 0] != -1).astype(np.int32)
        if T > 1:
            deg += (np.diff(big, axis=1) != 0).sum(axis=1).astype(np.int32)
        degrees[s:e] = deg
        n_edges += int(deg.sum())
        if verbose:
            print(f"    rows {e}/{N} | {time.time()-t0:6.1f}s")
    n_edges //= 2
    return degrees, n_edges


# ------------------------------------------------------------------ 主流程
KNOWN_AVG_DEG = {8: 3.60, 10: 4.92, 12: 6.15, 14: 7.39, 16: 8.43, 18: 9.50}
KNOWN_NODES = {8: 10, 10: 26, 12: 80, 14: 246, 16: 810, 18: 2704}


def run(T: int, eng: Engine, do_graph: bool = True):
    print(f"\n=== T={T} ===")
    reps, st = enumerate_necklaces(T, eng, verbose=False)
    ok = st["necklaces"] == st["necklaces_formula"]
    print(f"  necklaces = {st['necklaces']:,}  (formula {st['necklaces_formula']:,})  "
          f"{'OK' if ok else '*** MISMATCH ***'}  [{st['seconds']}s]")
    if T in KNOWN_NODES:
        exp = KNOWN_NODES[T]
        print(f"  vs 文档节点数 {exp}: {'OK' if exp == st['necklaces'] else '*** MISMATCH ***'}")
    out = {"T": T, **st}
    if do_graph:
        deg, ne = build_graph(T, reps, eng, verbose=False)
        avg = float(deg.mean())
        out.update({"avg_degree": round(avg, 4), "edges": int(ne),
                    "deg_min": int(deg.min()), "deg_max": int(deg.max()),
                    "deg_hist": np.bincount(deg).tolist()})
        tag = ""
        if T in KNOWN_AVG_DEG:
            tag = f"  vs 文档 {KNOWN_AVG_DEG[T]:.2f} -> {'OK' if abs(avg-KNOWN_AVG_DEG[T])<0.02 else '*** MISMATCH ***'}"
        print(f"  edges = {ne:,} | avg degree = {avg:.4f} (min {deg.min()}, max {deg.max()}){tag}")
    return reps, out


if __name__ == "__main__":
    eng = Engine()
    print("device:", eng.info())
    Ts = [int(a) for a in sys.argv[1:]] or [8, 10, 12, 14, 16, 18]
    results = []
    for T in Ts:
        _, o = run(T, eng)
        results.append(o)
    print("\n" + "=" * 78)
    print(f"{'T':>3} {'nodes':>14} {'edges':>16} {'avg deg':>9} {'deg_max':>8}")
    for o in results:
        print(f"{o['T']:>3} {o['necklaces']:>14,} {o.get('edges',0):>16,} "
              f"{o.get('avg_degree',0):>9.4f} {o.get('deg_max',0):>8}")
