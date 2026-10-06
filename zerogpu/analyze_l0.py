"""
analyze_l0.py --- L0 层的三项独立核验与扩展。

§A 双图对照
    语料里 "L0 内容图" 有两个**不同**的实现，此前未被并列比较：
      (i)  E1 词图  ：节点 = 长度 L 的全部平衡 ±1 **词**（C(L,L/2) 个）；
                      边 = **线性**相邻异号对换 i=0..L-2（不含绕回）。
                      出处：E1_NN_check.py:51-78（`moves()` 用 range(len(w)-1)）。
      (ii) 闭链图 ：节点 = 长度 T 的平衡 **项链**（旋转类，Burnside 计数）；
                      边 = **循环**相邻异号对换 (i, i+1 mod T)。
                      出处：simulations/zero_sum_geometry_probe.py:103-104。
    LAYER_LEDGER §3.1 把"1 个零模；λ1≈π²/L²；λmax→4L/3"列为 **L0** 事实，
    但未说明是哪张图。本节逐一实测，指出该谱律属于 (i)，而 (ii) 不同。

§B 类权重 q_L 与 DIM-COST-Q
    R87 H2： ω_C = o_C / C(L,L/2)（o_C = 类的轨道大小）， q_L = Σ_C ω_C²。
    文档只给到 L=4,6,8（5/9, 7/25, 133/1225）。本节用 GPU 枚举把 L 推到 32。

§C 首次通过桥 -> 循环类 的单射性
    L2 全分支动力学的闭合事件（见 l2_branch.py）推到长度 n 的项链类上。
    实测：在站点 '+' 上，n<=32 时**事件数 = 相异类数**（映射是单射）。
"""
from __future__ import annotations

import json
import os
import sys
import time
from math import comb, gcd
from collections import defaultdict

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from zcl import Engine
import l0_closure as L0


# =============================================================== §A 双图对照
def e1_word_graph(L: int):
    """E1 的词图：平衡词 + 线性相邻对换。返回 (n_nodes, rows, cols)。"""
    from itertools import combinations
    words = []
    for pos in combinations(range(L), L // 2):
        w = 0
        for i in pos:
            w |= 1 << i
        words.append(w)
    idx = {w: i for i, w in enumerate(words)}
    rows, cols = [], []
    for w in words:
        i = idx[w]
        for k in range(L - 1):                      # 注意：range(L-1)，无绕回
            bi = (w >> k) & 1
            bj = (w >> (k + 1)) & 1
            if bi != bj:
                v = w ^ ((1 << k) | (1 << (k + 1)))
                rows.append(i)
                cols.append(idx[v])
    return len(words), np.array(rows), np.array(cols)


def cyc_necklace_graph(T: int, reps: np.ndarray, eng: Engine, chunk: int = 1 << 19):
    """闭链图：项链 + 循环对换。返回 (n_nodes, rows, cols)。"""
    import pyopencl as cl
    N = len(reps)
    mask = np.uint64((1 << T) - 1)
    kf = eng.kernel("l0_fused", L0.FUSED, "neighbors_canon")
    rows_all, cols_all = [], []
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
        valid = (reps[idx] == can) & (idx != np.arange(s, e, dtype=np.int64)[:, None])
        rr = np.repeat(np.arange(s, e, dtype=np.int64), T).reshape(m, T)[valid]
        cc = idx[valid]
        rows_all.append(rr)
        cols_all.append(cc)
    return N, np.concatenate(rows_all), np.concatenate(cols_all)


def spectrum(n, rows, cols, k_extreme: int = 4):
    """返回 (零模数, lambda1, lambda_max)。用稀疏 Laplacian。"""
    data = np.ones(len(rows))
    A = sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    A = A.maximum(A.T)                      # 简单图
    deg = np.asarray(A.sum(axis=1)).ravel()
    Lap = sp.diags(deg) - A
    Lap = Lap.tocsr()
    # 最大特征值
    lmax = float(spla.eigsh(Lap, k=1, which="LA", return_eigenvectors=False)[0])
    # 最小特征值（0）与其上的第一个
    lo = spla.eigsh(Lap, k=min(k_extreme, n - 2), sigma=-1e-6, which="LM",
                    return_eigenvectors=False)
    lo = np.sort(lo)
    n_zero = int(np.sum(np.abs(lo) < 1e-6))
    l1 = float(lo[n_zero]) if len(lo) > n_zero else float("nan")
    return n_zero, l1, lmax, float(deg.mean())


def section_A(eng: Engine):
    print("\n" + "=" * 78)
    print("§A 双图对照：E1 词图（线性对换） vs 闭链图（循环对换）")
    print("=" * 78)
    out = {"e1_word_graph": [], "cyc_necklace_graph": []}

    print(f"\n{'L':>3} | {'E1词图 N':>10} {'边':>10} {'均度':>7} {'零模':>4} "
          f"{'λ1':>10} {'λ1·L²/π²':>10} {'λmax':>8} {'λmax/L':>7}")
    for L in range(4, 21, 2):
        n, r, c = e1_word_graph(L)
        nz, l1, lmax, dm = spectrum(n, r, c)
        ratio = l1 * L * L / (np.pi ** 2)
        out["e1_word_graph"].append(
            {"L": L, "N": n, "edges": len(r) // 2, "mean_deg": round(dm, 4),
             "zero_modes": nz, "lambda1": l1, "lambda1_over_pi2_L2": round(ratio, 6),
             "lambda_max": round(lmax, 6), "lambda_max_over_L": round(lmax / L, 6)})
        print(f"{L:>3} | {n:>10,} {len(r)//2:>10,} {dm:>7.4f} {nz:>4} "
              f"{l1:>10.6f} {ratio:>10.4f} {lmax:>8.4f} {lmax/L:>7.4f}")

    print(f"\n{'T':>3} | {'闭链图 N':>10} {'边':>10} {'均度':>7} {'零模':>4} "
          f"{'λ1':>10} {'λ1·T²/π²':>10} {'λmax':>8} {'λmax/T':>7}")
    for T in range(4, 25, 2):
        reps, _ = L0.enumerate_necklaces(T, eng, verbose=False)
        n, r, c = cyc_necklace_graph(T, reps, eng)
        if n < 6:
            continue
        nz, l1, lmax, dm = spectrum(n, r, c)
        ratio = l1 * T * T / (np.pi ** 2) if l1 == l1 else float("nan")
        out["cyc_necklace_graph"].append(
            {"T": T, "N": n, "edges": len(r) // 2, "mean_deg": round(dm, 4),
             "zero_modes": nz, "lambda1": l1,
             "lambda1_over_pi2_T2": round(float(ratio), 6) if ratio == ratio else None,
             "lambda_max": round(lmax, 6), "lambda_max_over_T": round(lmax / T, 6)})
        print(f"{T:>3} | {n:>10,} {len(r)//2:>10,} {dm:>7.4f} {nz:>4} "
              f"{l1:>10.6f} {ratio:>10.4f} {lmax:>8.4f} {lmax/T:>7.4f}")
    return out


# ====================================================== §B 类权重 q_L (大 L)
def minimal_period(rep: int, L: int) -> int:
    """项链规范形 rep 的最小周期 p（= 轨道大小 o_C）。"""
    mask = (1 << L) - 1
    for p in range(1, L + 1):
        if L % p:
            continue
        r = ((rep << p) | (rep >> (L - p))) & mask
        if r == rep:
            return p
    return L


def section_B(eng: Engine, Lmax: int = 32):
    print("\n" + "=" * 78)
    print("§B 类权重 q_L = Σ_C (o_C/N_L)²（R87 H2 的推广，文档只到 L=8）")
    print("=" * 78)
    print(f"{'L':>3} {'类数':>12} {'N_L=C(L,L/2)':>16} {'q_L':>14} {'L/N_L':>12} "
          f"{'比值':>7} {'最大类权':>9}")
    rows = []
    for L in range(2, Lmax + 1, 2):
        reps, _ = L0.enumerate_necklaces(L, eng, verbose=False)
        N = comb(L, L // 2)
        sizes = np.array([minimal_period(int(r), L) for r in reps], dtype=np.int64)
        assert sizes.sum() == N, f"L={L}: Σo_C={sizes.sum()} != N={N}"
        w = sizes / N
        q = float((w ** 2).sum())
        rows.append({"L": L, "classes": int(len(reps)), "N_L": N, "q_L": q,
                     "L_over_N_L": L / N, "ratio": q / (L / N) if L / N else None,
                     "max_class_weight": float(w.max())})
        print(f"{L:>3} {len(reps):>12,} {N:>16,} {q:>14.10g} {L/N:>12.6g} "
              f"{q/(L/N):>7.4f} {w.max():>9.6f}")
    return rows


# ================================================ §C 首次通过桥的单射性核验
def section_C(eng: Engine, nmax: int = 30):
    print("\n" + "=" * 78)
    print("§C 首次通过闭合词 -> 长度 n 项链类：是否单射？")
    print("=" * 78)
    import l2_branch as L2
    print(f"{'n':>3} {'事件数':>12} {'相异类数':>12} {'长度n项链总数':>14} {'单射?':>7} "
          f"{'最大重数':>8}")
    rows = []
    for n in range(4, nmax + 1, 2):
        hits = L2.site_closure_events(eng, "+", n)
        if len(hits) == 0:
            continue
        uniq, cnts = np.unique(hits, return_counts=True)
        nclass = L0.necklace_count(n, n // 2)
        inj = len(uniq) == len(hits)
        rows.append({"n": n, "events": int(len(hits)), "distinct": int(len(uniq)),
                     "necklaces_of_len_n": int(nclass), "injective": bool(inj),
                     "max_multiplicity": int(cnts.max())})
        print(f"{n:>3} {len(hits):>12,} {len(uniq):>12,} {nclass:>14,} "
              f"{'是' if inj else '否':>7} {cnts.max():>8}")
    return rows


def main():
    eng = Engine()
    print("device:", eng.info())
    Lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    res = {"device": eng.info()}
    res["A_dual_graph"] = section_A(eng)
    res["B_class_weight_q"] = section_B(eng, Lmax)
    res["C_injectivity"] = section_C(eng, min(Lmax, 30))
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "l0_analysis.json"), "w") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)
    print("\nsaved results/l0_analysis.json")


if __name__ == "__main__":
    main()
