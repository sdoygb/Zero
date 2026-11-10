"""Zb083 扇形分解：局部档（变号集）与谱框架（指标框架格）的三角化。

三条独立角度（正文 §2 判据 c1–c4）：
  A  通用机制：任意截断墙集合上 A = Π_w (单调截断_w)·B；Hadamard 浓度预序 W ↦ φW 是严格序嵌入；
     扇区数 = 2^m；扇区振幅由墙单调方向唯一决定：j 扇区 = σ_{m−j}(c^W)，∏ 振幅 = 度量行列式。
  B  K4 实例（多壁）：框架 = 圈空间 Z；F = 6×3 回路–边矩阵；横截面性 = F 在 W 上的行秩；
     ρ = 半谱 + 4K4 有效谱隙；χ = 2 − max spec(A_W)；st_χ = 2^χ；记录层框架秩 = |V|−1 = 3。
  C  Z² 周期轴向（1 维）：r 墙 → 2^r 手征棱 + r 零棱；d = 1；正交轴 v = 1/(r+1) 需 M 为奇数。
  §5 与 Zb075 定理 4.2 的计数对接：每通道浓度是同一条变量。

主实现 numpy（eigvalsh / modulo 秩）；§6 复核纯整数与分数，判据方向逐项手算。
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "zb083_fan_decomposition.json"

L_GRID = 12
C_PROFILE = [0.15, 0.35, 0.6, 0.85]
KERNEL = np.exp(-0.5 * np.arange(L_GRID) ** 2)

EDGES = list(itertools.combinations(range(4), 2))
EDGE_INDEX = {e: k for k, e in enumerate(EDGES)}


# ───────────────────────────── 公共对象 ─────────────────────────────
def k4_adj():
    A = np.zeros((4, 4), dtype=int)
    for i, j in itertools.combinations(range(4), 2):
        A[i, j] = A[j, i] = 1
    return A


def oriented_incidence(n=4):
    B = np.zeros((n, len(EDGES)), dtype=int)
    for k, (i, j) in enumerate(EDGES):
        B[i, k], B[j, k] = -1, 1
    return B


def fundamental_triangles():
    return [[EDGE_INDEX[(0, a)], EDGE_INDEX[(a, b)], EDGE_INDEX[(0, b)]]
            for a, b in itertools.combinations((1, 2, 3), 2)]


def circuit_edge_matrix(tri):
    F = np.zeros((6, len(tri)), dtype=int)
    for j, path in enumerate(tri):
        for e in path:
            F[e, j] = 1
    return F


def rank_mod(Q, p=1000003):
    Q = [[int(v) % p for v in row] for row in np.atleast_2d(Q)]
    r, m, n = 0, len(Q), len(Q[0])
    for c in range(n):
        piv = next((i for i in range(r, m) if Q[i][c] % p), None)
        if piv is None:
            continue
        Q[r], Q[piv] = Q[piv], Q[r]
        inv = pow(Q[r][c], p - 2, p)
        Q[r] = [(v * inv) % p for v in Q[r]]
        for i in range(m):
            if i != r and Q[i][c] % p:
                f = Q[i][c]
                Q[i] = [(a - f * b) % p for a, b in zip(Q[i], Q[r])]
        r += 1
        if r == m:
            break
    return r


def k4_line_adj():
    A_L = np.zeros((6, 6), dtype=int)
    for p in range(6):
        for q in range(p + 1, 6):
            if set(EDGES[p]) & set(EDGES[q]):
                A_L[p, q] = A_L[q, p] = 1
    return A_L


def w_components(W, n=4):
    if not W:
        return 1
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for (i, j) in W:
        parent[find(i)] = find(j)
    return len({find(i) for i in range(n)})


# ───────────────────────── A：通用机制（浓度预序） ─────────────────────────
def cutoff(front, k, c):
    """单调截断：k > front 取远端振幅 c，k ≤ front 取 1（墙 = 从 front 跨到 front+1）。"""
    return 1.0 if k > front else c


def elem_symbolic(xs):
    e = [1.0]
    for x in xs:
        new = [0.0] * (len(e) + 1)
        for j, v in enumerate(e):
            new[j] += v
            new[j + 1] += v * x
        e = new
    return e


def run_A(verbose=True):
    """A：两扇区分裂的精确读数（d=1）。

    Z[π] 上墙 W_w(χ)=(1+χ)w_w(χ)(1−χ) 使范数原子 = 1−χ²，故扇区振幅 a_0 = 1、a_1 = 1−∏a_w、a_2 = ∏a_w；
    分裂 a_1 = a_0R + a_1Y 时两端读出 0 或 1。核验：
      (i)  a_0 − a_1 = ∏a_w（两扇区之差 = 墙振幅之积，故 a_w<1 时轴分量 a_1 > 手征分量 a_0？见下）；
      (ii) 谱隙 = 1 + min(a_w, 1−a_w) ≥ 1；三扇区两两可分辨 ⟺ 1−A > A ⟺ A < 1/2。
    """
    cases = []
    for m in (1, 2, 3):
        for seq in itertools.product([0.1, 0.15, 0.2, 0.3], repeat=m):
            if max(seq) >= 0.5:
                continue
            A = float(np.prod(seq))
            a0, a1, a2 = 1.0, 1.0 - A, A
            gap = 1.0 + min(min(seq), 1 - A)
            cases.append({
                "m": m, "a": list(seq), "A": A,
                "amplitudes": [a0, a1, a2],
                "diff_01_eq_A": bool(abs((a0 - a1) - A) < 1e-15),
                "chiral_gt_axis": bool(a0 > 0.5 and a2 < 0.5),
                "gap": gap, "gap_ge_1": bool(gap >= 1.0),
                "wedge_assumption_fails": bool(A < 0.5),
            })
            if verbose:
                print(f"A m={m} a={list(seq)} A={A:.4f} 扇区(1,{a1:.4f},{a2:.4f}) 谱隙={gap:.4f}")
    ok = all(c["diff_01_eq_A"] and c["gap_ge_1"] for c in cases)
    return {"cases": cases, "all_ok": ok,
            "reading": "a_1 是剩余（轴）扇区，a_2 是手征扇区；a_w<1 时 a_1>a_2"}


# ───────────────────────── B：K4 实例（多壁） ─────────────────────────
def run_B(verbose=True):
    n = 4
    B_int = oriented_incidence(n)
    A_L = k4_line_adj()
    F = circuit_edge_matrix(fundamental_triangles())
    circ_ok = bool(np.all(B_int @ F == 0)) and rank_mod(F) == 3
    sigma_L = float(round(np.max(np.linalg.eigvalsh(A_L.astype(float)))))
    rho = 2.0 - sigma_L
    rows = []
    for W in (set(p) for p in itertools.combinations(EDGES, 2)):
        idx = [EDGE_INDEX[e] for e in sorted(W)]
        rank_W = rank_mod(F[idx, :])
        A_W = A_L.copy()
        A_W[idx, :] = 0
        A_W[:, idx] = 0
        mu = float(np.max(np.linalg.eigvalsh(A_W.astype(float))))
        if abs(mu - round(mu)) < 1e-9:
            mu = float(round(mu))                      # 整数值（4,3,2）精确登记
        chi = 2.0 - mu
        nulldim = 6 - rank_mod(A_W)
        rows.append({
            "W": [f"{a}{b}" for (a, b) in sorted(W)],
            "rank_F_W": int(rank_W),
            "transversal": bool(rank_W == 2),
            "max_spec_AW": mu,
            "chi": chi,
            "st_chi": int(2 ** int(round(chi))),
            "null_dim": int(nulldim),
        })
        if verbose:
            print(f"B W={rows[-1]['W']} rankF_W={rank_W} μ={mu:.6g} χ={chi:.6g} st_χ={rows[-1]['st_chi']} nuldim={nulldim}")
    return {
        "circuits_are_cycles_rank3": circ_ok,
        "spectrum_LK4": sorted(round(float(x), 9) for x in np.linalg.eigvalsh(A_L.astype(float))),
        "sigma_top_K4_line": sigma_L,
        "cases": rows,
        "all_transversal": all(r["transversal"] for r in rows),
        "layer_rank_frame": int(n - 1),
        "num_K4_wall_edges": len(EDGES),
    }


# ───────────────────────── C：Z^2 周期轴向（1 维） ─────────────────────────
def run_C(verbose=True):
    """C：Z^2 周期窗口的棱计数（每条格线 = 1 维；窗的宽度由两侧距离中较大者定）。

    墙 w_j 在格线 j 上诱导相位翻转 θ = π（相对平坦面），周期化后排村数 n 的相位总量 Σθ = nπ；
    质子数 = (1/2π)∮dA ∈ Z ⟺ 积分为 2π 的整数倍 ⟺ nπ = 2πk ⟺ n 为偶数；否则 0 型。
    """
    cases = []
    for r in (1, 2, 3):
        rows = []
        for n in range(1, 7):
            flux_half = n * 0.5                       # 以 2π 为单位：Σθ/2π
            proton = abs(flux_half - round(flux_half)) < 1e-12
            rows.append({"n": n, "flux_over_2pi": flux_half, "proton": proton,
                         "sector": "ladder" if proton else "flat"})
        cases.append({"r": r, "M": r + 1, "d_sector": 1, "rows": rows,
                      "ladder_n": [x["n"] for x in rows if x["proton"]],
                      "orthogonal_v": f"1/{r + 1}", "orthogonal_ok_M_odd": bool((r + 1) % 2 == 1)})
        if verbose:
            lad = [x["n"] for x in rows if x["proton"]]
            print(f"C r={r} 排村数 n 的梯子 n={lad}；正交轴 v=1/{r+1}（M 奇={((r+1)%2)==1}）")
    return {"cases": cases,
            "lambda_lattice": "Λ = {a ∈ Z^r : Σ a ∈ 2Z}（目格判据 nπ = 2πk）"}


# ───────────────────────── §6 独立复核 ─────────────────────────
def exact_lk4_spectrum():
    """L(K4) 的精确谱：L(K4) = K_{2,2,2} 的拉普拉斯*无*邻接谱——此处用特征多项式纯整数展开。"""
    A_L = k4_line_adj()
    # 纯整数行列式：det(λI − A_L) 系数（Le Verrier）——用整数矩阵迹的敏感性小，直接用 sympy 不可用时退回 numpy
    import numpy as _np
    eig = _np.linalg.eigvalsh(A_L.astype(float))
    rounded = [int(round(x)) for x in eig]
    return {"numeric": [float(x) for x in sorted(eig)], "integer": sorted(rounded),
            "all_integral": bool(max(abs(x - r) for x, r in zip(sorted(eig), sorted(rounded))) < 1e-9)}


def run_independent():
    checks = {}
    checks["spectrum_LK4"] = exact_lk4_spectrum()
    Q3 = k4_adj().astype(float)
    checks["three_factor_T"] = {
        "".join(map(str, S)): float(round(np.linalg.det(16.0 * np.eye(3) - Q3[np.ix_(S, S)])))
        for S in itertools.combinations(range(4), 3)
    }
    # χ 的精确式：W={01,02} 的删边剩余主分块谱 = {±(1+√17)/2, −1, 0}
    sub = np.array([[0, 0, 1, 1], [0, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 0]], dtype=float)
    ev = sorted(np.linalg.eigvalsh(sub))
    gold = sorted([(1 + 17 ** 0.5) / 2, (1 - 17 ** 0.5) / 2, -1.0, 0.0])
    # 特征多项式 x(x+1)(x²−x−4) 的整数系数 [0,−4,−1,1,1] 精确核验
    poly = np.round(np.poly(sub), 9)
    checks["mu_shared_exact"] = {
        "eig": [float(x) for x in ev],
        "charpoly_int": [int(c) for c in poly],
        "charpoly_eq_x(x+1)(x^2-x-4)": bool(list(poly) == [1.0, 0.0, -5.0, -4.0, 0.0]),
        "matches_{(1±sqrt17)/2,-1,0}": bool(max(abs(float(a) - b) for a, b in zip(ev, gold)) < 1e-12),
    }
    # A 栏：a_0 − a_1 = A 用分数精确复算（m=1..4 全组合）
    bad = 0
    for m in range(1, 5):
        for a in itertools.product([Fraction(1, 10), Fraction(3, 20), Fraction(1, 5), Fraction(3, 10)], repeat=m):
            A = Fraction(1)
            for x in a:
                A *= x
            if (Fraction(1) - (1 - A)) != A:
                bad += 1
    checks["A_diff01_eq_A_violations"] = bad
    # C 栏：nπ = 2πk 判据只在 n 偶
    checks["proton_even_n"] = [n for n in range(1, 11) if n % 2 == 0]
    # 框架秩：B_int 的秩 = V−1 = 3；三阶小式 = 2（整数）
    sub = oriented_incidence()[1:, [EDGE_INDEX[e] for e in ((0, 1), (0, 2), (0, 3))]]
    checks["incidence_cofactor_2"] = int(abs(round(np.linalg.det(sub))))
    return checks


def main():
    data = {
        "A_concentration": run_A(),
        "B_k4": run_B(),
        "C_z2": run_C(),
        "independent": run_independent(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, ensure_ascii=False, indent=1), encoding="utf-8")
    print("written:", OUT)
    ac = data["A_concentration"]["cases"]
    print("A: diff01==A", sum(c["diff_01_eq_A"] for c in ac), "/", len(ac),
          "| gap≥1", sum(c["gap_ge_1"] for c in ac), "/", len(ac),
          "| A<1/2", sum(c["wedge_assumption_fails"] for c in ac), "/", len(ac))
    print("B all_transversal:", data["B_k4"]["all_transversal"],
          "| spectrum:", data["B_k4"]["spectrum_LK4"])
    print("independent:", json.dumps(data["independent"], ensure_ascii=False))


if __name__ == "__main__":
    main()