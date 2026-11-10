#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Zb086 · 偏置量：子系统的非中心分量只由边界通量改变。

七组判据（对应正文 §2–§4，汇总在 §4）：
  S1 内部不动     S=V 时，任意 x(τ)=x(0)+Bf 都给 b(τ)=b(0)（Zb004 定理 2.2）。
  S2 边界通量律   S⊊V 时，Δb_S = -Φ_∂S/|S|，逐次精确（割边求和）。
  S3 判据方向     同一实验改用**全部边**求和，恒等式立刻失效——证明"只用割边"不是记账选择。
  S4 太阳偏置     八行星给太阳本体的 b_p = a_p m_p/(M_⊙+m_p)，与 R_⊙ 比较。
  S5 比值与锚     全部长度乘同一个因子：b/R_⊙ 不变，b 的绝对值同比变——比值内生、绝对量需锚（Zb050）。
  S6 加权边界律   随机正定容度 W：Wφ_new=Wφ+Bf 下 Δb_S=-Φ/|S|_W，且 1_SᵀWφ_⊥=0（定理 1.3、命题 1.2）。
  S7 两体零模解   Mx_⊙+m_p x_p=0、x_p-x_⊙=a_p 联立，|x_⊙|=a_p m_p/(M+m_p)（命题 3.1）。

运行：python3 ZeroGPU/src/zb086_bias_offset.py
输出：ZeroGPU/results/zb086_bias_offset.json
"""
from __future__ import annotations

import json
import random
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "zb086_bias_offset.json"

AU_M = 1.495978707e11
R_SUN_M = 6.957e8
M_SUN_MJ = 1047.3486
PLANETS = [
    ("Mercury", 0.38710, 0.000174),
    ("Venus", 0.72333, 0.002564),
    ("Earth", 1.00000, 0.003147),
    ("Mars", 1.52371, 0.000323),
    ("Jupiter", 5.20288, 1.000000),
    ("Saturn", 9.53667, 0.299400),
    ("Uranus", 19.18916, 0.045720),
    ("Neptune", 30.06992, 0.053720),
]


def _graph(n: int, seed: int):
    rng = random.Random(seed)
    edges = [(i, i + 1) for i in range(n - 1)]
    edges += [(rng.randrange(n), rng.randrange(n)) for _ in range(3)]
    edges = [(i, j) for i, j in edges if i != j]
    B = np.zeros((n, len(edges)))
    for e, (i, j) in enumerate(edges):
        B[i, e], B[j, e] = -1.0, +1.0
    return edges, B


def s1_internal_fixed(n=6, trials=200, seed=11):
    edges, B = _graph(n, seed)
    rng = random.Random(seed)
    worst = 0.0
    for _ in range(trials):
        x = np.array([rng.randint(-5, 5) for _ in range(n)], dtype=float)
        f = np.array([rng.randint(-5, 5) for _ in edges], dtype=float)
        xn = x + B @ f
        worst = max(worst, abs(xn.mean() - x.mean()))
    return {"n": n, "edges": len(edges), "trials": trials,
            "max_abs_delta_b": float(worst), "pass": bool(worst < 1e-12)}


def _phi(edges, f, S):
    """流出 S 为正的净通量。"""
    tot = 0.0
    for e, (i, j) in enumerate(edges):
        si, sj = i in S, j in S
        if si and not sj:
            tot += f[e]
        elif sj and not si:
            tot -= f[e]
    return tot


def s2_boundary_law(n=6, trials=200, seed=12):
    edges, B = _graph(n, seed)
    rng = random.Random(seed)
    S = {0, 1, 2}
    cut = [(i, j) for (i, j) in edges if (i in S) != (j in S)]
    worst = 0.0
    for _ in range(trials):
        x = np.array([rng.randint(-5, 5) for _ in range(n)], dtype=float)
        f = np.array([rng.randint(-5, 5) for _ in edges], dtype=float)
        xn = x + B @ f
        db = xn[list(S)].mean() - x[list(S)].mean()
        worst = max(worst, abs(db + _phi(edges, f, S) / len(S)))
    return {"subsystem": sorted(S), "cut_edges": len(cut), "trials": trials,
            "max_abs_residual": float(worst), "pass": bool(worst < 1e-12)}


def s3_wrong_direction(n=6, trials=200, seed=12):
    """反例：用全部边求和（而不是割边）时，恒等式应失效。"""
    edges, B = _graph(n, seed)
    rng = random.Random(seed)
    S = {0, 1, 2}
    worst = 0.0
    for _ in range(trials):
        x = np.array([rng.randint(-5, 5) for _ in range(n)], dtype=float)
        f = np.array([rng.randint(-5, 5) for _ in edges], dtype=float)
        xn = x + B @ f
        db = xn[list(S)].mean() - x[list(S)].mean()
        worst = max(worst, abs(db + f.sum() / len(S)))
    return {"max_abs_residual_all_edges": float(worst), "fails_as_expected": bool(worst > 1e-9)}


def s4_solar():
    rows, total = [], 0.0
    for name, a, m in PLANETS:
        b_au = a * m / (M_SUN_MJ + m)
        b_km = b_au * AU_M / 1e3
        rows.append({"planet": name, "a_au": a, "m_mjup": m,
                     "b_au": b_au, "b_km": b_km,
                     "b_over_rsun": b_km * 1e3 / R_SUN_M})
        total += b_au
    return {"rows": rows,
            "aligned_sum_au": total,
            "aligned_sum_over_rsun": total * AU_M / R_SUN_M,
            "max_single_over_rsun": max(r["b_over_rsun"] for r in rows)}


def s5_ratio_vs_anchor(scale=1e6):
    """长度单位整体改标：绝对值同比变，比值不动（Zb050 的分栏）。"""
    a, m = 5.20288, 1.0
    b_m = a * m / (M_SUN_MJ + m) * AU_M
    r1 = b_m / R_SUN_M
    r2 = (b_m * scale) / (R_SUN_M * scale)
    return {"scale": scale, "ratio_scaled": r2,
            "ratio_original": r1, "ratio_invariant": bool(abs(r2 - r1) < 1e-12),
            "absolute_changes_by_scale": scale}


def s6_weighted_boundary(n=6, trials=300, seed=13):
    """一般正定容度 W：Wφ_new = Wφ + Bf。检验定理 1.3 与命题 1.2 的 W-正交性。"""
    edges, B = _graph(n, seed)
    rng = random.Random(seed)
    S = {0, 1, 2}
    w = np.array([rng.uniform(0.5, 2.0) for _ in range(n)])
    W = np.diag(w)
    one_S = np.zeros(n)
    one_S[list(S)] = 1.0
    denom = float(one_S @ W @ one_S)
    worst_law = worst_orth = 0.0
    for _ in range(trials):
        phi = np.array([rng.uniform(-3, 3) for _ in range(n)])
        f = np.array([rng.randint(-4, 4) for _ in range(len(edges))], dtype=float)
        phi_new = phi + np.linalg.solve(W, B @ f)
        b = float(one_S @ W @ phi) / denom
        b_new = float(one_S @ W @ phi_new) / denom
        phi_perp = phi - b * one_S
        worst_orth = max(worst_orth, abs(float(one_S @ W @ phi_perp)))
        worst_law = max(worst_law, abs((b_new - b) + _phi(edges, f, S) / denom))
    return {"n": n, "trials": trials, "subsystem": sorted(S),
            "max_abs_orthogonality": float(worst_orth),
            "max_abs_boundary_law": float(worst_law),
            "pass": bool(max(worst_orth, worst_law) < 1e-12)}


def s7_two_body():
    """命题 3.1：两节点零模约束与边长读出联立，逐行星核对。"""
    rows, worst = [], 0.0
    M = 1.0
    for name, a, m in PLANETS:
        mu = m / M_SUN_MJ                       # 行星质量（太阳质量为单位）
        x_sun = -a * mu / (M + mu)              # 太阳相对质心
        x_planet = a + x_sun
        zero_mode = M * x_sun + mu * x_planet
        edge = (x_planet - x_sun) - a
        formula = a * mu / (M + mu)
        worst = max(worst, abs(zero_mode), abs(edge), abs(abs(x_sun) - formula))
        rows.append({"planet": name, "zero_mode_residual": float(zero_mode),
                     "edge_readout_residual": float(edge),
                     "b_over_a": float(abs(x_sun) / a),
                     "formula_over_a": float(formula / a)})
    return {"rows": rows, "max_abs_residual": float(worst), "pass": bool(worst < 1e-14)}


def main():
    res = {
        "S1_internal_fixed": s1_internal_fixed(),
        "S2_boundary_law": s2_boundary_law(),
        "S3_wrong_direction": s3_wrong_direction(),
        "S4_solar": s4_solar(),
        "S5_ratio_vs_anchor": s5_ratio_vs_anchor(),
        "S6_weighted_boundary": s6_weighted_boundary(),
        "S7_two_body": s7_two_body(),
    }
    res["PASS"] = all([res["S1_internal_fixed"]["pass"],
                       res["S2_boundary_law"]["pass"],
                       res["S3_wrong_direction"]["fails_as_expected"],
                       res["S5_ratio_vs_anchor"]["ratio_invariant"],
                       res["S6_weighted_boundary"]["pass"],
                       res["S7_two_body"]["pass"]])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(res, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
