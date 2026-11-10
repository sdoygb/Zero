#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Zb089 · 层间交换：潜热是层间通量的能量当量。

五组判据（对应正文 §2–§3）：
  S1 分层通量律   层 A（干空气）与层 W（水汽）各自的计数变化 = ∓ν；总分守恒。
  S2 单值当量     往返闭环 A→W→A 上，若两个方向的当量相等，总能量变化恒为 0（Zb045 定理 1.2 的实例）。
  S3 判据方向     当量随方向变（ε↑≠ε↓）时，同一个闭环净增 (ε↑-ε↓)ν——能量凭空出现。
  S4 通量任意性   ν 取遍整数，S3 的漂移随之线性变（Zb004 定理 3.1：零和分布间总有整数流）。
  S5 细致平衡     三态环：细致平衡时净环流为 0；乘 2.5／0.6 破坏后净环流非零（定理 3.3）。

运行：python3 ZeroGPU/src/zb089_layer_exchange.py
输出：ZeroGPU/results/zb089_layer_exchange.json
"""
from __future__ import annotations

import json
import random
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "zb089_layer_exchange.json"


def build(nA=4, nW=3):
    """层 A：0..nA-1；层 W：nA..nA+nW-1；层间两条边。"""
    edges = [(i, i + 1) for i in range(nA - 1)]
    edges += [(nA + i, nA + i + 1) for i in range(nW - 1)]
    edges += [(0, nA), (nA - 1, nA + 1)]           # 层间
    n = nA + nW
    B = np.zeros((n, len(edges)))
    for e, (i, j) in enumerate(edges):
        B[i, e], B[j, e] = -1.0, +1.0
    inter = [e for e, (i, j) in enumerate(edges) if (i < nA) != (j < nA)]
    return edges, B, inter, nA, nW


def s1_law(trials=300, seed=31):
    edges, B, inter, nA, nW = build()
    rng = random.Random(seed)
    worstA = worstW = worstTot = 0.0
    for _ in range(trials):
        x = np.array([rng.randint(-5, 5) for _ in range(nA + nW)], dtype=float)
        f = np.array([rng.randint(-5, 5) for _ in edges], dtype=float)
        xn = x + B @ f
        dA = xn[:nA].sum() - x[:nA].sum()
        dW = xn[nA:].sum() - x[nA:].sum()
        nu = 0.0
        for e in inter:
            i, j = edges[e]
            if i < nA:      # A → W
                nu += f[e]
            else:           # W → A
                nu -= f[e]
        worstA = max(worstA, abs(dA + nu))
        worstW = max(worstW, abs(dW - nu))
        worstTot = max(worstTot, abs(dA + dW))
    return {"trials": trials, "max_abs_residual_A": float(worstA),
            "max_abs_residual_W": float(worstW), "max_abs_total_change": float(worstTot),
            "pass": bool(max(worstA, worstW, worstTot) < 1e-12)}


def cycle(nu=3, eps_up=2.0, eps_down=2.0):
    """一次 A→W→A 的往返：先转 nu 单位上去，再转 nu 单位下来。"""
    dE = eps_up * nu - eps_down * nu
    return {"nu": nu, "eps_up": eps_up, "eps_down": eps_down, "delta_E": float(dE)}


def s5_detailed_balance(beta=1.0, seed=32):
    """定理 3.3：三态环上的细致平衡与净环流。

    速率 k_ij = exp(-β(E_j-E_i)/2) 使 π ∝ e^{-βE} 满足细致平衡。
    """
    E = np.array([0.0, 1.0, -0.4])

    def rates(mult=None):
        k = np.exp(-beta * (E[None, :] - E[:, None]) / 2.0)
        np.fill_diagonal(k, 0.0)
        if mult:
            for (i, j), m in mult.items():
                k[i, j] *= m
        return k

    def stationary(k):
        n = k.shape[0]
        Q = k.copy()
        np.fill_diagonal(Q, -k.sum(1))
        A = np.vstack([Q.T, np.ones(n)])
        b = np.zeros(n + 1)
        b[-1] = 1.0
        return np.linalg.lstsq(A, b, rcond=None)[0]

    def cycle_current(k, pi):
        return float(pi[0] * k[0, 1] * k[1, 2] * k[2, 0]
                     - pi[0] * k[0, 2] * k[2, 1] * k[1, 0])

    k0 = rates()
    pi0 = stationary(k0)
    k1 = rates({(0, 1): 2.5, (1, 2): 0.6})
    pi1 = stationary(k1)
    return {"E": E.tolist(), "beta": beta,
            "pi_balanced": pi0.tolist(),
            "cycle_current_balanced": cycle_current(k0, pi0),
            "pi_broken": pi1.tolist(),
            "cycle_current_broken": cycle_current(k1, pi1),
            "pass": bool(abs(cycle_current(k0, pi0)) < 1e-14
                         and abs(cycle_current(k1, pi1)) > 1e-3)}


def main():
    r = {"S1_law": s1_law()}
    r["S2_single_valued"] = {"cases": [cycle(nu, 2.0, 2.0) for nu in (1, 3, 7)],
                             "all_zero": all(abs(cycle(n, 2.0, 2.0)["delta_E"]) < 1e-15 for n in (1, 3, 7))}
    r["S3_direction"] = {"cases": [cycle(nu, 2.0, 1.5) for nu in (1, 3, 7)],
                         "drift_nonzero": all(abs(cycle(n, 2.0, 1.5)["delta_E"]) > 0.1 for n in (1, 3, 7))}
    r["S4_linearity"] = {str(nu): cycle(nu, 2.0, 1.5)["delta_E"] for nu in (1, 2, 3, 4, 5)}
    r["S5_detailed_balance"] = s5_detailed_balance()
    r["PASS"] = bool(r["S1_law"]["pass"] and r["S2_single_valued"]["all_zero"]
                     and r["S3_direction"]["drift_nonzero"]
                     and r["S5_detailed_balance"]["pass"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(r, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(r, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
