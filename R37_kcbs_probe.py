#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R37_kcbs_probe.py -- KCBS/contextuality 检验（定版）：Zero 的单体统计是量子的还是经典的？

判据（von Neumann 迹不等式，精确）：
  标准五角星射线 v_0..v_4（相邻正交），A := sum_i |v_i><v_i|，S(rho;U) = sum_i <v_i|U† rho U|v_i>。
  对全部取向 U 取最大：  S_max(rho) = sum_k lambda_k(rho) * mu_k(A)   （两侧同序降序）
  非语境界 2；S_max > 2 <=> 该态是语境的（单体量子性成立）。

输出：R37_kcbs_results.json
"""

from __future__ import annotations

import io
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NCTX = 2.0
QM_MAX = math.sqrt(5)

c45 = math.cos(4 * math.pi / 5)
th0 = math.atan(math.sqrt(-1.0 / c45))
vs = [np.array([math.cos(th0),
                math.sin(th0) * math.cos(4 * math.pi * i / 5),
                math.sin(th0) * math.sin(4 * math.pi * i / 5)], dtype=complex) for i in range(5)]
A = sum(np.outer(v, v.conj()) for v in vs)
mu = np.sort(np.linalg.eigvalsh(A))[::-1]

orth_err = max(abs(np.vdot(vs[i], vs[(i + 1) % 5])) for i in range(5))


def s_max(spec):
    s = np.sort(np.asarray(spec, dtype=float))[::-1]
    return float(np.dot(s, mu))


def s_at(U, rho):
    """给定取向的 S（用于交叉核对）"""
    return float(sum(np.real(np.vdot(vs[i], U.conj().T @ rho @ U @ vs[i])) for i in range(5)))


out = {}

# ---------------------------------------------------------------- P1 机制
out["P1_control"] = {
    "adjacent_orthogonality_error": float(orth_err),
    "mu_eigenvalues": [float(x) for x in mu],
    "mu_sum": float(mu.sum()),
    "mu1_equals_sqrt5": bool(abs(mu[0] - QM_MAX) < 1e-12),
    "S_maximally_mixed": s_max([1 / 3, 1 / 3, 1 / 3]),
    "S_optimal_pure": s_max([1, 0, 0]),
    "noncontextual_bound": NCTX,
}

# ---------------------------------------------------------------- P2 随机取向核对闭式
rng = np.random.default_rng(20261003)
worst = 0.0
for _ in range(3000):
    Z = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    Q, _ = np.linalg.qr(Z)
    lam = rng.dirichlet([1.0, 1.0, 1.0])
    rho = Q @ np.diag(lam) @ Q.conj().T
    worst = max(worst, s_at(np.eye(3), rho) - s_max(lam))
out["P2_closed_form_check"] = {
    "random_trials": 3000,
    "max_violation_of_bound_by_random_orientation": float(worst),
    "closed_form_upper_bounds_all_samples": bool(worst <= 1e-12),
}

# ---------------------------------------------------------------- P3 原生态谱
w_doc = np.array([2.0, 5.0, 20.0, 100.0])
w_doc = w_doc / w_doc.sum()
native = {
    "G29_4blocks_all": [float(x) for x in w_doc],
    "top3_normalized": [float(x) for x in (w_doc[:3] / w_doc[:3].sum())],
    "drop_smallest_3": [float(x) for x in ((w_doc[1:] / w_doc[1:].sum()))],
}
nat = {}
for k, sp in native.items():
    lam = np.sort(np.asarray(sp))[::-1]
    if len(lam) != 3:                      # 4 块谱只作数据登记，不代入 3 维判据
        nat[k] = {"spectrum": [float(x) for x in lam], "S_max": None, "contextual": None}
        continue
    S = s_max(lam)
    nat[k] = {"spectrum": [float(x) for x in lam], "S_max": S, "contextual": bool(S > NCTX)}
out["P3_native"] = nat
out["P3_native_source"] = "G72 §1 / G29 §2 的推前权重（块大小 2,5,20,100）"

# ---------------------------------------------------------------- P4 阈值
def lam_star():
    """族 (λ,(1-λ)/2,(1-λ)/2) 上 S_max = 2 的临界 λ"""
    return (NCTX - mu[1]) / (mu[0] - mu[1])


out["P4_threshold"] = {
    "family": "(lambda, (1-lambda)/2, (1-lambda)/2)",
    "formula": "S_max = mu1*lambda + mu2*(1-lambda)  =>  lambda* = (2-mu2)/(mu1-mu2)",
    "lambda_star": float(lam_star()),
    "native_lambda1_top3": float(np.max(native["top3_normalized"])),
    "native_margin": float(np.max(native["top3_normalized"]) - lam_star()),
}

out["verdict"] = {
    "criterion": "S_max(rho)=sum_k lambda_k mu_k > 2  <=>  语境（单体量子性）",
    "maximally_mixed": "非语境（5/3=1.667）",
    "native_state_contextual": bool(all(v["contextual"] for v in nat.values() if v["contextual"] is not None)),
    "native_S_values": {k: v["S_max"] for k, v in nat.items() if v["S_max"] is not None},
    "threshold_lambda_star": float(lam_star()),
    "caveats": [
        "state-dependent：测量集（pentagram）是按态选的，这是语境性检验的标准做法",
        "3 维归约的选择仍是 pi/E5 的缺口",
        "原生态谱用的是文档例（G72/G29 的推前权重），不是从 pi 显式算出的那一个",
    ],
}

with io.open(os.path.join(HERE, "R37_kcbs_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)

print(json.dumps(out, ensure_ascii=False, indent=1))
