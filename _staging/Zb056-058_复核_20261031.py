#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立复核 Zb056–Zb058（外部审校报告 + 文档自检）。"""
import json
import math

import numpy as np

out = {}

# ---------------------------------------------------------------- Zb056
# a = tau^2, c=1: true H = hc/a; naive t=tau*a gives H_naive = H_true/(1+tau h)
rows = []
for tau in [1.0, 2.0, 5.0]:
    a = tau ** 2
    h = 2 / tau
    H_true = h * 1.0 / a
    H_naive = h * 1.0 / (a * (1 + tau * h))
    rows.append({"tau": tau, "a": a, "h": h, "H_true": H_true,
                 "H_naive": H_naive, "ratio_true_naive": H_true / H_naive,
                 "1_plus_tau_h": 1 + tau * h})
out["Zb056"] = rows

# ---------------------------------------------------------------- Zb057
c = 1.0
p = 2.0
q = math.pi
omega = 2 * c * abs(math.sin(q / 2))
T = 2 * math.pi / omega
rows = []
for tau in [1.0, 3.0, 7.0]:
    a = tau ** p
    h = p / tau
    rows.append({
        "tau": tau, "a": a, "h": h,
        "omega": omega, "T": T,
        "h_times_T": h * T,
        "h_over_omega": h / omega,
        "h_bound_omega_over_2pi": omega / (2 * math.pi),
        "h_bound_omega": omega,
    })
out["Zb057"] = {
    "omega_max": 2 * c,
    "T_at_pi": T,
    "h_hard_limit_c_over_pi": c / math.pi,
    "tau_threshold_pi_p_over_c": math.pi * p / c,
    "rows": rows,
    "H_bound_c2_over_pi_a": "c^2/(pi a)",
}

# ---------------------------------------------------------------- Zb058
theta3 = math.acos(1 / 3)
theta4 = math.acos(1 / 4)
defects3 = {k: 2 * math.pi - k * theta3 for k in range(3, 8)}
defects4 = {k: 2 * math.pi - k * theta4 for k in range(3, 6)}

# Numerical check of regular simplex dihedral angles.
def simplex_dihedral(n):
    # standard n-simplex in R^(n+1): vertices e_i (i=0..n)
    V = np.eye(n + 1)
    # facets: omit vertex j; normal to facet is the vector pointing to omitted vertex
    # angle between outward normals of facets 0 and 1
    n0 = V[0] - V[1:].mean(axis=0)
    n1 = V[1] - np.vstack([V[0], V[2:]]).mean(axis=0)
    cosang = np.dot(n0, n1) / (np.linalg.norm(n0) * np.linalg.norm(n1))
    return math.pi - math.acos(max(-1.0, min(1.0, cosang)))

out["Zb058"] = {
    "theta3_deg": math.degrees(theta3),
    "theta3_2pi_over": 2 * math.pi / theta3,
    "defects3_deg": {k: math.degrees(v) for k, v in defects3.items()},
    "theta4_deg": math.degrees(theta4),
    "theta4_2pi_over": 2 * math.pi / theta4,
    "defects4_deg": {k: math.degrees(v) for k, v in defects4.items()},
    "simplex_dihedral_3d_deg": math.degrees(simplex_dihedral(3)),
    "simplex_dihedral_4d_deg": math.degrees(simplex_dihedral(4)),
    "arccos_1_3_deg": math.degrees(math.acos(1 / 3)),
    "arccos_1_4_deg": math.degrees(math.acos(1 / 4)),
}

print(json.dumps(out, ensure_ascii=False, indent=2))
with open("_staging/Zb056-058_复核_20261031.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
