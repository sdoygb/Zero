#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立复核 Zb051–Zb055（外部审校报告 + 文档自检）。"""
import json
import math

import numpy as np

out = {}

# ---------------------------------------------------------------- Zb051
# Binet/Kepler: for kappa=1, a=1, e=0.5, m=1, check T^2=4*pi^2*a^3/kappa
kappa = 1.0
a = 1.0
e = 0.5
p = a * (1 - e * e)
L = math.sqrt(kappa * a * (1 - e * e))  # m=1
T_formula = 2 * math.pi * a ** 1.5 / math.sqrt(kappa)
# area speed L/(2m); ellipse area pi*a*b
b = a * math.sqrt(1 - e * e)
T_area = math.pi * a * b / (L / 2)
out["Zb051"] = {
    "p": p, "L": L, "T_formula": T_formula, "T_area": T_area,
    "T2_over_4pi2_a3": T_formula ** 2 / (4 * math.pi ** 2 * a ** 3),
    "rel_diff": abs(T_formula - T_area) / T_formula,
}

# ---------------------------------------------------------------- Zb052
G = 6.67430e-11
M_E = 5.972168e24
M_M = 7.342e22
M_tot = M_E + M_M
mu = M_E * M_M / M_tot
kappa_tot = G * M_tot
R_E = 6.371e6
R_M = 1.7374e6
I_E = 0.3307 * M_E * R_E ** 2
I_M = 0.393 * M_M * R_M ** 2
I_tot = I_E + I_M
a0 = 3.844e8
omega_spin_E = 2 * math.pi / 86164.0905
omega_spin_M = 2 * math.pi / (27.321661 * 86400)
L_orb0 = mu * math.sqrt(kappa_tot * a0)
L_spin_E0 = I_E * omega_spin_E
L_spin_M0 = I_M * omega_spin_M
L_tot = L_orb0 + L_spin_E0 + L_spin_M0
a_star = math.sqrt(3 * I_tot / mu)

def F(a):
    return mu * math.sqrt(kappa_tot * a) + I_tot * math.sqrt(kappa_tot / a ** 3) - L_tot

# bisection on [a_star*1.01, 10*a0]
lo, hi = a_star * 1.01, 10 * a0
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if F(mid) > 0:
        hi = mid
    else:
        lo = mid
a_f = 0.5 * (lo + hi)
T_f_sec = 2 * math.pi * math.sqrt(a_f ** 3 / kappa_tot)
out["Zb052"] = {
    "a_star": a_star,
    "a_f": a_f,
    "a_f_over_a0": a_f / a0,
    "T_f_days": T_f_sec / 86400,
    "L_tot": L_tot,
    "F_at_a_star_minus_Ltot": F(a_star),
}

# ---------------------------------------------------------------- Zb053
markov = {}
for k in [0, 1, 2]:
    for eps in [1e-2, 1e-4, 1e-6]:
        markov[f"k={k},eps={eps}"] = 10 * eps ** (1 - k)
out["Zb053"] = markov

# ---------------------------------------------------------------- Zb054
rng = np.random.default_rng(20261031)
# random connected graph: tree + extra edges
n = 8
edges = set()
for i in range(1, n):
    j = int(rng.integers(0, i))
    edges.add(tuple(sorted((i, j))))
while len(edges) < 2 * n - 1:
    i, j = sorted(rng.choice(n, size=2, replace=False))
    edges.add((int(i), int(j)))
L = np.zeros((n, n))
for i, j in edges:
    L[i, i] += 1
    L[j, j] += 1
    L[i, j] -= 1
    L[j, i] -= 1
phi = rng.normal(size=n)
E = phi @ L @ phi
phi_perp = phi - phi.mean() * np.ones(n)
E_perp = phi_perp @ L @ phi_perp
out["Zb054"] = {
    "L1_norm": float(np.linalg.norm(L @ np.ones(n))),
    "ker_dim": int(np.sum(np.abs(np.linalg.eigvalsh(L)) < 1e-9)),
    "E_phi": float(E),
    "E_perp": float(E_perp),
    "energy_diff": float(abs(E - E_perp)),
    "cross": float(phi.mean() * np.ones(n) @ L @ phi_perp),
}

# ---------------------------------------------------------------- Zb055
# q=pi/2, c=1: omega=2 sin(pi/4)=sqrt2; period in steps=2pi/omega.
c = 1.0
q = math.pi / 2
omega = 2 * c * abs(math.sin(q / 2))
period_steps = 2 * math.pi / omega
for a_e, a_o in [(1.0, 2.0)]:
    f_e = c * omega / a_e
    f_o = c * omega / a_o
    out["Zb055"] = {
        "omega": omega,
        "period_steps": period_steps,
        "period_naive_1_over_omega": 1 / omega,
        "f_e": f_e, "f_o": f_o,
        "f_ratio": f_e / f_o,
        "a_ratio": a_o / a_e,
        "period_equals_2pi_over_omega": abs(period_steps - 2 * math.pi / omega) < 1e-15,
    }

print(json.dumps(out, ensure_ascii=False, indent=2))
with open("_staging/Zb051-055_复核_20261031.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
