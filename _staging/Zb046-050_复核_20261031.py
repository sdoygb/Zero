#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立复核 Zb046–Zb050（外部审校报告 + 文档自检）。
只做数值/代数核验，不修改正文。
"""
import itertools
import json
import math
from collections import defaultdict, deque

import numpy as np

out = {}

# ---------------------------------------------------------------- Zb046
rng = np.random.default_rng(20261031)
# (A1) block LU + determinant + Schur complement positive definiteness
max_lu = 0.0
max_det = 0.0
for _ in range(200):
    m, n = 3, 4
    X = rng.normal(size=(m + n, m + n))
    M = X @ X.T + np.eye(m + n)  # SPD
    A, B = M[:m, :m], M[:m, m:]
    C, D = M[m:, :m], M[m:, m:]
    S = A - B @ np.linalg.solve(D, C)
    L = np.block([[np.eye(m), B @ np.linalg.inv(D)], [np.zeros((n, m)), np.eye(n)]])
    R = np.block([[S, np.zeros((m, n))], [C, D]])
    max_lu = max(max_lu, float(np.max(np.abs(L @ R - M))))
    max_det = max(max_det, float(abs(np.linalg.det(M) - np.linalg.det(D) * np.linalg.det(S))))
    eigS = np.linalg.eigvalsh(S)
    if np.min(eigS) <= 0:
        raise RuntimeError("Schur complement not positive definite")
out["Zb046"] = {"max_block_LU_error": max_lu, "max_det_factor_error": max_det,
                "random_SPD_trials": 200, "min_eig_S": float(np.min(eigS))}

# ---------------------------------------------------------------- Zb047
def balanced_words(L):
    for p in itertools.combinations(range(L), L // 2):
        w = [-1] * L
        for idx in p:
            w[idx] = 1
        yield w

def layer_heights(L):
    b = defaultdict(int)
    for w in balanced_words(L):
        s = 0
        h = 0
        for x in w:
            s += x
            h = max(h, abs(s))
        b[h] += 1
    return dict(sorted(b.items()))

def fixed_counts(L):
    fs = fr = frs = 0
    for w in balanced_words(L):
        if w == [-x for x in w]:
            fs += 1
        if w == w[::-1]:
            fr += 1
        if w == [-x for x in w[::-1]]:
            frs += 1
    return fs, fr, frs

tab = {}
for L in [4, 6, 8, 10, 12, 14, 16]:
    fs, fr, frs = fixed_counts(L)
    tab[L] = {"fixed_s": fs, "fixed_r": fr, "fixed_rs": frs}
out["Zb047_fixed"] = tab

# DP distribution and checks
def dp_dist(L):
    dp = {(0, 0): 1}
    for step in range(L):
        ndp = defaultdict(int)
        for (s, m), cnt in dp.items():
            for dv in (1, -1):
                ns = s + dv
                if abs(ns) <= L - 1 - step:
                    nm = max(m, abs(ns))
                    ndp[(ns, nm)] += cnt
        dp = ndp
    b = defaultdict(int)
    for (s, m), cnt in dp.items():
        if s == 0:
            b[m] += cnt
    return dict(sorted(b.items()))

dp_tab = {}
for L in [4, 6, 8, 12, 16, 20, 24, 28, 32, 40, 48]:
    b = dp_dist(L)
    H = max(b)
    total = sum(b.values())
    dp_tab[L] = {
        "dist": b,
        "total": total,
        "peak_h": max(b, key=lambda h: b[h]),
        "w_peak": b[max(b, key=lambda h: b[h])] / total,
        "mean_h": sum(h * c for h, c in b.items()) / total,
        "mean_over_sqrtL": sum(h * c for h, c in b.items()) / total / math.sqrt(L),
        "b_H": b[H],
        "b_Hm1": b[H - 1],
        "2L-2": 2 * L - 2,
        "max_gap": max(math.log(b[h] / b[h + 1]) for h in range(1, H)),
        "ln_Lm1": math.log(L - 1),
    }
out["Zb047_dp"] = dp_tab

# ---------------------------------------------------------------- Zb048
def apply_perm(edge, p):
    u, v = edge
    return tuple(sorted((p[u - 1], p[v - 1])))

def apply_matching(m, p):
    return frozenset(apply_perm(e, p) for e in m)

matchings = [frozenset({(1, 2), (3, 4)}),
             frozenset({(1, 3), (2, 4)}),
             frozenset({(1, 4), (2, 3)})]
S4 = list(itertools.permutations([1, 2, 3, 4]))
stab1 = [p for p in S4 if apply_matching(matchings[0], p) == matchings[0]]
ker = [p for p in S4 if all(apply_matching(m, p) == m for m in matchings)]
preserve_edges = [p for p in S4 if set([p[0], p[1]]) == {1, 2} and set([p[2], p[3]]) == {3, 4}]
V4 = [(1, 2, 3, 4), (2, 1, 4, 3), (3, 4, 1, 2), (4, 3, 2, 1)]
out["Zb048"] = {
    "stab_M1_order": len(stab1),
    "stab_M1": stab1,
    "ker_order": len(ker),
    "ker": ker,
    "ker_equals_V4": set(ker) == set(V4),
    "preserve_edges_order": len(preserve_edges),
    "preserve_edges": preserve_edges,
    "preserve_edges_equals_V4": set(preserve_edges) == set(V4),
    "V4": V4,
}

# ---------------------------------------------------------------- Zb049
# (D1) finite memory does NOT imply linear: h_{n+1}=h_n^2
h = 0.7
h2 = h * h
out["Zb049_finite_memory_nonlinear"] = {
    "example": "h_{n+1}=h_n^2",
    "H0": h, "H1": h2,
    "linear_in_H0_would_give": 2 * h,  # if linear with slope 2
    "note": "nonlinear map, so Theorem 5.1's hypothesis must be linear recurrence",
}

# (D2) c=1 degeneracy
c = 1.0
h0, h1 = 0.0, 1.0
seq = [h0, h1]
for n in range(1, 8):
    seq.append((c + 1) * seq[-1] - c * seq[-2])
out["Zb049_c1_degeneracy"] = {
    "c": c, "h0": h0, "h1": h1, "seq": seq,
    "alpha_beta_cn_with_c1": [h1 for _ in seq],
    "note": "c=1 gives h_n=alpha+beta*n, not alpha+beta*1^n",
}

# (D3) c=1.7 numeric and S=3 roots
c = 1.7
h0, h1 = 1.0, 3.0
seq = [h0, h1]
for n in range(1, 12):
    seq.append((c + 1) * seq[-1] - c * seq[-2])
beta = (h1 - h0) / (c - 1)
alpha = h0 - beta
closed = [alpha + beta * c ** n for n in range(len(seq))]
out["Zb049_c17"] = {"max_abs_diff": float(max(abs(a - b) for a, b in zip(seq, closed))),
                    "alpha": alpha, "beta": beta}
S = 3.0
roots = [(S + math.sqrt(S * S + 4 * S)) / 2, (S - math.sqrt(S * S + 4 * S)) / 2]
out["Zb049_S3_roots"] = roots

# (D4) compensation graph degree and connectivity, stationarity
def neighbors(w):
    L = len(w)
    plus = [i for i, x in enumerate(w) if x == 1]
    minus = [i for i, x in enumerate(w) if x == -1]
    for i in plus:
        for j in minus:
            nw = list(w)
            nw[i], nw[j] = nw[j], nw[i]
            yield tuple(nw)

graph_tab = {}
for L in [4, 6, 8]:
    X = list(balanced_words(L))
    Xset = set(map(tuple, X))
    degs = {len(list(neighbors(w))) for w in Xset}
    # BFS
    start = next(iter(Xset))
    seen = {start}
    q = deque([start])
    while q:
        u = q.popleft()
        for v in neighbors(u):
            if v not in seen:
                seen.add(v)
                q.append(v)
    # stationarity check
    d = next(iter(degs))
    stat = all(abs(sum((1 / len(Xset)) * (1 / d) for u in Xset if w in set(neighbors(u))) - 1 / len(Xset)) < 1e-12
               for w in Xset)
    graph_tab[L] = {"n": len(Xset), "degrees": sorted(degs), "components": 1 if len(seen) == len(Xset) else "?",
                    "stationary_ok": stat}
out["Zb049_graph"] = graph_tab

# ---------------------------------------------------------------- Zb050
alpha_inv = 137.035999084
pred = 5 / alpha_inv
obs = 0.99726968 / 27.321661
out["Zb050_period"] = {
    "pred": pred,
    "obs_sidereal": obs,
    "rel_pred_minus_obs": (pred - obs) / obs,
    "rel_obs_minus_pred": (obs - pred) / obs,
    "solar_day_over_sidereal_month": 1.0 / 27.321661,
    "sidereal_day_over_synodic_month": 0.99726968 / 29.53059,
}
r_vals = [0.333, 0.234, 0.298, 0.309, 0.135, 0.081, 0.045, 0.036, 0.388, 0.005]
F_all = math.sqrt(sum(x * x for x in r_vals) / len(r_vals))
F_keep = math.sqrt(sum(x * x for x in r_vals[:5] + r_vals[8:]) / 7)
out["Zb050_F"] = {"F_all": F_all, "F_keep": F_keep, "ratio": F_keep / F_all}

Lambda, k0, LH = 3, 2, 150
exponent = math.sqrt(LH / (k0 * Lambda))
r0 = (Lambda / k0) * 137.036 * (1836.153 ** exponent) * 2.998e8 * 3.6165e-17
AU = 1.495978707e11
out["Zb050_r0"] = {"exponent": exponent, "r0_m": r0, "r0_AU": r0 / AU,
                    "rel_vs_0.3074": (r0 / AU - 0.3074) / 0.3074}

print(json.dumps(out, ensure_ascii=False, indent=2, default=str))
with open("_staging/Zb046-050_复核_20261031.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
