#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R44_survival_vs_contextuality_probe.py -- 选维（生存）与单体语境性互斥。

(S) 生存：F_D = C(D,2) q^D 的峰在 D=4  <=>  q in (1/2, 3/5)
(C) 语境性：S_max(q) = mu2 + (mu1-mu2)(1+sqrt(2q-1))/2 > 2  <=>  q > 3/5
输出：R44_survival_vs_contextuality_results.json
"""
from __future__ import annotations
import io, json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MU1 = math.sqrt(5.0)
MU2 = (5.0 - math.sqrt(5.0)) / 2.0


def peak_D(q, Dmax=40):
    F = [math.comb(D, 2) * q ** D for D in range(2, Dmax + 1)]
    return 2 + int(np.argmax(F)), F


def s_max_bound(q):
    """给定 q=sum lambda^2，S 的上确界（Cauchy-Schwarz 饱和）"""
    if 2 * q - 1 < 0:
        return None
    return MU2 + (MU1 - MU2) * (1.0 + math.sqrt(2 * q - 1)) / 2.0


def s_max_spectrum(lam):
    lam = np.sort(np.asarray(lam, float))[::-1]
    if lam.size < 3:
        lam = np.concatenate([lam, np.zeros(3 - lam.size)])
    return float(np.dot(lam[:3], [MU1, MU2, MU2]))


rows = []
for q in (0.40, 0.45, 0.50, 0.53, 5 / 9, 0.58, 0.60, 0.62, 0.65, 0.70, 0.80):
    pd, F = peak_D(q)
    sb = s_max_bound(q)
    rows.append({"q": q, "peak_D": pd, "F4_over_F3": F[2] / F[1], "F5_over_F4": F[3] / F[2],
                 "S_max_bound": sb,
                 "survival_peak_at_4": bool(pd == 4),
                 "contextual": bool(sb is not None and sb > 2 + 1e-12)})

# 解析临界
q_star = 0.6
lam1_star = (1 + math.sqrt(2 * q_star - 1)) / 2
two_block_a = (2 - MU2) / (MU1 - MU2)

# 两条既有正面结论的载体
routeA_L4 = [4 / 6, 2 / 6]
doc = [0.0157, 0.0394, 0.1575, 0.7874]

out = {
    "mu": [MU1, MU2, MU2],
    "survival_window_for_D4": [0.5, 0.6],
    "contextuality_threshold_q": q_star,
    "analytic": {
        "lambda1_max_of_q": "lambda1_max = (1+sqrt(2q-1))/2",
        "S_max_of_q": "S_max = mu2 + (mu1-mu2)*(1+sqrt(2q-1))/2",
        "lambda1_at_q_star": lam1_star,
        "two_block_critical_a": two_block_a,
        "S_max_at_q_star": s_max_bound(q_star),
    },
    "rows": rows,
    "carriers": {
        "routeA_L4": {"q": sum(w * w for w in routeA_L4), "S_max": s_max_spectrum(routeA_L4),
                      "peak_D": peak_D(sum(w * w for w in routeA_L4))[0]},
        "documented_2_5_20_100": {"q": sum(w * w for w in doc), "S_max": s_max_spectrum(doc),
                                  "peak_D": peak_D(sum(w * w for w in doc))[0]},
    },
    "verdict": {
        "windows_are_complementary": True,
        "no_go": "q<3/5 => 峰在 D=4 但 S_max<2；q>3/5 => S_max>2 但峰在 D>=5；二者在 q=3/5 相接不重叠",
        "routeA_L4": "生存 ✅ / 语境 ✗",
        "documented_case": "语境 ✅ / 生存 ✗",
        "escape": "只能改账本形式 F_D（R31/R32 框架内的 no-go）",
    },
}
with io.open(os.path.join(HERE, "R44_survival_vs_contextuality_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)

print("q       峰D  F4/F3   F5/F4   S_max上界  (S)生存  (C)语境")
for r in rows:
    print("%-7.4f %-4d %-7.4f %-7.4f %-9.4f %-8s %s" % (
        r["q"], r["peak_D"], r["F4_over_F3"], r["F5_over_F4"],
        (-1.0 if r["S_max_bound"] is None else r["S_max_bound"]),
        "✅" if r["survival_peak_at_4"] else "✗", "✅" if r["contextual"] else "✗"))
print()
print("临界 q*=3/5: lambda1*=%.6f  S_max=%.6f" % (lam1_star, s_max_bound(0.6)))
print(json.dumps(out["carriers"], ensure_ascii=False, indent=1))
