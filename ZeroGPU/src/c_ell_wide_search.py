#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""c_ell_wide_search.py —— c_ℓ 的结构式广搜（带 look-elsewhere 计数）。

目标：c_ℓ=μ_ℓ²/v=1.27462e-3（ln=6.663）。
形式：a₁^{e₁}·a₂^{e₂}·a₃^{e₃}（原子取 Zero 已导出的纯数；指数 ∈{−2,−1,1,2}）。
报告：最小相对误差；按「匹配数/组合数」给 look-elsewhere 期望，避免把巧合当结构。

输出：ZeroGPU/results/c_ell_wide_search.json
用法：python3 ZeroGPU/src/c_ell_wide_search.py
"""
from __future__ import annotations

import itertools
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)

V, M_TAU = 246.22, 1.77686
A_L = math.sqrt(2)
delta0 = math.pi / 12 - 2 / 9
delta_e = 3 * math.pi / 4 - delta0
f_tau = max(1 + A_L * math.cos(delta_e + 2 * math.pi * k / 3) for k in range(3))
c_ell = M_TAU / (V * f_tau**2)

atoms = {
    "2": 2.0, "3": 3.0, "pi": math.pi, "sqrt2": math.sqrt(2), "sqrt3": math.sqrt(3),
    "1/2": 0.5, "1/3": 1 / 3, "2/3": 2 / 3, "2/9": 2 / 9, "1/8": 1 / 8, "3/8": 3 / 8,
    "pi/12": math.pi / 12, "delta0": delta0, "ln2": math.log(2),
    "0.405": 0.405, "0.693": 0.693, "0.981": 0.981,
    "14": 14.0, "16": 16.0, "21": 21.0, "28": 28.0,
}
names = list(atoms)
exps = (-2, -1, 1, 2)
best = []
n_lt_1e4 = n_lt_1e3 = 0
total = 0
for a, b, c in itertools.permutations(names, 3):
    for ea in exps:
        va = atoms[a] ** ea
        for eb in exps:
            vab = va * atoms[b] ** eb
            for ec in exps:
                v = vab * atoms[c] ** ec
                if not (1e-8 < v < 1e8):
                    continue
                total += 1
                rel = abs(v / c_ell - 1)
                if rel < 1e-4:
                    n_lt_1e4 += 1
                    if rel < 1e-3:
                        n_lt_1e3 += 1
                if len(best) < 8 or rel < best[-1][0]:
                    best.append((rel, f"{a}^{ea}·{b}^{eb}·{c}^{ec}", v))
                    best.sort()
                    best = best[:8]
print(f"   c_ℓ={c_ell:.9e}；组合数 {total}")
print("   最小误差的 8 个：")
for rel, s, v in best:
    print(f"     {s:32s} = {v:.6e}  相对差 {100*rel:+.4f}%")
print(f"   <1e-4 的匹配：{n_lt_1e4} 个（期望 {total*1e-4/ (math.log(1e8)-math.log(1e-8)) *1:.1f} 量级）")
print(f"   <1e-3 的匹配：{n_lt_1e3} 个")

# 结论判据：若最小误差 >1e-4，则「三因子结构式」里没有精确候选
min_rel = best[0][0]
res = {"pass": [], "fail": [], "c_ell": c_ell, "n_combos": total,
       "n_lt_1e4": n_lt_1e4, "n_lt_1e3": n_lt_1e3, "min_rel": min_rel,
       "best": [{"form": s, "value": v, "rel": r} for r, s, v in best]}
with open(os.path.join(OUT, "c_ell_wide_search.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/c_ell_wide_search.json")
sys.exit(0 if min_rel > 1e-4 else 1)
