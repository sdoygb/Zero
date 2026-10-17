#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""octagon_shape_check.py —— 中微子形状的原生候选：八分圆正切（推导稿 §9）。

模型：中微子三代在质量零点附近的小偏移
    x_k = tan(k·π/8)·d,   k = 0,1,2
（tan 来自闭包 2π/8 的八分圆；√2−1 = tan(π/8)）
质量：m_k = (A²/2)·x_k²，A² = v/784.5（同一个原生质量单元）。
则 m1=0、R := m2/m3 = tan²(π/8) = 3−2√2、K_ν = 2−√2。
标度 d 由 m3 定（1/d≈5.6×10⁴，原生来源未定）。

输出：ZeroGPU/results/octagon_shape_check.json
用法：python3 ZeroGPU/src/octagon_shape_check.py
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


V = 246.22
A2 = V / 784.5
t = math.tan(math.pi / 8)
R = t * t
K = (1 + R) / (1 + math.sqrt(R)) ** 2

R_obs = math.sqrt(7.42e-5 / 2.510e-3)
K_obs = (1 + R_obs) / (1 + math.sqrt(R_obs)) ** 2
print(f"   tan²(π/8) = {R:.9f} = 3−2√2 = {3-2*math.sqrt(2):.9f}")
print(f"   R_obs = {R_obs:.9f}   偏差 {100*(R/R_obs-1):+.3f}%")
print(f"   K = 2−√2 = {K:.9f}   K_obs = {K_obs:.9f}   偏差 {100*(K/K_obs-1):+.3f}%")

# 标度：由 m3=50.10 meV 定
m3 = 0.05010e-9  # GeV
d = math.sqrt(2 * m3 / A2)
ms = [(A2 / 2) * (d * math.tan(k * math.pi / 8)) ** 2 * 1e12 for k in range(3)]
print(f"   d = {d:.6e}  1/d = {1/d:.1f}（原生来源未定）")
print(f"   m = {[round(x,3) for x in ms]} meV；obs m2=8.614, m3=50.10")
X_dev = (R * R) / (7.42e-5 / 2.510e-3) - 1
print(f"   X 偏差 {100*X_dev:+.3f}%（相位律模型为 +5.47%）")

check("R 与实测差 <0.5%", abs(R / R_obs - 1) < 0.005, f"{100*(R/R_obs-1):+.3f}%")
check("K 与实测差 <0.1%", abs(K / K_obs - 1) < 1e-3, f"{100*(K/K_obs-1):+.3f}%")
check("m2 与实测差 <1%", abs(ms[1] / 8.614 - 1) < 0.01, f"{100*(ms[1]/8.614-1):+.3f}%")
check("m1=0（正序）", ms[0] == 0)

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "R": R, "K": K, "R_obs": R_obs, "K_obs": K_obs,
       "d": d, "m_meV": ms, "X_dev": X_dev}
with open(os.path.join(OUT, "octagon_shape_check.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/octagon_shape_check.json")
sys.exit(0 if not FAIL else 1)
