#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_amplitude_winding.py —— 中微子振幅的卷绕公式 A=24·log2/(2π)（推导稿 §15）。

振幅 = 24 格累积的模流相位（以圈为单位）：
  A_nu = 24 * ln2 / (2π) = 2.6476272
其中 24 = 相位格数（Z8×Z3，Lean 已证）、ln2 = 非中心态内禀频率（Z087）、2π = 一圈。

由 A_nu 定零点 theta0 = arccos(-1/A)，三代在 theta0 + 2πk/3，
得 R = m2/m3、X = R^2，与观测和八分圆模型对照。

输出：ZeroGPU/results/neutrino_amplitude_winding.json
用法：python3 ZeroGPU/src/neutrino_amplitude_winding.py
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


X_obs = 7.42e-5 / 2.510e-3
R_obs = math.sqrt(X_obs)

A = 24 * math.log(2) / (2 * math.pi)
th0 = math.acos(-1 / A)
M = [abs(1 + A * math.cos(th0 + 2 * math.pi * k / 3)) for k in range(3)]
R = M[1] / M[2]
X = R * R
print(f"   A_nu = 24 ln2 / 2π = {A:.7f}   A^2 = {A*A:.6f}")
print(f"   零点 θ0 = {math.degrees(th0):.4f}°   M = {[round(x,6) for x in M]}")
print(f"   R = {R:.6f}（观测 {R_obs:.6f}，偏差 {100*(R/R_obs-1):+.3f}%）")
print(f"   X = {X:.6f}（观测 {X_obs:.6f}，偏差 {100*(X/X_obs-1):+.3f}%）")

R_oct = 3 - 2 * math.sqrt(2)
print(f"   八分圆 R = 3−2√2 = {R_oct:.6f}（偏差 {100*(R_oct/R_obs-1):+.3f}%）")

check("A_nu^2 与 7 差 <0.5%", abs(A * A / 7 - 1) < 0.005, f"{100*(A*A/7-1):+.3f}%")
check("X 与观测差 <0.2%", abs(X / X_obs - 1) < 0.002, f"{100*(X/X_obs-1):+.3f}%")
check("R 优于八分圆（|偏差| 更小）", abs(R / R_obs - 1) < abs(R_oct / R_obs - 1),
      f"{100*(R/R_obs-1):+.3f}% vs {100*(R_oct/R_obs-1):+.3f}%")
check("零点在 (0,π) 内", 0 < th0 < math.pi, f"{math.degrees(th0):.2f}°")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "A": A, "theta0_deg": math.degrees(th0),
       "M": M, "R": R, "X": X, "R_obs": R_obs, "X_obs": X_obs, "R_octagon": R_oct}
with open(os.path.join(OUT, "neutrino_amplitude_winding.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_amplitude_winding.json")
sys.exit(0 if not FAIL else 1)
