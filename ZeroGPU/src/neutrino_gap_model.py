#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_gap_model.py —— "16 步圆、走 15 步、差一步闭合"的统一模型（推导稿 §17）。

图像：
  · 16 槽圆等分 16 步，每步 2π/16 = π/8 = 22.5°；
  · 路径走 15 步（16 槽的 15 条边），第 16 步缺失 —— 缺口 = 一步 = π/8；
  · 缺口角给形状：R = tan^2(π/8) = 3−2√2；
  · 15 步每步带联合权重 2/9 = (2/3)(1/3)（非中心态两支，Z123）——标度 (2/9)^15；
  · m_k = A^2 (2/9)^15 tan^2(kπ/8),  k=0,1,2.

输出：ZeroGPU/results/neutrino_gap_model.json
用法：python3 ZeroGPU/src/neutrino_gap_model.py
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


A2 = 0.313856  # GeV，轻子质量单元 A^2
X_obs = 7.42e-5 / 2.510e-3
R_obs = math.sqrt(X_obs)
Dm21 = 7.42e-5
Dm31 = 2.510e-3

step = 2 * math.pi / 16
gap = 2 * math.pi - 15 * step
print(f"  每步 = {math.degrees(step):.4f}°   走 15 步 = {math.degrees(15*step):.4f}°   缺口 = {math.degrees(gap):.4f}°")

R_shape = math.tan(math.pi / 8) ** 2
print(f"  tan^2(π/8) = {R_shape:.9f}   3−2√2 = {3-2*math.sqrt(2):.9f}")
scale = (2 / 9) ** 15
m = [A2 * scale * math.tan(k * math.pi / 8) ** 2 for k in range(3)]
m_meV = [x * 1e12 for x in m]
print(f"  (2/9)^15 = {scale:.6e}")
print(f"  m = {[round(x,3) for x in m_meV]} meV")

m3 = 50.1  # meV 观测
R_model = math.sqrt(m[1] / m[2]) if m[2] > 0 else 0
X_model = m[1] / m[2]
print(f"  R = {X_model:.6f}（观测 {R_obs:.6f}，偏差 {100*(X_model/R_obs-1):+.3f}%）")
print(f"  Σ = {sum(m_meV):.2f} meV")
mbb = abs(sum(m_meV[k] * (0.0 if k == 0 else 0.0) for k in range(3)))
# 用 |U_ei|^2 近似：mββ=|Σ U_ei^2 m_i|，取 θ12=33.4°, θ13=8.6°
t12, t13 = math.radians(33.4), math.radians(8.6)
Ue = [math.cos(t13) * math.cos(t12), math.cos(t13) * math.sin(t12), math.sin(t13)]
mbb = abs(sum(Ue[i] ** 2 * m_meV[i] for i in range(3)))

check("每步 = π/8（22.5°）", abs(step - math.pi / 8) < 1e-12, f"{math.degrees(step):.4f}°")
check("缺口 = 一步 = π/8", abs(gap - step) < 1e-12, f"{math.degrees(gap):.4f}°")
check("tan²(π/8)=3−2√2", abs(R_shape - (3 - 2 * math.sqrt(2))) < 1e-12, f"{R_shape:.9f}")
check("m1=0（tan 0）", abs(m_meV[0]) < 1e-12, "0")
check("m3 与 50.1 meV 差 <1%", abs(m_meV[2] / 50.1 - 1) < 0.01, f"{m_meV[2]:.3f} meV")
check("Δm²21/Δm²31 与观测差 <1%", abs((m_meV[1] ** 2 / m_meV[2] ** 2) / X_obs - 1) < 0.01,
      f"{100*((m_meV[1]**2/m_meV[2]**2)/X_obs-1):+.3f}%")
print(f"  mββ ≈ {mbb:.2f} meV")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "step_deg": math.degrees(step), "gap_deg": math.degrees(gap),
       "tan2_pi8": R_shape, "scale_2_9_15": scale, "m_meV": m_meV, "mbb_meV": mbb,
       "X_model": m_meV[1] ** 2 / m_meV[2] ** 2, "X_obs": X_obs}
with open(os.path.join(OUT, "neutrino_gap_model.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_gap_model.json")
sys.exit(0 if not FAIL else 1)
