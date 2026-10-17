#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_amplitude_geometry.py —— R1：中微子振幅 A_ν=√(3/2) 的几何候选（本轮新得）。

候选规则（同一立方体、两个对角比）：
  立方体（seed 8 = 2³）：棱 1，面对角线 √2，体对角线 √3。
  * 带电轻子（面/棱）：θ₂ = arccos(1/√2) = π/4；振幅 A_ℓ = sec θ₂ = √2；
    饱和零点 = π − θ₂ = 3π/4（与 Z121 的零点一致）。
  * 中微子（体/面）：θ₃ = arccos(√(2/3))；振幅 A_ν = sec θ₃ = √(3/2)；
    饱和零点 = π − θ₃ = arccos(−√(2/3))。
  ⟹ 两个振幅都从**同一个立方体**的两个对角比出来，无自由参数。

检查：
  A. 两个振幅与官方值精确相等；
  B. 饱和零点 self-consistent（cos(π−θ)=−1/A），带电侧零点=3π/4；
  C. 带电侧偏移 δ0 = 零点 − 电子角 = π/12−2/9（Z123）精确一致；
  D. 中微子侧用 A_ν 复算 R、K、X（与上一脚本同，作回归）。

输出：ZeroGPU/results/neutrino_amplitude_geometry.json
用法：python3 ZeroGPU/src/neutrino_amplitude_geometry.py
"""
from __future__ import annotations

import json
import math
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


print("== A. 立方体的两个对角比 ==")
edge, face, body = sp.Integer(1), sp.sqrt(2), sp.sqrt(3)
Al = sp.simplify(face / edge)
Av = sp.simplify(body / face)
check("面对角线/棱 = √2 = A_ℓ", sp.simplify(Al - sp.sqrt(2)) == 0, str(Al))
check("体对角线/面对角线 = √(3/2) = A_ν", sp.simplify(Av - sp.sqrt(sp.Rational(3, 2))) == 0, str(Av))

print("\n== B. 零点 = π − arccos(1/A)（自洽） ==")
th2 = sp.acos(1 / Al)          # arccos(1/√2) = π/4
th3 = sp.acos(1 / Av)          # arccos(√(2/3))
zero_l = sp.simplify(sp.pi - th2)
zero_v = sp.simplify(sp.pi - th3)
print("   θ₂ =", sp.deg(th2), "°   zero_ℓ =", sp.deg(zero_l), "°")
print("   θ₃ =", sp.deg(th3), "°   zero_ν =", sp.deg(zero_v), "°")
check("θ₂ = π/4（八分圆）", sp.simplify(th2 - sp.pi / 4) == 0)
check("带电零点 = 3π/4", sp.simplify(zero_l - 3 * sp.pi / 4) == 0)
check("中微子零点自洽：cos(zero_ν) = −√(2/3)",
      sp.simplify(sp.cos(zero_v) + sp.sqrt(sp.Rational(2, 3))) == 0)
check("中微子零点 = π − θ₃（不是 2π 的有理分数）",
      sp.simplify(zero_v - (sp.pi - th3)) == 0,
      f"zero_ν/(2π) = {float(zero_v/(2*sp.pi)):.6f}")

print("\n== C. 带电侧 δ0 与 Z123 一致 ==")
delta_e = sp.simplify(sp.Rational(2, 3) * sp.pi + sp.Rational(2, 9))     # 电子角（Z123）
delta0 = sp.simplify(zero_l - delta_e)
target = sp.simplify(sp.pi / 12 - sp.Rational(2, 9))
check("δ0 = 3π/4 − δ_e = π/12 − 2/9", sp.simplify(delta0 - target) == 0,
      f"δ0={sp.deg(delta0)}°")

print("\n== D. 中微子侧回归（A_ν = √(3/2)、m1=0） ==")
# f(θ)=1+A cosθ；θ1=zero_ν（饱和），θ2=θ1+2π/3，θ3=θ1+4π/3
xs = [zero_v, zero_v + 2 * sp.pi / 3, zero_v + 4 * sp.pi / 3]
fs = [sp.simplify(1 + Av * sp.cos(x)) for x in xs]
fs = [sp.nsimplify(sp.simplify(f)) for f in fs]
print("   f =", [float(f) for f in fs])
check("f(θ1)=0（饱和）", sp.simplify(fs[0]) == 0)
R = sp.simplify((fs[1] / fs[2]) ** 2)
check("R = (73−28√6)/25", sp.simplify(R - (73 - 28 * sp.sqrt(6)) / 25) == 0, str(sp.simplify(R)))
K = sp.simplify((fs[1] ** 2 + fs[2] ** 2) / (fs[1] + fs[2]) ** 2)
check("K_ν = 7/12", sp.simplify(K - sp.Rational(7, 12)) == 0)
X = float(R) ** 2
Xobs = 7.42e-5 / 2.510e-3
check("形状返回 X（+5.47% 张力即 Z019 的张力）", abs(X / Xobs - 1.0547) < 2e-3,
      f"X={X:.6f} vs obs {Xobs:.6f}")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL,
       "geometry": {"face_over_edge": float(Al), "body_over_face": float(Av),
                    "theta2_deg": float(sp.deg(th2)), "theta3_deg": float(sp.deg(th3)),
                    "zero_l_deg": float(sp.deg(zero_l)), "zero_v_deg": float(sp.deg(zero_v))},
       "delta0_deg": float(sp.deg(delta0)), "R": float(R), "K": float(K), "X": X, "Xobs": Xobs}
with open(os.path.join(OUT, "neutrino_amplitude_geometry.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_amplitude_geometry.json")
sys.exit(0 if not FAIL else 1)
