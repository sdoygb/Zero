#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_lambda_amplitude.py —— 中微子形状／标度：把旧库 Z_Λ 律的结论在现库复算一遍。

问题（目标）：只锚定 v，推导中微子质量。本轮只做三件事（全部可复跑）：
  A. 精确核验形状式：Z_Λ 律 √m=μ(1+A cos(θ+2πn/3))，m1=0 饱和，A_ν=√(3/2)
     ⟹ R=(73−28√6)/25 精确；K_ν=7/12；并给出与实测 X 的张力。
  B. 对照现库机制（A=√2 + δ0=π/12−2/9）：复现 8.36 倍失配，说明**差别只在振幅 A**。
  C. 标度侧：列出「seesaw m_ν=m_D²/M_R」用库内已导出的标度（M_I、M_X^unif、M_X^E6、几何均值）
     与三个 Dirac 候选（m_τ、m_τ/2、m_μ）得到的 m_ν3，和实测 √Δm²31 对照（只作线索，不作结论）。

输出：ZeroGPU/results/neutrino_lambda_amplitude.json
用法：python3 ZeroGPU/src/neutrino_lambda_amplitude.py
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


# ---------- 数据（NuFIT 正序，与 Z019 同源） ----------
D21 = 7.42e-5      # eV^2
D31 = 2.510e-3     # eV^2
X_OBS = D21 / D31
M3_OBS = math.sqrt(D31)          # 0.0501 eV
V_EW = 246.22                    # GeV
M_Z = 91.1876
M_TAU = 1.77686                  # GeV
M_MU = 0.1056584
M_E = 0.000510999

# ---------- 库内标度 ----------
M_I = M_Z * math.exp(4 * math.pi)          # Z154：M_I = M_Z e^{4π}
M_XU = 1.009e13                            # Z155：SM-only 统一标度
M_XE6 = 9.821e16                           # Z153/Z156：M_X = v e^{33.6197}
M_GM = math.sqrt(M_I * M_XU)               # 几何均值（候选）

print("== A. 形状：Z_Λ 律 + A_ν=√(3/2)，m1=0 饱和 ==")
A = sp.sqrt(sp.Rational(3, 2))
# m1=0: 1 + A cos(x) = 0  ⇒  cos x = -1/A,  sin x = +sqrt(1-1/A^2)
cosx = -1 / A
sinx = sp.sqrt(1 - 1 / A**2)
sp.simplify(sinx)
# f(θ) = 1 + A cos θ；两代相角 x+2π/3、x+4π/3
c120, s120 = sp.Rational(-1, 2), sp.sqrt(3) / 2
# cos(x+2π/3) = cos x cos120 − sin x sin120
c2 = cosx * c120 - sinx * s120
c3 = cosx * c120 + sinx * s120
f2 = sp.simplify(1 + A * c2)
f3 = sp.simplify(1 + A * c3)
R_sym = sp.simplify((f2 / f3) ** 2)
print("   f2 =", sp.simplify(f2), "=", float(f2))
print("   f3 =", sp.simplify(f3), "=", float(f3))
print("   f2/f3 =", sp.simplify(f2 / f3), "=", float(f2 / f3))
print("   R =", sp.simplify(R_sym), "=", float(R_sym))

R_target = (73 - 28 * sp.sqrt(6)) / 25
check("R 与 (73−28√6)/25 精确相等", sp.simplify(R_sym - R_target) == 0,
      f"差={sp.N(R_sym - R_target, 30)}")

R = float(R_sym)
K = sp.simplify((f2**2 + f3**2) / (f2 + f3) ** 2)   # m=f², √m=f
print("   K_ν =", sp.simplify(K), "=", float(K))
check("K_ν = 7/12（精确）", sp.simplify(K - sp.Rational(7, 12)) == 0)
X_model = R**2
print(f"   X_model = R² = {X_model:.9f}；X_obs = {X_OBS:.6f}；相对差 {100*(X_model-X_OBS)/X_OBS:+.3f}%")
check("形状预测与实测同向（X_model 偏高；张力即 Z019 的 5.47%）",
      0.02 < X_model < 0.04, f"X_model={X_model:.6f}")

print("\n== B. 现库机制（A=√2、δ0=π/12−2/9）的失配复现 ==")
A2 = sp.sqrt(2)
delta0 = sp.pi / 12 - sp.Rational(2, 9)
delta_e = sp.Rational(3, 4) * sp.pi - delta0     # 电子相角 = 3π/4 − δ0
f = [sp.simplify(1 + A2 * sp.cos(delta_e + 2 * sp.pi * k / 3)) for k in range(3)]
ms = sorted([float(x)**2 for x in f])          # 质量 m = f²（f=√m 的振幅）
X_mech = (ms[1]**2 - ms[0]**2) / (ms[2]**2 - ms[0]**2)
print(f"   X_mech = {X_mech:.9f}；X_obs = {X_OBS:.6f}；比 {X_OBS/X_mech:.3f} 倍")
check("8.36 倍失配复现（A=√2 的机制）", abs(X_OBS / X_mech - 8.36) < 0.05,
      f"比值={X_OBS/X_mech:.3f}")

# 若保持 m1=0 与 2π/3 间隔、只改振幅，使 X 精确命中所需的 A（数值解）
def X_of_A2(A2):
    A_ = math.sqrt(A2)
    c = -1 / A_
    s = math.sqrt(1 - 1 / A2)
    f2 = 1 + A_ * (c * -0.5 - s * (math.sqrt(3) / 2))
    f3 = 1 + A_ * (c * -0.5 + s * (math.sqrt(3) / 2))
    Rr = (f2 / f3) ** 2
    return Rr ** 2
from scipy.optimize import brentq
A2_needed = brentq(lambda x: X_of_A2(x) - X_OBS, 1.05, 4.0)
print(f"   若只改振幅并要求 X 精确：A²_required ≈ {A2_needed:.5f}（现库 2；Z_Λ 律 1.5）")
check("只改振幅即可精确命中 X，所需 A²≈1.514（与 Z_Λ 律的 3/2 差 0.9%）",
      abs(A2_needed - 1.5) / 1.5 < 0.02, f"A²_required={A2_needed:.4f}")

print("\n== C. 标度线索：seesaw m_ν3 = m_D²/M_R（不是结论，只列数） ==")
print(f"   目标 m3 = √Δm²31 = {M3_OBS*1e3:.3f} meV = {M3_OBS:.4f} eV")
rows = []
for mD_name, mD in [("m_τ", M_TAU), ("m_τ/2", M_TAU / 2), ("m_μ", M_MU), ("m_e", M_E)]:
    for MR_name, MR in [("M_I", M_I), ("√(M_I·M_X^unif)", M_GM),
                        ("M_X^unif", M_XU), ("M_X^E6", M_XE6),
                        ("M_I·e^{2π}", M_I * math.exp(2 * math.pi))]:
        mnu = mD**2 / MR            # GeV
        mnu_eV = mnu * 1e9
        rows.append({"m_D": mD_name, "M_R": MR_name, "m_nu3_eV": mnu_eV,
                     "ratio_to_target": mnu_eV / M3_OBS})
        print(f"   m_D={mD_name:5s} M_R={MR_name:18s} → m_ν3 = {mnu_eV:.3e} eV"
              f"   (×{mnu_eV/M3_OBS:.3f} 目标)")
near = [r for r in rows if 0.5 < r["ratio_to_target"] < 2]
check("存在 0.5–2 倍量级命中的候选（只作线索）", len(near) > 0,
      f"{len(near)} 个：" + "; ".join(f"{r['m_D']}/{r['M_R']}" for r in near))
for r in near:
    print(f"   ★ 候选：m_D={r['m_D']}，M_R={r['M_R']} → {r['m_nu3_eV']:.4f} eV")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "R": R, "K": float(K), "X_model": X_model,
       "X_obs": X_OBS, "X_mech": X_mech, "A2_required": A2_needed,
       "scales": {"M_I": M_I, "M_X_unif": M_XU, "M_X_E6": M_XE6, "geo_mean": M_GM},
       "scale_rows": rows, "target_m3_eV": M3_OBS}
with open(os.path.join(OUT, "neutrino_lambda_amplitude.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_lambda_amplitude.json")
sys.exit(0 if not FAIL else 1)
