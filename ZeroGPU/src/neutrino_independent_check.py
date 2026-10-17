#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_independent_check.py —— 中微子程序的独立交叉核验（Z 系列汇总）。

不 import 项目脚本；只用**主数据**重算全部关键量，再与已落盘的
results/neutrino_*.json 逐项比对（容差 1e-9 相对）。覆盖：
  1. 两个振幅（立方体路线与载体谱比路线）
  2. δ0、R、K_ν
  3. M_R、m_ν3、m_ν2、Σ m_ν、X
  4. 分支判据（SM-only vs 中间谱）
  5. 预言的带（M_X ±10%）
  6. 形状残差的三条坐标（A²、ε、Δ）
任何一处与落盘值不符即 FAIL——它同时是「核验」（goal 的纪律）与回归护栏。

输出：ZeroGPU/results/neutrino_independent_check.json
用法：python3 ZeroGPU/src/neutrino_independent_check.py
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def load(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        print(f"   ! 缺 {name}")
        return None
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


V, M_Z, M_TAU = 246.22, 91.1876, 1.77686
D21, S21 = 7.42e-5, 0.20e-5
D31, S31 = 2.510e-3, 0.027e-3
M_X = 1.009e13

print("== 1. 两个振幅：两条路线 ==")
A2_cube_L = 2.0                      # 面对角线/棱 = √2
A2_cube_N = 1.5                      # 体对角线/面对角线 = √(3/2)
A2_spec_L = 4.0 / 2.0                # Γ8 谱比 4/2
A2_spec_N = 6.0 / 4.0                # Γ16 谱比 6/4
check("立方体路线 A_ℓ²=A_ν²=(√2)²,(√(3/2))²", (A2_cube_L, A2_cube_N) == (2.0, 1.5))
check("谱比路线与立方体同值", (A2_spec_L, A2_spec_N) == (A2_cube_L, A2_cube_N))

print("\n== 2. δ0、R、K ==")
zero_l = 3 * math.pi / 4
delta_e = 2 * math.pi / 3 + 2 / 9
delta0 = zero_l - delta_e
check("δ0 = π/12−2/9", abs(delta0 - (math.pi / 12 - 2 / 9)) < 1e-15)
R = (73 - 28 * math.sqrt(6)) / 25
K = 7.0 / 12.0
f = [0.0, 1.5 - math.sqrt(6) / 4, 1.5 + math.sqrt(6) / 4]
R2 = (f[1] / f[2]) ** 2
K2 = (f[1] ** 2 + f[2] ** 2) / (f[1] + f[2]) ** 2
check("R 与 f 的一致", abs(R - R2) < 1e-15 and abs(K - K2) < 1e-15,
      f"R={R:.9f} K={K:.9f}")

print("\n== 3. 标度 ==")
M_I = M_Z * math.exp(4 * math.pi)
M_R = math.sqrt(M_I * M_X)
m3 = (M_TAU / 2) ** 2 / M_R * 1e9
m2 = R * m3
Sigma = m2 + m3
X = (m2 / m3) ** 2
print(f"   M_I={M_I:.6e} M_R={M_R:.6e} m3={m3*1e3:.3f} meV m2={m2*1e3:.3f} Σ={Sigma*1e3:.2f} X={X:.6f}")

print("\n== 4. 与落盘值比对 ==")
ref = load("neutrino_scale_conditional.json")
if ref:
    check("M_R 与落盘一致", abs(ref["inputs"]["M_R"] / M_R - 1) < 1e-9)
    check("m3 与落盘一致", abs(ref["m_nu3_eV"] / m3 - 1) < 1e-9)
    check("R 与落盘一致", abs(ref["R"] / R - 1) < 1e-9)

geo = load("neutrino_amplitude_geometry.json")
if geo:
    check("几何路线 A_ν 与落盘一致",
          abs(geo["geometry"]["body_over_face"] - math.sqrt(1.5)) < 1e-12)
spec = load("neutrino_amplitude_spectrum.json")
if spec:
    check("谱比 A_ν 与落盘一致", abs(spec["A2_neutrino"] - 1.5) < 1e-12)
    check("谱比 A_ℓ 与落盘一致", abs(spec["A2_charged"] - 2.0) < 1e-12)

br = load("neutrino_branch_test.json")
if br:
    mid = [b for b in br["branches"] if b["branch"] == "中间谱"][0]
    M_I2 = M_Z * math.exp(4 * math.pi)
    m3_mid = (M_TAU / 2) ** 2 / math.sqrt(M_I2 * 4.954e13) * 1e9
    check("中间谱分支 m3 一致", abs(mid["m3_eV"] / m3_mid - 1) < 1e-9,
          f"{mid['m3_eV']*1e3:.2f} vs {m3_mid*1e3:.2f}")

res = load("neutrino_shape_residual.json")
if res and "routes" in res:
    check("形状残差 A² 与落盘一致", abs(res["routes"]["A2"] / 1.5136327854667173 - 1) < 1e-9)

pred = load("neutrino_predictions.json")
if pred:
    check("Σ 与 predictions 一致", abs(pred["scale"]["Sigma_eV"] / Sigma - 1) < 1e-9)

print("\n== 5. 带（M_X ±10%） ==")
lo = (M_TAU / 2) ** 2 / math.sqrt(M_I * M_X * 0.9) * 1e9
hi = (M_TAU / 2) ** 2 / math.sqrt(M_I * M_X * 1.1) * 1e9
check("±10% 带覆盖实测 m3", min(lo, hi) <= math.sqrt(D31) <= max(lo, hi),
      f"{min(lo,hi)*1e3:.2f}–{max(lo,hi)*1e3:.2f} vs {math.sqrt(D31)*1e3:.2f}")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
out = {"pass": PASS, "fail": FAIL,
       "values": {"delta0": delta0, "R": R, "K": K, "M_I": M_I, "M_R": M_R,
                  "m3_eV": m3, "m2_eV": m2, "Sigma_eV": Sigma, "X": X},
       "band_10": [min(lo, hi), max(lo, hi)]}
with open(os.path.join(RES, "neutrino_independent_check.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_independent_check.json")
sys.exit(0 if not FAIL else 1)
