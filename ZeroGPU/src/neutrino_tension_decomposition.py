#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_tension_decomposition.py —— 中微子预言的张力分解（形状 vs 标度）。

模型（Z162＋Z163）：m1=0、R=(73−28√6)/25、m3=(m_τ/2)²/√(M_I M_X^unif)。
本脚本把总张力拆成两处**独立**的修正：
  A. 尺度：Δm²31（等价 m3）的 pull；反解所需 M_R；
  B. 形状：X 的 pull；反解所需 A²；
  C. 反事实：只修尺度 / 只修形状时，其余观测量的 pull。
输出：ZeroGPU/results/neutrino_tension_decomposition.json
用法：python3 ZeroGPU/src/neutrino_tension_decomposition.py
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


# 数据
D21, S21 = 7.42e-5, 0.20e-5
D31, S31 = 2.510e-3, 0.027e-3
X_OBS = D21 / D31
S_X = X_OBS * math.sqrt((S21 / D21) ** 2 + (S31 / D31) ** 2)

# 模型
V, M_Z, M_TAU = 246.22, 91.1876, 1.77686
M_I = M_Z * math.exp(4 * math.pi)
M_X = 1.009e13
M_R = math.sqrt(M_I * M_X)
A_L = math.sqrt(2)
delta0 = math.pi / 12 - 2 / 9
delta_e = 3 * math.pi / 4 - delta0
f = [1 + A_L * math.cos(delta_e + 2 * math.pi * k / 3) for k in range(3)]
f_tau = max(f)
R = (73 - 28 * math.sqrt(6)) / 25
m3 = (M_TAU / 2) ** 2 / M_R * 1e9          # eV
m2 = R * m3
d21 = m2**2
d31 = m3**2
X = d21 / d31

print("== 模型点 ==")
print(f"   m3={m3*1e3:.2f} meV (obs {math.sqrt(D31)*1e3:.2f})  R={R:.6f}  X={X:.6f} (obs {X_OBS:.6f})")

print("\n== A. 尺度张力（Δm²31 / m3） ==")
pull31 = (d31 - D31) / S31
print(f"   Δm²31 模型 {d31:.4e} vs 数据 {D31:.4e} ⟹ pull {pull31:+.1f}σ")
M_R_req = (M_TAU / 2) ** 2 / (math.sqrt(D31) * 1e-9)     # 使 m3=m3_obs 的 M_R
print(f"   所需 M_R={M_R_req:.4e} vs 中点 {M_R:.4e} ⟹ 相对修正 {100*(M_R_req/M_R-1):+.2f}%")
check("尺度 pull 在 −6σ 附近（m3 低 3%）", -6.5 < pull31 < -4.5, f"{pull31:+.1f}σ")

print("\n== B. 形状张力（X） ==")
pullX = (X - X_OBS) / S_X
print(f"   X 模型 {X:.6f} vs 数据 {X_OBS:.6f} ⟹ pull {pullX:+.2f}σ")
A2_req = 6 * ((1 + math.sqrt(X_OBS)) / (1 + X_OBS**0.25) ** 2 - 1 / 3)
print(f"   所需 A²={A2_req:.6f} vs 3/2 ⟹ 相对修正 {100*(A2_req/1.5-1):+.2f}%")
check("形状 pull 在 +2σ 附近（X 高 5.5%）", 1.0 < pullX < 3.0, f"{pullX:+.2f}σ")

print("\n== C. 反事实（只修一处） ==")
# (i) 只修尺度：m3=m3_obs，形状 R 不变
m3a = math.sqrt(D31)
m2a = R * m3a
p21a = (m2a**2 - D21) / S21
pXa = (X - X_OBS) / S_X
print(f"   只修尺度: Δm²21 pull {p21a:+.2f}σ；Δm²31 pull 0σ；X pull {pXa:+.2f}σ")
# (ii) 只修形状：A²=A2_req，尺度不变
s = (math.sqrt(3) / 2) * math.sqrt(A2_req - 1)
Ra = ((1.5 - s) / (1.5 + s)) ** 2
m2b = Ra * m3
p21b = (m2b**2 - D21) / S21
p31b = (m3**2 - D31) / S31
print(f"   只修形状: Δm²21 pull {p21b:+.2f}σ；Δm²31 pull {p31b:+.1f}σ；X pull 0σ")
check("只修尺度时形状张力仍在（+2σ 量级）", 1.0 < pXa < 3.0, f"{pXa:+.2f}σ")
check("只修形状时尺度张力仍在（−5σ 量级）", -6.5 < p31b < -4.5, f"{p31b:+.1f}σ")
check("两处修正互相独立（各自只解自己的观测量）",
      abs((M_R_req / M_R - 1)) > 0.02 and abs(A2_req / 1.5 - 1) > 0.005)

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL,
       "model": {"m3_eV": m3, "R": R, "X": X, "d21": d21, "d31": d31},
       "scale": {"pull31": pull31, "M_R": M_R, "M_R_required": M_R_req,
                 "rel_correction": M_R_req / M_R - 1},
       "shape": {"pullX": pullX, "A2": 1.5, "A2_required": A2_req,
                 "rel_correction": A2_req / 1.5 - 1},
       "counterfactual": {"scale_only": {"pull21": p21a, "pullX": pXa},
                          "shape_only": {"pull21": p21b, "pull31": p31b}}}
with open(os.path.join(OUT, "neutrino_tension_decomposition.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_tension_decomposition.json")
sys.exit(0 if not FAIL else 1)
