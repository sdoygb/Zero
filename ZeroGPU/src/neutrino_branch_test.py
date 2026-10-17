#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_branch_test.py —— Z163 条件式的分支判据：SM-only vs 中间谱。

条件式 m_ν3 = (m_τ/2)²/√(M_I M_X^unif) 对 Z155 的两个中间谱分支给两个预言：
  * SM-only（M_X = 1.009e13 GeV）  ⟹ m_ν3 ≈ 48.6 meV
  * 中间谱（M_X = 4.954e13 GeV）   ⟹ m_ν3 ≈ 21.9 meV
用实测 Δm²31 与 Σmν 上界判哪一个分支与数据相容；并给出 σ 距离（按 Δm² 的 1.1% 误差传播）。

输出：ZeroGPU/results/neutrino_branch_test.json
用法：python3 ZeroGPU/src/neutrino_branch_test.py
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
M_Z = 91.1876
M_TAU = 1.77686
D21 = 7.42e-5
D31 = 2.510e-3
SIG_D31 = 0.027e-3        # NuFIT 5.2 量级（1.1%）
M3_OBS = math.sqrt(D31)
SIG_M3 = 0.5 * SIG_D31 / D31 * M3_OBS   # m3 ∝ sqrt(D31)
SIGMA_COSMO = 0.064       # eV，DESI DR2 联合（Z019）

A_L = math.sqrt(2)
delta0 = math.pi / 12 - 2 / 9
delta_e = 3 * math.pi / 4 - delta0
f_tau = max(1 + A_L * math.cos(delta_e + 2 * math.pi * k / 3) for k in range(3))
c_ell = M_TAU / (V * f_tau**2)

M_I = M_Z * math.exp(4 * math.pi)
r3 = math.log(1.009e13 / M_I)     # SM-only 的环③
BRANCHES = [
    ("SM-only", 1.009e13, r3),
    ("中间谱", 4.954e13, math.log(4.954e13 / M_I)),
]

R = (73 - 28 * math.sqrt(6)) / 25
rows = []
for name, MX, r in BRANCHES:
    MR = math.sqrt(M_I * MX)                     # = M_X e^{-r/2}
    m3 = (M_TAU / 2) ** 2 / MR * 1e9             # eV
    m2 = R * m3
    sigma = m2 + m3
    d21 = m2**2
    d31 = m3**2
    z = (m3 - M3_OBS) / SIG_M3
    rows.append({"branch": name, "M_X": MX, "ring3": r, "M_R": MR,
                 "m3_eV": m3, "m2_eV": m2, "Sigma_eV": sigma,
                 "d21": d21, "d31": d31, "sigma_m3": SIG_M3, "z_m3": z,
                 "sigma_vs_cosmo_ratio": sigma / SIGMA_COSMO})
    print(f"   {name:6s}: M_R={MR:.4e} 环③={r:.4f} ⟹ "
          f"m3={m3*1e3:6.2f} meV  m2={m2*1e3:5.2f}  Σ={sigma*1e3:6.2f} meV  "
          f"(Δm²31={d31:.3e} vs {D31:.3e})  z={z:+.1f}")

sm = rows[0]
mid = rows[1]
print(f"\n   实测 √Δm²31 = {M3_OBS*1e3:.2f} meV（σ≈{SIG_M3*1e3:.2f} meV）")

check("SM-only 分支在实测 5% 内（与 Z163 一致）", abs(sm["m3_eV"] / M3_OBS - 1) < 0.05,
      f"偏差 {100*(sm['m3_eV']/M3_OBS-1):+.2f}%")
check("中间谱分支与实测差 > 50%", abs(mid["m3_eV"] / M3_OBS - 1) > 0.5,
      f"偏差 {100*(mid['m3_eV']/M3_OBS-1):+.1f}%")
check("两分支的 Σmν 都在 DESI DR2 上界内（现数据不能只靠 Σ 区分）",
      max(r["Sigma_eV"] for r in rows) < SIGMA_COSMO,
      f"Σ: {sm['Sigma_eV']*1e3:.1f} vs {mid['Sigma_eV']*1e3:.1f} meV")
check("条件式在 SM-only 分支上的张力为 2–6σ（形状张力的同一处）",
      2 < abs(sm["z_m3"]) < 6, f"z={sm['z_m3']:+.1f}")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "target_eV": M3_OBS, "sigma_m3_eV": SIG_M3,
       "R": R, "c_ell": c_ell, "f_tau": f_tau, "M_I": M_I, "branches": rows,
       "cosmo_bound_eV": SIGMA_COSMO}
with open(os.path.join(OUT, "neutrino_branch_test.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_branch_test.json")
sys.exit(0 if not FAIL else 1)
