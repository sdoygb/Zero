#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""mass_unit_options.py —— 质量单元的三种选项及其后果（推导稿 §4）。

选项：
  O0（现方案）：c_ℓ = m_τ/(v f_τ²)（用测量 m_τ 定标）
  O1（计数公设）：c_ℓ = 1/C(8,2)² = 1/784
  O2（环组合近邻）：c_ℓ = e^{-(2r1+2r2-2r3+r4)}（624 组合里最近，0.93%）
对每个选项算出：带电轻子三代质量、中微子 m3/Σ/mββ、与数据的偏差。
用途：让"定理级还缺什么"变成可比较的数，而不是口号。

输出：ZeroGPU/results/mass_unit_options.json
用法：python3 ZeroGPU/src/mass_unit_options.py
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

r1 = math.log(M_Z / V)
r2 = 4 * math.pi
r3 = math.log(M_X / M_I)
r4 = math.log(9.821e16 / M_X)

c_meas = M_TAU / (V * f_tau**2)
options = {
    "O0 测量定标": c_meas,
    "O1 计数 1/784": 1 / math.comb(8, 2) ** 2,
    "O2 环组合近邻": math.exp(-(2 * r1 + 2 * r2 - 2 * r3 + r4)),
}
print(f"   c_ℓ 测量值 = {c_meas:.9e}")
for name, c in options.items():
    m = [V * c * fk**2 for fk in f]
    # 中微子：m_ν3 = (m_τ/2)²/M_R
    m_nu3 = (m[2] / 2) ** 2 / M_R
    m_nu2 = R * m_nu3
    Sigma = m_nu2 + m_nu3
    print(f"   {name:14s}: c={c:.6e}  m_τ={m[2]:.5f} GeV（测量 1.77686，{100*(m[2]/M_TAU-1):+.3f}%）"
          f"  m_ν3={m_nu3*1e12:.2f} meV  Σ={Sigma*1e12:.2f} meV")

check("O0 的 m_τ 与测量一致（自洽）",
      abs(options["O0 测量定标"] * V * f_tau**2 / M_TAU - 1) < 1e-12)
check("O1 的 m_τ 与测量差 >0.05%（近邻非等式）",
      abs(options["O1 计数 1/784"] * V * f_tau**2 / M_TAU - 1) > 5e-4,
      f"{100*(options['O1 计数 1/784']*V*f_tau**2/M_TAU-1):+.3f}%")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "c_meas": c_meas,
       "options": {k: v for k, v in options.items()}}
with open(os.path.join(OUT, "mass_unit_options.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/mass_unit_options.json")
sys.exit(0 if not FAIL else 1)
