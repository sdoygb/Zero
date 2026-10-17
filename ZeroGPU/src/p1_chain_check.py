#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""p1_chain_check.py —— 在单一计数公设 P1 下的完整链（推导稿 §5）。

P1（无量纲、计数）：质量单元 A² = v/C(8,2)² = v/784（28 = 8 槽载体的槽对子数）。
本脚本从 P1 出发算到底：
  带电轻子 m_i = A² f_i²（A_ℓ=√2、δ0 内生）；
  中微子 M_R=√(M_I M_X)、m_D=m_τ/2、m_ν3=m_D²/M_R、形状 R、Σ、mββ；
  与数据的全部偏差。
用途：把"接受 P1 之后会得到什么"变成可复跑的数。

输出：ZeroGPU/results/p1_chain_check.json
用法：python3 ZeroGPU/src/p1_chain_check.py
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


V, M_Z = 246.22, 91.1876
D21, D31 = 7.42e-5, 2.510e-3
Ue2, Ue3 = 0.5443, 0.1492
M_I = M_Z * math.exp(4 * math.pi)
M_X = 1.009e13
M_R = math.sqrt(M_I * M_X)

A_L = math.sqrt(2)
delta0 = math.pi / 12 - 2 / 9
delta_e = 3 * math.pi / 4 - delta0
f = [1 + A_L * math.cos(delta_e + 2 * math.pi * k / 3) for k in range(3)]
f_tau = max(f)

P1 = 1 / math.comb(8, 2) ** 2
A2 = V * P1
m = [A2 * fk**2 for fk in f]
m_tau = max(m)
print(f"   P1: c_ℓ=1/784={P1:.9e}  A²={A2:.6f} GeV  m_τ={m_tau:.6f} GeV")
print(f"   带电轻子（GeV）: {[round(x,6) for x in sorted(m)]}")

m_nu3 = (m_tau / 2) ** 2 / M_R * 1e12          # meV
R = (73 - 28 * math.sqrt(6)) / 25
m_nu2 = R * m_nu3
Sigma = m_nu2 + m_nu3
X = R**2
d21 = (m_nu2 * 1e-3) ** 2
d31 = (m_nu3 * 1e-3) ** 2
m_bb = (Ue2**2 * m_nu2 + Ue3**2 * m_nu3)       # meV

print(f"   m_ν3={m_nu3:.2f} meV  m_ν2={m_nu2:.2f}  Σ={Sigma:.2f} meV  mββ={m_bb:.2f} meV")
print(f"   Δm²21={d21:.3e}  Δm²31={d31:.3e}  X={X:.6f}")

check("m_τ 与测量差 <0.1%", abs(m_tau / 1.77686 - 1) < 1e-3,
      f"{100*(m_tau/1.77686-1):+.4f}%")
check("Σmν 在 DESI 上界（64 meV）内", Sigma < 64, f"{Sigma:.2f} meV")
check("mββ 低于下一代灵敏度（10 meV）", m_bb < 10, f"{m_bb:.2f} meV")
check("X 与数据同量级（差 <10%）", abs(X / (D21 / D31) - 1) < 0.1,
      f"{100*(X/(D21/D31)-1):+.1f}%")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "P1": P1, "A2_GeV": A2, "m_lep_GeV": sorted(m),
       "m_tau_GeV": m_tau, "m_nu3_meV": m_nu3, "m_nu2_meV": m_nu2,
       "Sigma_meV": Sigma, "m_bb_meV": m_bb, "X": X, "R": R}
with open(os.path.join(OUT, "p1_chain_check.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/p1_chain_check.json")
sys.exit(0 if not FAIL else 1)
