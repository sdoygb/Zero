#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""half_slot_chain_check.py —— 两条半槽修正下的完整链（推导稿 §2 更新）。

D1 : c_ℓ = 1/(C(8,2)² + 1/2) = 1/784.5                 （8 槽：槽对子² + 半槽）
D3': M_R = √(M_I M_X)·(1 − 1/(2·16)) = √(M_I M_X)·31/32 （16 槽：半槽/16 = 1/32）
D2 : m_D = m_τ/2                                        （半满轨道占据；Lean 已证）

检查：三代带电轻子、m_ν3、Σ、mββ、Δm²21/31、X 与数据的偏差（含 pull）。

输出：ZeroGPU/results/half_slot_chain_check.json
用法：python3 ZeroGPU/src/half_slot_chain_check.py
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
M_I = M_Z * math.exp(4 * math.pi)
M_X = 1.009e13
M_R0 = math.sqrt(M_I * M_X)
ALPHA = 1 - 1 / (2 * 16)          # 31/32
M_R = M_R0 * ALPHA

A = math.sqrt(2)
d0 = math.pi / 12 - 2 / 9
de = 3 * math.pi / 4 - d0
f = [1 + A * math.cos(de + 2 * math.pi * k / 3) for k in range(3)]
c_ell = 1 / (math.comb(8, 2) ** 2 + 0.5)
m_lep = [V * c_ell * fk**2 for fk in f]
m_tau = m_lep[2]
print(f"   D1 c_ℓ={c_ell:.9e}  三代={[round(x,8) for x in m_lep]}")
print(f"   D3' M_R=√(M_I M_X)·31/32 = {M_R:.6e}（未修正 {M_R0:.6e}，修正 {100*(ALPHA-1):+.3f}%）")

R = (73 - 28 * math.sqrt(6)) / 25
m_nu3 = (m_tau / 2) ** 2 / M_R * 1e12
m_nu2 = R * m_nu3
Sigma = m_nu2 + m_nu3
m_bb = (0.5443**2 * m_nu2 + 0.1492**2 * m_nu3)
X = R**2
d21 = (m_nu2 * 1e-3) ** 2
d31 = (m_nu3 * 1e-3) ** 2
print(f"   m_ν3={m_nu3:.2f} meV  m_ν2={m_nu2:.2f}  Σ={Sigma:.2f}  mββ={m_bb:.2f}  X={X:.6f}")
print(f"   Δm²21={d21:.3e}  Δm²31={d31:.3e}")

m_data = [0.000510999, 0.1056584, 1.77686]
devs = [100 * (a / b - 1) for a, b in zip(m_lep, m_data)]
print(f"   三代偏差% {[round(x,4) for x in devs]}")
check("三代带电轻子偏差 <0.01%", max(abs(x) for x in devs) < 0.01,
      f"max={max(abs(x) for x in devs):.4f}%")
check("m_ν3 与 √Δm²31 差 <0.5%", abs(m_nu3 / 50.10 - 1) < 0.005,
      f"{100*(m_nu3/50.10-1):+.3f}%")
check("Δm²31 的 pull 在 ±1σ 内", abs(d31 - 2.510e-3) / 0.027e-3 < 1,
      f"{(d31-2.510e-3)/0.027e-3:+.2f}σ")
check("Σmν 在 DESI 上界（64 meV）内", Sigma < 64, f"{Sigma:.2f}")
check("mββ 低于下一代灵敏度（10 meV）", m_bb < 10, f"{m_bb:.2f}")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "c_ell": c_ell, "m_lep_GeV": m_lep,
       "m_tau_GeV": m_tau, "M_R0": M_R0, "alpha": ALPHA, "M_R": M_R,
       "m_nu3_meV": m_nu3, "m_nu2_meV": m_nu2, "Sigma_meV": Sigma,
       "m_bb_meV": m_bb, "X": X, "d21": d21, "d31": d31, "dev_lep_pct": devs}
with open(os.path.join(OUT, "half_slot_chain_check.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/half_slot_chain_check.json")
sys.exit(0 if not FAIL else 1)
