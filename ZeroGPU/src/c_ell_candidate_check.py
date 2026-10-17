#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""c_ell_candidate_check.py —— c_ℓ 的计数候选：1/C(8,2)² = 1/784（Z166 §2.1）。

候选：c_ℓ = 1 / C(8,2)² = 1/28² = 1/784，其中 8 是 8 槽载体 Γ8=K4□K2 的槽数、
28=|E(K8)|（槽对子的数目）。等价写法 1/784 = 1/(16·7²)（16 槽 × 7 个满轨道²）。

检查：
  A. 候选值与 c_ℓ=m_τ/(v f_τ²) 的相对差；
  B. 若候选为精确式，反推的 m_τ 与测量值差多少（σ）；
  C. 邻近的其它计数候选（1/783、1/785、1/(8·98)、1/(C(8,2)·C(8,3))）——对照。

输出：ZeroGPU/results/c_ell_candidate_check.json
用法：python3 ZeroGPU/src/c_ell_candidate_check.py
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
M_TAU = 1.77686
SIG_MTAU = 0.00012
A_L = math.sqrt(2)
delta0 = math.pi / 12 - 2 / 9
delta_e = 3 * math.pi / 4 - delta0
f = [1 + A_L * math.cos(delta_e + 2 * math.pi * k / 3) for k in range(3)]
f_tau = max(f)
c_ell = M_TAU / (V * f_tau**2)

cand = 1 / math.comb(8, 2) ** 2
print(f"   c_ℓ = {c_ell:.9e}   候选 1/784 = {cand:.9e}   相对差 {100*(cand/c_ell-1):+.4f}%")
m_tau_if = cand * V * f_tau**2
print(f"   若候选精确：m_τ = {m_tau_if:.6f} GeV（测量 {M_TAU} ± {SIG_MTAU}）")
z = (m_tau_if - M_TAU) / SIG_MTAU
print(f"   ⟹ 偏差 {100*(m_tau_if/M_TAU-1):+.4f}% = {z:+.1f}σ")
check("候选相对差 <0.1%", abs(cand / c_ell - 1) < 1e-3, f"{100*(cand/c_ell-1):+.4f}%")
check("候选若精确，m_τ 偏差 <0.1%（虽 >5σ）", abs(m_tau_if / M_TAU - 1) < 1e-3,
      f"{100*(m_tau_if/M_TAU-1):+.4f}%")

print("\n   邻近候选对照：")
for name, val in [("1/783", 1 / 783), ("1/785", 1 / 785),
                  ("1/(8·98)", 1 / 784), ("1/(C(8,2)·C(8,3))", 1 / (28 * 56)),
                  ("1/(16·49)", 1 / (16 * 49)), ("1/28²", 1 / 28**2)]:
    print(f"     {name:20s} = {val:.9e}  差 {100*(val/c_ell-1):+.4f}%")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "c_ell": c_ell, "candidate": cand,
       "rel": cand / c_ell - 1, "m_tau_if_candidate": m_tau_if,
       "m_tau_sigma": z, "f_tau": f_tau}
with open(os.path.join(OUT, "c_ell_candidate_check.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/c_ell_candidate_check.json")
sys.exit(0 if not FAIL else 1)
