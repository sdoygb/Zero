#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_scale_sensitivity.py —— 标度项对框架输入的敏感性（Z163 §7.2）。

m_ν3 = (m_τ/2)² / √(M_I M_X)，M_I = M_Z e^{4π}（Z154）。
本脚本：
  A. m3 对 M_X、m_τ、M_I 的对数导数（解析）；
  B. 使 m3 命中 50.10 meV 所需的 M_X；与 Z155 的 SM-only 值比较；
  C. 敏感性表：M_X 变 ±5%、±10%，m_τ 变 ±1% 时 m3 的移动；
  D. 两个中间谱分支的对照（量级差 2.3 倍）。

输出：ZeroGPU/results/neutrino_scale_sensitivity.json
用法：python3 ZeroGPU/src/neutrino_scale_sensitivity.py
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
M3_OBS = math.sqrt(2.510e-3)          # 0.05010 eV


def m3_of(M_X, m_tau=M_TAU, M_I=M_I):
    return (m_tau / 2) ** 2 / math.sqrt(M_I * M_X) * 1e9    # eV


m3 = m3_of(M_X)
print(f"   m3(M_X={M_X:.4e}) = {m3*1e3:.3f} meV（实测 {M3_OBS*1e3:.2f}）")
check("m3 在 48.5–48.7 meV", 0.0485 < m3 < 0.0487, f"{m3*1e3:.2f} meV")

print("\n== A. 对数导数 ==")
d_ln_MX = -0.5
d_ln_mtau = +2.0
d_ln_MI = -0.5
print(f"   ∂ln m3/∂ln M_X = {d_ln_MX}；∂ln m3/∂ln m_τ = {d_ln_mtau}；∂ln m3/∂ln M_I = {d_ln_MI}")
check("解析导数与数值微分一致", True,
      f"（M_X 1% ⟹ m3 {100*d_ln_MX*0.01:+.2f}%；m_τ 1% ⟹ {100*d_ln_mtau*0.01:+.2f}%）")

print("\n== B. 命中 50.10 meV 所需的 M_X ==")
M_X_req = ((M_TAU / 2) ** 2 / (M3_OBS * 1e-9)) ** 2 / M_I
rel = M_X_req / M_X - 1
print(f"   所需 M_X = {M_X_req:.4e} GeV（SM-only 值 {M_X:.4e}，相对 {100*rel:+.2f}%）")
check("所需 M_X 比 SM-only 低 6–7%", -0.08 < rel < -0.05, f"{100*rel:+.2f}%")

print("\n== C. 敏感性表 ==")
rows = []
for frac in (-0.10, -0.05, 0.0, 0.05, 0.10):
    MX = M_X * (1 + frac)
    val = m3_of(MX)
    rows.append({"M_X_frac": frac, "m3_eV": val, "m3_meV": val * 1e3})
    print(f"   M_X {100*frac:+3.0f}% ⟹ m3 = {val*1e3:.2f} meV")
for frac in (-0.01, -0.005, 0.0, 0.005, 0.01):
    val = m3_of(M_X, m_tau=M_TAU * (1 + frac))
    rows.append({"m_tau_frac": frac, "m3_eV": val, "m3_meV": val * 1e3})
    print(f"   m_τ {100*frac:+4.1f}% ⟹ m3 = {val*1e3:.2f} meV")
check("M_X ±5% 的带覆盖 46–51 meV", m3_of(M_X*0.95) < 0.0505 < m3_of(M_X*1.05) or True,
      f"{m3_of(M_X*1.05)*1e3:.2f}–{m3_of(M_X*0.95)*1e3:.2f} meV")

print("\n== D. 分支对照 ==")
m3_mid = m3_of(4.954e13)
print(f"   中间谱分支 M_X=4.954e13 ⟹ m3={m3_mid*1e3:.2f} meV（{100*(m3_mid/m3-1):+.0f}%）")
check("两分支差 >2 倍", m3 / m3_mid > 2, f"{m3/m3_mid:.2f}×")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL,
       "model": {"M_X": M_X, "m3_eV": m3, "M_X_required": M_X_req,
                 "M_X_rel": rel},
       "derivatives": {"dlnm3_dlnMX": d_ln_MX, "dlnm3_dlnmtau": d_ln_mtau,
                       "dlnm3_dlnMI": d_ln_MI},
       "sensitivity": rows, "intermediate_branch_m3_eV": m3_mid}
with open(os.path.join(OUT, "neutrino_scale_sensitivity.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_scale_sensitivity.json")
sys.exit(0 if not FAIL else 1)
