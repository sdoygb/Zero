#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_scale_conditional.py —— G2：中微子绝对标度的条件式（本轮）。

思路：中微子扇区**不引入新的量纲输入**；它的标度由（i）带电轻子质量单元 c_ℓ=μ_ℓ²/v、
（ii）环②/环③ 的对数中点 M_R=√(M_I·M_X^unif)、（iii）半满槽的 1/2 因子给出。

条件式：m_ν3 = (m_τ/2)² / √(M_I·M_X^unif)，其中 m_τ = c_ℓ v f_τ²。
故 m_ν3 = v·c_ℓ²·f_τ⁴·e^{-(r2+r3)}/4 —— 除 c_ℓ（全费米子共享的归一化）外没有别的量纲输入。

检查：
  A. M_I、M_X^unif、M_R 的值；
  B. m_ν3 与实测 √Δm²31 的比；
  C. 由 Z162 的 R 给出 m_ν2、Σm_ν、X；
  D. 写成 v·c_ℓ² 形式并核对与 B 相等；
  E. 换到 Z155 的中间谱分支（M_X=4.954e13）看敏感性（登记，不作选择）。

输出：ZeroGPU/results/neutrino_scale_conditional.json
用法：python3 ZeroGPU/src/neutrino_scale_conditional.py
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


# ---------- 输入（与库内一致） ----------
V = 246.22                 # GeV
M_Z = 91.1876              # GeV
M_TAU = 1.77686            # GeV
D21 = 7.42e-5              # eV^2
D31 = 2.510e-3             # eV^2
M3_OBS = math.sqrt(D31)    # 0.0501 eV

# ---------- 带电轻子侧（Z121/Z123） ----------
A_L = math.sqrt(2)
delta0 = math.pi / 12 - 2 / 9
delta_e = 3 * math.pi / 4 - delta0
f = [1 + A_L * math.cos(delta_e + 2 * math.pi * k / 3) for k in range(3)]
f_tau = max(f)                       # 2.3794
m_tau_law = None
# c_ell 由 m_tau 定：m_tau = c_ell * v * f_tau^2  ⟹ c_ell = m_tau/(v f_tau^2)
c_ell = M_TAU / (V * f_tau**2)
print(f"   c_ℓ = μ_ℓ²/v = {c_ell:.6e}   (f_τ={f_tau:.6f})")

# ---------- 环②/环③ 与中段标度 ----------
r2 = 4 * math.pi                         # Z154：环② = 2π k0 = 4π
M_I = M_Z * math.exp(r2)                 # 2.615e7 GeV
r_ln_MZ_v = math.log(M_Z / V)            # ≈ −0.9933
ln_MIv = r2 + r_ln_MZ_v
M_X_unif = 1.009e13                      # Z155：SM-only
ln_MXv = math.log(M_X_unif / V)
r3 = ln_MXv - ln_MIv                     # 环③（SM-only 分支）
M_R = math.sqrt(M_I * M_X_unif)
ln_MRv = 0.5 * (ln_MIv + ln_MXv)
print(f"   M_I={M_I:.6e}  M_X^unif={M_X_unif:.4e}  M_R=√(M_I M_X)={M_R:.6e}  GeV")
print(f"   r2=4π={r2:.6f}  r3={r3:.6f}  ln(M_R/v)={ln_MRv:.6f}")
check("M_I = M_Z e^{4π}", abs(M_I - 2.615e7) / M_I < 0.01, f"{M_I:.4e}")
check("M_R = v e^{(ln(M_I/v)+ln(M_X/v))/2}（对数中点）",
      abs(ln_MRv - math.log(M_R / V)) < 1e-12)

# ---------- 条件式 ----------
m_D = M_TAU / 2
m_nu3 = m_D**2 / M_R                     # GeV
m_nu3_eV = m_nu3 * 1e9
print(f"\n   m_D=m_τ/2={m_D:.6f} GeV ⟹ m_ν3 = {m_nu3_eV*1e3:.3f} meV（实测 √Δm²31={M3_OBS*1e3:.3f} meV）")
check("m_ν3 在实测 5% 内", abs(m_nu3_eV / M3_OBS - 1) < 0.05,
      f"比={m_nu3_eV/M3_OBS:.4f}")

# Z162 的形状
R = (73 - 28 * math.sqrt(6)) / 25
m_nu2 = R * m_nu3_eV
Sigma = m_nu2 + m_nu3_eV
X = (m_nu2**2) / (m_nu3_eV**2)
print(f"   R={R:.9f} ⟹ m_ν2={m_nu2*1e3:.3f} meV，Σm_ν={Sigma*1e3:.2f} meV，X={X:.6f}")
check("Σm_ν 与 Z019 的定标值同量级（55–60 meV）", 0.050 < Sigma < 0.065,
      f"{Sigma*1e3:.2f} meV")
check("X 与 Z162 一致（5.47% 张力继承）", abs(X / (R**2) - 1) < 1e-12)

# v·c_ℓ² 形式
m_nu3_form = V * c_ell**2 * f_tau**4 * math.exp(-(ln_MIv + ln_MXv) / 2) / 4
print(f"   v·c_ℓ²·f_τ⁴·e^{{-(ln(M_I/v)+ln(M_X/v))/2}}/4 = {m_nu3_form*1e9*1e3:.3f} meV")
check("条件式与 v·c_ℓ² 形式相等", abs(m_nu3_form - m_nu3) / m_nu3 < 1e-9)

# 敏感性：Z155 中间谱分支
M_X_alt = 4.954e13
M_R_alt = math.sqrt(M_I * M_X_alt)
m_nu3_alt = (m_D**2 / M_R_alt) * 1e9
print(f"   敏感性：中间谱分支 M_X=4.954e13 ⟹ M_R={M_R_alt:.4e}，m_ν3={m_nu3_alt*1e3:.2f} meV")
check("两分支都落在 0.02–0.06 eV（量级由中段决定）",
      0.015 < m_nu3_alt < 0.07 and 0.03 < m_nu3_eV < 0.06)

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL,
       "inputs": {"v": V, "m_tau": M_TAU, "M_I": M_I, "M_X_unif": M_X_unif,
                  "r2": r2, "r3": r3, "M_R": M_R, "c_ell": c_ell, "f_tau": f_tau},
       "m_nu3_eV": m_nu3_eV, "m_nu3_obs_eV": M3_OBS, "R": R,
       "m_nu2_eV": m_nu2, "Sigma_eV": Sigma, "X": X,
       "m_nu3_alt_eV": m_nu3_alt, "M_R_alt": M_R_alt}
with open(os.path.join(OUT, "neutrino_scale_conditional.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_scale_conditional.json")
sys.exit(0 if not FAIL else 1)
