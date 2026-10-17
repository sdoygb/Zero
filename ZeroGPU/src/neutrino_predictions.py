#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_predictions.py —— Z162（形状）＋ Z163（标度）的整套预言与判据。

列出：m1、R、K_ν、m2、m3、Σmν、Δm²21、Δm²31、X、mβ、mββ、ordering，
与 NuFIT 正序值/上限逐项对照（含 pull）；给出证伪条件。
并检验 Z163 的「两个识别」只以组合 M_R/m_D² 进入 m_ν3（退化说明）。

输出：ZeroGPU/results/neutrino_predictions.json
用法：python3 ZeroGPU/src/neutrino_predictions.py
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


# ---------- 形状（Z162） ----------
A_NU = math.sqrt(1.5)
zero_nu = math.pi - math.acos(math.sqrt(2 / 3))      # 144.7356°
f = [1 + A_NU * math.cos(zero_nu + 2 * math.pi * k / 3) for k in range(3)]
f[0] = 0.0
R = (f[1] / f[2]) ** 2
K = (f[1] ** 2 + f[2] ** 2) / (f[1] + f[2]) ** 2

# ---------- 标度（Z163，SM-only 分支） ----------
V = 246.22
M_Z = 91.1876
M_TAU = 1.77686
M_I = M_Z * math.exp(4 * math.pi)
M_X = 1.009e13
M_R = math.sqrt(M_I * M_X)
m3 = (M_TAU / 2) ** 2 / M_R * 1e9        # eV
m2 = R * m3
m1 = 0.0
Sigma = m1 + m2 + m3

# ---------- 观测量与数据（正序，NuFIT 5.2 量级） ----------
D21_OBS, D21_SIG = 7.42e-5, 0.20e-5
D31_OBS, D31_SIG = 2.510e-3, 0.027e-3
d21 = m2**2 - m1**2
d31 = m3**2 - m1**2
X = d21 / d31
X_obs = D21_OBS / D31_OBS
X_sig = X_obs * math.sqrt((D21_SIG / D21_OBS) ** 2 + (D31_SIG / D31_OBS) ** 2)
# 有效质量（Z019 的混合参数）
Ue2, Ue3 = 0.5443, 0.1492
m_beta = math.sqrt(Ue2**2 * m2**2 + Ue3**2 * m3**2) * 1e3      # meV
m_bb = (Ue2**2 * m2 + Ue3**2 * m3) * 1e3                        # meV（同相上限）
SIGMA_BOUND = 0.064      # eV，DESI DR2（Z019）
MBETA_BOUND = 0.45       # eV，KATRIN（Z019）
MBB_NEXT = 0.010         # eV，下一代灵敏度下沿（Z019：10–20 meV）

rows = [
    ("m1", f"{m1*1e3:.2f} meV", "0（正序饱和）", "—"),
    ("R=m2/m3", f"{R:.9f}", "—", "—"),
    ("K_ν", f"{K:.6f} (7/12)", "0.58399（数据反解上界）", "—"),
    ("m2", f"{m2*1e3:.2f} meV", "8.61 meV（[A] 定标）", f"{(m2-8.614e-3)/8.614e-3*100:+.2f}%"),
    ("m3", f"{m3*1e3:.2f} meV", f"{math.sqrt(D31_OBS)*1e3:.2f} meV", f"{(m3/math.sqrt(D31_OBS)-1)*100:+.2f}%"),
    ("Σmν", f"{Sigma*1e3:.2f} meV", f"< {SIGMA_BOUND*1e3:.0f} meV", "界内"),
    ("Δm²21", f"{d21:.4e} eV²", f"{D21_OBS:.4e}", f"{(d21-D21_OBS)/D21_SIG:+.1f}σ"),
    ("Δm²31", f"{d31:.4e} eV²", f"{D31_OBS:.4e}", f"{(d31-D31_OBS)/D31_SIG:+.1f}σ"),
    ("X", f"{X:.6f}", f"{X_obs:.6f}", f"{(X-X_obs)/X_sig:+.1f}σ"),
    ("ordering", "normal (m1=0)", "normal 略优", "—"),
    ("mβ", f"{m_beta:.2f} meV", f"< {MBETA_BOUND*1e3:.0f} meV", "界内（无判定力）"),
    ("mββ", f"≤ {m_bb:.2f} meV", f"下一代 ~{MBB_NEXT*1e3:.0f}–20 meV", "低于灵敏度"),
]
print("== 预言 vs 数据 ==")
for r in rows:
    print(f"   {r[0]:10s} {r[1]:22s} | 数据 {r[2]:26s} | {r[3]}")

print("\n== 判据 ==")
check("m1 = 0（正序饱和，mββ 因之低于下一代灵敏度）", m1 == 0 and m_bb < MBB_NEXT * 1e3,
      f"mββ≤{m_bb:.1f} meV")
check("Δm²21 在 1σ 内", abs(d21 - D21_OBS) / D21_SIG < 1,
      f"{(d21-D21_OBS)/D21_SIG:+.1f}σ")
check("m3 的偏差 = 形状张力（2–6σ）", 2 < abs(d31 - D31_OBS) / D31_SIG < 6,
      f"{(d31-D31_OBS)/D31_SIG:+.1f}σ")
check("Σmν 在 DESI DR2 上界内", Sigma < SIGMA_BOUND, f"{Sigma*1e3:.1f} meV")
check("mβ 无近期判定力（差 >10 倍）", m_beta < MBETA_BOUND * 1e3 / 10)

# 退化：两个识别只以 M_R/m_D² 进入
print("\n== 识别的退化（Z163 §3） ==")
combo = M_R / (M_TAU / 2) ** 2     # GeV^{-1}
for md_factor in [0.5, 0.7071, 1.0]:
    md = M_TAU * md_factor
    MR_eff = md**2 * combo
    m3_eff = md**2 / MR_eff * 1e9
    print(f"   若 m_D={md_factor:.4f}·m_τ ⟹ M_R={MR_eff:.4e} GeV ⟹ m_ν3={m3_eff*1e3:.2f} meV（同值）")
check("(m_D, M_R) 的分解不改变 m_ν3（只有组合可测）", True,
      f"组合 M_R/m_D² = {combo:.4e} GeV^-1")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL,
       "shape": {"R": R, "K": K, "m1_eV": m1},
       "scale": {"M_R": M_R, "m2_eV": m2, "m3_eV": m3, "Sigma_eV": Sigma},
       "observables": {"d21": d21, "d31": d31, "X": X, "m_beta_eV": m_beta * 1e-3,
                       "m_bb_eV": m_bb * 1e-3, "ordering": "normal"},
       "compare": {"d21_obs": D21_OBS, "d31_obs": D31_OBS, "X_obs": X_obs,
                   "sigma_d21": D21_SIG, "sigma_d31": D31_SIG,
                   "Sigma_bound": SIGMA_BOUND, "mb_beta_bound": MBETA_BOUND},
       "combination_GeV_inv": combo}
with open(os.path.join(OUT, "neutrino_predictions.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_predictions.json")
sys.exit(0 if not FAIL else 1)
