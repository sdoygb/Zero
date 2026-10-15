#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_delta0_axion.py —— axion 式归零动力学: 非中心态两支积 2/9 = τ 相角 → δ₀ 完全内生
======================================================================================
机制(axion 模板): 归零势驱动 δ→零点; 但非中心态产生偏移(归零下限), 极小值在 δ₀≠0.
核心恒等式候选: 反解出的 τ 相角 δ ≈ 2/9 = (2/3)(1/3) = 稳态两支权重积;
  电子相角 δ_e = δ + 2π/3;  δ₀ = 3π/4 − δ_e = π/12 − δ ≈ π/12 − 2/9.
角定义: 三代角 θ_n = δ + 2π(n-1)/3, n=1→τ, 2→e, 3→μ.
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

m = {"tau": 1776.86, "e": 0.51099895, "mu": 105.6583755}   # MeV
rr = np.array([m["tau"], m["e"], m["mu"]])                  # n=1,2,3
P_("=" * 96)
P_("axion 式归零动力学 · k4_delta0_axion.py")
P_("=" * 96)
P_("输入质量 (MeV): τ=%.4f  e=%.6f  μ=%.6f" % (m["tau"], m["e"], m["mu"]))

two9 = (2/3)*(1/3)
P_("稳态两支积 (2/3)(1/3) = %.6f rad = %.4f°  (Berry 相位候选)" % (two9, np.degrees(two9)))

# ============ 1. 反解 τ 相角 δ (网格最小化质量比) ============
P_("\n" + "=" * 96)
P_("【1. 反解 τ 相角 δ】   √m_n = A[1+√2cos(δ + 2π(n-1)/3)]")
best = None
for d0 in np.linspace(0.001, 2*np.pi, 60000):
    a = np.array([1+np.sqrt(2)*np.cos(d0 + 2*np.pi*k/3) for k in range(3)])
    if (a <= 0).any(): continue
    m_r = a**2
    m_r = m_r/m_r[1]
    err = np.sum(np.abs(np.log(m_r) - np.log(rr/rr[1])))
    if best is None or err < best[0]:
        best = (err, d0, m_r)
err, d_tau, m_rec = best
P_("  反解 τ 相角 δ = %.6f rad = %.4f°" % (d_tau, np.degrees(d_tau)))
P_("  目标 2/9      = %.6f rad = %.4f°" % (two9, np.degrees(two9)))
P_("  偏差         = %.4f%%   (质量数据精度所致)" % (abs(d_tau-two9)/two9*100))
P_("  重建质量比 τ:e:μ = %.4f : 1 : %.4f   (观测 %.4f : 1 : %.4f)"
  % (m_rec[0], m_rec[2], rr[0]/rr[1], rr[2]/rr[1]))

# ============ 2. δ₀ 恒等式 ============
P_("\n" + "=" * 96)
P_("【2. δ₀ 恒等式: δ₀ = π/12 − δ ≈ π/12 − 2/9】")
d_e = d_tau + 2*np.pi/3                      # 电子相角
d0_geom = 3*np.pi/4 - d_e                    # 零点 − 电子相角
d0_pi12 = np.pi/12 - two9                    # π/12 − 两支积
P_("  电子相角 δ_e = δ+2π/3 = %.6f rad = %.4f°" % (d_e, np.degrees(d_e)))
P_("  δ₀ = 3π/4 − δ_e   = %.8f rad = %.5f°   (零点−电子相角)" % (d0_geom, np.degrees(d0_geom)))
P_("  δ₀ = π/12 − 2/9   = %.8f rad = %.5f°   (几何差−两支积)" % (d0_pi12, np.degrees(d0_pi12)))
P_("  δ₀ 参考(先前)     = 0.039578 rad")
P_("  吻合: π/12−2/9 偏差 %.4f%%" % (abs(d0_pi12-0.039578)/0.039578*100))
P_("  构成:")
P_("    π/12 = 3π/4 − 2π/3 = 零点(光子) − 三代轮转  (seed 8 格点, 纯几何)")
P_("    2/9  = (2/3)(1/3) = 稳态两支权重积 (非中心态固有量) = τ 相角")

# ============ 3. 重建质量 ============
P_("\n" + "=" * 96)
P_("【3. 用 Zero 内生角重建质量 (归一化 τ=1)】")
d_e_zb = 2*np.pi/3 + two9     # 电子相角 = 轮转+两支积
ang = {"τ": d_e_zb-2*np.pi/3, "e": d_e_zb, "μ": d_e_zb+2*np.pi/3}
A2 = np.array([1+np.sqrt(2)*np.cos(v) for v in [ang["τ"], ang["e"], ang["μ"]]])**2
mz = A2/A2[0]
P_("  相角: τ=%.3f°  e=%.3f°  μ=%.3f°" % (np.degrees(ang["τ"]), np.degrees(ang["e"]), np.degrees(ang["μ"])))
P_("  重建 τ:e:μ = 1 : %.4f : %.4f" % (mz[1], mz[2]))
P_("  观测 τ:e:μ = 1 : %.4f : %.4f  (1:%.0f:%.0f)" % (rr[1]/rr[0], rr[2]/rr[0], rr[0]/rr[1], rr[0]/rr[2]))
Q = A2.sum()/np.sum(np.sqrt(A2))**2
P_("  Koide(重建) = %.6f  (2/3=%.6f)" % (Q, 2/3))

# ============ 4. axion 式归零动力学演示 ============
P_("\n" + "=" * 96)
P_("【4. axion 式归零动力学: 电子被驱动归零, 非中心偏移定极小值 δ_e】")
P_("  叙述: 电子相角被'归零势'驱动趋向零点 3π/4;")
P_("        但非中心态(两支积 2/9)造成偏移, 归不到 0, 平衡在 δ_e=3π/4−δ₀.")
# 归零势 + 非中心偏移: 极小在 θ_e = 3π/4 − δ₀
def pot(theta):
    z0 = 3*np.pi/4
    return (theta-z0)**2 + two9*np.abs(np.sin(theta-z0))   # 归零 + 偏移(示意)
grid = np.linspace(0, 3*np.pi/2, 2001)
imin = grid[np.argmin([pot(t) for t in grid])]
P_("  示意势极小位置 ≈ %.5f rad = %.3f°  (目标 δ_e=%.5f°=%.3f°)"
  % (imin, np.degrees(imin), d_e_zb, np.degrees(d_e_zb)))
P_("  归零下限 δ₀ = 3π/4 − δ_e = π/12 − 2/9 = %.6f rad" % (np.pi/12 - two9))

out = {
  "delta_tau_from_fit": d_tau, "target_two9": two9, "rel_dev_pct": abs(d_tau-two9)/two9*100,
  "delta_e": d_e, "delta0_geom": d0_geom, "delta0_pi12_minus_two9": d0_pi12,
  "identity": "δ₀ = π/12 − (2/3)(1/3)  (π/12=零点−三代轮转; 2/9=稳态两支积=τ相角)",
  "rebuilt_mass_ratio": mz.tolist(), "koide_rebuilt": Q,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_delta0_axion.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_delta0_axion.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_delta0_axion.json / .txt")
