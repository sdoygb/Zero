#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_proton_mass.py —— 质子单粒子路径: 质子=三色闭合结构, 直接套质量机制 m=A²M(θ)²
======================================================================================
问题: 若质子共享从电子定的标度 A, 反解质子角 θ_p, 看是否落在 Zero 内生格点上.
机制(Z121): m=A²M(θ)², M(θ)=1+√2cosθ; 电子 θ_e=2π/3+2/9, δ₀=π/12−2/9.
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

me, m_p = 0.51099895, 938.272088    # MeV
P_("=" * 96)
P_("质子单粒子路径 · k4_proton_mass.py")
P_("=" * 96)

# ============ 1. 从电子定标度 A ============
theta_e = 2*np.pi/3 + 2/9            # 电子角 = 轮转 + 两支积
Me = 1 + np.sqrt(2)*np.cos(theta_e)
A = np.sqrt(me)/Me
P_("电子: θ_e=%.6f rad=%.3f°  M(θ_e)=%.6f  √m_e=%.6f" % (theta_e, np.degrees(theta_e), Me, np.sqrt(me)))
P_("标度 A = √m_e/M(θ_e) = %.6f MeV" % A)
P_("A² = %.6f" % A**2)

# ============ 2. 反解质子角 ============
Mp_target = np.sqrt(m_p)/A
# 1+√2cosθ = Mp_target → cosθ = (Mp_target-1)/√2
cos_p = (Mp_target - 1)/np.sqrt(2)
theta_p = np.arccos(np.clip(cos_p, -1, 1))
P_("\n质子: 需 M(θ_p)=%.6f (=√938.272/A)" % Mp_target)
P_("cosθ_p = %.6f → θ_p = %.6f rad = %.4f°" % (cos_p, theta_p, np.degrees(theta_p)))
P_("验证: 1+√2cosθ_p=%.6f, A²M²=%.3f MeV (m_p=%.3f)" % (1+np.sqrt(2)*np.cos(theta_p), A**2*(1+np.sqrt(2)*np.cos(theta_p))**2, m_p))

# 另一解 (cos 对称, 2π−θ)
P_("  镜像解 θ_p' = 2π−θ_p = %.4f°" % np.degrees(2*np.pi-theta_p))

# ============ 3. θ_p vs Zero 内生候选角 ============
P_("\n" + "=" * 96)
P_("【3. θ_p=%.4f° 与 Zero 内生候选角对比】" % np.degrees(theta_p))
tau = 2*np.pi; two9 = 2/9; d0 = np.pi/12 - two9
cands = {
    "格点 45°  (π/4, seed8)": np.pi/4,
    "格点 60°  (π/3, 1/6圈)": np.pi/3,
    "格点 90°  (π/2)": np.pi/2,
    "三代间隔 120°": 2*np.pi/3,
    "零点 135° (3π/4)": 3*np.pi/4,
    "电子角 132.7°": theta_e,
    "5×(2/9)=63.66°": 5*two9,
    "π/3 + 2δ₀": np.pi/3 + 2*d0,
    "π/4 + 2δ₀": np.pi/4 + 2*d0,
    "2π/9 + 2π/9=80°": 2*2*np.pi/9,
}
res = []
for nm, val in cands.items():
    diff = abs(np.degrees(val) - np.degrees(theta_p))
    # 也允许镜像/模 2π
    diff2 = min(diff, abs(np.degrees(2*np.pi-val) - np.degrees(theta_p)))
    res.append((min(diff,diff2), nm, val))
res.sort()
P_("  %-24s %12s %8s" % ("候选", "角", "偏差°"))
for diff, nm, val in res[:8]:
    P_("  %-24s %8.3f° %8.3f" % (nm, np.degrees(val), diff))

# ============ 4. 质子质量比 (若用某候选角) ============
P_("\n" + "=" * 96)
P_("【4. 候选角给出的质子质量 (共享标度 A)】")
best_nm, best_val = res[0][1], res[0][2]
M_best = 1 + np.sqrt(2)*np.cos(best_val)
m_pred = A**2*M_best**2
P_("  最佳候选 %s: M=%.5f → m=%.1f MeV (观测 938.3)" % (best_nm, M_best, m_pred))

# 若质子角 = π/3 (60°)
for nm, val in cands.items():
    Mv = 1 + np.sqrt(2)*np.cos(val)
    if Mv > 0:
        P_("  %-24s M=%.5f m=%.1f MeV (偏差 %.1f%%)" % (nm, Mv, A**2*Mv**2, abs(A**2*Mv**2-m_p)/m_p*100))

out = {"theta_p_rad": theta_p, "theta_p_deg": np.degrees(theta_p),
       "Mp_target": Mp_target, "A": A, "candidates": [{"name":n,"deg":np.degrees(v),"diff_deg":round(d,3)} for d,n,v in res[:8]]}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_proton_mass.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_proton_mass.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_proton_mass.json / .txt")
