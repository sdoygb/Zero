#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_proton_ellipse.py —— 质子三点结合 → 圆变形为椭圆(偏心), 质量中心偏移 0.35% 修正
======================================================================================
用户图像: 三个色荷结合成质子, 结合能导致圆→椭圆(偏心), 质量中心偏离零点一点点 (0.35%).
目标: 找椭圆变形量(偏心e)内生来源, 使合成从 3 → 2.9895 (0.35% 修正).
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

me, m_p = 0.51099895, 938.272088
theta_e = 2*np.pi/3 + 2/9
A = np.sqrt(me)/(1+np.sqrt(2)*np.cos(theta_e)); A2 = A**2
S0 = m_p/A2   # 2.989468
D = 3 - S0    # 0.010532
P_("=" * 96)
P_("质子椭圆变形修正 · k4_proton_ellipse.py")
P_("=" * 96)
P_("目标合成 S0=%.6f (修正 Δ=3−S0=%.6f, 0.35%%)" % (S0, D))

two9=2/9; d0=np.pi/12-two9
# 椭圆偏心率候选的内生量
e_cands = {
    "δ₀(rad)": d0, "2/9(rad)": two9, "√(2δ₀)": np.sqrt(2*d0),
    "2δ₀": 2*d0, "δ₀/2": d0/2, "2/9·δ₀": two9*d0,
    "1/√8": 1/np.sqrt(8), "π/8": np.pi/8, "π/12": np.pi/12,
}

P_("\n" + "="*96)
P_("【A. 椭圆几何量给修正: 需要的偏心 e】")
# 模型1: 周长修正 C=2πa(1−e²/4), 若合成 ∝ 周长比 → 需要 (1−e²/4)=S0/3
for model, need_factor in [("周长比(合成∝周长)", S0/3), ("M∝√周长", np.sqrt(S0/3))]:
    e2 = 4*(1-need_factor)
    if e2 >= 0:
        e = np.sqrt(e2)
        P_("  模型 %s: 需 e²=%.5f → e=%.4f" % (model, e2, e))
# 模型2: 面积比 A=πab, a=1+ε,b=1−ε, 面积不变 → 周长变
# 模型3: 形心偏移

P_("\n" + "="*96)
P_("【B. 三色荷在椭圆上: 合成量随偏心 e 变化】")
def ellipse_sum(e, theta):
    # 椭圆 x=a cosθ, y=b sinθ, a=1/c? 用参数 a=√(1+e), b=√(1-e) 保持面积 πab=π
    a = np.sqrt(1+e); b = np.sqrt(1-e)
    r = np.sqrt((a*np.cos(theta))**2 + (b*np.sin(theta))**2)  # 径向
    M = 1 + np.sqrt(2)*np.cos(theta)
    return r, M
# 三个色荷参数角 0,120,240, 沿椭圆某主轴方向
for ph0 in [0, np.pi/4, np.pi/8, 3*np.pi/8]:
    th = np.array([ph0, ph0+2*np.pi/3, ph0+4*np.pi/3])
    for emod, egrid in [("e=δ₀", [d0]), ("e=2δ₀", [2*d0]), ("e=2/9", [two9]), ("e扫描", np.linspace(0,0.3,301))]:
        for e in (egrid if isinstance(egrid,list) else egrid):
            r, M = ellipse_sum(e, th)
            # 合成: Σ(r·M) 或 Σr 或 ΣM(不变)
            SrM = np.sum(r*M); Sr = np.sum(r)
            rel = abs(SrM-S0)/S0*100
            if emod=="e扫描":
                if rel<1.0:
                    P_("    θ_c0=%.1f° e=%.4f: Σ(r·M)=%.5f (偏差 %.3f%%)" % (np.degrees(ph0), e, SrM, rel))
                continue
            P_("    θ_c0=%.1f° %s: Σ(r·M)=%.5f (偏差 %.3f%%)" % (np.degrees(ph0), emod, SrM, rel))

P_("\n" + "="*96)
P_("【C. 需要的偏心 e(若 Σ(r·M) 模型)】")
for ph0 in [0, np.pi/4, np.pi/8]:
    th = np.array([ph0, ph0+2*np.pi/3, ph0+4*np.pi/3])
    # 扫 e 找 Σ(r·M)=S0
    from scipy.optimize import brentq
    def f(e):
        r, M = ellipse_sum(e, th)
        return np.sum(r*M) - S0
    try:
        e_sol = brentq(f, 0, 0.5)
        r, M = ellipse_sum(e_sol, th)
        P_("  θ_c0=%.1f°: 需 e=%.5f (%.4f°)  对照 δ₀=%.4f°  2/9=%.4f°" % (np.degrees(ph0), e_sol, np.degrees(e_sol), np.degrees(d0), np.degrees(two9)))
    except Exception as ex:
        P_("  θ_c0=%.1f°: 无解 (%s)" % (np.degrees(ph0), ex))

out = {"delta_needed": D, "s0": S0}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_proton_ellipse.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_proton_ellipse.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_proton_ellipse.json / .txt")
