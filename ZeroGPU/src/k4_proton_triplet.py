#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_proton_triplet.py —— 质子=三个色荷点, 用三点偏差合成质量 (用户的猜想)
======================================================================================
猜想: 质子是圆上有三个点(色荷), 三个点跟零点一起偏差, 偏差有大有小, 用三点偏差合成.
模型: 三个色荷等间隔 120° (色闭合), 中心角 θ_c 扫描; 每个点偏离零点 3π/4.
质量机制(Z121): 偏离振幅 M(θ)=1+√2cosθ, m∝M². 标度 A 从电子定.
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
Me = 1 + np.sqrt(2)*np.cos(theta_e)
A = np.sqrt(me)/Me
A2 = A**2
P_("=" * 96)
P_("质子=三色荷点 · k4_proton_triplet.py (三点偏差合成)")
P_("=" * 96)
P_("标度 A=%.6f  A²=%.6f (从电子 θ_e=%.3f° 定)" % (A, A2, np.degrees(theta_e)))
P_("零点 = 3π/4 (135°, 光子站那); 目标 m_p=%.3f MeV" % m_p)

z0 = 3*np.pi/4
# 三个色荷等间隔 120°, 中心角 θ_c
def triplet(θc):
    return np.array([θc, θc+2*np.pi/3, θc+4*np.pi/3])
def M_of(t):
    return 1 + np.sqrt(2)*np.cos(t)

P_("\n" + "=" * 96)
P_("【A. 等间隔三色荷 + 各种合成 (中心角 θ_c 扫描)】")
theta_grid = np.linspace(0, 2*np.pi, 7201)
best_hits = []
for θc in theta_grid:
    th = triplet(θc)
    M = M_of(th)
    D = th - z0                      # 相对零点偏差角
    synth = {
        "ΣM²        ": (M**2).sum(),      # 恒 6
        "ΣM         ": M.sum(),           # 恒 3
        "ΣM²/2      ": (M**2).sum()/2,
        "ΣM²/3      ": (M**2).sum()/3,
        "Σ|M|       ": np.abs(M).sum(),
        "ΣD²        ": (D**2).sum(),
        "Σ|D|       ": np.abs(D).sum(),
        "RMS(M²)    ": np.sqrt((M**2).mean()),
        "Π|M|^(1/3) ": np.prod(np.abs(M))**(1/3),
    }
    for nm, s in synth.items():
        # m = A²·S  (S 视作平方偏离和的量纲)
        m_pred = A2 * s
        rel = abs(m_pred - m_p)/m_p*100
        if rel < 2.0:
            best_hits.append((rel, nm, np.degrees(θc), s, m_pred))

best_hits.sort(key=lambda x: x[0])
P_("  偏差<2%% 的命中 (前12):")
for rel, nm, θcd, s, mp in best_hits[:12]:
    P_("    %-10s θ_c=%.2f°  S=%.4f  m=%.1f MeV (偏差 %.2f%%)" % (nm, θcd, s, mp, rel))

# 恒定量验证
P_("\n  恒定量: ΣM²=%.4f (等间隔三色荷), ΣM=%.4f" % ((M_of(triplet(0))**2).sum(), M_of(triplet(0)).sum()))
P_("  ⟹ m=A²·ΣM²=%.0f,  m=A²·ΣM=%.1f,  m=A²·ΣM²/2=%.1f  (均与 θ_c 无关)" % (A2*6, A2*3, A2*3))

# ============ B. 三种"三点合成"的物理叙述 ============
P_("\n" + "=" * 96)
P_("【B. 三候选合成 vs 观测】")
for nm, s, label in [("ΣM²(三平方和)", 6.0, "三个点偏差平方和"),
                     ("ΣM (三和)", 3.0, "三个点偏差之和(=3)"),
                     ("ΣM²/2 (=3)", 3.0, "三平方和/2, 数值同ΣM")]:
    m_pred = A2*s
    P_("  %-16s S=%g → m=%.1f MeV (观测 938.3, 偏差 %.2f%%)  %s" % (nm, s, m_pred, abs(m_pred-m_p)/m_p*100, label))

# ============ C. 找让 m_p 精确的合成 S0 ============
P_("\n" + "=" * 96)
P_("【C. 观测反推: m_p 需要合成值 S0 = m_p/A²】")
S0 = m_p/A2
P_("  S0 = 938.272/%.2f = %.6f" % (A2, S0))
P_("  若合成=ΣM²/2=3: S0/S=%.6f (需 1)" % (S0/3))
P_("  若合成=ΣM=3:     S0/3=%.6f" % (S0/3))
P_("  S0=%.6f → 与 3 偏差 %.3f%%" % (S0, abs(S0-3)/3*100))

out = {"A": A, "S0_needed": S0, "sum_M2_over2": 3.0, "m_from_S0": A2*3.0,
       "sumM2_const": 6.0, "hits": [{"synth":n,"theta_c":t,"s":s,"m":m,"rel":r} for r,n,t,s,m in best_hits[:8]]}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_proton_triplet.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_proton_triplet.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_proton_triplet.json / .txt")
