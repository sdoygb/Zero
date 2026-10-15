#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_proton_correction.py —— 质子质量 0.35% 修正的来源搜索 (闭合到 m_p=938.272)
======================================================================================
现状: 等间隔三色荷 ΣM=3 (恒定) → m=A²·3=941.6 (偏差 0.35%). 需合成 S0=2.989468.
目标: 从 Zero 内生结构逼出 3−S0=0.010532 的修正, 让质子质量无参数闭合.
方向: ①0.0105 的 Zero 内生组合  ②非等间隔色荷(真实结构微扰)  ③合成加修正项.
"""
from __future__ import annotations
import os, json, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

me, m_p = 0.51099895, 938.272088
theta_e = 2*np.pi/3 + 2/9
A = np.sqrt(me)/(1 + np.sqrt(2)*np.cos(theta_e)); A2 = A**2
S0 = m_p/A2
P_("=" * 96)
P_("质子质量修正来源 · k4_proton_correction.py")
P_("=" * 96)
P_("A=%.6f  A²=%.6f  目标合成 S0=m_p/A²=%.6f  需修正 Δ=3−S0=%.6f" % (A, A2, S0, 3-S0))

two9 = 2/9; d0 = np.pi/12 - two9; pi= np.pi
tau = 2*pi
def M(t): return 1 + np.sqrt(2)*np.cos(t)

# ============ ① 0.0105 的 Zero 内生组合 ============
P_("\n" + "="*96)
P_("【① Δ=0.010532 的 Zero 内生量组合】")
D = 3 - S0
base = {"2/9":two9, "δ₀":d0, "π/12":np.pi/12, "π/8":np.pi/8, "π/4":np.pi/4,
        "π/3":np.pi/3, "π/2":np.pi/2, "π":pi, "2π/3":2*pi/3, "π/6":pi/6,
        "(2/3)(1/3)":two9, "1/3":1/3, "1/6":1/6, "1/12":1/12, "√2-1":np.sqrt(2)-1}
hits = []
for nm, v in base.items():
    rel = abs(v-D)/D*100
    hits.append((rel, nm, v))
for k in [2,3,4,5,6,8,9,12,16,24,25]:
    hits.append((abs(1/k-D)/D*100, f"1/{k}", 1/k))
    hits.append((abs(d0/k-D)/D*100, f"δ₀/{k}", d0/k))
    hits.append((abs(two9/k-D)/D*100, f"2/9/{k}", two9/k))
# 组合 a·δ₀ ± b·2/9, a·2/9·δ₀, 等
for a in range(1,6):
    for b in range(1,6):
        hits.append((abs(a*d0-b*two9-D)/D*100 if a*d0-b*two9>0 else 1e3, f"{a}δ₀−{b}·2/9", a*d0-b*two9))
        hits.append((abs(a*two9-d0-D)/D*100 if a*two9-d0>0 else 1e3, f"{a}·2/9−δ₀", a*two9-d0))
hits.append((abs(two9*d0-D)/D*100, "2/9·δ₀", two9*d0))
hits.sort(key=lambda x:x[0])
P_("  Top10 最接近 Δ=0.010532:")
for rel, nm, v in hits[:10]:
    P_("    %-14s = %.7f  (偏差 %.3f%%)" % (nm, v, rel))

# ============ ② 非等间隔色荷微扰 ============
P_("\n" + "="*96)
P_("【② 色荷间隔微扰: 三色荷非严格 120°, 微扰 ε 取内生量】")
# 色闭合微扰: 间隔 2π/3 + ε_i. 找 ε 使 ΣM=S0.
def sumM(eps1, eps2, θc):
    th = np.array([θc, θc+2*pi/3+eps1, θc+4*pi/3+eps2])
    return np.sum(M(th))
# 扫描 θc 和 ε, 找 ΣM=2.989468
best = None
for θc in np.linspace(0, 2*pi, 3601):
    for e1 in np.linspace(-0.05, 0.05, 401):
        e2 = -e1  # 对称微扰(保持中心)
        sm = sumM(e1, e2, θc)
        rel = abs(sm-S0)/S0*100
        if best is None or rel < best[0]:
            best = (rel, θc, e1)
P_("  非对称微扰(e2=-e1)最优: θ_c=%.2f° ε1=%.4f rad=%.3f°  ΣM=%.5f (偏差 %.3f%%)"
  % (np.degrees(best[1]), best[2], np.degrees(best[2]), sumM(best[2],-best[2],best[1]), best[0]))
P_("  所需 ε1 量级: %.4f rad = %.3f°  (对照 δ₀=%.3f°, 2/9=%.3f°)" % (best[2], np.degrees(best[2]), np.degrees(d0), np.degrees(two9)))

# ③ 合成加修正项
P_("\n" + "="*96)
P_("【③ 合成 = ΣM − c·(内生量): 找自然 c】")
for nm, v in [("δ₀",d0), ("2/9",two9), ("(2/9)²",two9**2), ("δ₀²",d0**2), ("(2/9)(1/3)",two9*(1/3))]:
    c = D/v
    P_("    c = Δ/%-5s = %.3f   (自然? %s)" % (nm, c, "整数/半" if abs(c-round(c))<0.01 or abs(2*c-round(2*c))<0.01 else "非自然"))

out = {"delta_needed": D, "s0": S0, "best_simple": hits[0], "best_micro": best,
       "note": "ΣM=3 恒定(等间隔), 修正需非等间隔微扰或合成修正项"}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_proton_correction.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_proton_correction.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_proton_correction.json / .txt")
