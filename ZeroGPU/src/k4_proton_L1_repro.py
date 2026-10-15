#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_proton_L1_repro.py —— 复现 Z088 稳态(geom seed, last, path16), 再用 L1 记录 ω 试三色荷合成
======================================================================================
诚实目标: (1) 确认 Z088 5/9 二分支稳态可复现; (2) 看二分支 ω 能否忠实映射到质子三色荷,
        而非自由拟合.
"""
from __future__ import annotations
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

from z0_engine import path_adj
from z0_vectorized import ZeroCoreVec

me, m_p = 0.51099895, 938.272088
theta_e = 2*np.pi/3 + 2/9
A = np.sqrt(me)/(1+np.sqrt(2)*np.cos(theta_e)); A2 = A**2
S0 = m_p/A2
P_("=" * 96)
P_("复现 Z088 稳态 + 三色荷映射 · k4_proton_L1_repro.py")
P_("=" * 96)
P_("目标合成 S0=%.6f (Δ=%.6f)" % (S0, 3-S0))

# ── (1) 复现 Z088: path16 geom last L4 ──
V = 16
seed = [float(2**(V-1-v)) for v in range(V)]
vec = ZeroCoreVec(path_adj(V), seed=seed, retention="last", L=4)
GENS = 60000
for g in range(GENS):
    vec.run_generation()
rec = vec.record / vec.record.sum()
qh = vec.q_hist(); ql = vec.q_live()
P_("\n【① 复现 Z088 稳态】")
P_("  q_hist=%.4f  q_live=%.4f  (Z088: 5/9=%.4f)" % (qh, ql, 5/9))
P_("  record 比例 ω(前8): %s" % np.round(rec[:8],4))
P_("  record 比例 ω(后8): %s" % np.round(rec[8:],4))
top = np.argsort(rec)[-3:][::-1]
P_("  ω 最大 3 顶点: %s 份额 %s" % (top, np.round(rec[top],4)))

# ── (2) 三色荷映射尝试 ──
P_("\n【② 二分支 ω → 质子三色荷合成】(诚实评估)")
def M(t): return 1 + np.sqrt(2)*np.cos(t)
# Z088 稳态: 两强顶点(top2) ω=(2/3,1/3), 第三色荷用邻居(ω小)
t2 = top[:2]
for neighbor in [t2[0]-1, t2[1]+1]:
    if 0 <= neighbor < V and neighbor not in t2:
        idx3 = [t2[0], t2[1], neighbor]
        w3 = rec[idx3]; w3n = w3/w3.sum()
        P_("  三色荷顶点 %s ω=%s(归一化)" % (idx3, np.round(w3n,4)))
        # 色荷角: θ_c 扫描, M 按顶点顺序
        for tc in [3*np.pi/4, np.pi/3, 2*np.pi/3, np.pi/4, np.pi/6, 5*np.pi/12, 7*np.pi/12]:
            th = np.array([tc, tc+2*np.pi/3, tc+4*np.pi/3])
            Mv = np.array([M(x) for x in th])
            # 分配: M 值排序匹配 ω 排序(强顶点拿大 M)
            Msort = np.sort(Mv)[::-1]   # 大到小
            wsort = np.sort(w3n)[::-1]  # 大到小
            S = 3*np.sum(wsort*Msort)/np.sum(wsort)
            if abs(S-S0)/S0 < 0.05:
                P_("    θ_c=%.1f° M降序%s ω降序%s → S'=%.4f (偏差 %.3f%%)" %
                   (np.degrees(tc), np.round(Msort,3), np.round(wsort,3), S, abs(S-S0)/S0*100))
        P_("")

out = {"q_hist": qh, "q_live": ql, "top3": top.tolist(), "rec_top3": rec[top].tolist()}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_proton_L1_repro.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_proton_L1_repro.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_proton_L1_repro.json / .txt")
