#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_proton_L1_dyn.py —— 质子=三色荷共享L1的动力学(第一版): L1记录稳态ω加权合成
======================================================================================
机制: 三个色荷顶点在图上(L2)演化→闭合→L1记录累积; 质子的有效合成用 L1 记录稳态ω
      加权(重数归一化 Z088 §3.1). 若 ω 内生不均匀 → 合成从 ΣM=3 偏移, 看是否 0.35%.
忠实: 用 Z088 配置(path, L4, 六边闭环) 跑稳态, 不调参.
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

from z0_engine import path_adj, torus_adj
from layer_closure import LayerClosure

me, m_p = 0.51099895, 938.272088
theta_e = 2*np.pi/3 + 2/9
A = np.sqrt(me)/(1+np.sqrt(2)*np.cos(theta_e)); A2 = A**2
S0 = m_p/A2      # 2.989468
D  = 3 - S0      # 0.010532
P_("=" * 96)
P_("质子=三色荷共享L1 · 动力学第一版 · k4_proton_L1_dyn.py")
P_("=" * 96)
P_("目标: 有效合成 S' = 3×Σ(ω·M)/Σω 从 3 内生偏移到 %.6f (Δ=%.6f)" % (S0, D))

def M(t): return 1 + np.sqrt(2)*np.cos(t)

# ── 三个色荷角候选 (等间隔120°, 中心角 θ_c) ──
cands = {
    "θ_c=3π/4(绕零点)": 3*np.pi/4,
    "θ_c=π/3(=1/6圈)": np.pi/3,
    "θ_c=2π/3(三代)": 2*np.pi/3,
    "θ_c=π/4": np.pi/4,
}
threes = {}
for nm, tc in cands.items():
    th = np.array([tc, tc+2*np.pi/3, tc+4*np.pi/3])
    threes[nm] = (th, np.array([M(x) for x in th]))

# ── 跑 layer_closure 稳态, 取 L1 记录稳态 ω ──
def run_steady(adj, L, gens, label, retention="last"):
    lc = LayerClosure(adj, L=L, cap=10**5, kac=True, cap_deg=True, retention=retention)
    lc.run(gens, log_every=0)
    # 稳态 ω: 用最后 window 代平均 (L1 记录份额)
    window = min(50, len(lc.w_hist_trace))
    ws = lc.w_hist_trace[-window:]
    wavg = np.mean([np.array(x) for x in ws if x is not None], axis=0)
    qh = lc.q_hist_trace[-1]
    P_("  [%s] q_hist=%.4f  ω_稳态: %s" % (label, qh if qh else float('nan'),
        np.round(wavg,4)))
    return lc, wavg

# path16 L4 (Z088 稳态 5/9, 二分支 2:1)
P_("\n" + "="*96)
P_("【① path16 L4: Z088 六边闭环稳态】")
lc16, w16 = run_steady(path_adj(16), 4, 30000, "path16 L4 30000代(last)")

# 三个色荷顶点: 稳态 ω 最大的三个 (内生选择) + 固定中间三顶点
top_idx = np.argsort(w16)[-3:]
P_("  稳态 ω 最大的 3 顶点: %s (份额 %s)" % (top_idx, np.round(w16[top_idx],4)))

P_("\n" + "="*96)
P_("【② 质子有效合成 S' = 3×Σ(ω·M)/Σω (L1记录加权)】")
for nm, (th, Mv) in threes.items():
    # a) 用 top3 顶点的 ω 加权
    #    M 值按大小与色荷角排序对应? 诚实做法: 三色荷角顺序与顶点顺序(按ω降序)对应
    wtop = w16[top_idx]; wtop = wtop/wtop.sum()
    S_a = 3 * np.sum(wtop * Mv) / np.sum(wtop)
    # b) 均匀加权(对照): 无 L1 共享 → ΣM=3
    S_uni = np.sum(Mv)
    # c) 用 ω 的集中度修正: 若 ω 均匀则=3
    rel_a = abs(S_a - S0)/S0*100
    P_("  %-22s M=%s" % (nm, np.round(Mv,4)))
    P_("    均匀(无L1共享): ΣM=%.4f | L1加权 top3: S'=%.4f (偏差 %.3f%%)" % (S_uni, S_a, rel_a))

# ── 也试: 三个色荷 = 稳态二分支(2:1)的两个强顶点 + 邻居 ──
P_("\n" + "="*96)
P_("【③ 三色荷 = 图上连续三顶点 (质子在中间的局域结构)】")
for center in range(7, 9):  # path16 中间
    idx3 = [center-1, center, center+1]
    for nm, (th, Mv) in threes.items():
        w3 = w16[idx3]; w3 = w3/w3.sum()
        S = 3*np.sum(w3*Mv)/np.sum(w3)
        P_("  中心v=%d %-20s ω=%s → S'=%.4f (偏差 %.3f%%)" %
           (center, nm, np.round(w3,3), S, abs(S-S0)/S0*100))

out = {"target_s0": S0, "top3_idx": top_idx.tolist(), "w_top3": w16[top_idx].tolist(),
       "note": "第一版 L1 记录加权合成机制"}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_proton_L1_dyn.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_proton_L1_dyn.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_proton_L1_dyn.json / .txt")
