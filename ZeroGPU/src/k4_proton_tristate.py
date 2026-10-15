#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_proton_tristate.py —— 三顶点稳态: 质子三色荷 = 三 hub 三角形图, L1 记录承载三份额
======================================================================================
设计: 三个 hub(色荷) 连成三角形(S3 对称), 各带叶子制造度不均(Z0: 结构唯一来源=度不均).
跑稳态, 取三 hub 的 L1 record 份额 ω, 看是否三分(非二分支), 再作三色荷合成.
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
P_("三顶点稳态 · 三hub三角形图 · k4_proton_tristate.py")
P_("=" * 96)
P_("目标合成 S0=%.6f (Δ=%.6f, 0.35%%)" % (S0, 3-S0))

def tristate_adj(L):
    """三 hub(0,1,2) 连三角形, hub i 带 L[i] 叶子. 返回 (adj, hub_idx)."""
    leaves = sum(L)
    V = 3 + leaves
    adj = [[] for _ in range(V)]
    for a in (0,1,2):
        for b in (0,1,2):
            if a < b: adj[a].append(b); adj[b].append(a)
    nxt = 3
    for i, nl in enumerate(L):
        for _ in range(nl):
            adj[i].append(nxt); adj[nxt].append(i); nxt += 1
    return [sorted(a) for a in adj], [0,1,2]

def run_tristate(L, gens, seed_kind, label):
    adj, hubs = tristate_adj(L)
    V = len(adj)
    if seed_kind == "geom":
        seed = [float(2**(V-1-v)) for v in range(V)]
    elif seed_kind == "hub3":
        seed = [1.0]*V
        for h in hubs: seed[h] = 100.0
    else:
        seed = [1.0]*V
    vec = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
    for _ in range(gens): vec.run_generation()
    rec = vec.record/vec.record.sum()
    w_hub = rec[hubs]
    qh = vec.q_hist()
    P_("\n【%s】V=%d q_hist=%.4f" % (label, V, qh))
    P_("  三hub L1记录份额 ω=%s (归一化 %s)" % (np.round(w_hub,4), np.round(w_hub/w_hub.sum(),4)))
    P_("  Σω²(三hub集中)=%.4f  度=%s" % (np.sum(w_hub**2), [len(adj[h]) for h in hubs]))
    return w_hub/w_hub.sum()

# 度不均的三中心图
for L, nm in [([2,1,0],"叶2,1,0 度4,3,2"), ([3,1,0],"叶3,1,0 度5,3,2"),
              ([2,2,0],"叶2,2,0 度4,4,2"), ([3,2,1],"叶3,2,1 度5,4,3"),
              ([1,1,1],"叶1,1,1 度3,3,3(对称)")]:
    try:
        run_tristate(L, 30000, "geom", nm)
    except Exception as e:
        P_("  %s 失败: %s" % (nm, e))

# 三色荷合成: 用给定三份额 ω, 加权三色荷 M(θ_c 绕零点等)
P_("\n" + "="*96)
P_("【三色荷合成 S'=3×Σ(ω·M)/Σω, 三份额内生于稳态】")
def M(t): return 1 + np.sqrt(2)*np.cos(t)

# 用第一个度不均图的稳态 ω
adj0, hubs0 = tristate_adj([2,1,0])
V0 = len(adj0)
seed0 = [float(2**(V0-1-v)) for v in range(V0)]
vec0 = ZeroCoreVec(adj0, seed=seed0, retention="last", L=4)
for _ in range(30000): vec0.run_generation()
rec0 = vec0.record/vec0.record.sum()
w0 = rec0[hubs0]; w0 = w0/w0.sum()
P_("  稳态三份额 ω=%s (叶2,1,0 图)" % np.round(w0,4))
for tc, tcnm in [(3*np.pi/4,"绕零点"), (np.pi/3,"1/6圈"), (np.pi/4,"π/4"), (5*np.pi/12,"5/12圈")]:
    th = np.array([tc, tc+2*np.pi/3, tc+4*np.pi/3])
    Mv = np.array([M(x) for x in th])
    for ord_desc in [True, False]:
        Ms = np.sort(Mv)[::-1] if ord_desc else np.sort(Mv)
        S = 3*np.sum(w0*Ms)/np.sum(w0)
        P_("    θ_c=%-7s M%s → S'=%.4f (偏差 %.3f%%)" %
           (tcnm, np.round(Ms,3), S, abs(S-S0)/S0*100))
    P_("")

out = {"target_s0": S0, "note": "三hub三角形图三份额稳态"}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_proton_tristate.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_proton_tristate.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_proton_tristate.json / .txt")
