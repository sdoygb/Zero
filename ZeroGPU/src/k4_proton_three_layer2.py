#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_proton_three_layer2.py —— 作用原理: 非last(retention) 下 L1→L2 是否出现变化
======================================================================================
last 下三层完全自洽(4:2:1), 无 L1→L2 变化. 试 top2 / all: 历史累积作种子,
看 L1 record vs L2 gen_closures/seed 是否出现差异, 差异能否给 0.35%.
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

from z0_vectorized import ZeroCoreVec
me, m_p = 0.51099895, 938.272088
theta_e = 2*np.pi/3 + 2/9
A = np.sqrt(me)/(1+np.sqrt(2)*np.cos(theta_e)); A2 = A**2
S0 = m_p/A2
P_("=" * 96)
P_("非last保留 → L1→L2 变化 · k4_proton_three_layer2.py")
P_("=" * 96)

def tristate_adj(L):
    leaves=sum(L); V=3+leaves
    adj=[[] for _ in range(V)]
    for a in (0,1,2):
        for b in (0,1,2):
            if a<b: adj[a].append(b); adj[b].append(a)
    nxt=3
    for i,nl in enumerate(L):
        for _ in range(nl): adj[i].append(nxt); adj[nxt].append(i); nxt+=1
    return [sorted(a) for a in adj], [0,1,2]
def M(t): return 1 + np.sqrt(2)*np.cos(t)

adj, hubs = tristate_adj([1,1,1]); V=len(adj)
seed0 = [float(2**(V-1-v)) for v in range(V)]

for retention in ["top2", "all", "last"]:
    vec = ZeroCoreVec(adj, seed=seed0, retention=retention, L=4)
    GENS = 100000
    rec_trace=[]; gc_trace=[]
    for g in range(GENS):
        vec.run_generation()
        if (g+1)%5000==0:
            rr=vec.record/vec.record.sum(); rec_trace.append(rr[hubs]/rr[hubs].sum())
            gg=vec.gen_closures/vec.gen_closures.sum(); gc_trace.append(gg[hubs]/gg[hubs].sum())
    rec=rec_trace[-1]; gc=gc_trace[-1]
    P_("\n【%s】" % retention)
    P_("  L1 record  : %s" % np.round(rec,4))
    P_("  L2 gen_cl  : %s" % np.round(gc,4))
    diff = gc-rec
    P_("  L1→L2 变化 : %s  最大|Δ|=%.4f" % (np.round(diff,4), np.max(np.abs(diff))))
    P_("  变化相对 record: %s" % np.round(np.abs(diff)/np.maximum(rec,1e-9),4))
    # 三色荷合成: 用 L2 gen_closures(本代, 质量读出) vs L1 record
    from itertools import permutations
    for lname, w in [("L2 gen_closures", gc), ("L1 record", rec)]:
        w=np.array(w)
        for tc, tcnm in [(3*np.pi/4,"绕零点"), (np.pi/3,"1/6圈"), (np.pi/4,"π/4"), (2*np.pi/3,"三代")]:
            th=np.array([tc,tc+2*np.pi/3,tc+4*np.pi/3]); Mv=np.array([M(x) for x in th])
            for perm in permutations([0,1,2]):
                S=3*np.sum(np.array(w)[list(perm)]*Mv)
                if abs(S-S0)/S0 < 0.03:
                    P_("    %-13s θ=%-6s ω序%s → S'=%.5f (偏差 %.3f%%)" %
                       (lname, tcnm, [int(p) for p in perm], S, abs(S-S0)/S0*100))

out={"s0":S0}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES,"k4_proton_three_layer2.json"),"w",encoding="utf-8") as f_:
    json.dump(out,f_,ensure_ascii=False,indent=1)
with open(os.path.join(RES,"k4_proton_three_layer2.txt"),"w",encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_proton_three_layer2.json / .txt")
