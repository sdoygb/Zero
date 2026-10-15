#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_proton_tristate2.py —— 三份额稳态 4:2:1 确认 + 三色荷合成
======================================================================================
对称三中心图(度3,3,3) 内生出三份额 ω=[4/7,2/7,1/7]. 这是真三色荷稳态.
用三份额加权三色荷 M(θ), 找能否内生合成 3→2.9895 (0.35%).
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
P_("三份额 4:2:1 确认 + 三色荷合成 · k4_proton_tristate2.py")
P_("=" * 96)
P_("目标合成 S0=%.6f (Δ=%.6f)" % (S0, 3-S0))

def tristate_adj(L):
    leaves = sum(L); V = 3+leaves
    adj = [[] for _ in range(V)]
    for a in (0,1,2):
        for b in (0,1,2):
            if a<b: adj[a].append(b); adj[b].append(a)
    nxt=3
    for i,nl in enumerate(L):
        for _ in range(nl): adj[i].append(nxt); adj[nxt].append(i); nxt+=1
    return [sorted(a) for a in adj], [0,1,2]

def M(t): return 1 + np.sqrt(2)*np.cos(t)

# ── ① 确认对称图稳态 4:2:1 (度3,3,3, 多代/多seed) ──
P_("\n【① 对称三中心图(度3,3,3) 稳态三份额】")
for seed_kind in ["geom", "hub3", "uniform"]:
    adj, hubs = tristate_adj([1,1,1]); V=len(adj)
    if seed_kind=="geom": seed=[float(2**(V-1-v)) for v in range(V)]
    elif seed_kind=="hub3":
        seed=[1.0]*V
        for h in hubs: seed[h]=100.0
    else: seed=[1.0]*V
    vec = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
    for _ in range(60000): vec.run_generation()
    rec = vec.record/vec.record.sum()
    w = rec[hubs]; w=w/w.sum()
    if w[2] > 1e-12:
        P_("  %-8s: q_hist=%.4f  三份额ω=%s (比 %d:%d:%d)" %
           (seed_kind, vec.q_hist(), np.round(w,4),
            round(w[0]/w[2]), round(w[1]/w[2]), 1))
    else:
        P_("  %-8s: q_hist=%.4f  ω=%s (三份额退化)" % (seed_kind, vec.q_hist(), np.round(w,4)))

# ── ② 用 4:2:1 三份额加权三色荷 M ──
P_("\n" + "="*96)
P_("【② 三色荷合成 S'=3×Σ(ω·M), ω=4:2:1 内生】")
w = np.array([4,2,1])/7.0
P_("  内生三份额 ω=%s (4:2:1)" % np.round(w,4))
best=[]
for tc, tcnm in [(3*np.pi/4,"绕零点"), (np.pi/3,"1/6圈"), (np.pi/4,"π/4"),
                 (2*np.pi/3,"三代"), (5*np.pi/12,"5/12"), (np.pi/6,"1/12圈"),
                 (7*np.pi/12,"7/12"), (np.pi/2,"1/4圈")]:
    th = np.array([tc, tc+2*np.pi/3, tc+4*np.pi/3])
    Mv = np.array([M(x) for x in th])
    # ω 与 M 的各种对应(排序组合 3!=6)
    from itertools import permutations
    for perm in permutations([0,1,2]):
        S = 3*np.sum(np.array(w)[list(perm)]*Mv)
        rel = abs(S-S0)/S0*100
        best.append((rel, tcnm, perm, Mv, S))
best.sort(key=lambda x:x[0])
P_("  Top6 最接近 S0=%.4f:" % S0)
for rel, tcnm, perm, Mv, S in best[:6]:
    P_("    θ_c=%-8s M%s ω序%s → S'=%.4f (偏差 %.3f%%)" %
       (tcnm, np.round(Mv,3), np.round(w[list(perm)],3), S, rel))

# ── ③ 其他合成: Σ(ω M²), Σω·|M|, (Σ√(ωM))² 等 ──
P_("\n" + "="*96)
P_("【③ 其他三份额合成形式(诚实探测)】")
adj, hubs = tristate_adj([1,1,1]); V=len(adj)
seed=[float(2**(V-1-v)) for v in range(V)]
vec = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
for _ in range(60000): vec.run_generation()
rec = vec.record/vec.record.sum()
w = rec[hubs]; w=w/w.sum()
for tc, tcnm in [(3*np.pi/4,"绕零点"), (np.pi/3,"1/6圈")]:
    th = np.array([tc, tc+2*np.pi/3, tc+4*np.pi/3]); Mv=np.array([M(x) for x in th])
    Ms = np.sort(Mv)[::-1]; ws = np.sort(w)[::-1]
    S1 = 3*np.sum(ws*Ms)                       # 加权和
    S2 = 3*np.sum(ws*Ms**2)/3 if np.sum(ws)>0 else 0
    S3 = (np.sqrt(np.sum(ws*Ms**2))*3)
    P_("  θ_c=%s: Σ(ωM)=%.4f | (1/3)Σ(ωM²)=%.4f | √Σ(ωM²)=%.4f" % (tcnm, S1, np.sum(ws*Ms**2)/3, np.sqrt(np.sum(ws*Ms**2))))

out = {"target_s0": S0, "w_421": [float(x) for x in w],
       "best": (best[0][0], best[0][1], [int(x) for x in best[0][2]], [float(x) for x in best[0][3]], best[0][4]) if best else None}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_proton_tristate2.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_proton_tristate2.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_proton_tristate2.json / .txt")
