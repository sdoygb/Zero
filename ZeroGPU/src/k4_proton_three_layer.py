#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_proton_three_layer.py —— 作用原理建模: L0-L1-L2 三层互动, 找质子的质量读出
======================================================================================
机制: 三色荷 = 对称三中心图(三hub三角形, 度3,3,3). 跑六边闭环稳态.
取三层的三色荷份额: L1 record(累积) / L2 gen_closures(本代闭合) / seed(L2演化态).
"L1→L2 的变化" = 这些量的差异. 质量读出 = 用哪层加权三色荷 M 给 0.35%?
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
P_("三层互动作用原理 · 质子质量读出 · k4_proton_three_layer.py")
P_("=" * 96)
P_("目标合成 S0=%.6f (Δ=%.6f)" % (S0, 3-S0))

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
seed = [float(2**(V-1-v)) for v in range(V)]
P_("图: 三hub三角形+各1叶, V=%d, seed=geom" % V)

# ── 跑稳态, 采样各层三色荷份额 ──
vec = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
GENS=60000
traces = {"rec":[], "gc":[], "seed":[]}
for g in range(GENS):
    vec.run_generation()
    if (g+1)%2000==0:
        rr = vec.record/vec.record.sum(); traces["rec"].append(rr[hubs]/rr[hubs].sum())
        gg = vec.gen_closures/vec.gen_closures.sum(); traces["gc"].append(gg[hubs]/gg[hubs].sum())
        ss = vec.seed/vec.seed.sum(); traces["seed"].append(ss[hubs]/ss[hubs].sum())

P_("\n【① 三层三色荷份额随代演化(每2000代)】")
P_("  gen    L1 record(4:2:1?)      L2 gen_closures         L2 seed")
for i in range(0, len(traces["rec"]), max(1,len(traces["rec"])//6)):
    g = (i+1)*2000
    P_("  %6d  [%.3f %.3f %.3f]   [%.3f %.3f %.3f]   [%.3f %.3f %.3f]" %
       (g, *traces["rec"][i], *traces["gc"][i], *traces["seed"][i]))

# 末态
rec_last = traces["rec"][-1]; gc_last = traces["gc"][-1]; seed_last = traces["seed"][-1]
P_("\n【② 稳态末态三份额(归一化)】")
P_("  L1 record  : %s  Σω²=%.4f" % (np.round(rec_last,4), np.sum(rec_last**2)))
P_("  L2 gen_cl  : %s  Σω²=%.4f" % (np.round(gc_last,4), np.sum(gc_last**2)))
P_("  L2 seed    : %s  Σω²=%.4f" % (np.round(seed_last,4), np.sum(seed_last**2)))
P_("  L1→L2 变化 (record→gen_closures): %s" % np.round(gc_last-rec_last,4))

# ── ③ 用各层量加权三色荷 M ──
P_("\n" + "="*96)
P_("【③ 质量读出候选: 各层份额加权三色荷 M → 找 2.9895】")
layers = {"L1 record": rec_last, "L2 gen_closures": gc_last, "L2 seed": seed_last}
from itertools import permutations
best_all = []
for lname, w in layers.items():
    w = np.array(w)
    for tc, tcnm in [(3*np.pi/4,"绕零点"), (np.pi/3,"1/6圈"), (np.pi/4,"π/4"),
                     (2*np.pi/3,"三代"), (5*np.pi/12,"5/12"), (np.pi/6,"1/12")]:
        th = np.array([tc, tc+2*np.pi/3, tc+4*np.pi/3]); Mv=np.array([M(x) for x in th])
        for perm in permutations([0,1,2]):
            S = 3*np.sum(np.array(w)[list(perm)]*Mv)
            best_all.append((abs(S-S0)/S0*100, lname, tcnm, perm, S))
best_all.sort(key=lambda x:x[0])
P_("  Top8 全局最接近 S0:")
for rel, lname, tcnm, perm, S in best_all[:8]:
    P_("    %-14s θ=%-7s ω序%s → S'=%.5f (偏差 %.3f%%)" % (lname, tcnm, [int(p) for p in perm], S, rel))

# ④ 若 L1→L2 变化量 = 0.35%?
P_("\n" + "="*96)
P_("【④ L1→L2 的变化量(record→seed 相对差)】")
rel_change = np.max(np.abs(seed_last/rec_last - 1))
P_("  每三色荷 record→seed 比: %s" % np.round(seed_last/rec_last,4))
P_("  最大相对变化 = %.4f (0.35%%?)" % rel_change)

out = {"s0": S0, "rec": rec_last.tolist(), "gc": gc_last.tolist(), "seed": seed_last.tolist()}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_proton_three_layer.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_proton_three_layer.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_proton_three_layer.json / .txt")
