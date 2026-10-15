# -*- coding: utf-8 -*-
"""Path E：分支+湮灭（"全分支被实例化"）⟹ 计数分布是否非均匀？

前几轮都假设"荷守恒"（补偿搬运的总量不变）。但 Z001 §3.2 说的是
**每个允许的运动都被实例化**（全分支），重数是**计数**——所以离散度
逐代**增加**，不是守恒。引擎实测绝对值每代爆涨（Z087 §1.2：1.2e98→3.0e107）
正是这个结构的直接后果，不是数值故障。

于是关键问题变了：在"分支创造 + 湮灭删除"的计数动力学下，
  平稳的**占比**分布 ω_v 是否非均匀？是否入 D4 窗？
本文件检验三类零参数的计数规则：
  (R1) 全分支 + 硬容量（每顶点至多 K 个走者）
  (R2) 全分支 + 自饱和（增长率 ∝ (1 − m_v/K)）
  (R3) 全分支 + 线性湮灭（死亡率 ∝ m_v）
"""
from __future__ import annotations
import sys, os
import numpy as np
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "..", "ZeroGPU", "src")))
from z0_engine import path_adj, torus_adj, star_adj

CHECKS=[]
def check(n,c,d=""):
    CHECKS.append((n,bool(c))); print(f"[{'OK' if c else 'FAIL'}] {n}  {d}")

def branch(adj, m0, rule, K=1.0, steps=40000, kappa=0.2):
    """m_v ← m_v·g(m_v) + κ·Σ_{u~v} m_u·g(m_u)/deg(u)   —— 分支(创造)+扩散+规则调整"""
    V=len(adj); deg=np.array([len(a) for a in adj],float); deg[deg==0]=1
    m=np.array(m0,float)
    for t in range(steps):
        if   rule=="R1-cap":    g=np.maximum(0.0, 1.0-m/K)
        elif rule=="R2-sat":    g=np.maximum(0.0, 1.0-m/K)
        elif rule=="R3-death":  g=np.ones(V)-0.02*m
        produced=m*g
        inflow=np.zeros(V)
        for u,nb in enumerate(adj):
            for v in nb: inflow[v]+=produced[u]/deg[u]
        m=produced+kappa*inflow
        if not np.isfinite(m).all() or m.max()>1e12: return m,t,False
        if t>1000 and np.abs(m-m.mean()).max()<1e-13*max(1,m.mean()): return m,t,True
    return m,steps,False

V=16; rng=np.random.RandomState(3)
print("="*92)
print("Path E：分支+湮灭动力学（不守恒；对应'全分支被实例化'）")
print("="*92)
print(f"{'图':10s} {'规则':11s} {'K':>6s} {'末态var/mean²':>13s} {'q(m)':>9s} {'入窗?':>6s} {'稳定':>5s}")
for gname,adj in [("path16",path_adj(16)),("cycle16",torus_adj(1,16)),("star16",star_adj(16))]:
    n=len(adj)
    for rule in ("R1-cap","R2-sat","R3-death"):
        for K in (0.5,1.0,2.0):
            m0=0.3+0.4*rng.rand(n)
            m,t,conv=branch(adj,m0,rule,K=K)
            w=m/m.sum(); q=float((w*w).sum())
            cv=float(np.var(m)/max(1e-30,m.mean()**2))
            print(f"{gname:10s} {rule:11s} {K:6.2f} {cv:13.5f} {q:9.6f} "
                  f"{str(0.5<q<0.6):>6s} {str(conv):>5s}")
    print()
print("d=1 参考（均匀 q=1/16=0.0625；完全集中 q=1）")
print()
print("="*92); print("判据"); print("="*92)
hits=[]
for gname,adj in [("path16",path_adj(16)),("cycle16",torus_adj(1,16)),("star16",star_adj(16))]:
    n=len(adj)
    for rule in ("R1-cap","R2-sat","R3-death"):
        for K in (0.5,1.0,2.0):
            m,_,conv=branch(adj,0.3+0.4*np.random.RandomState(3).rand(n),rule,K=K)
            q=float((((m/m.sum()))**2).sum())
            if 0.5<q<0.6 and conv: hits.append((gname,rule,K,q))
print(f"入 D4 窗且稳定的组合: {len(hits)}")
for h in hits[:10]: print("   ",h)
check("分支-湮灭族存在入窗的非均匀平稳分布", len(hits)>0, f"{len(hits)} 个组合")
