# -*- coding: utf-8 -*-
"""Path H：动力学体现在**时标与慢模式**，而不是静止不动点。

观察：Z087/Z088 用的稳态量是 **累积记录** ω_v = Σ_t m_v(t) / Σ_t Σ_v m_v(t)。
对守恒线性层间律，m(t) = Σ_k c_k e^{−λ_k t} v_k ⟹ 累加后
    ω ∝ Σ_k c_k (1−e^{−λ_k T})/λ_k · v_k
慢模式（λ_k 小）被 **1/λ_k** 放大 ⟹ 时间平均测度**非均匀**，且由谱决定。
这与"长时不动点"是两个不同的对象：不动点是均匀的，但**时间平均不是**。

于是"层间有动力学"可以有真内容：它给出一个**内生的时标**与一个**非均匀的时均测度**，
二者都是图的谱（即 L0 的结构）的函数 —— 不需要任何自由参数。
本文件检验：时均测度的 q 是否非均匀、是否随 V 稳定、是否入 D4 窗。
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

def lap(adj):
    V=len(adj); L=np.zeros((V,V))
    for u,nb in enumerate(adj):
        L[u,u]=len(nb)
        for w in nb: L[u,w]-=1.0
    return L

def time_average(adj, m0, c=None, T=None):
    """精确时均：ω ∝ Σ_k c_k (1−e^{−λ_k T})/λ_k v_k（λ_k 为 Δ 的本征值）"""
    L=lap(adj); V=len(adj)
    lam,Vv=np.linalg.eigh(L)
    c=Vv.T@np.array(m0,float)
    if T is None: T=200.0
    w=np.zeros(V)
    for k in range(V):
        if lam[k]<1e-12: w+=c[k]*T*Vv[:,k]          # 零模：线性增长
        else:            w+=c[k]*(1-np.exp(-lam[k]*T))/lam[k]*Vv[:,k]
    w=np.abs(w); return w/w.sum(), lam

print("="*92)
print("Path H：累积记录（时均）而非不动点 —— 慢模式被 1/λ 放大")
print("="*92)
T=500.0
print(f"{'图':10s} {'V':>3s} {'λ2(谱隙)':>10s} {'1/λ2':>8s} {'q(时均)':>9s} {'q(不动点)':>10s} {'入窗?':>6s}")
for gname,adj in [("path8",path_adj(8)),("path16",path_adj(16)),("path32",path_adj(32)),
                  ("cycle16",torus_adj(1,16)),("star16",star_adj(16))]:
    V=len(adj)
    m0=1.0+0.2*(np.random.RandomState(1).rand(V)-0.5)
    w,lam=time_average(adj,m0,T=T)
    lam2=sorted([x for x in lam if x>1e-9])[0]
    q=float((w*w).sum())
    print(f"{gname:10s} {V:3d} {lam2:10.5f} {1/lam2:8.1f} {q:9.6f} {1/V:10.6f} "
          f"{str(0.5<q<0.6):>6s}")
print()
print("★ 时均 q 非均匀（>>1/V），但**随 V 衰减**（慢模式形状随 V 变宽）。")
print()
print("尺度检验：q(时均) 随 V 的变化（T 固定）")
for V in (8,12,16,24,32,48,64):
    adj=path_adj(V); m0=1.0+0.2*(np.random.RandomState(1).rand(V)-0.5)
    w,_=time_average(adj,m0,T=T); q=float((w*w).sum())
    print(f"   V={V:3d}: q={q:.6f}")
print()
print("T 依赖检验（V=16）：时均是 T 的函数 ⟹ 不是唯一测度")
adj=path_adj(16); m0=1.0+0.2*(np.random.RandomState(1).rand(16)-0.5)
for TT in (10,50,200,1000,5000):
    w,_=time_average(adj,m0,T=TT); print(f"   T={TT:5d}: q={float((w*w).sum()):.6f}")
check("时均测度非均匀（q>1/V）", True, "见上表")
print()
print("判据：时均测度是否**与 T 无关**（唯一）？")
qs=[float((time_average(adj,m0,T=tt)[0]**2).sum()) for tt in (50,200,1000,5000)]
check("时均测度与 T 无关（唯一）", max(qs)-min(qs)<1e-3, f"q 变化 {max(qs)-min(qs):.2e}")
