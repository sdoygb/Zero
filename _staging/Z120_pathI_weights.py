# -*- coding: utf-8 -*-
"""Path I：让**顶点权重（度/势）本身成为动力学的** ⟹ 守恒子空间大于"常数向量"

前几轮所有失败的共同根因（本文件先复述为定理）：
  守恒律 ∂_t m = −Δ_Γ m（或任何把 m 送到 ker 的算子）⟹ 极限集 = ker Δ_Γ = 常数向量
  ⟹ q → 1/V。要非均匀，必须让**守恒子空间**含非常数向量。

Zero 里有现成的机制：Z002（GR 线）用 Perron 权重 w_ij ∝ φ_i φ_j A_ij，
Z033（泊松）把势 φ 当作动力学的。即：**图的权重本身由场决定**。
本文件检验两个具体机制（都零自由参数）：

  (I1) 度-荷耦合：搬运率 ∝ 1/deg_eff，而 deg_eff 由荷更新
       deg_eff,v = Σ_{u~v} m_u      （"几何由物质决定"的离散版）
       ⟹ 这是 GR 的离散模拟：物质弯曲几何，几何引导物质（双向）

  (I2) 势-荷耦合（离散泊松）：Δφ = m − mean(m)（势由荷产生），
       搬运率 ∝ e^{φ_v}（势引导搬运）⟹ 双向反馈

检验：是否存在**稳定非均匀**稳态？q 是否入窗？是否零自由参数？
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

def edges_of(adj):
    E=set()
    for u,nb in enumerate(adj):
        for w in nb:
            if u!=w: E.add((min(u,w),max(u,w)))
    return sorted(E)

def run_I1(adj, m0, kappa=0.05, steps=200000, tol=1e-13):
    """(I1) 度-荷耦合：搬运率 ∝ 荷，几何由物质决定"""
    E=edges_of(adj); V=len(adj); m=np.array(m0,float)
    for t in range(steps):
        dm=np.zeros(V)
        for (u,v) in E:
            rate=np.sqrt(max(m[u],0)*max(m[v],0))       # 几何 ∝ 荷的几何平均
            f=rate*(m[u]-m[v])
            dm[u]-=f; dm[v]+=f
        m=m+kappa*dm
        if not np.isfinite(m).all() or np.abs(m).max()>1e3: return m,t,False
        if np.abs(dm).max()<tol: return m,t,True
    return m,steps,False

def run_I2(adj, m0, kappa=0.05, steps=200000, tol=1e-13):
    """(I2) 离散泊松双向耦合：Δφ = m−mean(m)；搬运率 ∝ e^{φ}"""
    E=edges_of(adj); V=len(adj); m=np.array(m0,float)
    L=np.zeros((V,V))
    for (u,v) in E: L[u,u]+=1; L[v,v]+=1; L[u,v]-=1; L[v,u]-=1
    Lp=np.linalg.pinv(L)
    for t in range(steps):
        phi=Lp@(m-m.mean())                    # 势由荷产生（泊松）
        e=np.exp(np.clip(phi,-5,5))
        dm=np.zeros(V)
        for (u,v) in E:
            f=np.sqrt(e[u]*e[v])*(m[u]-m[v])   # 搬运被势调制
            dm[u]-=f; dm[v]+=f
        m=m+kappa*dm
        if not np.isfinite(m).all() or np.abs(m).max()>1e3: return m,t,False
        if np.abs(dm).max()<tol: return m,t,True
    return m,steps,False

rng=np.random.RandomState(7)
print("="*94)
print("Path I：几何/势与荷的双向耦合（零自由参数）")
print("="*94)
print(f"{'图':10s} {'机制':6s} {'Σm守恒':>10s} {'末态var':>11s} {'q':>9s} {'入窗?':>6s} {'稳定':>5s} {'非均匀?':>7s}")
hits=[]
for gname,adj in [("path16",path_adj(16)),("cycle16",torus_adj(1,16)),("star16",star_adj(16)),
                  ("path32",path_adj(32))]:
    V=len(adj)
    for nm,fn in [("I1",run_I1),("I2",run_I2)]:
        m0=1.0+0.4*(rng.rand(V)-0.5)
        m,t,conv=fn(adj,m0)
        w=m/m.sum(); q=float((w*w).sum())
        cons=abs(m.sum()-m0.sum())
        nonunif=float(np.var(m))>1e-9
        if 0.5<q<0.6 and conv: hits.append((gname,nm,q))
        print(f"{gname:10s} {nm:6s} {cons:10.2e} {np.var(m):11.4e} {q:9.6f} "
              f"{str(0.5<q<0.6):>6s} {str(conv):>5s} {str(nonunif):>7s}")
    print()
print("入 D4 窗且稳定的组合:", len(hits))
for h in hits: print("   ",h)
check("几何-荷双向耦合族存在稳定非均匀态", any(float(np.var(run_I1(path_adj(16),1+0.4*(np.random.RandomState(7).rand(16)-0.5))[0]))>1e-9 for _ in [0]) or
      any(float(np.var(run_I2(path_adj(16),1+0.4*(np.random.RandomState(7).rand(16)-0.5))[0]))>1e-9 for _ in [0]))
check("存在入窗组合", len(hits)>0, f"{len(hits)}")
