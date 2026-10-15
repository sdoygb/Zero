# -*- coding: utf-8 -*-
"""Path B：L1 的"定向/循环次序" ⟹ 非对称输运算子 ⟹ 漂移 ⟹ 非均匀稳态

前轮的 no-go 有一个被忽略的前提：算子**对称**。Z024 定理 2 的证明用的是
"守恒量只能是常数"，但那是对**状态函数**说的，不是对输运算子说的。
Z001 §2 的双覆盖 iota（iota^2=id，无不动点）与 L1 的**循环次序**都给 L1
一个**定向**：每条边有两个方向，走者可以偏向其中一个方向。

若输运算子不对称 ⟹ 细致平衡被破 ⟹ 出现净漂移 ⟹ 平稳测度**非均匀**。
这正是 Z088 §5.2 声称"图类型选择 D4"却没有实现的那个机制。
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

def directed_move(adj, p):
    """定向输运：沿边 u→v 的搬运率 = p（顺向）/ (1−p)（逆向）。
       p=1/2 ⟹ 对称（回到前轮 no-go 的世界）；p≠1/2 ⟹ 有序。"""
    V=len(adj)
    W=np.zeros((V,V))
    for u,nb in enumerate(adj):
        for v in nb:
            if u<v: W[u,v]=p; W[v,u]=1-p
    return W

def evolve(adj, m0, p, kappa=0.1, steps=200000, tol=1e-14):
    V=len(adj); W=directed_move(adj,p)
    m=np.array(m0,float)
    for t in range(steps):
        # 补偿搬运：每个顶点按出边率搬走、按入边率收进（总量守恒）
        out=W.sum(axis=1)           # 出率
        inflow=W.T@m                # 入流
        dm=inflow-out*m
        m=m+kappa*dm
        if not np.isfinite(m).all() or np.abs(m).max()>1e6: return m,t,False
        if np.abs(dm).max()<tol: return m,t,True
    return m,steps,False

V=16
rng=np.random.RandomState(2)
m0=np.ones(V)
print("="*86)
print("Path B：定向输运（p = 顺向率）——对称 p=1/2 时回到前轮 no-go")
print("="*86)
print(f"{'图':10s} {'p':>6s} {'Σm守恒':>10s} {'末态var':>11s} {'q(m)':>9s} {'入窗?':>6s} {'收敛':>5s}  平稳分布形态")
for gname,adj in [("path16",path_adj(16)),("cycle16",torus_adj(1,16)),("star16",star_adj(16))]:
    for p in (0.5, 0.5001, 0.51, 0.7, 0.9):
        m,t,conv=evolve(adj,m0,p)
        w=m/m.sum(); q=float((w*w).sum())
        shape="均匀" if np.var(m)<1e-12 else ("单调漂移" if np.all(np.diff(m)>=-1e-12) else "非平凡")
        print(f"{gname:10s} {p:6.4f} {abs(m.sum()-m0.sum()):10.2e} {np.var(m):11.4e} "
              f"{q:9.6f} {str(0.5<q<0.6):>6s} {str(conv):>5s}  {shape}")
    print()

print("="*86)
print("判定")
print("="*86)
for gname,adj in [("path16",path_adj(16)),("cycle16",torus_adj(1,16))]:
    m,_,_=evolve(adj,m0,0.7)
    q=float((((m/m.sum()))**2).sum())
    check(f"{gname} p=0.7: 平稳分布非均匀", np.var(m)>1e-10, f"var={np.var(m):.4e}")
    check(f"{gname} p=0.7: Σm 守恒", abs(m.sum()-V)<1e-9, f"Σm={m.sum():.6f}")
    print(f"       q={q:.6f}  (D4 窗需 0.5<q<0.6)")
