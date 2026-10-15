# -*- coding: utf-8 -*-
"""Z126 复核 C4：float64 的"丢失"是实现路径敏感，不是精度必然。

把整代更新写成一条线性算子 M_step（作用在 (seed, W) 这个 272 维状态上），
于是可以在毫秒级跑数十万代，并能精确对比不同浮点路径。
  M_step = [[B, D], [C, E]]：seed' = B·seed + D·vec(W)；W' = C·seed + E·vec(W)
由单位输入逐列测出（每条列一次一代），再对大 N 迭代。
"""
from __future__ import annotations
import sys, os
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE,"..","ZeroGPU","src")))
from z0_engine import path_adj
V=16; adj=path_adj(V); N=V*V
A=np.zeros((V,V))
for u,nb in enumerate(adj):
    for w in nb: A[u,w]=1.0

def step_vec(x, dtype=np.float64):
    """x = [seed(V), vec(W)(N)] -> 一代之后的同型向量（用矩阵乘，模拟引擎路径）"""
    seed=x[:V].reshape(V,1); W=x[V:].reshape(V,V)
    Wn=W.copy()
    idx=np.arange(V); Wn[idx,idx]+=x[:V]
    for _ in range(4): Wn=Wn@A
    d=np.array([Wn[i,i] for i in range(V)])
    Wn[idx,idx]=0.0
    return np.concatenate([d, Wn.reshape(-1)])

DIM=V+N
M=np.zeros((DIM,DIM))
for j in range(DIM):
    e=np.zeros(DIM); e[j]=1.0
    M[:,j]=step_vec(e)
print(f"状态维数 {DIM}（seed {V} + W {N}）；算子 M 已构造")
print(f"M 的谱半径 = {max(abs(np.linalg.eigvals(M))):.6f}")
print()

def run(G, dtype=np.float64, normalize=True):
    x=np.zeros(DIM,dtype=dtype)
    x[:V]=np.array([2**(V-1-i) for i in range(V)],dtype=dtype)
    Md=M.astype(dtype)
    for _ in range(G):
        x=Md@x
        if normalize:
            tot=x[:V].sum()
            if tot>0: x=x/tot          # 同一公因子 ⟹ 比例严格不变
    seed=x[:V]; w=seed/seed.sum()
    return float((w*w).sum()), int((w>1e-12).sum())

print("="*84)
print("路径 Q（矩阵乘 = 引擎口径）在不同代数下的结果")
print("="*84)
print(f"{'gen':>7s} {'q':>16s} {'活顶点':>6s}")
for G in (100,1000,10000,30000,60000,100000,200000,500000):
    q,nz=run(G)
    print(f"{G:7d} {q:16.10f} {nz:6d}")
print()
print("="*84)
print("与精确值对照")
print("="*84)
print("  精确（整数/Fraction 与 mpmath 50 位，双独立实现）：q=0.156399403095  16 活")
q1,n1=run(1000); q2,n2=run(500000)
print(f"  float64 算子化， 1000 代：q={q1:.10f}  {n1} 活")
print(f"  float64 算子化，5×10^5 代：q={q2:.10f}  {n2} 活")
print()
print("  引擎默认（不逐代归一化，靠内部 1e250 阈值）：6 万代 q=0.5555555  2 活")
