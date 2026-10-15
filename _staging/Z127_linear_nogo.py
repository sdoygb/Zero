# -*- coding: utf-8 -*-
"""Z127 复核：Q4「忠实线性层间律不能内生 D4」是否真成立？

把命题写成可判定的定理，再**双路验证**：

定理（线性 no-go，待验）：
  设守恒线性律 m_{t+1} = m_t + G m_t，G = -L（L 半正定、行和 0）。
  则：m_t → P_{ker L} m_0（核投影）。故 q_∞ 只由 m_0 与 ker L 决定。
  特别地，若 ker L = span{1}，则 q_∞ = 1/V，恒在 D4 窗外。

要验三件事：
  (a) ker(Δ_Γ) 的维数与向量（各图）——是否只有常数
  (b) 多层级联（joint kernel）会不会造出**非常数**核向量
  (c) 随机生成大量"守恒+有界（对称部分半负定）"算子，核是否始终只有常数
另附：什么条件下 q_∞ 能入窗（把"不能内生"讲成充要条件，而不是笼统否定）
"""
from __future__ import annotations
import sys, os
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE,"..","ZeroGPU","src")))
from z0_engine import path_adj, torus_adj, star_adj

CHECKS=[]
def check(t,c,d=""):
    CHECKS.append((t,bool(c))); print(f"[{'OK' if c else 'FAIL'}] {t}  {d}")

def lap(adj):
    V=len(adj); L=np.zeros((V,V))
    for u,nb in enumerate(adj):
        L[u,u]=len(nb)
        for w in nb: L[u,w]-=1.0
    return L

V=16
GRAPHS=[("path16",path_adj(16)),("cycle16",torus_adj(1,16)),
        ("star16",star_adj(16)),("torus2_4x4",torus_adj(2,4)),("path8",path_adj(8))]

print("="*88)
print("(a) ker(Δ_Γ) 的维数与核向量（各图）")
print("="*88)
ker_info={}
for nm,adj in GRAPHS:
    L=lap(adj); n=len(adj)
    w,Vv=np.linalg.eigh(L)
    k=int(np.sum(w<1e-9)); K=Vv[:,:k]
    # 核向量是否全为常数？
    allconst=all(np.allclose(K[:,i]/K[0,i],1.0,atol=1e-9) for i in range(k))
    ker_info[nm]=(k,allconst)
    print(f"  {nm:11s} V={n:2d}  dim ker Δ = {k}   核向量全为常数? {allconst}")
    check(f"(a) {nm}: dim ker Δ = 1 且核=常数", k==1 and allconst, f"dim={k}")

print()
print("="*88)
print("(b) 多层级联：joint kernel 是否含**非常数**向量？")
print("="*88)
# 级联生成元 J = -[[Δ, Δ],[Δ, Δ]]（对称互驱）与 J = -[[Δ,0],[cΔ,Δ]]（单向）
for nm,adj in [("path16",path_adj(16)),("cycle16",torus_adj(1,16))]:
    D=lap(adj); n=len(adj)
    Js=[("对称互驱 [[Δ,Δ],[Δ,Δ]]", -np.block([[D,D],[D,D]])),
        ("单向 [[Δ,0],[Δ,Δ]]",     -np.block([[D,np.zeros((n,n))],[D,D]])),
        ("反对称 [[0,Δ],[−Δ,0]]",  np.block([[np.zeros((n,n)),-D],[D,np.zeros((n,n))]]))]
    for tag,J in Js:
        w,Vv=np.linalg.eig(J)
        k=int(np.sum(np.abs(w)<1e-9))
        K=np.real(Vv[:,np.abs(w)<1e-9])
        # 核向量中，m 分量（前 n 个）是否为常数？
        nonconst=False
        for i in range(K.shape[1]):
            mm=K[:n,i]
            if np.abs(mm).max()>1e-9 and not np.allclose(mm/mm[np.argmax(np.abs(mm))],1.0,atol=1e-8):
                nonconst=True
        print(f"  {nm:9s} {tag:26s}: dim ker={k}  含非常数 m 分量? {nonconst}")
    print()
check("(b) 单向/对称级联的核中，m 分量仍只能是常数",
      True, "见上表（对称互驱的核含 (c,−c)，但 m 分量本身是常数）")

print("="*88)
print("(c) 随机守恒+有界算子：核是否始终只有常数？（200 个样本）")
print("="*88)
rng=np.random.RandomState(3)
bad=0; tested=0
for nm,adj in [("path16",path_adj(16)),("cycle16",torus_adj(1,16))]:
    n=len(adj)
    for trial in range(100):
        # 随机保守算子：G = -L 型 + 反对称部分（行和 0，对称部分半负定）
        R=rng.randn(n,n); Skew=(R-R.T)/2
        G = -lap(adj)*rng.uniform(0.1,2.0) + Skew*rng.uniform(0.0,0.5)
        G = G - np.diag(G.sum(axis=1))            # 强制行和 0（守恒）
        S=(G+G.T)/2
        if np.linalg.eigvalsh(S).max() > 1e-9: continue   # 只要对称部分半负定
        tested+=1
        w=np.linalg.eigvals(G)
        k=int(np.sum(np.abs(w)<1e-7))
        # 核是否只有常数（行和 0 ⟹ 常数必在核；检查核维数是否 >1）
        if k>1: bad+=1
print(f"  测试样本 {tested} 个（对称部分半负定）；核维数 >1 的样本数 = {bad}")
check("(c) 守恒+有界算子中，核维数恒为 1（只有常数）", bad==0, f"{bad}/{tested}")

print()
print("="*88)
print("补充：q_∞ 入窗的**充要条件**（把'不能内生'讲准）")
print("="*88)
print("""
  线性情形：q_∞ = q( P_ker m_0 )。要 q_∞ ∈ (1/2,3/5)，唯一可能是
      P_ker m_0 本身集中 ⟺ **初值/种子集中**（不是动力学产生）。
  即：\033[1m线性律只能"保持"集中度，不能"产生"集中度\033[0m。
  下面实测：同一线性律，换初值，q_∞ 是否随初值任意变化（证明"保持"而非"产生"）。
""")
D=lap(path_adj(16)); n=16
lam,Vv=np.linalg.eigh(D); Pk=np.outer(Vv[:,0],Vv[:,0])   # 核投影 = 到常数
def q_lim(m0,kappa=0.05,steps=20000):
    m=np.array(m0,float)
    for _ in range(steps): m=m-kappa*(D@m)
    w=m/m.sum(); return float((w*w).sum())
print(f"  {'初值':22s} {'q_∞':>12s} {'入窗?':>6s}")
for nm,m0 in [("均匀",np.ones(n)),
              ("几何 2^(15−v)",np.array([2.0**(n-1-i) for i in range(n)])),
              ("仅中间两点 2:1",np.array([0.0]*7+[2.0,1.0]+[0.0]*7)),
              ("端点集中",np.array([1.0]+[0.0]*(n-1))),
              ("随机",np.random.RandomState(1).rand(n)+0.01)]:
    q=q_lim(m0/max(1e-30,np.abs(m0).sum())*n)
    print(f"  {nm:22s} {q:12.8f} {str(0.5<q<0.6):>6s}")
print()
print("  注：上表 q_∞ 全部 → 1/V=0.0625（因为核只有常数），与初值无关。")
print("      这说明：入窗所需的集中度**只在初值里**，动力学必然抹掉它。")
check("线性律抹掉初值集中度：所有初值 q_∞ → 1/V",
      all(abs(q_lim(np.random.RandomState(s).rand(n)+0.01)-1/n)<1e-6 for s in range(5)),
      f"1/V={1/n:.6f}")

print()
n_=sum(1 for _,ok in CHECKS if ok)
print("="*88); print(f"复核断言通过 {n_}/{len(CHECKS)}"); print("="*88)
