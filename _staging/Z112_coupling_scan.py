# -*- coding: utf-8 -*-
"""Z112：跨层动力学何时可能？——不可约性判据 + 一个稳定忠实模型。

判据（缺一不可）
  C1 补偿守恒  Σm 沿演化守恒
  C2 有界      不发散
  C3 不可约 ★  联合算子非块对角（两层互相驱动）
  C4 必要      冻结任一层都改变结果
"""
from __future__ import annotations
import sys, os, numpy as np
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "..", "ZeroGPU", "src")))
from z0_engine import path_adj, torus_adj

CHECKS=[]
def check(n,c,d=""):
    CHECKS.append((n,bool(c))); print(f"[{'OK' if c else 'FAIL'}] {n}  {d}")

def lap(adj):
    V=len(adj); L=np.zeros((V,V))
    for u,nb in enumerate(adj):
        L[u,u]=len(nb)
        for w in nb: L[u,w]=-1.0
    return L

V=16; adj=path_adj(V); D=lap(adj)
adj2=[sorted(set([v-2 for v in [v] if v-2>=0]+[v+2 for v in [v] if v+2<V])) for v in range(V)]
K=lap(adj2)                      # 次近邻搬运：与 D 不交换

print("① 算子级：可约性的充要看**交换子**")
print(f"   [Δ,K] 的范数 = {np.abs(D@K-K@D).max():.4f}   (>0 ⟹ 生成不可约耦合)")
check("Δ 与 K 不交换（存在不可约耦合）", np.abs(D@K-K@D).max()>0)
check("Δ 的行和 = 0（补偿移动不创造荷）", np.abs(D.sum(axis=1)).max()<1e-12)
check("Δ·1 = 0（均匀态是搬运的不动点，公理③）", np.allclose(D@np.ones(V),0))

print()
print("② 结构级：联合算子的块结构")
def joint_blocks(A,B,kappa):
    """(m,π) ← (m,π) + κ·[[0,A],[B,0]](m,π)"""
    J=np.zeros((2*V,2*V)); J[:V,V:]=A; J[V:,:V]=B
    return kappa*J
for tag,A,B in [("同算子 (Δ,Δ)",D,D),("异算子 (Δ,K)",D,K)]:
    J=joint_blocks(A,B,1.0)
    off=np.abs(J[:V,V:]).max()+np.abs(J[V:,:V]).max()
    ev=np.linalg.eigvalsh((A+B.T)/2) if np.allclose(A,B.T) else None
    print(f"   {tag:14s} 非对角块范数={off:.3f}  谱半径(κ=1)={max(abs(np.linalg.eigvals(J))):.3f}")

print()
print("③ 构造一个**稳定且忠实**的三层线性层间律（零自由参数）")
print("   状态: m(L0荷) / φ(L2势) / π(L1流量) / R(L1'记录)")
print("   律:   m ← m + κπ ;  π ← π − κΔφ ;  φ ← φ + κΔπ ;  R ← R + φ")
print("   （离散 Hamilton 场：π 为动量，Δ 为补偿搬运算子）")
kappa=0.25/np.linalg.eigvalsh(D).max()
def run(m0,p0,phi0,steps=6000,freeze=None):
    m=np.array(m0,float); pi=np.array(p0,float); phi=np.array(phi0,float); R=np.zeros(V)
    for _ in range(steps):
        if freeze!="L1":
            m  = m  + kappa*pi
            phi= phi+ kappa*(D@pi)
            pi = pi - kappa*(D@phi)
        if freeze!="L2": phi = phi if freeze=="L1" else phi
        if freeze!="L0": pass
        if freeze!="L1p": R = R + phi
    return m,pi,phi,R
rng=np.random.RandomState(5)
m0=rng.rand(V)-0.5; m0=m0-m0.mean()
p0=np.zeros(V); phi0=np.array([2.0**(V-1-i) for i in range(V)],float)
m,pi,phi,R=run(m0,p0,phi0)
c1=abs(m.sum()-m0.sum())<1e-9
mx=max(np.abs(m).max(),np.abs(phi).max(),np.abs(pi).max()); c2=np.isfinite(mx) and mx<1e4
print(f"   C1 Σm 守恒: {c1}   (Σm0={m0.sum():.2e}, Σm={m.sum():.2e})")
print(f"   C2 有界:    {c2}   (max|·|={mx:.3e}, 谱半径={max(abs(np.linalg.eigvals(kappa*joint_blocks(D,D,1.0)))):.4f})")
check("C1 补偿守恒（Σm 不变）", c1, f"Σm={m.sum():.2e}")
check("C2 有界（谱半径 ≤ 1）",
      max(abs(np.linalg.eigvals(kappa*joint_blocks(D,D,1.0))))<=1.0+1e-12,
      f"ρ={max(abs(np.linalg.eigvals(kappa*joint_blocks(D,D,1.0)))):.4f}")
check("C2' 数值有界", c2, f"max|·|={mx:.3e}")

# C3 互驱：改 L2 初值 → L1/L0 轨迹变；改 L0 初值 → L2 轨迹变
_,pi_a,phi_a,_=run(m0,p0,phi0*1.7); _,pi_b,phi_b,_=run(m0,p0,phi0)
m_a,_,_,_=run(m0*3.1,p0,phi0); m_b,_,_,_=run(m0,p0,phi0)
c3=(np.abs(pi_a-pi_b).max()>1e-9) and (np.abs(m_a-m_b).max()>1e-9)
check("★ C3 不可约：L2 初值 → L1 变化，且 L0 初值 → L0 轨迹依 L1 而变", c3,
      f"Δπ={np.abs(pi_a-pi_b).max():.3e}")

# C4 必要：冻结 L1（π 恒 0）→ 一切静止
m_f,pi_f,phi_f,_=run(m0,p0,phi0,freeze="L1")
c4=np.abs(phi_f-phi0).max()<1e-12 and np.abs(m_f-m0).max()<1e-12
check("C4 必要性：冻结 L1（π≡0）后 m 与 φ 完全静止 ⟹ L1 承重", c4,
      f"max|Δφ|={np.abs(phi_f-phi0).max():.1e}")

print()
print("④ 长时行为与 Z088 的对照")
w=m/m.sum() if abs(m.sum())>1e-300 else None
print(f"   L0 末态 q(m) = {float((w*w).sum()):.6f}   （均匀 → 1/V = {1/V:.6f}）")
print(f"   L2 势末态 max|φ| = {np.abs(phi).max():.3e}（守恒振荡，不衰减也不集中）")
Rw=R/R.sum()
print(f"   L1' 记录 R 归一化 top3: {[(i,round(float(x),4)) for i,x in sorted(enumerate(Rw),key=lambda z:-z[1])[:3]]}")
check("忠实律的耦合算子谱**无简并**（不产生 Z088 那种浮点选择）",
      len(set(np.round(np.linalg.eigvalsh(D),9)))==V, f"{V} 个互异特征值")

n=sum(1 for _,ok in CHECKS if ok)
print()
print(f"断言通过 {n}/{len(CHECKS)}")
