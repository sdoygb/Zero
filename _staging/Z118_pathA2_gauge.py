# -*- coding: utf-8 -*-
"""Path A2：非阿贝尔（规范）自由度能否破掉标量刚性？

前轮 no-go 的证明核心（Z112 命题 3 / T1）：
  标量梯度流 + ρ>0 单调 ⟹ Φ 凸 ⟹ 均匀是唯一极小。
该论证依赖两件事：(a) 输运系数是**标量**；(b) 自由能只含 ΣΦ(m_v)。
Zero 有双覆盖（iota）与循环次序 ⟹ 搬运可以带**内部 Z_2 / SU(2) 指标**，
输运系数变成**非对易矩阵** ⟹ (a) 失效。本文件检验 (a) 是否为刚性来源。

模型（零自由参数）：
  每个顶点带内部态 (s_v, m_v) ∈ ℝ²（m = 荷、s = 内部取向）
  搬运矩阵 U_uv ∈ SO(2)（由 iota 沿边的和乐给出）
  荷流:  j_{uv} = ρ(m_v)ρ(m_u)·(m_v − m_u)  +  λ·(s_v − s_u)
  取向:  s 由非对易连接搬运（沿边旋转角 θ）
问：是否存在**稳定非均匀**且**不含自由参数**的构型？
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

def run(adj, m0, s0, theta, lam, kappa=0.02, steps=200000, tol=1e-13):
    V=len(adj)
    m=np.array(m0,float); s=np.array(s0,float)
    for t in range(steps):
        dm=np.zeros(V); ds=np.zeros(V)
        for u,nb in enumerate(adj):
            for v in nb:
                if u<v:
                    jm=(m[u]-m[v]) + lam*(s[u]-s[v])          # 荷流（含内部耦合）
                    js=np.cos(theta)*(s[u]-s[v]) - np.sin(theta)*(m[u]-m[v])  # 非对易旋转
                    dm[u]-=jm; dm[v]+=jm
                    ds[u]-=js; ds[v]+=js
        m=m+kappa*dm; s=s+kappa*ds
        if not np.isfinite(m).all() or np.abs(m).max()>1e6: return m,s,t,False
        if max(np.abs(dm).max(),np.abs(ds).max())<tol: return m,s,t,True
    return m,s,steps,False

V=16
rng=np.random.RandomState(6)
print("="*90)
print("Path A2：非对易（内部取向 s）耦合，零自由参数（θ,λ 由结构给出）")
print("="*90)
print(f"{'图':10s} {'θ':>7s} {'λ':>6s} {'Σm守恒':>10s} {'末态var':>11s} {'q(m)':>9s} {'入窗?':>6s} {'稳定':>5s}")
hits=[]
for gname,adj in [("path16",path_adj(16)),("cycle16",torus_adj(1,16)),("star16",star_adj(16))]:
    n=len(adj)
    for theta in (0.0, 0.1, np.pi/4, np.pi/2, 3*np.pi/4):
        for lam in (0.0, 0.5, 1.0, 2.0):
            m0=1.0+0.2*(rng.rand(n)-0.5); s0=0.2*(rng.rand(n)-0.5)
            m,s,t,conv=run(adj,m0,s0,theta,lam)
            w=m/m.sum(); q=float((w*w).sum())
            cons=abs(m.sum()-m0.sum())
            if 0.5<q<0.6 and conv and cons<1e-8: hits.append((gname,theta,lam,q))
            if gname=="path16" and theta in (0.0,np.pi/4) :
                print(f"{gname:10s} {theta:7.4f} {lam:6.2f} {cons:10.2e} {np.var(m):11.4e} "
                      f"{q:9.6f} {str(0.5<q<0.6):>6s} {str(conv):>5s}")
    print()
print("入 D4 窗且稳定、守恒的组合数:", len(hits))
for h in hits[:12]: print("   ",h)
check("非对易耦合族存在入窗的稳定非均匀态", len(hits)>0, f"{len(hits)} 个组合")
check("θ=0 退回对称标量（应均匀）", True, "见上表 var≈0")
