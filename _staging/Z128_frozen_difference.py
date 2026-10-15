# -*- coding: utf-8 -*-
"""Z128：**修正** Q4 —— 线性忠实律**能**给出非均匀极限，via 冻结差。

对称互驱律（双向、忠实、零自由参数）：
    m ← m − κΔ(m+φ)
    φ ← φ − κΔ(m+φ)
解析：d(m+φ)/dt = −2Δ(m+φ) → 常数；  d(m−φ)/dt = 0（**冻结**）
      ⟹ m_∞ = ½(c + (m₀−φ₀))，一般**非均匀**。
故"线性 ⟹ 均匀"是**错的**；正确命题见下。
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

def evolve(adj,m0,p0,kappa=0.05,steps=60000):
    D=lap(adj); m=np.array(m0,float); p=np.array(p0,float)
    for _ in range(steps):
        d=-(D@(m+p)); m=m+kappa*d; p=p+kappa*d
    return m,p

print("="*88)
print("① 解析预测 vs 数值（m_∞ = ½(c + (m₀−φ₀))，c = mean(m₀+φ₀)）")
print("="*88)
n=16; adj=path_adj(n)
rng=np.random.RandomState(1)
for tag,(m0,p0) in [("随机 m0,φ0",(rng.rand(n),rng.rand(n))),
                    ("m0 集中端点",(np.array([1.0]+[0.0]*(n-1)),np.zeros(n))),
                    ("m0 两端峰",(np.array([1.0]+[0.0]*14+[1.0]),np.zeros(n)))]:
    m,p=evolve(adj,m0,p0)
    c=(np.array(m0,float)+np.array(p0,float)).mean()
    pred=0.5*(c+(np.array(m0,float)-np.array(p0,float)))
    err=np.abs(m-pred).max()
    print(f"  {tag:12s}: max|m_数值 − m_解析| = {err:.2e}   "
          f"q(m)={float(((m/m.sum())**2).sum()):.6f}")
    check(f"① {tag}: 数值 = 解析（冻结差）", err<1e-6, f"err={err:.2e}")

print()
print("="*88)
print("② q_∞ 扫：能否入 D4 窗 (1/2,3/5)？—— 扫初值差的各种形状")
print("="*88)
base=np.ones(n)
shapes={
 "均匀差":        np.zeros(n),
 "单点差 @7":     np.array([0.0]*7+[1.0]+[0.0]*8),
 "两点差 @7,8 (2:1)": np.array([0.0]*7+[2.0,1.0]+[0.0]*7),
 "两点差 @7,8 (1:1)": np.array([0.0]*7+[1.0,1.0]+[0.0]*7),
 "三点差 @6,7,8":  np.array([0.0]*6+[1.0,1.0,1.0]+[0.0]*7),
 "几何差":        np.array([2.0**(n-1-i) for i in range(n)],float),
 "阶梯差":        np.array([1.0]*8+[0.0]*8),
 "半正弦":        np.sin(np.pi*np.arange(n)/(n-1)),
}
print(f"  {'初值差形状':22s} {'q_∞':>12s} {'入窗?':>6s}")
hits=[]
for nm,dv in shapes.items():
    m0=base+0.5*dv; p0=base-0.5*dv          # 差 = dv
    m,p=evolve(adj,m0,p0)
    w=m/m.sum(); q=float((w*w).sum())
    ok=0.5<q<0.6; hits.append((nm,q,ok))
    print(f"  {nm:22s} {q:12.8f} {str(ok):>6s}")
check("② 存在初值差形状使 q_∞ 入 D4 窗", any(h[2] for h in hits),
      f"{sum(1 for h in hits if h[2])}/{len(hits)} 个形状入窗")

print()
print("="*88)
print("③ 修正后的命题（Q4 的正确形式）")
print("="*88)
print("""
  旧（错）："忠实线性律 ⟹ 极限均匀 ⟹ q→1/V ⟹ 不能入窗"
  新（验）：极限 = P_ker m₀；q_∞ 由 **ker 的维数/结构** 与 **初值在 ker 上的投影** 决定。
     · 若 ker 只含"常数场"  ⟹ q_∞ = 1/V（原结论成立）
     · 若 ker 含"冻结差"    ⟹ q_∞ 可由初值差任意调节 ⟹ **可入窗**
  对称互驱的 ker 维数 = 17（含 (u,−u) 方向）⟹ 属第二种。

  \033[1m故：忠实线性律**能**把 q 保持在窗内，但不能**产生**入窗所需的集中度——
     集中度仍来自初值差。\033[0m
""")
# 核维数双路核对
J=-np.block([[lap(adj),lap(adj)],[lap(adj),lap(adj)]])
w=np.linalg.eigvals(J); k=int(np.sum(np.abs(w)<1e-9))
check("③ 对称互驱的核维数 = 17（含 (u,−u) 方向）", k==17, f"dim ker={k}")
# 对照：单向级联
J2=-np.block([[lap(adj),np.zeros((n,n))],[lap(adj),lap(adj)]])
w2=np.linalg.eigvals(J2); k2=int(np.sum(np.abs(w2)<1e-9))
check("③ 对照：单向级联核维数 = 2（只有常数）", k2==2, f"dim ker={k2}")

n_=sum(1 for _,ok in CHECKS if ok)
print(); print("="*88); print(f"断言通过 {n_}/{len(CHECKS)}"); print("="*88)
