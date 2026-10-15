# -*- coding: utf-8 -*-
"""Path K：手性 p 的内生来源 —— 把走者放在 (k,h) 网格上。

关键纠正：单条空间边上的偏向在路径图上**自动可逆**（Path F 已证），
所以"手性"不能定义在空间边上。但 Zero 的走者实际走在

        (k, h)   k = 第几步 (0..n),  h = 部分和 Σ w  (平衡词要求 h_n = 0)

这个网格是**有边界的**（h 不能越过结束条件），而边界正是
"左零空间 ≠ 常数"的来源。本文件检验三条**零自由参数**的规则：

  (K1) 均匀：每步 ±1 等率            → 平稳分布 = ？
  (K2) 计数序：每步的走向偏向"剩余可闭合路径更多"的一侧
       （即：优先选择使完成方式数更大的方向）—— 由组合数定义，无参数
  (K3) 双覆盖序：iota 把 ±走向配对，配对后的次序给出唯一偏向

对每条规则求平稳分布 ω（= 左零空间），算 q、算 D*，与 D4 窗比对。
"""
from __future__ import annotations
import sys, os
import numpy as np
from math import comb

CHECKS=[]
def check(n,c,d=""):
    CHECKS.append((n,bool(c))); print(f"[{'OK' if c else 'FAIL'}] {n}  {d}")

def build_grid(n):
    """状态 (k,h)：k=0..n 步，h 为部分和。h 的范围 [-n,n]，但只保留能回到 0 的。"""
    states=[]; idx={}
    # 从 (0,0) 出发，只保留能在剩余步数内回到 0 的状态
    for k in range(n+1):
        for h in range(-n, n+1):
            if (n-k-h) % 2 == 0 and abs(h) <= n-k:
                if (k,h) not in idx:
                    idx[(k,h)]=len(states); states.append((k,h))
    return states, idx

def completions(n,k,h):
    """从 (k,h) 出发、n 步内回到 h=0 的路径数"""
    r=n-k
    if (r-h) % 2: return 0
    up=(r+h)//2
    if up<0 or up>r: return 0
    return comb(r, up)

def stationary_of(rate_up, rate_dn, states, idx, n):
    """构造生成元（列和 0），求左零空间（平稳分布）"""
    S=len(states); G=np.zeros((S,S))
    for i,(k,h) in enumerate(states):
        if k==n: continue
        for dh, r in ((+1, rate_up(k,h)), (-1, rate_dn(k,h))):
            j=idx.get((k+1,h+dh))
            if j is None or r<=0: continue
            G[j,i]+=r; G[i,i]-=r
    w,v=np.linalg.eig(G.T)
    m=np.argmin(np.abs(w))
    p=np.real(v[:,m]); p=np.abs(p); s=p.sum()
    return p/s if s>1e-300 else p

def q_of(p):
    return float((p*p).sum())

n=16
states,idx=build_grid(n)
print(f"(k,h) 网格：n={n}，状态数={len(states)}（含边界约束）")
print()
RULES={
 "K1 均匀（±等率）":      (lambda k,h: 1.0, lambda k,h: 1.0),
 "K2 计数序（偏向可完成侧）": (lambda k,h: completions(n,k+1,h+1),
                          lambda k,h: completions(n,k+1,h-1)),
 "K3 反计数序":            (lambda k,h: completions(n,k+1,h-1),
                          lambda k,h: completions(n,k+1,h+1)),
}
print(f"{'规则':26s} {'平稳q':>10s} {'D*':>4s} {'入窗?':>6s}  边缘 h 分布 (|h|)")
import math
for nm,(ru,rd) in RULES.items():
    p=stationary_of(ru,rd,states,idx,n)
    q=q_of(p)
    # 归并到 |h|（对 k 求和）
    hm={}
    for i,(k,h) in enumerate(states):
        hm[abs(h)]=hm.get(abs(h),0.0)+p[i]
    tot=sum(hm.values()); hm={k:v/tot for k,v in hm.items()}
    qh=sum(v*v for v in hm.values())
    # 账本峰
    Dbest=max(range(2,30), key=lambda D: math.comb(D,2)*q**D)
    top=[(h,round(hm.get(h,0),4)) for h in range(0,8)]
    print(f"{nm:26s} {q:10.6f} {Dbest:4d} {str(0.5<q<0.6):>6s}  {top}")

print()
print("="*90)
print("关键检验：|h| 边缘分布（这才是与 Z001 层高划分同型的对象）")
print("="*90)
bh_target=[256,4118,4858,2518,880,208,30,2]
tot_t=sum(bh_target); bt=[x/tot_t for x in bh_target]
print(f"{'h':>3s} " + " ".join(f"{nm[:10]:>10s}" for nm in RULES) + f" {'Z001 b_h':>10s}")
hm_all={}
for nm,(ru,rd) in RULES.items():
    p=stationary_of(ru,rd,states,idx,n)
    hm={}
    for i,(k,h) in enumerate(states): hm[abs(h)]=hm.get(abs(h),0.0)+p[i]
    s=sum(hm.values()); hm_all[nm]={k:v/s for k,v in hm.items()}
for h in range(1,9):
    row=" ".join(f"{hm_all[nm].get(h,0):10.5f}" for nm in RULES)
    print(f"{h:3d} {row} {bt[h-1]:10.5f}")
print()
print("★ 判读：K2（计数序）是否给出与 Z001 的 b_h 同型的非均匀 |h| 分布？")
for nm in RULES:
    hm=hm_all[nm]; qh=sum(v*v for v in hm.values())
    print(f"   {nm:26s} q(|h|)={qh:.6f}  入窗={0.5<qh<0.6}")
