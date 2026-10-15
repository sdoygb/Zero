# -*- coding: utf-8 -*-
"""Path F：把不对称性交给**高度**（词自己的结构），而不是交给自由参数 p。

Path B 已证：定向输运 p≠1/2 给出非均匀平稳分布，q(p) 可入 D4 窗。
但 p 是自由参数 ⟹ 违反公理③（"不设选择规则"）。出路：
让 p 由**零和约束本身**决定 —— 走者的高度 h 就是词的部分和 Σw。
零和 ⟹ 每次闭合必须回到 h=0 ⟹ 高度是内生的。

约定：移动 ±1，但**相对高度的搬运率不对称**，由高度定：
   rate(h→h+1)/rate(h→h-1) = r(h)，r 只能依赖 h（不能依赖外部取舍）。
本文件扫描 r(h) 的若干**零参数**形式，看平稳分布能否非均匀、且 q 入窗。
"""
from __future__ import annotations
import sys, os, itertools
import numpy as np
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "..", "ZeroGPU", "src")))
from z0_engine import path_adj

CHECKS=[]
def check(n,c,d=""):
    CHECKS.append((n,bool(c))); print(f"[{'OK' if c else 'FAIL'}] {n}  {d}")

def stationary(generator, V):
    """求 Markov 生成元（列和 0）的平稳分布：左零空间，归一化。"""
    A=generator.T
    w,v=np.linalg.eig(A)
    i=np.argmin(np.abs(w))
    p=np.real(v[:,i]); p=p/p.sum()
    return np.abs(p)

def gen_height(H, r):
    """高度 h∈[0,H]（0 为闭合点）。rate(h→h+1)=r(h), rate(h→h−1)=1。
       边界 h=0: 只能向上；h=H: 只能向下（反射，对应"最长词"）。"""
    n=H+1; G=np.zeros((n,n))
    for h in range(n):
        up = r(h) if h<H else 0.0
        dn = 1.0  if h>0 else 0.0
        if h<H: G[h+1,h]+=up; G[h,h]-=up
        if h>0: G[h-1,h]+=dn; G[h,h]-=dn
    return G

H=20
FORMS={
 "r=1（对称）":            lambda h: 1.0,
 "r=(h+1)/(h+2)":          lambda h: (h+1)/(h+2),
 "r=1/(h+1)":              lambda h: 1.0/(h+1),
 "r=(H−h)/(h+1)":          lambda h: (H-h)/(h+1),
 "r=(h+1)/(H−h+1)":        lambda h: (h+1)/(H-h+1),
 "r=1+1/(h+1)":            lambda h: 1.0+1.0/(h+1),
 "r=1−1/(h+2)":            lambda h: 1.0-1.0/(h+2),
}
print("="*92)
print(f"Path F：高度依赖的搬运率，H={H}（高度格点），平稳分布与 q")
print("="*92)
print(f"{'r(h) 形式':20s} {'均匀?':>6s} {'q':>10s} {'入D4窗?':>8s}   平稳分布在 h 上的形态")
best=[]
for nm,r in FORMS.items():
    G=gen_height(H,r); p=stationary(G,H+1)
    q=float((p*p).sum())
    unif=float(np.var(p))
    best.append((abs(q-5/9),nm,q))
    top=[(h,round(float(x),5)) for h,x in enumerate(p) if x>0.02][:5]
    print(f"{nm:20s} {str(unif<1e-12):>6s} {q:10.6f} {str(0.5<q<0.6):>8s}   {top}")
print()
best.sort()
print(f"最接近 5/9 的形式: {best[0][1]}  q={best[0][2]:.6f}  偏差={best[0][0]:.2e}")
check("存在高度依赖的零参数形式给出非均匀平稳分布",
      any(float(np.var(stationary(gen_height(H,r),H+1)))>1e-9 for r in FORMS.values()))
print()
print("注：H=20 是自由参数（词长上界）。下面检验 q 是否随 H 稳定（尺度无关性）：")
for nm,r in [("r=1/(h+1)",FORMS["r=1/(h+1)"]),("r=(h+1)/(h+2)",FORMS["r=(h+1)/(h+2)"])]:
    qs=[]
    for HH in (10,20,40,80,160):
        G=gen_height(HH,r); p=stationary(G,HH+1); qs.append(float((p*p).sum()))
    print(f"  {nm:16s} q(H) = "+", ".join(f"{x:.5f}" for x in qs))
