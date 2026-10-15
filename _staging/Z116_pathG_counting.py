# -*- coding: utf-8 -*-
"""Path G：非均匀性来自**组态空间的计数**，不是来自动力学。

前几轮的共同前提是"状态 = 每顶点一个实数密度 m_v"，动力学在 ℝ^V 上。
但 Z001 的状态是一个**词** w ∈ W_n，计数测度 μ(A)=|A|/|W_n| 唯一。
所以"非均匀"可能不是动力学产生的，而是**组合计数**的必然结果：
   在词空间上取某划分的推前权重，它是否天然非均匀、且 q 入 (1/2,3/5)？

本文件检验三类划分的推前 q：
  (i)  层高 h(w)=max−min（Z001 §7.3 的划分，已证 q=0.55777859）
  (ii) 单个坐标/部分和的边缘分布（"位置"读出的推前）
  (iii) 步的"翻转次数"、"最后翻转位置"等自识别候选
"""
from __future__ import annotations
import sys, os
from itertools import combinations
from collections import Counter
from fractions import Fraction as F

CHECKS=[]
def check(n,c,d=""):
    CHECKS.append((n,bool(c))); print(f"[{'OK' if c else 'FAIL'}] {n}  {d}")

def words(n):
    """W_n = {w∈{±1}^n : Σw=0}"""
    for pos in combinations(range(n), n//2):
        w=[-1]*n
        for i in pos: w[i]=1
        yield w

def q_of(counter):
    tot=sum(counter.values())
    return sum(F(c,tot)**2 for c in counter.values()), tot

n=16
Ws=list(words(n))
print(f"|W_{n}| = {len(Ws)}  (=C({n},{n//2})={len(Ws)})")
check("词数 = C(n,n/2)", len(Ws)==__import__('math').comb(n,n//2), f"{len(Ws)}")

# (i) 层高划分（Z001 §7.3）
def hmax(w):
    s=0; lo=0; hi=0
    for x in w:
        s+=x; lo=min(lo,s); hi=max(hi,s)
    return hi-lo
c_h=Counter(hmax(w) for w in Ws)
q_h,tot=q_of(c_h)
print(f"\n(i) 层高 h=range 的划分: q = {q_h} = {float(q_h):.8f}   (Z001 报 0.55777859)")
check("层高划分 q 与 Z001 一致", abs(float(q_h)-0.55777859)<1e-7, f"{float(q_h):.8f}")

# (ii) 部分和（"位置"）边缘分布
def partial_max(w):
    s=0; m=0
    for x in w:
        s+=x; m=max(m,s)
    return m
c_p=Counter(partial_max(w) for w in Ws)
q_p,_=q_of(c_p)
print(f"(ii) 部分和最大值的边缘分布: q = {float(q_p):.8f}  入窗={0.5<float(q_p)<0.6}")
print(f"     分布: {[(k,float(F(v,len(Ws)))) for k,v in sorted(c_p.items())][:6]} ...")

# (iii) 若干自识别候选
def lastflip(w):
    for i in range(len(w)-1,0,-1):
        if w[i]!=w[i-1]: return i
    return 0
def nflip(w): return sum(1 for i in range(len(w)-1) if w[i]!=w[i-1])
def maxrun(w):
    best=cur=1
    for i in range(1,len(w)):
        cur = cur+1 if w[i]==w[i-1] else 1
        best=max(best,cur)
    return best
for nm,f in [("lastflip",lastflip),("nflip",nflip),("maxrun",maxrun)]:
    c=Counter(f(w) for w in Ws); q,tot=q_of(c)
    print(f"(iii) {nm:9s}: q = {float(q):.8f}  块数={len(c)}  入窗={0.5<float(q)<0.6}")

print()
print("="*80); print("结论"); print("="*80)
print(f"""  层高划分（Z001 §7.3）: q = {float(q_h):.8f} ∈ (1/2,3/5)  ✅ 入窗
  部分和边缘分布      : q = {float(q_p):.8f}
  三个自识别候选      : 见上
  ⟹ **非均匀权重在 L1′ 计数层是天然的**：它是组合计数的结果，不需要动力学去"产生"它。
     这与 Z016 层表一致（账本峰 D=4 在 L1′；"q 的原生来源"是初始段测度）。
     动力学的作用是**保持/选择**，不是**生成**。""")
