# -*- coding: utf-8 -*-
"""Path L：决定性检验 —— 分区（层高）在"忠实动力学"下是否被保持？

前十一轮的正面对撞给出一个二难。本文件把它变成可判定命题：
  L0 的状态空间 = 平衡词集 W_n（大小 12870），计数测度 μ = 均匀。
  Z001 的划分 π = 层高（max|部分和|），块权 ω_b = |B_b|/|W_n|。
  q = Σω_b² = 0.55777859 ∈ D4 窗  —— 这是**静态**事实（本轮已复现）。

问题：**任何**忠实于"全分支/计数"的动力学（即 W_n 上的双射，保持 μ）
      是否保持 π？若不保持，块权随时间漂移，q 就守不住。

本文件用三种 W_n 上的双射检验层高标签的保持性：
  (D1) 相邻位置交换（若仍平衡）
  (D2) 交换一个 +1 与一个 −1（补偿移动的置换实现）
  (D3) 循环移位
"""
from __future__ import annotations
from itertools import combinations
from collections import Counter
import numpy as np

def hlabel(w):
    s=0;m=0
    for x in w: s+=x; m=max(m,abs(s))
    return m

def words(n):
    for pos in combinations(range(n),n//2):
        w=[-1]*n
        for i in pos: w[i]=1
        yield tuple(w)

n=16
Ws=list(words(n))
c0=Counter(hlabel(w) for w in Ws)
b0=[c0[k] for k in sorted(c0)]
print(f"|W_{n}| = {len(Ws)}")
print(f"初始块大小 b_h = {b0}")
from fractions import Fraction as F
def q_of(b): 
    t=sum(b); return float(sum(F(x,t)**2 for x in b))
q0=q_of(b0)
print(f"初始 q(不合并) = {q0:.6f}   （合并 h≤3 后 = 0.55777859）")
print()

rng=np.random.RandomState(11)
def apply_and_measure(name, T, steps=20000):
    """T: W_n 上的双射（元组→元组）。度量块分布如何漂移。"""
    w=Ws[rng.randint(len(Ws))]
    counts=Counter(); 
    qs=[]
    for t in range(steps):
        w=T(w, rng)
        counts[hlabel(w)]+=1
        if (t+1)%2000==0:
            b=[counts[k] for k in sorted(counts)]
            qs.append(q_of(b))
    return qs

def T_adjswap(w, rng):
    """随机相邻交换（保持平衡当且仅当 w_i ≠ w_{i+1}）"""
    w=list(w)
    for _ in range(20):
        i=rng.randint(n-1)
        if w[i]!=w[i+1]:
            w[i],w[i+1]=w[i+1],w[i]; return tuple(w)
    return tuple(w)

def T_compensate(w, rng):
    """交换一个 +1 与一个 −1 的位置（补偿移动的置换实现）：仍是平衡词"""
    w=list(w)
    P=[i for i,x in enumerate(w) if x==1]; M=[i for i,x in enumerate(w) if x==-1]
    i=P[rng.randint(len(P))]; j=M[rng.randint(len(M))]
    w[i],w[j]=w[j],w[i]; return tuple(w)

def T_rotate(w, rng):
    k=rng.randint(n)
    return w[k:]+w[:k]

print(f"{'双射':22s} {'q 轨迹（每 2000 步）':<60s} {'漂移'}")
for nm,T in [("相邻交换",T_adjswap),("补偿置换",T_compensate),("循环移位",T_rotate)]:
    qs=apply_and_measure(nm,T,steps=20000)
    print(f"{nm:22s} " + " ".join(f"{x:.4f}" for x in qs))
print()
print("判读：若 q 轨迹**漂移离开** 0.5578，则块权不守恒 ⟹ 该动力学破坏 Z001 的分区结构。")
print()
# 更精确：单步下块标签的转移是否闭块
print("="*88)
print("决定性检验：层高标签在单步双射下的**平均逗留**")
print("="*88)
for nm,T in [("相邻交换",T_adjswap),("补偿置换",T_compensate),("循环移位",T_rotate)]:
    same=0; tot=0
    for _ in range(20000):
        w=Ws[rng.randint(len(Ws))]
        w2=T(w,rng)
        same += (hlabel(w)==hlabel(w2)); tot+=1
    print(f"  {nm:12s}: P(标签不变) = {same/tot:.4f}   (1.0 ⟹ 完全保持)")
print()
print("★ 结论：若 P(标签不变) < 1，则忠实动力学**不保持**分区 ⟹ 块权是瞬态量，")
print("   不能作为动力学的守恒读出。")
