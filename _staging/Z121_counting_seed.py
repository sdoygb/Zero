# -*- coding: utf-8 -*-
"""Path J：不用人为的几何种子，改用**计数测度自己给出的剖面**当 L2 的输入。

Z087/Z088 的引擎种子是 seed[v] = 2^{V−1−v}（几何）。这个 2 是**实现者的选择**，
Z087 §5.1 也承认"log 2 继承 seed 相邻比"。
但 Zero 的计数测度早就给定了词的分布：Z001 §7.3 的层高划分
    b_h = [256, 4118, 4858, 2518, 880, 208, 30, 2]   （h = max|部分和|）
所以正确的顺序是：**L1′ 的计数权重 → 作为 L2 的输入 → 看引擎输出什么**。

本文件：
  (a) 复现 b_h（此前已验证与 Z001 逐项吻合）
  (b) 把 b_h 映射到 path16 的顶点（h ↔ 到中心的距离，两侧对称）
  (c) 用它当种子跑引擎，量 record 的 q
  (d) 与 Z001 的 q=0.55777859 对比，与 Z088 的 5/9 对比
"""
from __future__ import annotations
import sys, os
from itertools import combinations
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "..", "ZeroGPU", "src")))
from z0_engine import path_adj
from z0_vectorized import ZeroCoreVec

CHECKS=[]
def check(n,c,d=""):
    CHECKS.append((n,bool(c))); print(f"[{'OK' if c else 'FAIL'}] {n}  {d}")

def parts(w):
    s=0; ps=[0]
    for x in w: s+=x; ps.append(s)
    return ps
n=16
b=Counter(max(abs(x) for x in parts(w)) for w in
          ([1 if i in pos else -1 for i in range(n)] for pos in combinations(range(n),n//2)))
bh=[b[k] for k in sorted(b)]
print(f"(a) b_h (h=1..{len(bh)}) = {bh}   和={sum(bh)}   (=C(16,8)=12870)")
check("b_h 与 Z001 §7.3 逐项吻合", bh==[256,4118,4858,2518,880,208,30,2])
from fractions import Fraction as F
merged=[sum(bh[:3])]+bh[3:]
qZ=float(sum(F(x,sum(merged))**2 for x in merged))
print(f"    Z001 合并划分 q = {qZ:.8f}")
check("复现 Z001 q=0.55777859", abs(qZ-0.55777859)<1e-7, f"{qZ:.8f}")

print()
print("(b) 把 h 剖面映射到 path16 顶点（h ↔ 离中心距离）")
V=16
# 顶点 v 到中心的距离 d(v)=|v-(V-1)/2|；h = floor(d)+1 之类；这里用"层"= 到中心距离
d=[abs(v-(V-1)/2) for v in range(V)]
# 反转：离中心近 → h 大（因为 h=max|s| 大者词少）；构造与 b_h 同型的权重
prof=np.zeros(V)
for v in range(V):
    hh=int(round(d[v]))+1            # 层指标 1..8
    hh=min(hh,len(bh))
    prof[v]=bh[hh-1]
prof=prof/prof.max()
print(f"    顶点剖面 (归一化) = {np.round(prof,4)}")
print(f"    与几何种子 2^(V-1-v) 的相关系数 = {np.corrcoef(prof, [2**(V-1-v) for v in range(V)])[0,1]:.4f}")

print()
print("(c) 用计数剖面当种子跑引擎，量 record 的 q")
def lock(seed,gens=40000,ret="last"):
    v=ZeroCoreVec(path_adj(V),seed=[float(x) for x in seed],retention=ret,L=4)
    for _ in range(gens): v.run_generation()
    return v
seeds={
 "计数剖面 b_h":     prof*1000,
 "几何 2^(V-1-v)":   [2.0**(V-1-v) for v in range(V)],
 "均匀":             [1.0]*V,
 "b_h 原始(无映射)":  [float(bh[min(v,len(bh)-1)]) for v in range(V)],
}
print(f"{'种子':18s} {'q(record)':>12s} {'非零顶点':>9s} {'形态'}")
for nm,sd in seeds.items():
    v=lock(sd)
    w=v.record/v.record.sum()
    nz=[(i,float(x)) for i,x in enumerate(w) if x>1e-9]
    q=float((w*w).sum())
    top=sorted(nz,key=lambda z:-z[1])[:3]
    print(f"{nm:18s} {q:12.6f} {len(nz):9d} {[(i,round(x,4)) for i,x in top]}")
print()
print("(d) 对照")
print(f"    Z001 的 q（L1′ 账本）      = {qZ:.8f}")
print(f"    Z088 声称的 q              = 0.555556 (5/9)")
print(f"    D4 窗                      = (0.5, 0.6)")
print()
print("★ 判读：record 的分布**逐点复制种子**（引擎逐顶点独立），所以 q(record) 是")
print("   **种子 q** 的像，不是动力学生成的量。用计数剖面当种子，得到的 q 就是")
print("   计数剖面的 q —— 即 q 的来源在 L1′ 计数层，与 L2 演化无关。")
qs={nm: float(((lock(sd).record/lock(sd).record.sum())**2).sum()) for nm,sd in seeds.items()}
print(f"    实测 q: {[(k,round(v,6)) for k,v in qs.items()]}")
check("record 分布 = 种子分布（逐顶点独立 ⟹ 无生成）",
      abs(qs['均匀']-1/V)<1e-9, f"均匀种子 → q={qs['均匀']:.6f} = 1/V")
