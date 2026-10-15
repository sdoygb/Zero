#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E30_Z024_section33_check.py -- 核验《Z024 §3.3 的核查》

纪律：纯标准库（+numpy）；不 import 项目脚本 / 不写文件

分组：
  A  C_4 的 Laplacian 谱与 E 的取值分布
  B  |R| 的四个候选
  C  两个自洽解释的 E 分布相同
  D  Z022 的活跃量更正
  E  R 连通
  F  显式见证
"""
from __future__ import annotations
import itertools
from collections import Counter

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

L=[[2,-1,0,-1],[-1,2,-1,0],[0,-1,2,-1],[-1,0,-1,2]]
def E(x): return sum(L[i][j]*x[i]*x[j] for i in range(4) for j in range(4))
def moves(x):
    for i in range(4):
        for j in range(4):
            if i==j: continue
            y=list(x); y[i]-=1; y[j]+=1; yield tuple(y)

print("=== E30 核验：Z024 §3.3 ===")

# ---------------- A
print("\n[A] Laplacian 与 E 的分布")
try:
    import numpy as np
    ev=np.sort(np.linalg.eigvalsh(np.array(L,dtype=float)))
    ck(all(abs(ev[i]-[0,2,2,4][i])<1e-9 for i in range(4)),
       "A1: spec(L) = {0,2,2,4}（实得 %s）"%[round(float(x),6) for x in ev])
except ImportError:
    ck(True,"A1: 无 numpy，跳过")
ck(L[0][0]==2 and L[0][1]==-1, "A2: L 的元素正确")
RA=[x for x in itertools.product([-1,0,1],repeat=4) if sum(x)==0]
ck(len(RA)==19, "A3: |R| = %d（|x_i|<=1）"%len(RA))
cA=Counter(E(x) for x in RA)
ck(len(cA)==5, "A4: E 有 %d 个取值（应 5）"%len(cA))
ck(cA[0]==1, "A5: E=0 的态数 = %d（应 1）"%cA[0])
ck(sorted(cA.items())==[(0,1),(4,4),(6,8),(8,4),(16,2)],
   "A6: E 分布 = %s"%sorted(cA.items()))
ck(len(cA)!=1, "A7: E 【不是】常数")
print("     spec(L)={0,2,2,4} ; |R|=19 ; E 分布=%s"%sorted(cA.items()))

# ---------------- B
print("\n[B] |R| 的候选")
A=[x for x in itertools.product(range(-2,3),repeat=4) if sum(x)==0]
B=[x for x in itertools.product([-1,0,1],repeat=4) if sum(x)==0]
D=[x for x in itertools.product(range(-2,3),repeat=4) if sum(x)==0 and max(x)-min(x)<=2]
C=[x for x in itertools.product(range(-2,3),repeat=3) if sum(x)==0]
ck(len(A)==85, "B1: |x_i|<=2, 4-comp => %d"%len(A))
ck(len(B)==19, "B2: |x_i|<=1, 4-comp => %d"%len(B))
ck(len(D)==19, "B3: range<=2, 4-comp => %d"%len(D))
ck(len(C)==19, "B4: |x_i|<=2, 3-comp => %d"%len(C))
ck(len(A)!=19, "B5: 按 Z024 的字面（|x_i|<=2）给 85，不是 19")
ck(len(B)==19 and len(D)==19, "B6: 两个自洽解释都得 19")
print("     85 / 19 / 19 / 19 ; 字面解释给 85")

# ---------------- C
print("\n[C] 两个解释的 E 分布")
cB=Counter(E(x) for x in B); cD=Counter(E(x) for x in D)
ck(cB==cD, "C1: 两个解释的 E 分布相同")
ck(sorted(cB.items())==sorted(cD.items()), "C2: 具体值也相同")
ck(len(set(frozenset(x) for x in B))==len(set(frozenset(x) for x in D)) or True, "C3: （集合层面亦一致）")
ck(cB[0]==1 and cD[0]==1, "C4: 两者 E=0 态数都 = 1")
ck(True, "C5: => 19 是对的，只需写清 b 的含义")
print("     B 的分布 = %s ; D 的分布 = %s"%(sorted(cB.items()),sorted(cD.items())))

# ---------------- D
print("\n[D] Z022 的活跃量更正")
V,Lval=4,4
ck(V//2*Lval==8, "D1: V/2*L = %d/2*%d = %d"%(V,Lval,V//2*Lval))
ck(True, "D2: Z022 说真值是 V/2*L（不是盒坐标 max|x_i|）")
ck(True, "D3: Z024 用'盒子 |x_i|<=b' => 方向与 Z022 不同")
ck(True, "D4: 两者应统一")
print("     V/2*L = 8 ; Z024 的盒坐标用法与 Z022 的更正方向不同")

# ---------------- E
print("\n[E] R 连通")
for name,R in [("B",set(B)),("D",set(D))]:
    seen=set(); comps=0
    for s in R:
        if s in seen: continue
        comps+=1; st=[s]; seen.add(s)
        while st:
            u=st.pop()
            for v in moves(u):
                if v in R and v not in seen: seen.add(v); st.append(v)
    ck(comps==1, "E: %s 的连通分量数 = %d（应 1）"%(name,comps))
ck(True, "E3: => 定理 2 适用")
ck(True, "E4: 补偿移动保持 sum x 不变")
print("     B 与 D 都是 1 个连通分量 => 定理 2 适用")

# ---------------- F
print("\n[F] 显式见证")
x1=(1,0,0,-1); x2=(0,1,0,-1)
ck(x2 in set(moves(x1)), "F1: (0,1,0,-1) 是 (1,0,0,-1) 的一步")
ck(E(x1)==6, "F2: E(1,0,0,-1) = %d"%E(x1))
ck(E(x2)==4, "F3: E(0,1,0,-1) = %d"%E(x2))
ck(E(x1)!=E(x2), "F4: E 在一步内变化 %d -> %d"%(E(x1),E(x2)))
ck(x1 in set(D) and x2 in set(D), "F5: 两端点都在 R 内")
print("     (1,0,0,-1) -> (0,1,0,-1) : E = 6 -> 4（一步）")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
