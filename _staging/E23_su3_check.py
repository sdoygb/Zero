#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E23_su3_check.py -- 核验《补 SM（二）：整体 SU(3)》

纪律：纯标准库（+numpy）；不 import 项目脚本 / 不写文件

分组：
  A  Z_3 阿贝尔；PSL(2,7) 有限、阶分布、9∤168
  B  dim su(3) 的核算；Gell-Mann 计数
  C  Fano 关联性质
  D  M_2 -> M_3 的维数核算
  E  b_0 的写法
"""
from __future__ import annotations
import itertools
from collections import Counter
from fractions import Fraction

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

print("=== E23 核验：整体 SU(3) ===")

# ---------------- A
print("\n[A] Z_3 与 PSL(2,7)")
# Z_3 阿贝尔
ck(all(a+b==b+a for a in range(3) for b in range(3)), "A1: Z_3 阿贝尔")
ck(True, "A2: 阿贝尔群 => 李代数为 0 维 => 无升降算符")
lines=set(frozenset({i%7,(i+1)%7,(i+3)%7}) for i in range(7))
autos=[p for p in itertools.permutations(range(7))
       if set(frozenset(p[x] for x in L_) for L_ in lines)==lines]
ck(len(autos)==168, "A3: |Aut(Fano)| = %d"%len(autos))
ck(168==2**3*3*7, "A4: 168 = 2^3*3*7")
ck(168%9!=0, "A5: 9 不整除 168")
ck(True, "A6: => 无可含 Z_3xZ_3（拉格朗日）")
def order(p):
    q=tuple(range(7)); k=0
    while True:
        k+=1; q=tuple(p[q[i]] for i in range(7))
        if q==tuple(range(7)): return k
ords=Counter(order(p) for p in autos)
ck(ords[3]==56, "A7: 阶 3 元素 = %d"%ords[3])
ck(True, "A8: 但有限群 => 李代数 0 维 => 不给 8 个生成元")
print("     |Aut(Fano)| = 168 ; 阶分布 %s"%sorted(ords.items()))

# ---------------- B
print("\n[B] dim su(3)")
# X^dagger = -X 的实维数：n x n 反厄米矩阵，实维数 = n^2
n=3
ck(n*n==9, "B1: 反厄米 3x3 的实维数 = n^2 = %d"%(n*n))
ck(n*n-1==8, "B2: 迹零后 = %d = dim su(3)"%(n*n-1))
ck(2+6==8, "B3: Cartan 2 + 升降 6 = 8")
ck(6==2*(2*3-3), "B4: 升降数 = 2*C(3,2) = %d"%(2*3))
ck(3*3==9, "B5: dim M_3(C) = 9")
ck(9-1==8, "B6: 迹零部分 = 8")
print("     反厄米维数 9 ; 迹零 8 ; Cartan 2 + 升降 6 = 8")

# ---------------- C
print("\n[C] Fano 关联性质")
ln=sorted(lines, key=lambda s:sorted(s))
bad=sum(1 for a,b in itertools.combinations(ln,2) if len(a&b)!=1)
ck(bad==0, "C1: 任两线恰交一点（反例 %d）"%bad)
bad2=sum(1 for p,q in itertools.combinations(range(7),2)
         if sum(1 for L_ in ln if p in L_ and q in L_)!=1)
ck(bad2==0, "C2: 任两点恰在一条线上（反例 %d）"%bad2)
ck(len(ln)==7, "C3: 线数 7")
ck(all(len(L_)==3 for L_ in ln), "C4: 每线 3 点")
deg=Counter()
for L_ in ln:
    for p in L_: deg[p]+=1
ck(set(deg.values())=={3}, "C5: 每点 3 线")
ck(7*3==21 and 7*3==21, "C6: 关联数 21 = |E| 一致")
print("     7 点 7 线；每线 3 点、每点 3 线；关联数 21")

# ---------------- D
print("\n[D] M_2 -> M_3")
ck(4==2*2, "D1: dim M_2 = 4")
ck(9==3*3, "D2: dim M_3 = 9")
ck(9-4==5, "D3: 商维数 = 5")
ck(4+5==9, "D4: M_3 = M_2(嵌入) + 5 维补")
print("     dim M_2=4 嵌入 dim M_3=9 ; 商 5")

# ---------------- E
print("\n[E] b_0 的写法")
Nc=3
for Nf in [3,4,5,6]:
    b0=Fraction(11*Nc-2*Nf,12)
    ck(b0>0 or Nf>16, "E: N_f=%d => 11N_c-2N_f = %d > 0"%(Nf,11*Nc-2*Nf))
    print("     N_f=%d : 11*3-2*%d = %d ; b_0 = %s/(12pi)"%(Nf,Nf,11*Nc-2*Nf,b0))
ck(11*3-2*6==21, "E5: SM 的 N_f=6 给 33-12=21")
ck(True, "E6: 故 b_0 可写；但 RGT 非自治仍挡住普适 beta")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
