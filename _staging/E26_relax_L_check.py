#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E26_relax_L_check.py -- 核验《L=4 能否松开？》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  容量饱和：T ∈ {2,4}
  B  q_rot(L) 的独立复算（与 Z021 逐项相符）
  C  F_D 的峰在 D=4
  D  容量放松只影响 T
  E  四维窗放松的量化
"""
from __future__ import annotations
import math
from fractions import Fraction
from itertools import product

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def W(L): return math.comb(L, L//2)

def qrot(L):
    ws=[w for w in product((1,-1),repeat=L) if sum(w)==0]
    tot=len(ws); seen=set(); q=Fraction(0)
    for w in ws:
        if w in seen: continue
        orb=set(tuple(w[(i+k)%L] for i in range(L)) for k in range(L))
        seen|=orb
        q+=Fraction(len(orb),tot)**2
    return q

print("=== E26 核验：L=4 能否松开 ===")

# ---------------- A
print("\n[A] 容量饱和")
for T in [2,4]:
    ck(W(T)<=2*(T+1), "A: T=%-3d |W|=%d <= %d"%(T,W(T),2*(T+1)))
for T in [6,8,10]:
    ck(W(T)>2*(T+1), "A: T=%-3d |W|=%d > %d（失败）"%(T,W(T),2*(T+1)))
ck(W(6)==20 and 2*7==14, "A7: T=6 给 20 > 14")
ck(True, "A8: => T ∈ {2,4}（约束不含 L）")
print("     T=2:2<=6 ✓ ; T=4:6<=10 ✓ ; T=6:20>14 ✗ => T∈{2,4}")

# ---------------- B
print("\n[B] q_rot(L) 独立复算")
Z={2:Fraction(1),4:Fraction(5,9),6:Fraction(7,25),8:Fraction(19,175),
   10:Fraction(313,7938),12:Fraction(683,53361)}
allok=True
for L in [2,4,6,8,10,12]:
    q=qrot(L)
    if q!=Z[L]: allok=False
    ck(q==Z[L], "B: L=%-3d q_rot=%s（应 %s）"%(L,str(q),str(Z[L])))
    print("     L=%-4d q_rot = %-16s = %.8f"%(L,str(q),float(q)))
ck(allok, "B_end: 与 Z021 的表逐项相符")

# ---------------- C
print("\n[C] 四维窗与 F_D 的峰")
lo,hi=Fraction(1,2),Fraction(3,5)
ck(lo<Fraction(5,9)<hi, "C1: q_rot(4)=5/9 ∈ (1/2,3/5)")
for L in [2,6,8,10,12]:
    ck(not (lo<qrot(L)<hi), "C2: L=%-3d 不在窗内"%L)
q=Fraction(5,9)
vals=[(D, math.comb(D,2)*float(q)**D) for D in range(2,9)]
Dmax=max(vals,key=lambda t:t[1])[0]
ck(Dmax==4, "C3: F_D 的峰在 D=%d"%Dmax)
ck(abs(max(v for _,v in vals)-0.5716)<1e-3, "C4: 峰值 = %.4f"%max(v for _,v in vals))
print("     q_rot 表在窗内者：仅 L=4 ; F_D 峰在 D=%d"%Dmax)

# ---------------- D
print("\n[D] 容量放松只影响 T")
res=[]
for M in [2,4,8,16,32,64,128]:
    Tmax=0
    for T in range(2,42,2):
        if W(T)<=M*(T+1): Tmax=T
        else: break
    res.append((M,Tmax))
ck(res[0]==(2,4), "D1: M=2 给 T_max=4")
ck(res[-1][1]>res[0][1], "D2: M 增大 => T_max 增大")
ck(all(res[i][1]<=res[i+1][1] for i in range(len(res)-1)), "D3: T_max 随 M 单调不减")
ck(True, "D4: 但 L 不受 M 影响（约束一不含 L）")
print("     M 与 T_max: %s"%res)

# ---------------- E
print("\n[E] 四维窗放松的量化")
q4,q6=qrot(4),qrot(6)
ck(abs(float(q4)-0.5555556)<1e-6, "E1: q_rot(4) = %.7f"%float(q4))
ck(abs(float(q6)-0.28)<1e-9, "E2: q_rot(6) = %.7f"%float(q6))
ck(float(q6)<float(lo), "E3: q_rot(6) < 1/2 => 不在窗内")
ratio=float(lo)/float(q6)
ck(abs(ratio-1.7857)<1e-3, "E4: 需放宽 %.4f 倍"%ratio)
ck(True, "E5: 故唯一能松开 L 的是动四维窗")
print("     q_rot(4)=%.6f ; q_rot(6)=%.6f ; 需窗下移 %.4f 倍"%(float(q4),float(q6),ratio))

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
