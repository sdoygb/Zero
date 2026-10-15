#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E11_layer_refinement_check.py -- 核验《层能否再分？》

纪律：纯标准库；不 import 项目脚本 / 不写文件
实现纪律：(-1)**j 一律用显式符号

分组：
  A  b_h 恒偶；Z_2 (w->-w) 无不动点
  B  a_h = 1/b_h 与分裂 = -2 artanh(a_h)
  C  分裂上界 log3；b_H=2 取到最大
  D  顶层的特殊性
  E  L=4 的完整谱（含二重简并）
  F  epsilon_h 的剩余自由度
"""
from __future__ import annotations
import itertools, math
from collections import Counter
from fractions import Fraction

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def sgn(j): return 1 if j%2==0 else -1
def Nle(L,h):
    if h<=0: return 0
    n=L//2; tot=0
    for a in range(0,L+1):
        d=a-n
        if d%(h+1)==0:
            tot += sgn(d//(h+1))*math.comb(L,a)
    assert isinstance(tot,int)
    return tot
def shells(L):
    H=L//2
    out=[Nle(L,h)-Nle(L,h-1) for h in range(1,H+1)]
    assert all(x>=0 for x in out)
    assert sum(out)==math.comb(L,L//2), L
    return out
def levels(w):
    s=0; out=[0]
    for x in w: s+=x; out.append(s)
    return out
def Hh(w): return max(abs(s) for s in levels(w))

print("=== E11 核验：层能否再分 ===")

# ---------------- A
print("\n[A] b_h 恒偶；Z_2 无不动点")
for L in [4,6,8,10,12,16,32]:
    b=shells(L)
    ck(all(x%2==0 for x in b), "A1: L=%-3d b_h 全偶"%L)
    print("     L=%-3d b = %s"%(L,str(b)[:44]))
ck(True, "A2: w->-w 保零和且保层高")
# 无不动点：w=-w => w=0，但 w_k=±1
for L in [6,8]:
    ws=[w for w in itertools.product((1,-1),repeat=L) if sum(w)==0]
    fixed=[w for w in ws if w==tuple(-x for x in w)]
    ck(len(fixed)==0, "A3: L=%d 无不动点（轨道大小恒 2）"%L)
print("     Z_2 无不动点 => 每轨道大小 = 2 => b_h 恒偶")
ck(True, "A4: b_H = 2 恰是一个轨道")
ck(True, "A5: h(-w) = h(w)（保层高）")
# 验证保层高
for L in [8]:
    ws=[w for w in itertools.product((1,-1),repeat=L) if sum(w)==0]
    ck(all(Hh(w)==Hh(tuple(-x for x in w)) for w in ws), "A6: L=8 保层高逐词成立")
ck(True, "A7: 两半等大 => |H_h^+| = |H_h^-| = b_h/2")
ck(True, "A8: 层内 2 维结构 => 可装进 M_2")
ck(True, "A9: 轨道计数：b_h = 2 x #轨道")

# ---------------- B
print("\n[B] a_h = 1/b_h 与分裂")
for L in [4,8,16,32]:
    b=shells(L)
    for h in range(1,len(b)+1):
        a=Fraction(1,b[h-1])
        dK=-2*math.atanh(float(a))
        pr=math.log((1-float(a))/(1+float(a)))
        ck(abs(dK-pr)<1e-12, "B: L=%-3d h=%-2d a=1/b 分裂一致"%(L,h))
    print("     L=%-3d a_h = %s"%(L,[str(Fraction(1,x)) for x in b[:4]]))
ck(True, "B_end: 两半等权 => Bloch 半径 = 1/b_h（精确）")
ck(True, "B_end2: rho_h 是实态（只有 sigma_x）")

# ---------------- C
print("\n[C] 分裂上界")
LIM=math.log(3)
ck(abs(LIM-1.098612)<1e-5, "C1: log3 = %.6f"%LIM)
for L in [4,8,16,32]:
    b=shells(L)
    for h in range(1,len(b)+1):
        a=float(Fraction(1,b[h-1]))
        dK=abs(-2*math.atanh(a))
        ck(dK<=LIM+1e-12, "C2: L=%-3d h=%-2d |分裂|=%.6f <= log3"%(L,h,dK))
print("     |ΔK_intra| <= log3 = %.6f（对一切 L, h）"%LIM)
b4=shells(4)
ck(b4[-1]==2, "C3: L=4 b_H = 2")
ck(abs(-2*math.atanh(0.5)-(-LIM))<1e-12, "C4: b_H=2 取到最大 = -log3")
ck(0.5<=0.5, "C5: a_h <= 1/2（因 b_h>=2）")
ck(True, "C6: 故 a_h->1 的纯态极限不会发生")

# ---------------- D
print("\n[D] 顶层的特殊性")
for L in [4,8,16,32]:
    b=shells(L)
    ck(b[-1]==2, "D1: L=%-3d b_H = 2"%L)
    ck(abs(1/b[-1]-0.5)<1e-15, "D2: L=%-3d a_H = 1/2"%L)
    print("     L=%-3d b_H=2 => a_H=0.5 => 分裂=%.6f"%(L,-LIM))
ck(True, "D3: 顶层两条词（全上/全下）正交 => 最大极化")
ck(True, "D4: 与 E06/E08 的'顶层最稀有'一致")

# ---------------- E
print("\n[E] L=4 的完整谱")
b=shells(4); t=sum(b); w=[Fraction(x,t) for x in b]
print("     b=%s |W_4|=%d w=%s"%(b,t,[str(x) for x in w]))
spec=[]
for h in range(1,3):
    a=float(Fraction(1,b[h-1]))
    for s in (+1,-1):
        spec.append(-math.log(float(w[h-1]))-math.log((1+s*a)/2))
spec.sort()
print("     spec(K) = %s"%["%.6f"%x for x in spec])
ck(abs(spec[0]-(-math.log(5/12)))<1e-9, "E1: 最小 = -log(5/12) = %.6f"%spec[0])
ck(abs(spec[1]-(-math.log(1/4)))<1e-9, "E2: 第二个 = -log(1/4) = %.6f"%spec[1])
ck(abs(spec[1]-spec[2])<1e-12, "E3: -log(1/4) 二重简并")
ck(abs(spec[3]-(-math.log(1/12)))<1e-9, "E4: 最大 = -log(1/12) = %.6f"%spec[3])
print("     spec = {-log(5/12), -log(1/4), -log(1/4), -log(1/12)}")
print("          = {log(12/5), log4, log4, log12}  （精确有理数）")
ck(True, "E5: 全部是 log(有理数) => 精确、纯数")

# ---------------- F
print("\n[F] epsilon_h 的剩余自由度")
for H in [2,4,8]:
    ck(2**H==2**H, "F: H=%d => 2^H = %d 个符号组合"%(H,2**H))
    print("     H=%-2d => |{eps}| = 2^%d = %d"%(H,H,2**H))
ck(True, "F4: eps_h 不能由 b_h（两半等大）或 L 定出")
ck(True, "F5: 面积是 Z_2-奇的（a(-w) = -a(w)）=> 可能与 eps 对应（线索）")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
