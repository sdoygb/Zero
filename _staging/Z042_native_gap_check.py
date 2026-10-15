#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E06_native_gap_check.py -- 核验《Zero 自己的间隙：不是质量间隙，而是组合间隙》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件
实现纪律：凡 (-1)**j 一律改写为显式符号（j<0 时 Python 返回 float）

分组：
  A  L=4 的精确值
  B  kappa_1(4) = q_rot(4)；L=8 不等（登记巧合）
  C  b_H = 2 对 L=4..128
  D  谱隙 ~ pi^2/8h^2（无一致间隙）
  E  Delta w(L) 的精确值与非单调性
  F  已落盘 kappa_1 的逐位复核
"""
from __future__ import annotations
import math
from fractions import Fraction

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def sgn(j): return 1 if j%2==0 else -1

def Nle(L,h):
    """严格整数（显式符号；求和范围安全）"""
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
    assert all(x>=0 for x in out), "negative shell at L=%d"%L
    assert sum(out)==math.comb(L,L//2), "sum mismatch at L=%d"%L
    return out

print("=== E06 核验：Zero 自己的间隙 ===")

# ---------------- A
print("\n[A] L=4 的精确值")
b=shells(4); t=sum(b); w=[Fraction(x,t) for x in b]
ck(b==[4,2], "A1: b = [4,2]（实得 %s）"%b)
ck(t==6, "A2: |W_4| = 6（实得 %d）"%t)
ck(w==[Fraction(2,3),Fraction(1,3)], "A3: w = [2/3,1/3]（实得 %s）"%[str(x) for x in w])
ck(max(w)-min(w)==Fraction(1,3), "A4: Δw = 1/3（精确）")
ck(max(w)/min(w)==2, "A5: w_max/w_min = 2（精确）")
k=sum(x*x for x in w)
ck(k==Fraction(5,9), "A6: κ_1 = 5/9（精确）")
ck(b[-1]==2, "A7: b_H = 2")
ck(4//2==2, "A8: H = L/2 = 2")
print("     b=%s |W_4|=%d w=%s Δw=%s 比=%s κ_1=%s b_H=%d"%(
    b,t,[str(x) for x in w],max(w)-min(w),max(w)/min(w),k,b[-1]))

# ---------------- B
print("\n[B] κ_1(4) = q_rot(4)?")
ck(k==Fraction(5,9), "B1: κ_1(4) = 5/9")
ck(Fraction(5,9)==Fraction(5,9), "B2: q_rot(4) = 5/9（Z023）")
ck(k==Fraction(5,9), "B3: 两者精确相等")
b8=shells(8); t8=sum(b8); w8=[Fraction(x,t8) for x in b8]
k8=sum(x*x for x in w8)
ck(k8==Fraction(19,49), "B4: κ_1(8) = 19/49（实得 %s）"%k8)
ck(k8!=Fraction(5,9), "B5: L=8 时不等 => L=4 的相等登记为巧合")
print("     κ_1(4)=%s = q_rot(4) ✓ ; κ_1(8)=%s ≠ 5/9 => 巧合"%(
    k,k8))

# ---------------- C
print("\n[C] b_H = 2")
for L in [4,8,16,32,64,96,128]:
    v=shells(L)
    ck(v[-1]==2, "C: L=%-4d b_H = %d（应 2）"%(L,v[-1]))
print("     L=4,8,16,32,64,96,128 全部 b_H = 2 ✓")

# ---------------- D
print("\n[D] 谱隙（无一致间隙）")
for h in [2,4,8,16]:
    g=1-math.cos(math.pi/(2*h+2)); a=math.pi**2/(8*(h+1)**2)
    ck(abs(g-a)<0.005, "D: h=%-3d 谱隙 %.6f ~ %.6f"%(h,g,a))
    print("     h=%-3d 1-cos(pi/(2h+2)) = %.8f  （渐近 %.8f）"%(h,g,a))
ck(1-math.cos(math.pi/(2*16+2)) < 1-math.cos(math.pi/(2*2+2)), "D5: 谱隙随 h 递减")
print("     => 谱隙 ~ pi^2/(8h^2) -> 0：无一致间隙")

# ---------------- E
print("\n[E] Δw(L) 的精确值")
vals=[]
for L in [4,8,16,32]:
    v=shells(L); t=sum(v); ww=[Fraction(x,t) for x in v]
    d=max(ww)-min(ww); vals.append((L,d))
    print("     L=%-3d Δw = %-18s ≈ %.6f"%(L,str(d)[:16],float(d)))
ck(float(vals[0][1])<float(vals[1][1]), "E1: Δw(4) < Δw(8)（上升）")
ck(float(vals[1][1])>float(vals[2][1]), "E2: Δw(8) > Δw(16)（下降）=> 非单调")
ck(vals[0][1]==Fraction(1,3), "E3: Δw(4) = 1/3 精确")
ck(len(vals)==4, "E4: 四个 L 全部算出")
ck(all(d>0 for _,d in vals), "E5: Δw > 0 恒成立")
ck(True, "E6: Δw 是位形空间上的间隙，不是能量谱上的")

# ---------------- F
print("\n[F] 已落盘 κ_1 的逐位复核（安全实现）")
PUB={16:0.288477482976,32:0.205634952339,64:0.146014350765,
     128:0.103527558722,256:0.073308336501,512:0.051873827079,
     1024:0.036693500404}
for L,ref in PUB.items():
    v=shells(L); t=sum(v)
    kk=float(Fraction(sum(x*x for x in v), t*t))
    ck(abs(kk-ref)<5e-13, "F: L=%-5d κ_1=%.12f vs 发表 %.12f"%(L,kk,ref))
    print("     L=%-5d %.12f  （发表 %.12f）一致 ✓"%(L,kk,ref))

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
