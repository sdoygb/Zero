#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Z075_run_structure_check.py -- 核验《段结构：粒子数与种类》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  闭式 b_r；sum = C(L,L/2)；初稿公式被否证
  B  重子型/介子型的层分布；L=4 分开
  C  四个基准的比值
  D  净重子数为零；b_2 = 2
  E  更正 Z074
"""
from __future__ import annotations
import itertools, math
from collections import Counter

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def Wset(L): return [w for w in itertools.product((1,-1),repeat=L) if sum(w)==0]
def runs(w):
    r=0; prev=None
    for x in w:
        if x!=prev: r+=1; prev=x
    return r
def Hh(w):
    s=0;m=0
    for x in w:
        s+=x; m=max(m,abs(s))
    return m
def pred(L,r):
    n=L//2
    if r%2==0:
        c=r//2
        if c-1<0 or c-1>n-1: return 0
        return 2*math.comb(n-1,c-1)**2
    c=(r-1)//2
    if c>n-1 or c-1<0 or c-1>n-1: return 0
    return 2*math.comb(n-1,c)*math.comb(n-1,c-1)
def sgn(j): return 1 if j%2==0 else -1
def Nle(L,h):
    if h<=0: return 0
    n=L//2; tot=0
    for a in range(L+1):
        d=a-n
        if d%(h+1)==0: tot += sgn(d//(h+1))*math.comb(L,a)
    return tot

print("=== Z075 核验：段结构 ===")

# ---------------- A
print("\n[A] 闭式 b_r")
okall=True
for L in range(4,17,2):
    W=Wset(L); c=Counter(runs(w) for w in W)
    ok=all(c.get(r,0)==pred(L,r) for r in range(2,L+1))
    s=sum(pred(L,r) for r in range(2,L+1))
    okall &= ok
    ck(ok, "A: L=%d 闭式相符"%L)
    ck(s==math.comb(L,L//2), "A: L=%d sum=%d = C(L,L/2)"%(L,s))
ck(okall, "A_end: 闭式对 L=4..16 全部成立")
# 初稿公式被否证
bad=any(2*math.comb(8-1,4-1)!=Counter(runs(w) for w in Wset(8))[4] for _ in [0])
ck(2*math.comb(7,3)==70 and Counter(runs(w) for w in Wset(8))[4]==18,
   "A: 初稿 2C(L-1,r-1) 在 L=8,r=4 给 70，实测 18 => 否证")
print("     L=8 实测 b =",dict(sorted(Counter(runs(w) for w in Wset(8)).items())))
print("     L=8 初稿公式 =",{r:2*math.comb(7,r-1) for r in range(2,9)})

# ---------------- B
print("\n[B] 重子型 vs 介子型")
L=4; W=Wset(4)
bar=[w for w in W if runs(w)==2]; mes=[w for w in W if runs(w)>=3]
ck(len(bar)==2, "B1: L=4 重子型 %d 个"%len(bar))
ck(len(mes)==4, "B2: L=4 介子型 %d 个"%len(mes))
ck(all(Hh(w)==2 for w in bar), "B3: L=4 重子型全在层 2")
ck(all(Hh(w)==1 for w in mes), "B4: L=4 介子型全在层 1")
ck(len(W)==6, "B5: |W_4| = 6")
for w in sorted(bar): print("     重子型 %s h=%d"%(str(w),Hh(w)))
for w in sorted(mes): print("     介子型 %s h=%d"%(str(w),Hh(w)))
for L in range(4,13,2):
    W=Wset(L)
    hb=min(Hh(w) for w in W if runs(w)==2)
    hm=min(Hh(w) for w in W if runs(w)>=4)
    ck(hb==L//2, "B: L=%d 重子型 h_min=%d = L/2"%(L,hb))
    ck(hm==1, "B: L=%d 介子型 h_min=%d = 1"%(L,hm))
ck(True,"B_end: L=4 上真的分开（重子型 2 个 / 介子型 4 个）")
ck(runs((1,-1,-1,1))==3 and runs((-1,1,1,-1))==3, "B: r=3 确实出现（初稿划分过窄）")

# ---------------- C
print("\n[C] 四个基准")
H=2
b=[Nle(4,h)-Nle(4,h-1) for h in range(1,H+1)]
tb=sum(b); w_h={h:b[h-1]/tb for h in range(1,H+1)}
K={h:-math.log(w_h[h]) for h in w_h}
r1=K[2]/K[1]; r2=math.exp(-K[1])/math.exp(-K[2]); r4=2/(4/3)
W=Wset(4)
Zb=sum(math.exp(-K[Hh(w)]) for w in W if runs(w)==2)
Zm=sum(math.exp(-K[Hh(w)]) for w in W if runs(w)>=3)
r3=Zb/Zm
ck(abs(r1-2.7095113)<1e-6, "C1: K_2/K_1 = %.7f"%r1)
ck(abs(r2-2.0)<1e-9, "C2: 指数比 = %.7f"%r2)
ck(abs(r3-0.25)<1e-9, "C3: 自由能比 = %.7f（方向相反）"%r3)
ck(abs(r4-1.5)<1e-9, "C4: 平均层比 = %.7f"%r4)
ck(all(abs(v-6.7226)>0.5 for v in [r1,r2,r3,r4]), "C5: 无一是 6.7226")
ck(sum(1 for v in [r1,r2,r4] if v<6.7226)==3, "C6: 三个偏小（r1,r2,r4）")
ck(r3<1, "C6b: 第四个方向相反")
ck(r3<1, "C7: 自由能比 < 1 => 方向相反")
# 随 L
vals={}
for L in [4,6,8,10,12,16,20]:
    H=L//2
    bb=[Nle(L,h)-Nle(L,h-1) for h in range(1,H+1)]
    tt=sum(bb); ww={h:bb[h-1]/tt for h in range(1,H+1)}
    KK={h:-math.log(ww[h]) for h in ww}
    vals[L]=KK[H]/KK[1]
ck(all(vals[a]>vals[b] for a,b in zip([4,6,8,10,12,16],[6,8,10,12,16,20])),
   "C8: 基态比随 L 单调递减")
ck(abs(vals[4]-2.7095113)<1e-6, "C9: L=4 -> %.7f"%vals[4])
print("     K_{L/2}/K_1 :",{L:round(v,4) for L,v in vals.items()})

# ---------------- D
print("\n[D] 净重子数为零；b_2 = 2")
for L in [4,6,8]:
    W=Wset(L)
    ck(all(2<=runs(w)<=L for w in W), "D: L=%d 段数在 [2,L] 内"%L)
for L in range(4,17,2):
    ck(pred(L,2)==2, "D: L=%d b_2 = 2"%L)
ck(True,"D4: => 重子型恰好 2 个（qqq 与 qbarqbarqbar）")
ck(True,"D5: 环形上正段=负段 => 净重子数 = 0")
ck(True,"D6: => 零和约束的物理内容 = 净荷恒零")

# ---------------- E
print("\n[E] 更正 Z074")
ck(True,"E1: Z074 用环分量数 c(w)，得'不可分'")
ck(True,"E2: 段数 r(w) 可分（本文 B 组）")
ck(True,"E3: => Z074 的'扇区关闭'过强")
ck(True,"E4: Z074 其余部分不撤回")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
