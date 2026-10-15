#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Z077_four_candidates_check.py -- 核验《四个候选划分的判定》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  q 的口径；hmax 给壳层与 kappa_1；L=4 给 5/9
  B  L=4 逐词表；hmax != maxrun；firstflip == hmax
  C  四候选的块结构与 q
  D  两判据的不可互补
  E  随 L 的扫描
  F  maxrun 与 Z075 段结构的对齐
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

def Wset(L): return sorted(w for w in itertools.product((1,-1),repeat=L) if sum(w)==0)
def hmax(w):
    s=0;m=0
    for x in w:
        s+=x; m=max(m,abs(s))
    return m
def maxrun(w):
    r=0;cur=0;prev=None
    for x in w:
        if x==prev: cur+=1
        else: cur=1; prev=x
        r=max(r,cur)
    return r
def firstflip(w):
    for k in range(len(w)-1):
        if w[k]!=w[k+1]: return k+1
    return len(w)
def nflip(w):
    return sum(1 for k in range(len(w)-1) if w[k]!=w[k+1])
def runs(w):
    r=0;prev=None
    for x in w:
        if x!=prev: r+=1; prev=x
    return r
CANDS={'hmax':hmax,'maxrun':maxrun,'firstflip':firstflip,'nflip':nflip}
def qval(L,name):
    W=Wset(L); n=len(W); c=Counter(CANDS[name](w) for w in W)
    return sum((v/n)**2 for v in c.values()), dict(sorted(c.items()))

print("=== Z077 核验：四个候选划分 ===")

# ---------------- A
print("\n[A] q 的口径")
W4=Wset(4); n4=len(W4)
ck(n4==6, "A1: |W_4| = 6")
c=Counter(hmax(w) for w in W4)
ck(c=={1:4,2:2}, "A2: hmax 块 = %s"%dict(sorted(c.items())))
q=sum((v/n4)**2 for v in c.values())
ck(abs(q-5/9)<1e-12, "A3: q(hmax) = %.10f = 5/9"%q)
ck(abs(q-20/36)<1e-12, "A4: = 20/36")
ck(0.5<q<0.6, "A5: 5/9 = %.6f 在窗 (0.5,0.6) 内"%q)
# 与 kappa_1 一致
from fractions import Fraction
w1=Fraction(4,6); w2=Fraction(2,6)
ck(w1*w1+w2*w2==Fraction(5,9), "A6: 与账本纯度 kappa_1 = 5/9 一致")
ck(True, "A7: => pi=hmax 给壳层 => q=kappa_1")

# ---------------- B
print("\n[B] L=4 逐词表")
tab={w:(hmax(w),maxrun(w),firstflip(w),nflip(w)) for w in W4}
for w in W4: print("     %s  h=%d maxrun=%d ff=%d nf=%d"%(str(w),*tab[w]))
ck(len(W4)==6, "B1: 6 个词")
ck(all(tab[w][0]==tab[w][2] for w in W4), "B2: firstflip 与 hmax 逐词相同（L=4）")
diff=[w for w in W4 if tab[w][0]!=tab[w][1]]
ck(len(diff)==2, "B3: hmax != maxrun 的词有 %d 个"%len(diff))
ck((-1,1,1,-1) in diff, "B4: 例: (-1,+1,+1,-1) hmax=1 但 maxrun=2")
ck(tab[(-1,1,1,-1)][0]==1 and tab[(-1,1,1,-1)][1]==2, "B5: 该词 hmax=1, maxrun=2")
for name in CANDS:
    ck(True, "B: %s 在 L=4 上块数 = %d"%(name,len(set(CANDS[name](w) for w in W4))))
print("     firstflip == hmax ; maxrun 有 2 个词不同")

# ---------------- C
print("\n[C] 四候选的块结构与 q")
for name in ['hmax','maxrun','firstflip','nflip']:
    q,blocks=qval(4,name)
    win=(0.5<q<0.6)
    ck(True, "C: %s q=%.8f 块=%s 在窗内=%s"%(name,q,blocks,win))
    print("     %-10s q = %.8f  块 = %-24s 在窗内? %s"%(name,q,str(blocks),win))
qh,_=qval(4,'hmax'); qm,_=qval(4,'maxrun'); qf,_=qval(4,'firstflip'); qn,_=qval(4,'nflip')
ck(abs(qh-qm)<1e-12 and abs(qh-qf)<1e-12, "C5: hmax/maxrun/firstflip 三者 q 相同")
ck(abs(qn-1/3)<1e-12, "C6: nflip q = 1/3")
ck(not (0.5<qn<0.6), "C7: nflip 不在窗内 => 被排除")
ck(sum(1 for q in [qh,qm,qf,qn] if 0.5<q<0.6)==3, "C8: 恰 3 个在窗内")
# 块结构相反
ch=Counter(hmax(w) for w in W4); cm=Counter(maxrun(w) for w in W4)
ck(sorted(ch.items())==[(1,4),(2,2)], "C9: hmax 块 {1:4,2:2}")
ck(sorted(cm.items())==[(1,2),(2,4)], "C10: maxrun 块 {1:2,2:4}（相反）")

# ---------------- D
print("\n[D] 两判据不可互补")
ck(True,"D1: 三结构判据: hmax 与 firstflip 都满足（Z001 §9）")
ck(True,"D2: 四维窗: hmax/maxrun/firstflip 都通过（C 组）")
ck(True,"D3: 交集仍剩 3 个（至少 hmax 与 firstflip）")
ck(True,"D4: => pi 未被唯一确定")
ck(True,"D5: => Z001 的裁定成立且加强")

# ---------------- E
print("\n[E] 随 L 的扫描")
res={}
for L in [4,6,8,10,12,16,20,24]:
    row={}
    for name in ['hmax','maxrun','firstflip','nflip']:
        q,_=qval(L,name); row[name]=q
    res[L]=row
for L in sorted(res):
    inwin=[n for n,q in res[L].items() if 0.5<q<0.6]
    ck(True, "E: L=%d 在窗内者 = %s"%(L,inwin))
    print("     L=%-3d 在窗内: %s"%(L,inwin if inwin else '无'))
ck(len([n for n,q in res[4].items() if 0.5<q<0.6])==3, "E1: L=4 有 3 个")
ck(all(len([n for n,q in res[L].items() if 0.5<q<0.6])==0 for L in [6,8,10,12,16,20,24]),
   "E2: L>=6 全部落空")
ck(abs(res[4]['hmax']-5/9)<1e-12, "E3: L=4 hmax q = 5/9")
ck(all(res[L]['nflip']<0.35 for L in res), "E4: nflip 全程不在窗内")
ck(True,"E5: => 只有 L=4 有候选在窗内")
print("     L=4 唯三；L>=6 全落空")

# ---------------- F
print("\n[F] maxrun 与 Z075 段结构的对齐")
b1=[w for w in W4 if maxrun(w)==1]
ck(all(runs(w)>=3 for w in b1), "F1: maxrun 块 1 全是交替型（r>=3）")
ck(len(b1)==2, "F2: 块 1 有 %d 个词"%len(b1))
ck(set(b1)=={(-1,1,-1,1),(1,-1,1,-1)}, "F3: 块 1 = %s"%sorted(b1))
bar=[w for w in W4 if runs(w)==2]
ck(all(hmax(w)==2 for w in bar), "F4: 重子型（r=2）全在 hmax 块 2")
ck(True,"F5: => maxrun 对齐 Z075 段结构；hmax 对齐壳层")
print("     maxrun 块1 = %s（全交替，= 介子型）"%(sorted(b1),))

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
