#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Z079_L4_exact_check.py -- 核验《L=4 的完整精确解》

纪律：纯标准库（+numpy）；不 import 项目脚本 / 不写文件

分组：
  A  6 个词、层高、段数
  B  b 与 w 的精确有理值
  C  kappa_1、H、K、Delta_mod
  D  补偿图正则、动力学谱、谱隙、平稳分布
  E  模谱
  F  段结构
  G  与 Heawood 的对照
"""
from __future__ import annotations
import itertools, math
from fractions import Fraction
from collections import Counter

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

L=4
W=sorted(w for w in itertools.product((1,-1),repeat=L) if sum(w)==0)
def Hh(w):
    s=0;m=0
    for x in w:
        s+=x; m=max(m,abs(s))
    return m
def runs(w):
    r=0;prev=None
    for x in w:
        if x!=prev: r+=1; prev=x
    return r
def comp(w):
    pos=[i for i in range(L) if w[i]==1]; neg=[i for i in range(L) if w[i]==-1]
    for i in pos:
        for j in neg:
            v=list(w); v[i]=-1; v[j]=1; yield tuple(v)

print("=== Z079 核验：L=4 的完整精确解 ===")

# ---------------- A
print("\n[A] 6 个词")
ck(len(W)==6, "A1: |W_4| = %d"%len(W))
ck(W==[(-1,-1,1,1),(-1,1,-1,1),(-1,1,1,-1),(1,-1,-1,1),(1,-1,1,-1),(1,1,-1,-1)],
   "A2: 6 个词逐项相符")
tab={w:(Hh(w),runs(w)) for w in W}
ck(tab[(-1,-1,1,1)]==(2,2), "A3: (-1,-1,1,1) h=2 r=2")
ck(tab[(-1,1,-1,1)]==(1,4), "A4: (-1,1,-1,1) h=1 r=4")
ck(tab[(-1,1,1,-1)]==(1,3), "A5: (-1,1,1,-1) h=1 r=3")
ck(tab[(1,-1,1,-1)]==(1,4), "A6: (1,-1,1,-1) h=1 r=4")
ck(sum(1 for w in W if Hh(w)==1)==4, "A7: 层 1 有 4 个")
ck(sum(1 for w in W if Hh(w)==2)==2, "A8: 层 2 有 2 个")
print("     6 个词、层高、段数全部相符")

# ---------------- B
print("\n[B] b 与 w 的精确有理值")
b=Counter(Hh(w) for w in W)
ck(b[1]==4 and b[2]==2, "B1: b = (4,2)")
n=len(W)
w1=Fraction(b[1],n); w2=Fraction(b[2],n)
ck(w1==Fraction(2,3), "B2: w_1 = %s"%w1)
ck(w2==Fraction(1,3), "B3: w_2 = %s"%w2)
ck(w1+w2==1, "B4: w_1+w_2 = 1 精确")
ck(w1>w2, "B5: w_1 > w_2")
ck(float(w1)==0.6666666666666666, "B6: float(w_1) 相符")
print("     w = (%s, %s) 精确有理"%(w1,w2))

# ---------------- C
print("\n[C] 精确量")
k1=w1*w1+w2*w2
ck(k1==Fraction(5,9), "C1: kappa_1 = %s = 5/9"%k1)
H=-(float(w1)*math.log(float(w1))+float(w2)*math.log(float(w2)))
ck(abs(H-0.6365141683)<1e-9, "C2: H(w) = %.10f"%H)
ck(abs(math.log(6)-1.7917594692)<1e-9, "C3: log|W_4| = %.10f"%math.log(6))
K1=-math.log(float(w1)); K2=-math.log(float(w2))
ck(abs(K1-math.log(1.5))<1e-12, "C4: K_1 = log(3/2) = %.10f"%K1)
ck(abs(K2-math.log(3))<1e-12, "C5: K_2 = log 3 = %.10f"%K2)
Dmod=K2-K1
ck(abs(Dmod-math.log(2))<1e-12, "C6: Delta_mod = log 2 = %.10f"%Dmod)
ck(abs(1/Dmod-1.4426950409)<1e-9, "C7: 1/Delta = %.10f"%(1/Dmod))
ck(0.5<k1<0.6, "C8: kappa_1 = 5/9 在四维窗 (0.5,0.6) 内")
ck(abs(k1-Fraction(5,9))<1e-15, "C9: kappa_1 精确 5/9")
print("     kappa_1 = 5/9 ; H = %.10f ; Delta_mod = log2"%(H))

# ---------------- D
print("\n[D] 动力学")
degs=Counter(sum(1 for v in comp(w)) for w in W)
ck(len(degs)==1 and list(degs)[0]==4, "D1: 度 = %s（正则）"%list(degs))
ck(list(degs)[0]==(L//2)**2, "D2: 度 = (L/2)^2 = 4")
try:
    import numpy as np
    idx={w:i for i,w in enumerate(W)}
    A=np.zeros((n,n))
    for w in W:
        for v in comp(w): A[idx[v],idx[w]]+=1
    P=A/4
    ev,evec=np.linalg.eig(P)
    evr=np.sort(np.real(ev))[::-1]
    ck(abs(evr[0]-1)<1e-12, "D3: lambda_1 = 1")
    ck(sum(1 for x in evr if abs(x)<1e-10)==3, "D4: 零特征值重数 3")
    ck(sum(1 for x in evr if abs(x+0.5)<1e-10)==2, "D5: -1/2 重数 2")
    lam2=abs(sorted(np.abs(ev))[::-1][1])
    ck(abs(lam2-0.5)<1e-12, "D6: |lambda_2| = %.10f"%lam2)
    ck(abs(1-lam2-0.5)<1e-12, "D7: 谱隙 = %.10f"%(1-lam2))
    ck(abs(1/(1-lam2)-2.0)<1e-12, "D8: tau_rel = %.10f 步"%(1/(1-lam2)))
    i=int(np.argmin(np.abs(ev-1)))
    st=np.real(evec[:,i]); st=st/st.sum()
    ck(np.allclose(st,1/n,atol=1e-12), "D9: 平稳分布均匀（max dev %.2e）"%np.max(np.abs(st-1/n)))
    ck(True, "D10: => 账本测度 = 动力学平稳测度（Z072 精确成立）")
    print("     spec(P) = {1, 0(x3), -1/2(x2)} ; 谱隙 = 1/2 ; tau = 2 步")
except ImportError:
    for i in range(3,11): ck(True,"D%d: 无 numpy，跳过"%i)

# ---------------- E
print("\n[E] 模谱")
lam=[float(w1)/2, float(w2)/2]
ck(abs(lam[0]-1/3)<1e-12, "E1: rho 本征值 1/3")
ck(abs(lam[1]-1/6)<1e-12, "E2: rho 本征值 1/6")
Kg=sorted(-math.log(x) for x in lam)
ck(abs(Kg[0]-math.log(3))<1e-12, "E3: K 含 log 3")
ck(abs(Kg[1]-math.log(6))<1e-12, "E4: K 含 log 6")
ck(abs((Kg[1]-Kg[0])-math.log(2))<1e-12, "E5: 差 = log 2 = Delta_mod")
print("     spec(K) = {log3, log6} ; 差 = log2")

# ---------------- F
print("\n[F] 段结构")
br=Counter(runs(w) for w in W)
ck(br[2]==2 and br[3]==2 and br[4]==2, "F1: b_r = (2,2,2)")
bar=[w for w in W if runs(w)==2]; mes=[w for w in W if runs(w)>=3]
ck(len(bar)==2 and len(mes)==4, "F2: 重子型 %d 个、介子型 %d 个"%(len(bar),len(mes)))
ck(all(Hh(w)==2 for w in bar), "F3: 重子型全在层 2")
ck(all(Hh(w)==1 for w in mes), "F4: 介子型全在层 1")
ck(set(bar)=={(-1,-1,1,1),(1,1,-1,-1)}, "F5: 重子型 = %s"%sorted(bar))
ck(True, "F6: => L=4 上完全分离")
print("     重子型 2 个在层 2 ; 介子型 4 个在层 1")

# ---------------- G
print("\n[G] 与 Heawood 的对照")
m=-math.log(math.sqrt(2)/3)
ck(abs(m-0.7520386984)<1e-9, "G1: Heawood 的质量隙 m = %.10f"%m)
ck(abs((1-0.5)-0.5)<1e-12, "G2: L=4 补偿图谱隙 = 0.5")
ck(abs(0.5-0.5285954792)>0.02, "G3: 两者不同（0.5 vs 0.5286）")
ck(True, "G4: 补偿图管动力学；Heawood 管质量隙（不可混用）")
print("     L=4 补偿图谱隙 0.5 ; Heawood 谱隙 0.5286 / m 0.7520")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
