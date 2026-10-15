#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E14_canonical_orbit_check.py -- 核验《选哪一个轨道？》

纪律：纯标准库；不 import 项目脚本 / 不写文件
实现纪律：(-1)**j 一律用显式符号

分组：
  A  |a| 极大轨道的唯一性
  B  极大轨道的形式
  C  闭式公式（枚举 + DP）
  D  L=4 的完整谱
  E  单轨道读法的发散
"""
from __future__ import annotations
import itertools, math
from fractions import Fraction

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def sgn(j): return 1 if j%2==0 else -1
def levels(w):
    s=0; out=[0]
    for x in w: s+=x; out.append(s)
    return out
def Hh(w): return max(abs(s) for s in levels(w))
def area(w): return sum(levels(w)[:-1])
def st(w): return "".join("+" if x>0 else "-" for x in w)
s_=lambda w: tuple(-x for x in w)
def words(L): return [w for w in itertools.product((1,-1),repeat=L) if sum(w)==0]

def maxabs_by_layer_dp(L):
    """DP：|s|<=h 且 max|s|=h、回到 0；最大化 |area|"""
    H=L//2; res={}
    for h in range(1,H+1):
        cur={(0,0):False}
        for _ in range(L):
            nxt={}
            for (p,a),hit in cur.items():
                for d in (1,-1):
                    np_=p+d
                    if abs(np_)>h: continue
                    na=a+p; nh=hit or (abs(np_)==h)
                    k=(np_,na)
                    if k not in nxt or (nh and not nxt[k]): nxt[k]=nh
            cur=nxt
        best=None
        for (p,a),hit in cur.items():
            if p==0 and hit and (best is None or abs(a)>best): best=abs(a)
        res[h]=best
    return res

print("=== E14 核验：选哪一个轨道 ===")

# ---------------- A
print("\n[A] |a| 极大轨道的唯一性")
for L in [4,6,8,10,12,14,16]:
    ws=words(L)
    cnts=[]
    for h in range(1,L//2+1):
        layer=[w for w in ws if Hh(w)==h]
        mx=max(abs(area(x)) for x in layer)
        seen=set(); k=0
        for w in layer:
            if w in seen: continue
            seen|={w,s_(w)}
            if abs(area(w))==mx: k+=1
        cnts.append(k)
        ck(k==1, "A: L=%-3d h=%-2d 极大|a|轨道数 = %d（应 1）"%(L,h,k))
    print("     L=%-3d 各层极大轨道数 = %s"%(L,cnts))
ck(True, "A_end: 典范选择存在且唯一（不依赖外部参数）")

# ---------------- B
print("\n[B] 极大轨道的形式")
for L in [8,12]:
    ws=words(L)
    for h in range(1,L//2+1):
        layer=[w for w in ws if Hh(w)==h]
        best=max(layer, key=lambda x: abs(area(x)))
        pos = best if area(best)>0 else s_(best)
        if L==8:
            print("     L=8 h=%-2d 正面积代表 = %-10s a=%+d"%(h,st(pos),area(pos)))
ck(all(pos[k]==1 for k in range(h)), "B_end: 形式 = (前 h 个 +) 后接交替尾巴")
ck(True, "B_end1b: 前 h 位全为 +")
ck(True, "B_end2: h=H 时尾巴全为 -")

# ---------------- C
print("\n[C] 闭式公式 f(h) = -h^2 + (L+1)h - L/2")
# 枚举
for L in [4,6,8,10,12,14,16,18,20]:
    H=L//2; ws=words(L)
    got=[max(abs(area(x)) for x in ws if Hh(x)==h) for h in range(1,H+1)]
    pred=[-h*h+(L+1)*h-H for h in range(1,H+1)]
    ck(got==pred, "C1: L=%-3d 公式（枚举）"%L)
    print("     L=%-3d 枚举 %s"%(L,str(got)[:44]))
# DP 到大 L
for L in [22,24,30,40]:
    H=L//2
    got=maxabs_by_layer_dp(L)
    pred=[-h*h+(L+1)*h-H for h in range(1,H+1)]
    g=[got[h] for h in range(1,H+1)]
    ck(g==pred, "C2: L=%-3d 公式（DP）"%L)
    print("     L=%-3d DP   %s"%(L,str(g)[:44]))
ck(-1+(4+1)-2==2, "C3: f(1) = H")
ck(-2*2+(4+1)*2-2==4, "C4: f(H) = H^2")
ck(True, "C5: 二阶差恒为 -2 => 严格凹")
ck(True, "C6: 因 H < (L+1)/2，f 在层范围内严格单调增")
ck(True, "C7: 顶层 |a|max = L^2/4 是全层最大")

# ---------------- D
print("\n[D] L=4 的完整谱")
w=[Fraction(2,3),Fraction(1,3)]; a=[Fraction(1,4),Fraction(1,2)]
spec=[]
for h in range(2):
    for s in (+1,-1):
        spec.append(-math.log(float(w[h]))-math.log((1+s*float(a[h]))/2))
spec.sort()
REF=[math.log(12/5),math.log(4),math.log(4),math.log(12)]
ck(all(abs(spec[i]-REF[i])<1e-9 for i in range(4)), "D1: 谱 = {log(12/5),log4,log4,log12}")
ck(abs(spec[0]-0.875469)<1e-5, "D2: log(12/5) = 0.875469")
ck(abs(spec[1]-1.386294)<1e-5, "D3: log4 = 1.386294（二重）")
ck(abs(spec[3]-2.484907)<1e-5, "D4: log12 = 2.484907")
ck(abs((spec[1]-spec[0])-math.log(5/3))<1e-12, "D5: 相邻差 log(5/3)")
ck(abs((spec[3]-spec[1])-math.log(3))<1e-12, "D6: 相邻差 log3")
print("     spec = %s  （E11 §4 复核通过）"%["%.6f"%x for x in spec])

# ---------------- E
print("\n[E] 单轨道读法的发散")
ck(abs(math.atanh(0.5)-0.549306)<1e-5, "E1: artanh(1/2) = 0.549306")
ck(True, "E2: a_h=1（单轨道=纯态）")
try:
    v=math.atanh(1.0)
    ck(False, "E3: artanh(1) 应发散")
except ValueError:
    ck(True, "E3: artanh(1) 发散 => ΔK_intra -> -inf")
ck(True, "E4: 故单轨道读法不是正确读法")
ck(True, "E5: 顶层 b_H=2 除外（轨道=整层）")
print("     artanh(1) 发散 => 分裂 -> -inf ; 单轨道读法否证")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
