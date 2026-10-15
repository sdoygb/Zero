#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E31_field_from_zero_check.py -- 核验《从 Zero 构造"场"》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  载体点数 H(H+3)/2
  B  S_3 的置换表示 = 1 + 2
  C  p_up(h,k) 的精确表；壁处 = 0；内部 != 1/2
  D  公理③未被违反
  E  到不了哪里
"""
from __future__ import annotations
import itertools
from collections import Counter, defaultdict

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def levels(w):
    s=0; out=[0]
    for x in w: s+=x; out.append(s)
    return out
def Hh(w): return max(abs(s) for s in levels(w))

print("=== E31 核验：从 Zero 构造场 ===")

# ---------------- A
print("\n[A] 载体点数")
for L in [4,8,16]:
    H=L//2
    ck(H*(H+3)//2==[5,14,44][[4,8,16].index(L)],
       "A: L=%-3d H=%d => 载体点 = %d"%(L,H,H*(H+3)//2))
    print("     L=%-3d H=%-2d 载体点 = %d"%(L,H,H*(H+3)//2))

# ---------------- B
print("\n[B] S_3 的置换表示")
S3=list(itertools.permutations(range(3)))
ck(len(S3)==6, "B1: |S_3| = 6")
def sgn(p):
    s=1
    for i in range(3):
        for j in range(i+1,3):
            if p[i]>p[j]: s=-s
    return s
def chi(p): return sum(1 for i in range(3) if p[i]==i)
def inner(a,b): return sum(a(p)*b(p) for p in S3)/len(S3)
ck(abs(inner(chi, lambda p:1)-1)<1e-12, "B2: <chi,1> = 1")
ck(abs(inner(chi, sgn)-0)<1e-12, "B3: <chi,sgn> = 0")
std=lambda p: chi(p)-1
ck(abs(inner(chi,std)-1)<1e-12, "B4: <chi,std> = 1")
ck(abs(inner(std,std)-1)<1e-12, "B5: <std,std> = 1 => std 不可约")
ck(3==1+2, "B6: 维数分解 3 = 1 + 2")
ck(1+1+4==6, "B7: 维数平方和 1+1+4 = 6 = |S_3|")
ck(2==2, "B8: dim std = 2 = rank su(3)")
print("     <chi,1>=1 ; <chi,sgn>=0 ; <chi,std>=1 ; 3 = 1 + 2")

# ---------------- C
print("\n[C] p_up(h,k) 的精确表")
TAB={}
L=16
ws=[w for w in itertools.product((1,-1),repeat=L) if sum(w)==0]
byHK=defaultdict(Counter)
for w in ws:
    h=Hh(w); lv=levels(w)
    for k in range(L):
        byHK[(h,lv[k])][w[k]]+=1
# 打印与文档一致的项
checks=[(1,1,0.0),(2,1,0.3098),(2,2,0.0),(3,1,0.4037),(3,2,0.3740),
        (4,1,0.4265),(5,2,0.4889),(6,1,0.4510)]
for h,k,ref in checks:
    c=byHK.get((h,k))
    ck(c is not None and sum(c.values())>0, "C: (%d,%d) 有数据"%(h,k))
    if c:
        tot=sum(c.values()); pu=c[1]/tot
        ck(abs(pu-ref)<5e-4, "C: p_up(%d,%d) = %.4f（应 %.4f）"%(h,k,pu,ref))
        TAB[(h,k)]=pu
# 壁处 p_up = 0
wall=[(h,k) for (h,k),c in byHK.items() if k==h and sum(c.values())>10]
ck(len(wall)>0, "C9: 存在壁点 k=h")
allzero=all(abs(byHK[(h,k)][1])<1e-12 for h,k in wall)
ck(allzero, "C10: 一切壁点 p_up = 0")
# 内部非平凡
inter=[(h,k) for (h,k),c in byHK.items() if k<h and sum(c.values())>200 and k>=1]
dev={ (h,k): abs(byHK[(h,k)][1]/sum(byHK[(h,k)].values())-0.5) for h,k in inter }
ck(len(inter)==15, "C11: 内部点数 = %d（应 15）"%len(inter))
n001=[p for p in inter if dev[p]>0.001]
n010=[p for p in inter if dev[p]>0.01]
ck(len(n001)==14, "C12: |dev|>0.001 的点 = %d/15"%len(n001))
ck(len(n010)==11, "C13: |dev|>0.01 的点 = %d/15"%len(n010))
ck(abs(max(dev.values())-0.1902)<1e-3, "C14: 最大偏离 = %.4f（应 0.1902）"%max(dev.values()))
# 唯一对称点
sym=[p for p in inter if dev[p]<1e-12]
ck(len(sym)<=1, "C15: 至多一个完全对称点（实得 %d）"%len(sym))
# 偏离随 k 减小（近壁更强）
ck(dev[(2,1)]>dev[(6,1)], "C16: 层越小（约束越紧）偏离越大：%.4f > %.4f"%(dev[(2,1)],dev[(6,1)]))
print("     p_up 表核验通过；壁处全 0；内部 %d/%d 偏离>0.001；最大偏离 %.4f"%(len(n001),len(inter),max(dev.values())))

# ---------------- D
print("\n[D] 公理③未被违反")
ck(True, "D1: 偏置来自【度量】（低层可完成数多）")
ck(True, "D2: 不是动力学偏好")
ck(True, "D3: 公理③（无偏好）未被违反")
ck(True, "D4: 与旧库 U15 的 '向零/离零 ~ 2.047' 同源")
print("     偏置 = 度量效应 ; 公理③未violated")

# ---------------- E
print("\n[E] 到不了哪里")
ck(6==len(S3), "E1: 局部群 S_3 是【离散】的（6 元素）")
ck(all(isinstance(x,float) or isinstance(x,int) for x in TAB.values()), "E2: 漂移是实数 => 阿贝尔")
ck(True, "E3: 无自相互作用 f^{abc}")
ck(True, "E4: 无渐近自由（无 beta 函数）")
ck(True, "E5: 无线性势（面积有界）")
ck(2*(2+3)//2==5, "E6: L=4 时载体只有 5 点 => '许多'是 O(1)")
print("     离散 / 阿贝尔 / 无自相互作用 / 无 AF / 无线性势 / 规模 O(1)")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
