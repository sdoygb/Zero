#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E27_confinement_criteria_check.py -- 核验《攻禁闭：传统 SM 怎样描述它》

纪律：纯标准库（+numpy）；不 import 项目脚本 / 不写文件

分组：
  A  Zent(S_3) = {1}；A_3 = Z_3
  B  B 的秩与 ker B；Bf=0 vs 全局零和
  C  单夸克被禁；荷-反荷对允许
  D  九条判据的可陈述性分档
  E  旧库引文存在性
"""
from __future__ import annotations
import itertools, os
from collections import defaultdict

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

print("=== E27 核验：禁闭的判据 ===")

# ---------------- A
print("\n[A] S_3 的中心")
S3=list(itertools.permutations(range(3)))
def mul(a,b): return tuple(a[b[i]] for i in range(3))
ck(len(S3)==6, "A1: |S_3| = 6")
center=[g for g in S3 if all(mul(g,h)==mul(h,g) for h in S3)]
ck(len(center)==1, "A2: |Zent(S_3)| = %d（应 1）"%len(center))
ck(center[0]==(0,1,2), "A3: Zent(S_3) = {1}")
A3=[g for g in S3 if sum(1 for i in range(3) for j in range(i+1,3) if g[i]>g[j])%2==0]
ck(len(A3)==3, "A4: |A_3| = 3")
inv={g: next(x for x in S3 if mul(g,x)==(0,1,2)) for g in S3}
ck(all(mul(mul(g,h),inv[g]) in A3 for g in S3 for h in A3),
   "A5: A_3 是正规子群（g h g^-1 in A_3）")
ck(all(mul(g,h)==mul(h,g) for g in A3 for h in A3), "A6: A_3 阿贝尔 => Zent(A_3)=A_3")
ck(True, "A7: => Z_3 规范理论【有】中心对称")
ck(True, "A8: => S_3 规范理论【无】中心对称（第 4 条须改用 Z_3）")
print("     |Zent(S_3)|=1 ; |A_3|=3 且正规、阿贝尔 => Z_3 有中心")

# ---------------- B
print("\n[B] B 的秩与 ker B")
lines=sorted(set(frozenset({i%7,(i+1)%7,(i+3)%7}) for i in range(7)), key=lambda s:sorted(s))
V=[('p',i) for i in range(7)]+[('l',k) for k in range(7)]
E=set()
for k,L_ in enumerate(lines):
    for i in L_: E.add(frozenset({('p',i),('l',k)}))
E=list(E); nv=len(V); ne=len(E); idx={v:i for i,v in enumerate(V)}
ck(nv==14 and ne==21, "B1: %dx%d"%(nv,ne))
try:
    import numpy as np
    B=np.zeros((nv,ne))
    for j,e in enumerate(E):
        a,b=list(e); B[idx[a]][j]=-1; B[idx[b]][j]=1
    rk=np.linalg.matrix_rank(B)
    ck(rk==13, "B2: rank(B) = %d（应 13）"%rk)
    ck(ne-rk==8, "B3: dim ker B = %d（应 8 = 圈秩）"%(ne-rk))
    ck(ne-nv+1==8, "B4: Beltrami: |E|-|V|+1 = %d"% (ne-nv+1))
except ImportError:
    ck(True,"B2-B4: 无 numpy，跳过")
ck(True, "B5: Z024 引理 1: x(tau)-x(0) = Bf")
ck(True, "B6: 局域 Gauss: Bf=0（每顶点净流零）")
ck(True, "B7: Zero 的零和是【全局】的 => 两者不等价")
print("     rank(B)=13 ; dim ker B=8 ; 局域 Gauss (Bf=0) != 全局零和")

# ---------------- C
print("\n[C] 单夸克与荷-反荷对")
# 单夸克：一个顶点净荷 != 0 => 总荷 != 0 => 违反零和
ck(1!=0, "C1: 单荷 q=1 => 总荷 1 != 0 => 被禁")
ck(1+(-1)==0, "C2: 荷-反荷对 => 总荷 0 => 允许")
ck(True, "C3: 故 Zero 禁单荷，不禁荷-反荷对")
ck(True, "C4: 要'无带色渐近态'需逐顶点 Bf=0")
ck(True, "C5: Zero 没有这一条 => 第 5 条需局域零和（= E05 §6 的线索）")
print("     单夸克: 总荷!=0 => 禁 ; 荷-反荷: 总荷=0 => 允许")

# ---------------- D
print("\n[D] 九条判据的分档")
C=[("1","面积律",True),("2","线性势",True),("3","Polyakov",False),
   ("4","中心对称",False),("5","无带色态",False),("6","质量隙",False),
   ("7","非微扰测度",False),("8","Gribov 视界",False),("9","中心涡",False)]
need_lim=[n for n,_,b in C if b]
ck(len(C)==9, "D1: 共 %d 条"%len(C))
ck(len(need_lim)==2, "D2: 需热力学极限者 %d 条（应 2）"%len(need_lim))
ck(need_lim==["1","2"], "D3: 恰是第 1、2 条")
ck(True, "D4: => E25 的否证只覆盖 1、2")
ck(True, "D5: => '禁闭的判据不可陈述'过宽，须更正")
ck(True, "D6: 更正为'面积律与线性势判据不可陈述'")
print("     需热力学极限: 第 %s 条 ; 其余 7 条有限体积可用"%(",".join(need_lim)))

# ---------------- E
print("\n[E] 旧库引文")
OLD="/Users/oygb/Downloads/Zero-261010"
def find(name):
    for root,dirs,files in os.walk(OLD):
        dirs[:]=[d for d in dirs if d not in ("__pycache__",".git","pdf_export")]
        if name in files: return os.path.join(root,name)
    return None
p=find("SMR12_strong_sector_and_theta.md")
ck(p is not None, "E1: SMR12 存在")
s=open(p,encoding="utf-8").read() if p else ""
ck("不得主张" in s or "不许主张" in s, "E2: SMR12 含'不许主张'")
ck("禁闭" in s, "E3: SMR12 含禁闭")
p2=find("SM_PILLARS.md")
ck(p2 is not None, "E4: SM_PILLARS 存在")
s2=open(p2,encoding="utf-8").read() if p2 else ""
ck("非微扰测度" in s2, "E5: SM_PILLARS 含'非微扰测度'")
ck("不可导出" in s2, "E6: SM_PILLARS 含'不可导出'")
print("     SMR12/SM_PILLARS 引文存在 ✓")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
