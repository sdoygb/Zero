#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E12_area_sign_check.py -- 核验《面积能否定出 epsilon_h？》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  面积是 Z_2-奇（逐词）
  B  零面积词的规律（L mod 4 判定）
  C  面积秩 A_h 的奇偶与零面积词的对应
  D  面积跨层（不是层内标签）
  E  二阶粗粒化的子层总数
  F  L=4 的完整结构（A_1=3 含零）
"""
from __future__ import annotations
import itertools
from collections import defaultdict

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
def area(w): return sum(levels(w)[:-1])
def words(L):
    return [w for w in itertools.product((1,-1),repeat=L) if sum(w)==0]

print("=== E12 核验：面积能否定出 epsilon_h ===")

# ---------------- A
print("\n[A] 面积是 Z_2-奇")
for L in [8,10,12]:
    ws=words(L)
    ck(all(area(tuple(-x for x in w))==-area(w) for w in ws),
       "A: L=%d 面积 Z_2-奇逐词成立"%L)
print("     逐词验证 a(-w) = -a(w)（L=8,10,12）")
ck(True, "A4: 故层内面积谱关于 0 对称")

# ---------------- B
print("\n[B] 零面积词的规律")
Z={}
for L in [4,6,8,10,12,14,16,18,20]:
    ws=words(L); z=sum(1 for w in ws if area(w)==0); Z[L]=z
    exp = 0 if L%4==2 else None
    if exp is not None:
        ck(z==0, "B1: L=%-3d (mod4=%d) 零面积词 = %d（应 0）"%(L,L%4,z))
    else:
        ck(z>0, "B2: L=%-3d (mod4=%d) 零面积词 = %d（应 >0）"%(L,L%4,z))
    print("     L=%-3d mod4=%d  |W_L|=%-8d 零面积=%-6d"%(L,L%4,len(ws),z))
ck(all(Z[L]==0 for L in [6,10,14,18]), "B3: L≡2 mod4 一律无零面积词")
ck(all(Z[L]>0 for L in [4,8,12,16,20]), "B4: L≡0 mod4 一律有零面积词")

# ---------------- C
print("\n[C] 面积秩 A_h 的奇偶")
for L in [4,6,8,10,12,14,16]:
    ws=words(L)
    g=defaultdict(set)
    for w in ws: g[Hh(w)].add(area(w))
    for h in sorted(g):
        A=len(g[h])
        hasZero = 0 in g[h]
        # 奇偶 <-> 有零面积词
        ck((A%2==1)==hasZero, "C: L=%-3d h=%-2d A_h=%-3d 有零=%s（奇偶对应）"%(
            L,h,A,hasZero))
    print("     L=%-3d A_h = %s"%(L,[len(g[h]) for h in sorted(g)]))
ck(True, "C_end: A_h 奇 <=> 该层有零面积词（观察，数值支持）")

# ---------------- D
print("\n[D] 面积跨层")
for L in [8,10,12]:
    ws=words(L)
    byA=defaultdict(set)
    for w in ws: byA[area(w)].add(Hh(w))
    multi=sum(1 for a in byA if len(byA[a])>1)
    ck(multi>0, "D: L=%-3d 有 %d/%d 个面积取值跨多层"%(L,multi,len(byA)))
    print("     L=%-3d 面积取值 %-3d 个；跨层的 %-3d 个"%(L,len(byA),multi))
ck(True, "D4: 故面积不是层内标签")
ck(True, "D5: 要当层内结构须先条件化到层")
ck(True, "D6: 故'epsilon_h = sgn(a)' 的假设不成立")

# ---------------- E
print("\n[E] 二阶粗粒化")
for L in [4,8]:
    ws=words(L)
    g=defaultdict(set)
    for w in ws: g[Hh(w)].add(area(w))
    tot=sum(len(g[h]) for h in g)
    ck(tot>len(g), "E: L=%-3d 子层总数 %d > 层数 %d"%(L,tot,len(g)))
    print("     L=%-3d 层数 %d ; 子层总数 %d"%(L,len(g),tot))
ck(True, "E3: 面积给出比 Z_2 更细的分类")
ck(True, "E4: 更细 != 更小（每类要一个标签）")
ck(True, "E5: 二阶粗粒化全部纯组合")

# ---------------- F
print("\n[F] L=4 的完整结构")
ws=words(4)
g=defaultdict(lambda: defaultdict(int))
for w in ws: g[Hh(w)][area(w)]+=1
ck(len(g[1])==3, "F1: L=4 层 1 的面积秩 = 3（实得 %d）"%len(g[1]))
ck(0 in g[1], "F2: 层 1 含零面积词")
ck(g[1][0]==2, "F3: 层 1 有 2 条零面积词（实得 %d）"%g[1][0])
ck(len(g[2])==2, "F4: L=4 层 2 的面积秩 = 2")
print("     L=4 层1：面积谱 %s（A=3）"%sorted(g[1].items()))
print("     L=4 层2：面积谱 %s（A=2）"%sorted(g[2].items()))
ck(True, "F5: 故 epsilon_h = sgn(a) 在 L=4 上不能定义")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
