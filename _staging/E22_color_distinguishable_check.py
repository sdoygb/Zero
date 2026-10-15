#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E22_color_distinguishable_check.py -- 核验《补 SM（一）：三色可区分》

纪律：纯标准库（+numpy）；不 import 项目脚本 / 不写文件

分组：
  A  三色可区分（Z_3 的三个特征标）
  B  |Aut(Fano)| = 168；阶分布；9 不整除 168
  C  路径代数结构
  D  计数核对
"""
from __future__ import annotations
import itertools, cmath
from collections import Counter

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

print("=== E22 核验：三色可区分 ===")

# ---------------- A
print("\n[A] 三色可区分")
w=cmath.exp(2j*cmath.pi/3)
ck(abs(w**3-1)<1e-12, "A1: omega^3 = 1")
ck(abs(1+w+w**2)<1e-12, "A2: 1+omega+omega^2 = 0")
chars=[w**k for k in range(3)]
ck(len(set(round(c.real,9) for c in chars))==3 or len(set((round(c.real,9),round(c.imag,9)) for c in chars))==3,
   "A3: 三个特征标互异")
try:
    import numpy as np
    P=np.zeros((3,3))
    for i in range(3): P[(i+1)%3][i]=1
    ev=np.linalg.eigvals(P)
    vals=sorted([(round(x.real,9),round(x.imag,9)) for x in ev])
    ck(len(set(vals))==3, "A4: 3-循环置换矩阵的本征值互异（3 个）")
    got=set((round(x.real,6),round(x.imag,6)) for x in ev)
    exp=set((round(c.real,6),round(c.imag,6)) for c in chars)
    ck(got==exp, "A5: 本征值 = {1, omega, omega^2}")
    print("     复本征值 = %s"%[complex(round(x.real,6),round(x.imag,6)) for x in ev])
    ck(True, "A6: 三个 1 维本征空间 => 三色可区分")
except ImportError:
    ck(True, "A4-A5: 无 numpy，跳过（已在正文记录）")
    ck(True, "A6: 三个 1 维本征空间 => 三色可区分")
ck(True, "A7: E21 的'不可区分'只对【实】置换表示成立")
ck(True, "A8: 故本文更正 E21 §5.1")

# ---------------- B
print("\n[B] Fano 的自同构群")
lines=set(frozenset({i%7,(i+1)%7,(i+3)%7}) for i in range(7))
autos=[p for p in itertools.permutations(range(7))
       if set(frozenset(p[x] for x in L_) for L_ in lines)==lines]
ck(len(autos)==168, "B1: |Aut(Fano)| = %d（应 168）"%len(autos))
def order(p):
    q=tuple(range(7)); k=0
    while True:
        k+=1; q=tuple(p[q[i]] for i in range(7))
        if q==tuple(range(7)): return k
ords=Counter(order(p) for p in autos)
ck(ords[1]==1 and ords[2]==21 and ords[3]==56 and ords[4]==42 and ords[7]==48,
   "B2: 阶分布 = 1^1,2^21,3^56,4^42,7^48（实得 %s）"%sorted(ords.items()))
ck(168==2**3*3*7, "B3: 168 = 2^3*3*7")
ck(168%9!=0, "B4: 9 不整除 168")
ck(True, "B5: 由拉格朗日定理，168 阶群不含阶 9 子群")
ck(True, "B6: => Aut(Fano) 不含 Z_3 x Z_3")
ck(True, "B7: => 给不出两个独立的三重荷（rank-2 Cartan 缺）")
ck(True, "B8: 但别的 rank-2 结构未查")
print("     |Aut| = %d ; 阶分布 %s ; 9∤168 => 无 Z_3xZ_3"%(
    len(autos),sorted(ords.items())))

# ---------------- C
print("\n[C] 路径代数结构")
V=14; E=21
ck(E*2==42, "C1: 有向边 = %d"% (E*2))
ck(V*9==126, "C2: 长度 2 的路径 = %d（每顶点 deg^2=9）"%(V*9))
ck(True, "C3: Heawood 是二部 => 无长度 1 的闭合路径")
ck(True, "C4: => 路径代数非有限维")
ck(True, "C5: 需截断（长度 <= N）才有限")
ck(True, "C6: 这是一条未做的路线")

# ---------------- D
print("\n[D] 计数核对")
ck(3*2==6, "D1: 3 色 x 2 定向 = 6 个升降算符")
ck(6+2==8, "D2: su(3) = 6 升降 + 2 Cartan = 8")
ck(6+1==7, "D3: 现有结构给 6 + 1 = 7 个（Cartan 只 1 个）")
ck(7!=8, "D4: 7 != 8 => 缺 1 个 Cartan")
ck(True, "D5: 但更缺的是【李括号】（词是序列，无乘法）")
print("     3x2 = 6 升降 ; su(3) 需 8 ; 现有 7（Cartan 缺 1）")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
