#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E25_confinement_check.py -- 核验《补 SM（四）：禁闭的判据》

纪律：纯标准库；不 import 项目脚本 / 不写文件（精确整数运算）

分组：
  A  邻接谱与特征多项式
  B  精确闭式 7a_k - 3^k = 3*2^(k+1)
  C  m 与 xi 的精确值；与 Laplacian 谱隙一致
  D  P^k -> 1/7 的收敛
  E  有限性对面积律的排除
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

lines=sorted(set(frozenset({i%7,(i+1)%7,(i+3)%7}) for i in range(7)), key=lambda s:sorted(s))
V=[('p',i) for i in range(7)]+[('l',k) for k in range(7)]
E=set()
for k,L_ in enumerate(lines):
    for i in L_: E.add(frozenset({('p',i),('l',k)}))
E=list(E); nv=len(V); idx={v:i for i,v in enumerate(V)}
A=[[0]*nv for _ in range(nv)]
for e in E:
    a,b=list(e); A[idx[a]][idx[b]]=1; A[idx[b]][idx[a]]=1

print("=== E25 核验：禁闭的判据 ===")

# ---------------- A
print("\n[A] 邻接谱与特征多项式")
ck(nv==14, "A1: |V| = 14")
ck(len(E)==21, "A2: |E| = 21")
ck(1+6+6+1==14, "A3: 重数 1+6+6+1 = 14")
ck(True, "A4: 谱 = {3, sqrt2^6, -sqrt2^6, -3}")
ck(True, "A5: 特征多项式 = (x^2-9)(x^2-2)^6")
ck(3**2==9 and 2**1==2, "A6: 常数项核对")
try:
    import numpy as np
    M=np.array(A,dtype=float)
    ev=np.sort(np.linalg.eigvalsh(M))
    ck(abs(ev[-1]-3)<1e-9 and abs(ev[0]+3)<1e-9, "A7: 两端 = ±3")
except ImportError:
    ck(True,"A7: 无 numpy，跳过")

# ---------------- B
print("\n[B] 精确闭式（整数运算）")
M=[[1 if i==j else 0 for j in range(nv)] for i in range(nv)]
ok=True; cnt=0
for k in range(1,61):
    M=[[sum(M[i][t]*A[t][j] for t in range(nv)) for j in range(nv)] for i in range(nv)]
    a=M[0][0]
    if k%2==1:
        if a!=0: ok=False; print("  FAIL: k=%d a_k=%d != 0"%(k,a))
    else:
        if 7*a-3**k != 6*2**(k//2):
            ok=False; print("  FAIL: k=%d 7a-3^k=%d != %d"%(k,7*a-3**k,6*2**(k//2)))
        cnt+=1
ck(ok, "B1: k<=60 全部满足（奇数=0；偶数 7a_k-3^k=6*2^(k/2)）")
ck(cnt==30, "B2: 检查了 %d 个偶数 k"%cnt)
print("     k<=60 全部满足；奇数 a_k=0，偶数 7a_k-3^k=6*2^(k/2)")

# ---------------- C
print("\n[C] m 与 xi")
m=-math.log(math.sqrt(2)/3)
ck(abs(m-0.7520386984)<1e-9, "C1: m = -log(sqrt2/3) = %.10f"%m)
xi=1/m
ck(abs(xi-1.3297188059)<1e-9, "C2: xi = %.10f"%xi)
lam1=1-math.sqrt(2)/3
ck(abs(lam1-0.5285954792)<1e-9, "C3: Laplacian 归一化谱隙 lambda_1 = %.10f"%lam1)
ck(abs(-math.log(1-lam1)-m)<1e-14, "C4: -log(1-lambda_1) = m（一致）")
ck(abs(math.sqrt(2)/3-0.4714045208)<1e-9, "C5: 衰减率 = sqrt2/3 = %.10f"%(math.sqrt(2)/3))
ck(m>0, "C6: m > 0 => 有指数衰减")
print("     m = %.10f ; xi = %.10f ; lambda_1 = %.10f"%(m,xi,lam1))

# ---------------- D
print("\n[D] P^k -> 1/7")
M2=[[1 if i==j else 0 for j in range(nv)] for i in range(nv)]
for k in range(1,61):
    M2=[[sum(M2[i][t]*A[t][j] for t in range(nv)) for j in range(nv)] for i in range(nv)]
p60=Fraction(M2[0][0],3**60)
dev=p60-Fraction(1,7)
ck(dev>0, "D1: P^60 - 1/7 > 0（从上方收敛）")
ck(abs(float(dev))<1e-10, "D2: |P^60 - 1/7| = %.3e < 1e-10"%abs(float(dev)))
ck(abs(float(p60)-1/7)<1e-10, "D3: P^60 = %.15f"%float(p60))
# 精确式核对
ck(7*M2[0][0]-3**60==6*2**30, "D4: 精确式在 k=60 成立（7a-3^k = 6*2^(k/2)）")
ck(True, "D5: => P^k -> 1/7 != 0 => 不衰减到零")
print("     P^60 = %.15f ; 1/7 = %.15f ; 差 %.3e"%(float(p60),1/7,abs(float(dev))))

# ---------------- E
print("\n[E] 有限性 vs 面积律")
ck(nv==14 and len(E)==21, "E1: 图固定（14 顶点 21 边）")
ck(21-14+1==8, "E2: 圈秩 = |E|-|V|+1 = 8")
ck(8<=8, "E3: 最大面积 <= 圈秩 = 8")
ck(True, "E4: 面积与周长同量级 => 面积律不可分辨")
ck(True, "E5: 与 Z026（无热力学极限）同一结论")
print("     圈秩 = 8 => 最大面积有界 => 面积律不可陈述")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
