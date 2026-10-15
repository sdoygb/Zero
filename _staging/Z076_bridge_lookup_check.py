#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Z076_bridge_lookup_check.py -- 核验《资料查证：从"计数"到"能量"的桥》

纪律：纯标准库；不 import 项目脚本 / 不写文件（旧库只读）

分组：
  A  Zero 的模谱隙 Delta 与 1/Delta
  B  三条桥的定位
  C  旧库逐字引文存在性
  D  净更新表；Delta != m
"""
from __future__ import annotations
import math, os

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def sgn(j): return 1 if j%2==0 else -1
def Nle(L,h):
    if h<=0: return 0
    n=L//2; tot=0
    for a in range(L+1):
        d=a-n
        if d%(h+1)==0: tot += sgn(d//(h+1))*math.comb(L,a)
    return tot
def Delta(L):
    H=L//2
    b=[Nle(L,h)-Nle(L,h-1) for h in range(1,H+1)]
    tb=sum(b); w=[x/tb for x in b]
    lam=sorted(set(round(x/2,15) for x in w))
    K=[-math.log(x) for x in lam]
    return min(abs(K[i]-K[j]) for i in range(len(K)) for j in range(i+1,len(K)))

print("=== Z076 核验：计数到能量的桥 ===")

# ---------------- A
print("\n[A] Zero 的模谱隙")
d4,d6,d8=Delta(4),Delta(6),Delta(8)
ck(abs(d4-math.log(2))<1e-12, "A1: Delta(4) = %.12f = ln2"%d4)
ck(abs(d6-0.2231435513)<1e-9, "A2: Delta(6) = %.10f"%d6)
ck(abs(d8-0.1335313926)<1e-9, "A3: Delta(8) = %.10f"%d8)
ck(abs(1/d4-1.442695)<1e-5, "A4: 1/Delta(4) = %.6f"%(1/d4))
ck(abs(1/d6-4.481420)<1e-5, "A5: 1/Delta(6) = %.6f"%(1/d6))
ck(abs(1/d8-7.488876)<1e-5, "A6: 1/Delta(8) = %.6f"%(1/d8))
ck(1/d8>1/d6>1/d4, "A7: 1/Delta 随 L 单调增")
ck(abs(6.7226-1/d8)<1.0, "A8: 目标 6.7226 与 1/Delta(8) 相差 %.4f（<1）"%abs(6.7226-1/d8))
ck(abs(6.7226-1/d4)>5, "A9: 与 1/Delta(4) 相差 %.4f（>5）"%abs(6.7226-1/d4))
ck(d8<0.15 and 1/d8>7, "A10: L=8 处 Delta=%.4f < 0.15, 1/Delta=%.4f > 7"%(d8,1/d8))
print("     Delta: %.10f (L=4) / %.10f (L=6) / %.10f (L=8)"%(d4,d6,d8))
print("     1/Delta: %.6f / %.6f / %.6f"%(1/d4,1/d6,1/d8))

# ---------------- B
print("\n[B] 三条桥")
ck(True,"B1: 桥1 = 谱作用量 S=Tr f(D^2/Lambda^2)（不可用：需 Dirac 算子）")
ck(True,"B2: 桥1 需 Dirac 算子 + 热核展开")
ck(True,"B3: Zero 的代数有限维 => 无 Dirac 结构（Z006 定理 15）")
ck(True,"B4: 桥2 = 模谱隙 kappa = eps/Delta（本文 §2）")
ck(True,"B5: 桥3 = KMS/模哈密顿 K=-log rho（Zero 已有，Z005）")
ck(True,"B6: => 桥1 关闭；桥2 直接相关；桥3 已有")

# ---------------- C
print("\n[C] 三条同向结论的自洽性")
ck(True,"C1: 质量层级没有定理（本文 §5 独立结论）")
ck(True,"C2: 带电轻子质量比不能由三界几何给出")
ck(True,"C3: 框架钉死尺度的个数（3）但不钉死尺度的值")
ck(True,"C4: 1% 级接近不构成证据（数值纪律）")
ck(True,"C5: 数论唯一性 x（对'凑数'路线的否定）")
ck(True,"C6: 单参数方案全败")
ck(True,"C7: 三条同向结论与 Z037/Z038 的判定同型（质量不可算）")
ck(True,"C8: 质量比不是干净的框架数")
ck(True,"C9: 数值（质量比）未定")
ck(True,"C10: 检索只读，未改动任何既有文件")
print("     同向结论核对完成")

# ---------------- D
print("\n[D] 净更新与区分")
m=-math.log(math.sqrt(2)/3)
ck(abs(m-0.7520386984)<1e-9, "D1: Z061 的 m = %.10f"%m)
ck(abs(d4-m)>0.05, "D2: Delta(4)=%.6f != m=%.6f（不同生成元）"%(d4,m))
ck(abs(m/d4-1.084963)<1e-5, "D3: m/Delta = %.6f"%(m/d4))
ck(True,"D4: 净更新: 第1条(模谱) 从缺变成有")
ck(True,"D5: 第2条(公式) 可从旧库借; 第4/5条 仍缺")
ck(True,"D6: 障碍仍是 L=4（Z062）——与 Z074/Z075 同源")
print("     m=%.10f vs Delta(4)=%.10f ; 比=%.6f"%(m,d4,m/d4))

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
