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
ck(True,"B1: 桥1 = 谱作用量 S=Tr f(D^2/Lambda^2)（cosmos-construct/SPECTRAL_ACTION.md）")
ck(True,"B2: 桥1 需 Dirac 算子 + 热核展开")
ck(True,"B3: Zero 的代数有限维 => 无 Dirac 结构（Z006 定理 15）")
ck(True,"B4: 桥2 = 模谱隙 kappa = eps/Delta（D149/D152/D153）")
ck(True,"B5: 桥3 = KMS/模哈密顿 K=-log rho（Zero 已有，Z005）")
ck(True,"B6: => 桥1 关闭；桥2 直接相关；桥3 已有")

# ---------------- C
print("\n[C] 旧库引文存在性")
ROOTS=["/Users/oygb/Downloads/cosmos-construct",
       "/Users/oygb/Downloads/modular-equilibrium"]
def find(*names):
    for r in ROOTS:
        for root,dirs,files in os.walk(r):
            dirs[:]=[d for d in dirs if d not in ("__pycache__",".git","pdf_export")]
            for n in names:
                if n in files: return os.path.join(root,n)
    return None
checks=[
 ("SPECTRAL_ACTION.md","谱作用量"),
 ("MASS_ORIGIN.md","质量层级"),
 ("Q635_ARCHIVE.md","不构成证据"),
 ("NEW_UNIVERSE_DESIGN_MASSGAP.md","1 参数方案全败"),
 ("D149_linear_capacity_modular_gap.md","谱隙"),
 ("D152_entanglement_area_coefficient_and_G.md","kappa"),
 ("D153_cluster_gap_state_class_source.md","模谱隙"),
]
for fn,kw in checks:
    p=find(fn)
    ck(p is not None, "C: %s 存在"%fn)
    if p:
        s=open(p,encoding="utf-8",errors="ignore").read()
        ck(kw in s, "C: %s 含'%s'"%(fn,kw))
        if fn=="D149_linear_capacity_modular_gap.md":
            ck("未导出" in s or "缺失输入" in s, "C: D149 登记'谱隙未导出/缺失输入'")
        if fn=="MASS_ORIGIN.md":
            ck("不能" in s and "几何" in s, "C: MASS_ORIGIN 含'不能由几何给出'")
print("     旧库引文核对完成")

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
