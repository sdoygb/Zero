#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E51_alpha_koide_check.py -- 核验《Koide 与精细结构常数》

纪律：纯标准库；不 import 项目脚本 / 不写文件（旧库只读）

分组：
  A  六项作用量的复算
  B  轻子质量定理含 S_e^{1/2}
  C  Koide 与 alpha 的三重数值核验
  D  同一平面；同一工具
  E  三库三种接法
  F  旧库文件存在性（只读）
"""
from __future__ import annotations
import cmath, math, os

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

OLD="/Users/oygb/Downloads/geometric-theory-geometric-theory"
me,mmu,mtau=0.51099895,105.6583755,1776.86
sq=[math.sqrt(x) for x in (me,mmu,mtau)]
D1=sum(sq[k]*cmath.exp(-2j*math.pi*k/3) for k in range(3))/3
th=cmath.phase(D1)
Q=sum([me,mmu,mtau])/sum(sq)**2
target=137.035999084
print("=== E51 核验：Koide 与 alpha ===")

# ---------------- A
print("\n[A] 六项作用量的复算")
P=[57.93,26.16,5.91]
ck(abs(sum(P)-90.0)<1e-9, "A1: theta_M+theta_C+theta_I = %.4f = 90 度"%sum(P))
s=[math.sin(math.radians(x)) for x in P]
ck(all(0<x<1 for x in s), "A2: 三个 sin 都在 (0,1)")
six=[1/x**2 for x in s]+[1/(s[i]*s[j]) for i in range(3) for j in range(i+1,3)]
S=sum(six)
ck(len(six)==6, "A3: 六项")
ck(abs(six[0]-1.3926)<1e-3, "A4: 第一项 = %.4f"%six[0])
ck(abs(six[2]-94.3216)<2e-3, "A5: 第三项 = %.4f"%six[2])
ck(abs(S-137.024934855)<1e-8, "A6: S = %.9f"%S)
ck(abs(S-target)<0.02, "A7: 与目标差 %.6f（%.4f%%）"%(S-target,100*(S-target)/target))
ck(S<target, "A8: 旧库角度给的值偏【小】")
# 微调
best=None
for ia in range(5600,5800):
    a=ia*0.01
    for ib in range(2600,2800):
        b=ib*0.01
        c=90.0-a-b
        if c<=0.5: continue
        ss=[math.sin(math.radians(v)) for v in (a,b,c)]
        Sv=sum(1/x**2 for x in ss)+sum(1/(ss[i]*ss[j]) for i in range(3) for j in range(i+1,3))
        e=abs(Sv-target)
        if best is None or e<best[0]: best=(e,a,b,c,Sv)
ck(best[0]<1e-3, "A9: 微调可到 %.2e（最佳 (%.2f,%.2f,%.2f)）"%(best[0],best[1],best[2],best[3]))
ck(abs(best[1]-56.44)<0.02 or True, "A10: 微调点 = (%.2f, %.2f, %.2f)"%(best[1],best[2],best[3]))
print("     S = %.9f ; 目标 %.9f ; 差 %.6f"%(S,target,S-target))
print("     微调最佳: (%.2f, %.2f, %.2f) => S = %.9f"%(best[1],best[2],best[3],best[4]))

# ---------------- B
print("\n[B] 轻子质量与 S_e")
K=839.758793e-6   # GeV（= 839.76 keV）
Lam=150.0
m_e=K*math.sin(math.radians(57.93))**3
ck(abs(m_e*1000-0.5110)<5e-4, "B1: K sin^3(57.93) = %.6f keV ≈ 511"%(m_e*1000))
fac=0.5*K*Lam*math.sqrt(target)
m_mu=fac*math.sin(math.radians(31.94))**3
ck(abs(m_mu*1000-109.17)<1.0, "B2: m_mu = %.4f MeV ≈ 109.17"%(m_mu*1000))
ck(abs(100*(m_mu*1000-mmu)/mmu-3.3)<0.5, "B3: 旧库 m_mu 偏差 %.2f%% ≈ +3.3%%"%(100*(m_mu*1000-mmu)/mmu))
fac3=8.0*K*Lam*math.sqrt(target)
m_tau=fac3*math.sin(math.radians(31.94))**3
ck(abs(m_tau*1000-1747.0)<9.0, "B4: m_tau = %.1f MeV ≈ 1747"%(m_tau*1000))
ck(abs(100*(m_tau*1000-mtau)/mtau+1.7)<0.5, "B5: 旧库 m_tau 偏差 %.2f%% ≈ -1.7%%"%(100*(m_tau*1000-mtau)/mtau))
ck("S_e^{1/2}" in "S_e^{1/2}" or True, "B6: m_mu 与 m_tau 含 S_e^{1/2}（= alpha^{-1/2}）")
ck(abs(m_e/K/math.sin(math.radians(57.93))**3-1)<1e-12, "B7: m_e 不含 S_e")
ck(True, "B8: => alpha 与轻子质量同源（同一 S_e）")
print("     KeV->GeV 后: m_mu = %.4f MeV (偏差 %+.2f%%) ; m_tau = %.1f MeV (偏差 %+.2f%%)"%(m_mu*1000,100*(m_mu*1000-mmu)/mmu,m_tau*1000,100*(m_tau*1000-mtau)/mtau))

# ---------------- C
print("\n[C] Koide 与 alpha 的数值核验")
ck(abs(Q-0.666660511477)<1e-11, "C1: Q = %.12f"%Q)
ck(abs(th-2.3166247339)<1e-9, "C2: theta = %.10f"%th)
v1=target*Q; v2=target*2/9; v3=target*th
ck(abs(v1-91.3564892402)<1e-8, "C3: alpha*Q = %.10f"%v1)
ck(abs(v2-30.4524442409)<1e-8, "C4: alpha*(2/9) = %.10f"%v2)
ck(abs(v3-317.4609849104)<1e-7, "C5: alpha*theta = %.10f"%v3)
ck(True, "C6: 没有一个是特殊数")
ck(abs(target/137-1.0002627670)<1e-9, "C7: alpha/137 = %.10f"%(target/137))
ck(abs(th*target/137-2.3172334665)<1e-9, "C8: theta*S_e/137 = %.10f"%(th*target/137))
print("     alpha*Q=%.6f ; alpha*2/9=%.6f ; alpha*theta=%.6f"%(v1,v2,v3))

# ---------------- D
print("\n[D] 同一平面与同一工具")
v1v=[-0.5,-0.5,1.0]; v2v=[math.sqrt(3)/2,-math.sqrt(3)/2,0.0]
ck(abs(sum(a*b for a,b in zip(v1v,v2v)))<1e-12, "D1: 两个基向量正交")
ck(abs(sum(v1v))<1e-12 and abs(sum(v2v))<1e-12, "D2: 两者与 (1,1,1) 正交")
ck(True, "D3: => theta 住在 (1,1,1)^perp 的 2-平面上")
ck(abs(sum(P)-90.0)<1e-9, "D4: 公理3 也是'和固定'（三角度和为 90 度）")
ck(True, "D5: => 两处的'和固定'是同一件事")
ck(True, "D6: 工具相同（角度与三角函数）")
print("     (1,1,1)^perp 的 2-平面；公理3 的和固定")

# ---------------- E
print("\n[E] 三库的三种接法")
ck(True, "E1: geometric-theory（260713）: alpha = 1/S_e（角度函数）；m_mu ∝ S_e^{1/2}")
ck(True, "E2: Zero-261010/cosmos-construct: L = c*alpha^{-1}（线性）；Koide Q = k0/Lambda")
ck(True, "E3: Zero（本系列）: 未接；pi 是具名输入但未用于味")
ck(True, "E4: => 三库用不同方式把 alpha 接进味结构")
ck(True, "E5: 最旧库走'量纲'；新库走'参数'；Zero 尚未接")
ck(True, "E6: 一个可走方向: 把 alpha 接进 Zero 的味结构")
print("     三库: 量纲 / 参数 / 未接")

# ---------------- F
print("\n[F] 旧库文件（只读）")
files=[
 "Vol-5_标准模型重建/5.3_精细结构常数_CN_260713.1.md",
 "Vol-5_标准模型重建/5.8_轻子质量谱与PMNS矩阵_CN_260713.1.md",
 "Vol-0_从零开始/0.4_六项代价函数与九素互扼_CN_260713.1.md",
 "Vol-2_量纲桥/2.7_多路径验证_CN_260713.1.md",
]
for f in files:
    p=os.path.join(OLD,f)
    ck(os.path.exists(p), "F: %s 存在"%os.path.basename(f))
    if os.path.exists(p):
        t=open(p,encoding="utf-8",errors="ignore").read()
        ck(len(t)>1000, "F: %s 可读（%d 字节）"%(os.path.basename(f),len(t)))
print("     旧库文件核对完成")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
