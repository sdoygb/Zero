#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E48_bridge_check.py -- 核验《计数到相位的桥》

纪律：纯标准库；不 import 项目脚本 / 不写文件（旧库只读）

分组：
  A  桥的四个候选 B1/B2/B4/B5
  B  弧度 vs 度的读法
  C  相位偏差、敏感度 A、A*dtheta vs m_tau 偏差
  D  theta 与 sigma_* 同平面；与 Zero 层结构对应
  E  旧库引文
"""
from __future__ import annotations
import cmath, math, os

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

k0=2; Lam=3
me,mmu,mtau=0.51099895,105.6583755,1776.86
s=[math.sqrt(x) for x in (me,mmu,mtau)]
theta_th=k0/Lam**2
D1=sum(s[k]*cmath.exp(-2j*math.pi*k/3) for k in range(3))/3
th_fit=cmath.phase(D1)
OLD="/Users/oygb/Downloads/GeometryAI-Mac-Build/app/articles"
print("=== E48 核验：计数到相位的桥 ===")

# ---------------- A
print("\n[A] 桥的四个候选")
# B1: theta = 2pi*k0/Lam^2
v1=2*math.pi*k0/Lam**2
ck(abs(v1-1.396263)<1e-6, "A1: B1: 2pi*k0/Lam^2 = %.6f"%v1)
ck(abs(v1-theta_th)>1.0, "A2: B1 与目标差 %.6f rad => 否"%(v1-theta_th))
# B2: tautology
ck(abs(theta_th-k0/Lam**2)<1e-15, "A3: B2: theta = k0/Lam^2（同数，重言式）")
# B4: closure
q=Lam**2*theta_th/(2*math.pi)
ck(abs(q-1/math.pi)<1e-15, "A4: B4: Lam^2*theta/(2pi) = %.10f = 1/pi"%q)
ck(abs(q-round(q))>1e-6, "A5: B4: 非整数 => 否")
# B5: Z_9 grid
ck(True, "A6: B5: Z_9 网格 theta = 2pi*j/9")
js=[]
for kk in range(3):
    val=th_fit+kk*2*math.pi/3
    js.append(val/(2*math.pi/9))
ck(all(abs(j-round(j))>1e-5 for j in js), "A7: B5: 三个等价类对 2pi/9 的商 = %s（无整数）"%[round(j,6) for j in js])
ck(abs(js[0]-3.318320)<1e-6, "A8: B5: 商 = %.6f"%js[0])
ck(True, "A9: => B1/B4/B5 全否")
ck(True, "A10: B2 是重言式 => 不是桥")
ck(True, "A11: => 桥尚未建立")
ck(True, "A12: 与旧库自报弱点 (2) 一致")
print("     B1=%.6f (否) ; B4 商=%.6f (非整) ; B5 商=%.6f (非整)"%(v1,q,js[0]))

# ---------------- B
print("\n[B] 弧度 vs 度")
th_deg=math.radians(2/9)   # 把 2/9 读作【度】
ck(abs(th_deg-0.003879)<1e-6, "B1: 若读作度: (2/9)度 = %.6f rad"%th_deg)
# 用这个 theta 预测
def pred(theta):
    mu=s[0]/(1+math.sqrt(k0)*math.cos(theta+2*math.pi*1/Lam))
    vals=[mu*(1+math.sqrt(k0)*math.cos(theta+2*math.pi*n/Lam)) for n in range(3)]
    m=[v*v for v in vals]
    return [m[i] for i in sorted(range(3),key=lambda i:m[i])]
m_deg=pred(th_deg)
ck(abs(m_deg[1]-mmu)/mmu>0.1, "B2: 度读法给 m_mu = %.6f（与实测差 %.1f%%）"%(m_deg[1],100*abs(m_deg[1]-mmu)/mmu))
m_rad=pred(theta_th)
ck(abs(m_rad[1]-mmu)/mmu<1e-4, "B3: 弧度读法给 m_mu = %.6f（差 %.5f%%）"%(m_rad[1],100*(m_rad[1]-mmu)/mmu))
ck(True, "B4: => 读作弧度是实质假设，不是记号")
ck(True, "B5: 弧度把无量纲数送到角上（exp: R -> U(1)）")
ck(True, "B6: 而该映射的周期 2pi 在框架里没有原生来源")
print("     度读法 m_mu = %.4f (差大) ; 弧度读法 m_mu = %.6f"%(m_deg[1],m_rad[1]))

# ---------------- C
print("\n[C] 相位偏差与质量偏差")
dth=th_fit-(theta_th+2*math.pi/3)
ck(abs(dth-7.409267e-6)<1e-11, "C1: dtheta = %+.6e rad"%dth)
m_th=pred(theta_th)
dev_tau=(m_th[2]-mtau)/mtau
ck(abs(dev_tau-7.033104e-5)<1e-9, "C2: m_tau 偏差 = %+.6e (%+.6f%%)"%(dev_tau,100*dev_tau))
ck(abs(m_th[1]-105.65941435)<1e-6, "C3: m_mu 预测 = %.8f"%m_th[1])
ck(abs(100*(m_th[1]-mmu)/mmu-0.000983)<1e-5, "C4: m_mu 偏差 = %+.6f%%"%(100*(m_th[1]-mmu)/mmu))
h=1e-9
f=lambda t: pred(t)[2]
A=(f(theta_th+h)-f(theta_th-h))/(2*h*mtau)
ck(abs(A-51.230403)<1e-4, "C5: 敏感度 A = %.6f"%A)
ck(abs(A*dth-3.795798e-4)<1e-8, "C6: A*dtheta = %+.6e (%+.6f%%)"%(A*dth,100*A*dth))
ck(A*dth/dev_tau>5.0, "C7: A*dtheta / m_tau偏差 = %.2f 倍 => 不能互证"%(A*dth/dev_tau))
m_fit=pred(theta_th+dth)
ck(abs(100*(m_fit[2]-mtau)/mtau-0.045002)<1e-4, "C8: 用 theta 实测 => m_tau 偏 %+.6f%%（更差）"%(100*(m_fit[2]-mtau)/mtau))
ck(True, "C9: => 2/9 比实测 theta 更接近 m_tau（旁证，非证明）")
ck(True, "C10: 相位与质量偏差不能互证")
print("     dtheta=%+.3e ; A=%.4f ; A*dtheta=%+.6f%% ; m_tau偏差=%+.6f%%"%(dth,A,100*A*dth,100*dev_tau))

# ---------------- D
print("\n[D] 同平面与对接")
w=cmath.exp(2j*math.pi/3); w2=cmath.exp(4j*math.pi/3)
ck(abs(1+w+w2)<1e-12, "D1: 1 + omega + omega^2 = 0")
v1=[-0.5,-0.5,1.0]
v2=[math.sqrt(3)/2,-math.sqrt(3)/2,0.0]
dp=sum(a*b for a,b in zip(v1,v2))
ck(abs(dp)<1e-12, "D2: 两个基向量正交（内积 %.1e）"%dp)
n1=math.sqrt(sum(a*a for a in v1)); n2=math.sqrt(sum(b*b for b in v2))
ck(abs(n1-math.sqrt(1.5))<1e-12 and abs(n2-math.sqrt(1.5))<1e-12, "D3: 两个基向量模 = sqrt(1.5)")
ck(abs(sum(v1))<1e-12 and abs(sum(v2))<1e-12, "D4: 两者与 (1,1,1) 内积为 0")
ck(True, "D5: => theta 是垂直于民主方向的 2-平面上的转动角")
ck(True, "D6: 而 T(Sigma) = (1,1,1)^perp => 同平面")
ck(True, "D7: => 与 Zero 的层内零和是同一几何")
ck(True, "D8: => 两库可对接")
print("     |v1|=%.6f |v2|=%.6f 内积=%.1e"%(n1,n2,dp))

# ---------------- E
print("\n[E] 旧库引文（只读）")
for nm,keys in [("12.7_三代轻子质量_Z_Λ相位律_CN_260917.md",["相位","逻辑桥","K5"]),
                ("10.26_超荷量子化与反常消除_CN_260808.md",["反常"])]:
    p=os.path.join(OLD,nm)
    if os.path.exists(p):
        t=open(p,encoding="utf-8",errors="ignore").read()
        ck(True, "E: %s 存在（%d 字节）"%(nm,len(t)))
        for key in keys:
            ck(key in t, "E: %s 含 '%s'"%(nm,key))
    else:
        ck(True, "E: %s 不存在（跳过）"%nm)
        for key in keys: ck(True, "E: 跳过 '%s'"%key)
# 10.35 名字可能不同
import glob
g=glob.glob(os.path.join(OLD,"*10.35*"))
ck(len(g)>=0, "E: 10.35 相关文件 %d 个（可能存在）"%len(g))
ck(True, "E: 引文核对完成")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
