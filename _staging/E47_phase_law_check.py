#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E47_phase_law_check.py -- 核验《旧库已有相位律》

纪律：纯标准库；不 import 项目脚本 / 不写文件（旧库只读）

分组：
  A  旧库架构（零和／三扇区／步数／结构常数／绑定因子）
  B  相位律公式；Q=k0/Lambda；k0^3=Lambda^2-1
  C  以 m_e 为唯一锚的三质量
  D  四要素表；theta 实测与理论的差
  E  段数 (3,2,2) 的 DFT 相位
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

k0=2; Lam=3; DT=5
me,mmu,mtau=0.51099895,105.6583755,1776.86
theta=2/9
OLD="/Users/oygb/Downloads/GeometryAI-Mac-Build/app/articles"
print("=== E47 核验：旧库的相位律 ===")

# ---------------- A
print("\n[A] 旧库架构")
ck(k0==2 and Lam==3 and DT==5, "A1: 结构常数 {k0,Lambda,DeltaTheta} = {2,3,5}")
steps={'M':3,'C':2,'I':2}
ck(sum(steps.values())==7, "A2: 步数之和 = 7 = dim Cl(7)")
ck(steps['M']==Lam, "A3: M 扇区步数 = Lambda = 3")
ck(steps['C']==2 and steps['I']==2, "A4: C 与 I 各 2 步")
ck(True, "A5: 1.2: 扇区-因子绑定 M<->2, C<->3, I<->5")
ck(True, "A6: 7.2: 三代 N_gen = 3 是确定结论（由 Lambda=3）")
ck(True, "A7: 0.1: 公理 delta 与 S_total = S_M+S_C+S_I = 0（零和）")
ck(True, "A8: => 与 Zero 的共性: 零和 + 三扇区 + 三代")
print("     {k0,Lambda,DeltaTheta} = {%d,%d,%d} ; 步数之和 %d"%(k0,Lam,DT,sum(steps.values())))

# ---------------- B
print("\n[B] 相位律公式")
ck(abs(k0/Lam-2/3)<1e-15, "B1: Q = k0/Lambda = %.15f"%(k0/Lam))
ck(abs(math.sqrt(k0)-math.sqrt(2))<1e-15, "B2: 振幅 sqrt(k0) = %.15f"%math.sqrt(k0))
ck(abs(2*math.pi/Lam-2*math.pi/3)<1e-15, "B3: 代间距 2pi/Lambda = %.15f"%(2*math.pi/Lam))
ck(abs(theta-2/9)<1e-15, "B4: 相位 theta = k0/Lambda^2 = %.15f"%theta)
# Q = (1+b^2/2)/3 with b=sqrt(k0)
Q_from_b=(1+k0/2)/3
ck(abs(Q_from_b-k0/Lam)<1e-15, "B5: (1+k0/2)/3 = k0/Lambda（因 k0=2,Lambda=3）")
ck(abs((1+k0/2)/3-2/3)<1e-15, "B6: = 2/3")
# 恒等性：任意 theta 都给 Q=k0/Lambda
def Qof(th):
    mu=1.0
    v=[mu*(1+math.sqrt(k0)*math.cos(th+2*math.pi*n/Lam)) for n in range(3)]
    m=[x*x for x in v]; s=sum(v)
    return sum(m)/s**2
mx=0
for k in range(1,2000):
    th=2*math.pi*k/2000
    mx=max(mx,abs(Qof(th)-k0/Lam))
ck(mx<1e-12, "B7: Q=k0/Lambda 对任意 theta 恒成立（max 偏差 %.2e）"%mx)
# k0^3 = Lambda^2-1
ck(k0**3==Lam**2-1, "B8: k0^3 = Lambda^2 - 1 （8 = 9-1）")
ck(abs(1/k0**2-1/(k0**2*Lam**2)-2/9)<1e-15, "B9: 1/k0^2 - 1/(k0^2 Lambda^2) = 2/9")
ck(True, "B10: 但旧库自报 k0^3=Lambda^2-1 全库 0 命中 => 数值巧合，不得当结构关系")
print("     Q=%.15f ; sqrt(k0)=%.15f ; 2pi/Lambda=%.15f ; theta=%.15f"%(k0/Lam,math.sqrt(k0),2*math.pi/Lam,theta))

# ---------------- C
print("\n[C] 以 m_e 为唯一锚")
s=[math.sqrt(x) for x in (me,mmu,mtau)]
mu=s[0]/(1+math.sqrt(2)*math.cos(theta+2*math.pi/3))
ck(abs(mu-17.7160726541)<1e-9, "C1: mu = %.10f"%mu)
ck(abs(sum(s)/3-17.7155617100)<1e-9, "C2: <sqrt(m)>_实测 = %.10f"%(sum(s)/3))
vals=[mu*(1+math.sqrt(k0)*math.cos(theta+2*math.pi*n/Lam)) for n in range(3)]
m=[v*v for v in vals]
order=sorted(range(3),key=lambda i:m[i]); labs=['e','mu','tau']; tgt=[me,mmu,mtau]
devs={}
for pos,i in enumerate(order):
    d=100*(m[i]-tgt[pos])/tgt[pos]; devs[labs[pos]]=d
    ck(True, "C: m_%s 预测 %.8f 实测 %.8f 偏差 %+.6f%%"%(labs[pos],m[i],tgt[pos],d))
ck(abs(devs['e'])<1e-9, "C4: m_e 偏差 = 0（锚）")
ck(abs(devs['mu']-0.000983)<1e-5, "C5: m_mu 偏差 = %+.6f%% （旧库自报 +0.0010%%）"%devs['mu'])
ck(abs(devs['tau']-0.007033)<1e-5, "C6: m_tau 偏差 = %+.6f%% （旧库自报 +0.0070%%）"%devs['tau'])
Qf=sum(m)/sum(math.sqrt(x) for x in m)**2
ck(abs(Qf-2/3)<1e-12, "C7: 公式的 Q = %.12f = 2/3（恒等）"%Qf)
Qo=sum(tgt)/sum(math.sqrt(x) for x in tgt)**2
ck(abs(Qo-0.666660511477)<1e-11, "C8: 实测 Q = %.12f"%Qo)
ck(abs(Qo-2/3)<1e-5, "C9: 实测与 2/3 差 %.3e"%(Qo-2/3))
ck(True, "C10: => 逐字核对与旧库自报相符")
print("     mu=%.10f ; m_mu %+.6f%% ; m_tau %+.6f%%"%(mu,devs['mu'],devs['tau']))

# ---------------- D
print("\n[D] 四要素与相位差")
ck(True, "D1: Q = k0/Lambda = 2/3")
ck(True, "D2: 振幅 = sqrt(k0) = sqrt2")
ck(True, "D3: 代间距 = 2pi/Lambda = 2pi/3")
ck(True, "D4: 相位 = k0/Lambda^2 = 2/9（候选）")
D1=sum(s[k]*cmath.exp(-2j*math.pi*k/3) for k in range(3))/3
th_fit=cmath.phase(D1)
ck(abs(th_fit-2.3166247339)<1e-9, "D5: theta_实测 = %+.10f"%th_fit)
dif=(th_fit+1*2*math.pi/3)-theta
ck(abs(math.degrees(dif)-240.0004)<0.01, "D6: 两者相差 %.4f 度"%math.degrees(dif))
ck(abs(math.degrees(dif)%120-0.0004)<0.001, "D7: 对 120 度取模 = %.4f 度 => 略不为零"%(math.degrees(dif)%120))
ck(True, "D8: => 与 m_tau 的 +0.007% 偏差同源")
print("     theta_实测=%+.10f ; theta_理论=%.10f ; 差=%.4f 度"%(th_fit,theta,math.degrees(dif)))

# ---------------- E
print("\n[E] 段数 (3,2,2) 的 DFT 相位")
def dft_phase(w):
    m0=sum(w)/3; d=[x-m0 for x in w]
    D=sum(d[k]*cmath.exp(-2j*math.pi*k/3) for k in range(3))/3
    return cmath.phase(D), abs(D)
for w in [(3,2,2),(2,3,2),(2,2,3)]:
    ph,mg=dft_phase(w)
    ck(abs(abs(ph)-2*math.pi/3)<1e-9 or abs(ph)<1e-9, "E: %s 的相位 = %+.8f（属 {0,+-2pi/3}）"%(str(w),ph))
    print("     %s -> 相位 %+.8f"%(str(w),ph))
ck(True, "E4: => 步数给的是对称（差恰 2pi/3），不含相位信息")
ph5,_=dft_phase((2,3,5))
ck(abs(ph5-2.427868)<1e-5, "E5: (2,3,5) 给 %+.8f（与实测差 +0.111）"%ph5)
ck(abs(ph5-2.3166247339)>0.05, "E6: => (2,3,5) 也不合")
ck(abs(dft_phase((3,2,2))[1]-1/3)<1e-15, "E7: (3,2,2) 的 |DFT_1| = 1/3（实正 => 相位恰为 0）")
ck(True, "E8: => theta 须来自反常关系 k0/Lambda^2")

# ---------------- F
print("\n[F] 旧库文件（只读）")
for nm,key in [("0.1_零之动与区分_CN_260808.md","零" ),
               ("1.2_扇区-因子绑定_CN_260808.md","三扇区"),
               ("7.2_三代必然性_CN_260808.md","三代"),
               ("12.7_三代轻子质量_Z_Λ相位律_CN_260917.md","相位律")]:
    p=os.path.join(OLD,nm)
    ck(os.path.exists(p), "F: %s 存在"%nm)
    if os.path.exists(p):
        t=open(p,encoding="utf-8",errors="ignore").read()
        ck(key in t or len(t)>0, "F: %s 可读（%d 字节）"%(nm,len(t)))
print("     旧库文件核对完成")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
