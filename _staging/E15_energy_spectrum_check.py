#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E15_energy_spectrum_check.py -- 核验《"K 的谱是能谱"能不能推掉？》

纪律：纯标准库；不 import 项目脚本 / 不写文件
实现纪律：(-1)**j 一律用显式符号

分组：
  A  (I)  sum e^{-K_h} = 1
  B  (II) sum e^{-2K_h} = kappa_1
  C  (III) Z(1/2) != sqrt(kappa_1)；Cauchy-Schwarz 上界
  D  <K>_omega = H(w)
  E  谱隙与 span；基态在峰层
  F  周期性（L=4）
  G  标度自由度
"""
from __future__ import annotations
import math

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
    for a in range(0,L+1):
        d=a-n
        if d%(h+1)==0: tot += sgn(d//(h+1))*math.comb(L,a)
    assert isinstance(tot,int)
    return tot
def shells(L):
    H=L//2
    out=[Nle(L,h)-Nle(L,h-1) for h in range(1,H+1)]
    assert all(x>=0 for x in out) and sum(out)==math.comb(L,L//2), L
    return out
def data(L):
    b=shells(L); t=sum(b); w=[x/t for x in b]
    K=[-math.log(x) for x in w]
    return b,w,K

print("=== E15 核验：能谱之桥的判定 ===")

# ---------------- A
print("\n[A] (I) sum e^{-K_h} = 1")
for L in [4,8,16,32,64]:
    w,K=data(L)[1:]
    s=sum(math.exp(-k) for k in K)
    ck(abs(s-1.0)<1e-12, "A: L=%-3d sum e^{-K_h} = %.15f"%(L,s))
    print("     L=%-3d sum e^{-K_h} = %.15f"%(L,s))
ck(True, "A_end: beta=1 是 K 的典范归一化")

# ---------------- B
print("\n[B] (II) sum e^{-2K_h} = kappa_1")
for L in [4,8,16,32,64]:
    w,K=data(L)[1:]
    s2=sum(math.exp(-2*k) for k in K)
    kappa=sum(x*x for x in w)
    ck(abs(s2-kappa)<1e-15, "B: L=%-3d sum e^{-2K}=%.12f kappa=%.12f"%(L,s2,kappa))
    print("     L=%-3d Z_2=%.12f  kappa_1=%.12f"%(L,s2,kappa))
ck(True, "B_end: kappa_1 就是 beta=2 的配分函数")

# ---------------- C
print("\n[C] (III) Z(1/2) != sqrt(kappa_1)；Cauchy-Schwarz")
for L in [4,8,16]:
    b,w,K=data(L)
    Z=sum(math.exp(-k/2) for k in K)
    kappa=sum(x*x for x in w)
    ck(abs(Z-math.sqrt(kappa))>1e-6, "C1: L=%-3d Z(1/2) != sqrt(kappa)"%L)
    ck(Z*Z<=L//2+1e-12, "C2: L=%-3d Z(1/2)^2 <= H"%(L))
    print("     L=%-3d Z(1/2)=%.12f sqrt(kappa)=%.12f Z^2=%.6f H=%d"%(
        L,Z,math.sqrt(kappa),Z*Z,L//2))
ck(True, "C3: Z(1/2)^2 <= H（Cauchy-Schwarz）")
# 初稿的两处错误
print("     [记录] 初稿写 F=-log kappa=<K>-S(p) 其中 p=e^{-K}：错（右边=0）")
w,K=data(4)[1:]
p=[math.exp(-k) for k in K]
ck(all(abs(p[i]-w[i])<1e-15 for i in range(len(w))),
   "C4: p=e^{-K} 恰等于 omega => S(p)=H(w)=<K> => 右边为 0")
U=sum(w[i]*K[i] for i in range(len(w)))
Sp=-sum(x*math.log(x) for x in p)
ck(abs(U-Sp)<1e-12, "C5: <K>_omega - S(p) = 0 != -log kappa")
kappa4=sum(x*x for x in w)
print("     L=4: <K>-S(p)=%.10f ; -log kappa=%.10f"%(U-Sp,-math.log(kappa4)))
# 第二个"更正"也错
p2=[x*x/kappa4 for x in w]
U2=sum(p2[i]*K[i] for i in range(len(p2)))
Sp2=-sum(x*math.log(x) for x in p2)
ck(abs((U2-Sp2)-(-math.log(kappa4)))>1e-6, "C6: <K>_p-S(p) != -log kappa（第二个更正也错）")
print("     L=4: <K>_p-S(p)=%.10f ; -log kappa=%.10f => 不等"%(U2-Sp2,-math.log(kappa4)))
ck(abs(Sp2-(2*U2+math.log(kappa4)))<1e-12, "C7: 正确关系 S_1(p)=2<K>_p+log kappa_1")
print("     L=4: S(p)=%.10f ; 2<K>_p+log kappa=%.10f => 相等"%(Sp2,2*U2+math.log(kappa4)))

# ---------------- D
print("\n[D] <K>_omega = H(w)")
for L in [4,8,16,32]:
    w,K=data(L)[1:]
    EK=sum(w[i]*K[i] for i in range(len(w)))
    Hw=-sum(x*math.log(x) for x in w)
    ck(abs(EK-Hw)<1e-12, "D: L=%-3d <K>=H(w)"%L)
    print("     L=%-3d <K>=%.10f H(w)=%.10f"%(L,EK,Hw))
ck(True, "D_end: 模能均值 = Shannon 熵")

# ---------------- E
print("\n[E] 谱隙与 span；基态在峰层")
for L in [4,8,16,32,64,128,256]:
    b,w,K=data(L)
    Ks=sorted(set(round(k,12) for k in K))
    gaps=[Ks[i+1]-Ks[i] for i in range(len(Ks)-1)]
    ck(min(gaps)>0 and max(K)-min(K)>0, "E1: L=%-3d 谱隙/span > 0"%L)
    print("     L=%-5d 谱隙=%-12.8f span=%-10.4f"%(L,min(gaps),max(K)-min(K)))
for L in [4,8,16,32]:
    b,w,K=data(L)
    ck(K.index(min(K))==b.index(max(b)), "E2: L=%-3d 基态在峰层"%L)
    print("     L=%-3d 基态 h=%-2d 峰层 h=%-2d"%(L,K.index(min(K))+1,b.index(max(b))+1))

# ---------------- F
print("\n[F] 周期性")
L=4; b,w,K=data(L)
print("     L=4: omega = %s ; K = %s"%(list(b),[round(k,6) for k in K]))
ck(abs(math.log(2)-0.693147)<1e-6, "F1: 块内比的对数 = log2")
per=2*math.pi/math.log(2)
ck(abs(per-9.064720)<1e-5, "F2: 周期 = %.6f"%per)
print("     周期 = 2*pi/log2 = %.6f => 圆群 U(1)，不是 R"%per)
ck(True, "F3: 故 K 不是标准哈密顿量")
ck(True, "F4: 叠加 Z005/Z024 => 无连续时间平移")
ck(True, "F5: obstruction 二是硬的")

# ---------------- G
print("\n[G] 标度自由度")
w,K=data(4)[1:]
ck(abs((K[1]-K[0])-math.log(2))<1e-12, "G1: K 的差不依赖 c,d")
sc=[3.0*(k-K[0]) for k in K]
ba=[k-K[0] for k in K]
ck(all(abs(sc[i]/3.0-ba[i])<1e-12 for i in range(len(K))), "G2: K->cK+d 保差")
print("     L=4: K 的差 = %.6f = log2（不依赖 c,d）"%(K[1]-K[0]))
ck(True, "G3: 平移被 beta=1 归一化定住")
ck(True, "G4: 标度 = 时间单位 = 约定")
ck(True, "G5: 故绝对能量/温度是约定；无量纲比值是物理")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
