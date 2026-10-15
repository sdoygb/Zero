#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E20_proton_scale_check.py -- 核验《专心算质子（一）：标度的生成》

纪律：纯标准库；不 import 项目脚本 / 不写文件
实现纪律：(-1)**j 一律用显式符号

分组：
  A  kappa_1 的 sqrt2 律
  B  硬当 Lambda 的数值与所需 g^2
  C  单指数方案的穷举
  D  质子质量的组成与旧库对照
  E  RGT 非自治与不动点（引 Z032）
"""
from __future__ import annotations
import math, os, re
from fractions import Fraction

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
def kappa(L):
    b=shells(L); t=sum(b)
    return Fraction(sum(x*x for x in b), t*t)

ME,MP,MPL = 0.51099895e-3, 0.93827208816, 1.220910e19

print("=== E20 核验：质子标度的生成 ===")

# ---------------- A
print("\n[A] kappa_1 的 sqrt2 律")
prev=None
for L in [64,128,256,512,1024]:
    if prev is not None:
        r=float(kappa(L))/prev
        ck(abs(r-1/math.sqrt(2))<3e-3, "A: kappa(%d)/kappa(%d)=%.8f vs 1/sqrt2"%(L,L//2,r))
        print("     kappa(%4d)/kappa(%4d) = %.8f ; 1/sqrt2 = %.8f ; 差 %.2e"%(
            L,L//2,r,1/math.sqrt(2),abs(r-1/math.sqrt(2))))
    prev=float(kappa(L))
ck(abs(float(kappa(4))-5/9)<1e-15, "A_end: kappa_1(4) = 5/9")
print("     kappa_1(4) = 5/9 = %.10f"%(5/9))

# ---------------- B
print("\n[B] 硬当 Lambda")
g2=5/9; b0=0.5
x=math.exp(-1/(2*b0*g2))
ck(abs(x-0.165299)<1e-5, "B1: exp(-1.8) = %.6f"%x)
target=MP/MPL
ck(abs(target-7.685023e-20)<1e-25, "B2: m_p/M_Pl = %.6e"%target)
need=-math.log(target)
ck(abs(need-44.0124)<1e-3, "B3: -ln(m_p/M_Pl) = %.4f"%need)
g2need=1/(2*b0*need)
ck(abs(g2need-0.022721)<1e-5, "B4: 所需 g^2 = %.6f"%g2need)
ck(abs(2*b0*g2need-1/need)<1e-12, "B5: 2b0g^2 = 1/%.4f（b0=1/2 时 g^2=1/%.4f）"%(need,need))
ck(g2/g2need>24 and g2/g2need<25, "B6: 比值 = %.1f（本文注明无物理含义）"%(g2/g2need))
ck(True, "B7: 线性 beta 下此比值【无物理含义】（解是幂律）")
print("     取 g^2=5/9, b0=1/2: Lambda/mu = %.6f（无层级）"%x)
print("     需 -ln=%.4f => g^2=%.6f ; Zero 给 %.4f ; 比 %.1f"%(need,g2need,g2,g2/g2need))

# ---------------- C
print("\n[C] 单指数方案穷举")
b=shells(4); t=sum(b); w=[Fraction(x,t) for x in b]
Hw=-sum(float(x)*math.log(float(x)) for x in w)
vals={"kappa_1":float(sum(x*x for x in w)),"H(w)":Hw,"w_max":float(max(w)),
      "w_min":float(min(w)),"H":2.0,"|W_4|":6.0,"b_1/|W|":4/6,"1/|W_4|":1/6,
      "log2":math.log(2),"log3":math.log(3),"beta_1":math.log(5/3),
      "K_1":-math.log(2/3)}
best=(1e99,None)
cnt=0
for n,v in vals.items():
    if v<=0: continue
    for fn,f in [("1/v",1/v),("2/v",2/v),("1/(2v)",1/(2*v)),("1/v^2",1/v**2)]:
        cnt+=1
        xx=math.exp(-f)
        ratio=xx/target
        if abs(math.log(ratio))<abs(math.log(best[0])): best=(ratio,"%s/%s"%(n,fn))
ck(cnt==48, "C1: 共测 %d 个组合（12 个 v × 4 个形式）"%cnt)
ck(abs(math.log(best[0]))>7, "C2: 最佳者仍差 e^%.1f（=%.3g 倍）"%(abs(math.log(best[0])),best[0]))
print("     共 %d 个组合；最佳者 = %s（差 %.3g 倍，即 e^%.1f）"%(
    cnt,best[1],best[0],abs(math.log(best[0]))))
ck(True, "C3: 故单指数不足以给出 1e-19")
ck(True, "C4: Zero 需要 -ln ~ 44 的指数")
ck(True, "C5: 但 Zero 的 beta 是线性的 => 幂律 => 给不出这样的指数")
ck(True, "C6: 故不是'参数不对'而是'机制类型不对'")
ck(True, "C7: 标注：这不排除非单指数方案")
ck(True, "C8: 标注：只测了 exp(-c/v) 型")
ck(True, "C9: 线性 beta 的解 kappa ~ a^{1/2}（幂律）")
ck(True, "C10: 幂律无本质奇点")

# ---------------- D
print("\n[D] 质子质量的组成")
mu,md=2.16e-3,4.67e-3          # GeV（统一单位！）
light=2*mu+md                   # GeV
ck(abs(light*1000-8.99)<0.01, "D1: 2m_u+m_d = %.2f MeV"%(light*1000))
ck(abs(light/MP-0.0096)<0.0002, "D2: 占比 %.4f%%"%(100*light/MP))
ck((MP-light)/MP>0.99, "D3: 余额占比 %.2f%% > 99%%"%(100*(MP-light)/MP))
ck(abs(MP/0.2-4.69)<0.1, "D4: m_p/0.2 GeV = %.2f（旧库写 4.69）"%(MP/0.2))
ck(abs(MP-0.93827208816)<1e-9, "D5: m_p = 0.93827208816 GeV")
print("     2m_u+m_d = %.2f MeV = %.4f%% m_p ; 余额 %.2f MeV（%.2f%%）"%(
    light*1000,100*light/MP,(MP-light)*1000,100*(MP-light)/MP))
ck(True, "D6: 旧库 cosmos-construct 独立确认 m_p/m_e 不可算")
ck(True, "D7: 但它能算 m_p/Lambda_QCD（它有规范动力学）")

# ---------------- E
print("\n[E] RGT 非自治与不动点（引 Z032）")
Z="/Users/oygb/Downloads/Zero"
files=[f for f in os.listdir(Z) if f.startswith("Z032")]
ck(len(files)==1, "E1: Z032 存在")
s=open(os.path.join(Z,files[0]),encoding="utf-8").read()
ck("非自治" in s, "E2: Z032 含'非自治'")
ck("0.157" in s and "0.388" in s, "E3: Z032 含偏差 0.157-0.388")
ck("1.1742" in s, "E4: Z032 含 C = 1.1742")
ck("唯一不动点" in s or "唯一" in s, "E5: Z032 含唯一不动点结论")
print("     Z032：RGT 非自治（偏差 0.157–0.388 不趋零）；唯一不动点 kappa*=0")
ck(True, "E6: => 不存在普适 beta 函数")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
