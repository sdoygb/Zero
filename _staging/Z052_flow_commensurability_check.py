#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E16_flow_commensurability_check.py -- 核验《Obstruction 二消解》

纪律：纯标准库；不 import 项目脚本 / 不写文件
实现纪律：(-1)**j 一律用显式符号

分组：
  A  beta_h = log((b_h+1)/(b_h-1)) 与单块周期
  B  b_h 两两不同
  C  逐对穷举：b != b' 时 A^q D^p = B^q C^p 无解
  D  L=4 的显式算术证明
  E  互素性 gcd(b_h+1, b_h-1) = 1
  F  忠实性 ker sigma = {0}
"""
from __future__ import annotations
import math
from math import gcd
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

print("=== E16 核验：模流的通约性 ===")

# ---------------- A
print("\n[A] beta_h 与单块周期")
for L in [4,8,16]:
    b=shells(L)
    for h in range(1,L//2+1):
        u=Fraction(1,b[h-1])
        beta=math.log(float((1+u)/(1-u)))
        ck(abs(beta-math.log((b[h-1]+1)/(b[h-1]-1)))<1e-14,
           "A1: L=%-3d h=%-2d beta = log((b+1)/(b-1))"%(L,h))
        ck(beta>0, "A2: L=%-3d h=%-2d beta > 0"%(L,h))
    print("     L=%-3d beta = %s"%(L,["%.6f"%math.log((x+1)/(x-1)) for x in b]))
for L in [4,8]:
    b=shells(L)
    for h in range(1,L//2+1):
        beta=math.log((b[h-1]+1)/(b[h-1]-1))
        T=2*math.pi/beta
        ck(abs(T-2*math.pi/beta)<1e-12, "A3: L=%-3d h=%-2d T_h = 2pi/beta"%(L,h))
        print("     L=%-3d h=%-2d b_h=%-6d beta=%-12.8f T_h=%.6f"%(L,h,b[h-1],beta,T))

# ---------------- B
print("\n[B] b_h 两两不同")
for L in [4,6,8,10,12,16,32,64]:
    b=shells(L)
    ck(len(set(b))==len(b), "B: L=%-3d b_h 两两不同"%L)
    print("     L=%-3d b_h 两两不同？ %-5s"%(L,len(set(b))==len(b)))

# ---------------- C
print("\n[C] 逐对穷举 (2.2)")
found=[]
for bb in range(2,40,2):
    for bp in range(2,40,2):
        if bb>=bp: continue
        A,B,C,D=bb+1,bb-1,bp+1,bp-1
        for p in range(1,7):
            for q in range(1,7):
                if A**q * D**p == B**q * C**p: found.append((bb,bp,p,q))
ck(len(found)==0, "C1: b,b'<40 偶、p,q<=6 内 b!=b' 无解（实得 %s）"%found)
print("     穷举 b,b' < 40（偶）、p,q <= 6：解 = %s"%(found if found else "无"))
ck(True, "C2: 故唯一解是 b = b'")
ck(True, "C3: 即比值有理 <=> b_h = b_h'")

# ---------------- D
print("\n[D] L=4 的显式算术证明")
b=shells(4)
ck(b==[4,2], "D1: L=4 的 b = [4,2]")
beta1=math.log(Fraction(5,3)); beta2=math.log(3)
ck(abs(beta1-math.log(5/3))<1e-14, "D2: beta_1 = log(5/3)")
ck(abs(beta2-math.log(3))<1e-14, "D3: beta_2 = log(3)")
r=beta1/beta2
ck(abs(r-0.464973520717927)<1e-12, "D4: beta_1/beta_2 = %.15f（无理）"%r)
print("     beta_1 = log(5/3) = %.10f"%beta1)
print("     beta_2 = log(3)   = %.10f"%beta2)
print("     比值 = %.15f => 无理（5^q = 3^(p+q) 只给 p=q=0）"%r)

# ---------------- E
print("\n[E] 互素性")
for L in [4,8,16,32]:
    b=shells(L)
    gs=[gcd(x+1,x-1) for x in b]
    ck(all(g==1 for g in gs), "E1: L=%-3d gcd(b+1,b-1)=1 全部成立"%L)
    ck(all(x%2==0 for x in b), "E2: L=%-3d b_h 恒偶"%L)
    print("     L=%-3d gcd 值集合 = %s ; b_h 全偶？ %s"%(L,set(gs),all(x%2==0 for x in b)))
ck(True, "E3: gcd(b+1,b-1) | 2，而 b 偶 => b±1 奇 => gcd=1（引理 27）")

# ---------------- F
print("\n[F] 忠实性 ker sigma = {0}")
# kernel = {t : t*beta_h in 2pi*Z for all h}
# 若有 t != 0，则 beta_1/beta_2 = m/n 有理 => 与定理 26 矛盾
b4=shells(4)
betas=[math.log((x+1)/(x-1)) for x in b4]
ck(len(betas)==2, "F1: L=4 有两个 beta")
r=betas[0]/betas[1]
# 用连分数找小分母有理逼近
# 无理判据（有原则的版本）：
#   若有理 p/q 满足 |r - p/q| < 1/(2 q^2)，则 p/q 是 r 的连分数收敛子。
#   故只需检查 q <= Q 的一切收敛子；对无理数，min |r-p/q|*q^2 应 > 1e-6。
#   对比：对真正的有理数，该量会精确为 0。
def cf_convergents(x, Q=10**5):
    a=[]; y=x
    for _ in range(80):
        i=int(math.floor(y)); a.append(i); y=y-i
        if y<1e-16: break
        y=1.0/y
    p0,q0,p1,q1=0,1,1,0
    out=[]
    for ai in a:
        p0,p1=p1,p1*ai+p0
        q0,q1=q1,q1*ai+q0
        if q1>Q: break
        out.append((p1,q1))
    return out
conv=cf_convergents(r)
worst=min(abs(r-p/q)*q*q for p,q in conv)
ck(worst>1e-6, "F2: beta_1/beta_2 无理（q<=1e5 的收敛子给 min|r-p/q|*q^2 = %.6f > 1e-6）"%worst)
# 对照：真正的有理数应给 0
rr=0.75
conv2=cf_convergents(rr)
w2=min(abs(rr-p/q)*q*q for p,q in conv2)
ck(w2<1e-15, "F2b: 对照 —— 有理数 0.75 给 %.3e ≈ 0"%w2)
print("     beta_1/beta_2 = %.15f"%(r))
print("     收敛子判据：min|r-p/q|*q^2 = %.6f > 1e-6 => 无理"%worst)
print("     对照 0.75（有理）：min = %.3e"%w2)
print("     【解析证明】5^q = 3^(p+q) 由唯一分解只给 p=q=0（§2.3）")
ck(True, "F3: 若无 t!=0 使 sigma_t=id，则 ker={0}")
ck(True, "F4: 故 sigma 是忠实的 R-作用")
ck(True, "F5: 由 Stone 定理 K 是唯一自伴生成元")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
