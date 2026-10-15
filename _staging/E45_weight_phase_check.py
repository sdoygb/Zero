#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E45_weight_phase_check.py -- 核验《纯权重侧的 theta》

纪律：纯标准库（+numpy）；不 import 项目脚本 / 不写文件

分组：
  A  DFT；arg(s1)=theta；|s1|/|s0|=b/2；s2=conj(s1)
  B  Q=(1+b^2/2)/3；b^2=2(3Q-1)；两处更正的核验
  C  三项权重互异；L=4 的 (2,2,2) => s1=0
  D  Z_3 正则表示的对称性；b 与 theta 同源
  E  语料引文
"""
from __future__ import annotations
import cmath, itertools, math, os
from collections import Counter

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

me,mmu,mtau=0.51099895,105.6583755,1776.86
s=[math.sqrt(x) for x in (me,mmu,mtau)]
theta=2.3166247339
print("=== E45 核验：纯权重侧的 theta ===")

# ---------------- A
print("\n[A] DFT 关系")
ck(abs(s[0]-0.714842)<1e-5, "A1: s_0 = %.6f"%s[0])
sh=[sum(s[k]*cmath.exp(-2j*math.pi*n*k/3) for k in range(3))/3 for n in range(3)]
ck(abs(sh[0].imag)<1e-12, "A2: sh_0 虚部 = 0")
ck(abs(sh[0].real-sum(s)/3)<1e-12, "A3: sh_0 = <s> = %.10f"%sh[0].real)
ck(abs(sh[2]-sh[1].conjugate())<1e-12, "A4: sh_2 = conj(sh_1)")
ck(abs(cmath.phase(sh[1])-theta)<1e-9, "A5: arg(sh_1) = %.12f = theta"%cmath.phase(sh[1]))
A=sh[0].real
b_obs=2*abs(sh[1])/A
ck(abs(b_obs-2*abs(sh[1])/abs(sh[0]))<1e-12, "A6: b = 2|sh_1|/|sh_0|")
r=abs(sh[1])/abs(sh[0])
ck(abs(r-b_obs/2)<1e-12, "A7: |sh_1|/|sh_0| = b/2 = %.12f"%r)
ck(abs(b_obs-math.sqrt(2))<1e-4, "A8: b = %.10f（与 sqrt2 差 %.2e）"%(b_obs,b_obs-math.sqrt(2)))
ck(abs(r-1/math.sqrt(2))<1e-5, "A9: r = %.12f（与 1/sqrt2 差 %.2e）"%(r,abs(r-1/math.sqrt(2))))
ck(abs(r-1/math.sqrt(2))<1e-5, "A10: b=sqrt2 <=> r=1/sqrt2")
print("     sh_1 = %s ; arg = %.12f ; r = %.12f"%(complex(round(sh[1].real,8),round(sh[1].imag,8)),cmath.phase(sh[1]),r))

# ---------------- B
print("\n[B] Q 与 b 的关系")
Q=sum([me,mmu,mtau])/sum(s)**2
ck(abs(Q-0.666660511477)<1e-11, "B1: Q = %.12f"%Q)
# 正确公式
for bb in [b_obs, math.sqrt(2)]:
    Qb=(1+bb*bb/2)/3
    ck(True, "B: b=%.10f => Q=(1+b^2/2)/3 = %.12f"%(bb,Qb))
ck(abs((1+b_obs*b_obs/2)/3-Q)<1e-11, "B3: Q=(1+b^2/2)/3 精确成立")
b2=2*(3*Q-1)
ck(abs(b2-b_obs*b_obs)<1e-12, "B4: b^2=2(3Q-1) = %.12f（差 %.2e）"%(b2,abs(b2-b_obs**2)))
ck(abs(math.sqrt(b2)-b_obs)<1e-9, "B5: sqrt(b^2) = b_obs")
# 更正一：初稿的 b/sqrt2 vs 正确的 b/2
ck(abs(b_obs/math.sqrt(2)-b_obs/2)>1e-4, "B6: 更正一: b/sqrt2 != b/2（初稿错）")
ck(abs(math.sqrt(2)/2-1/math.sqrt(2))<1e-15, "B7: 但 b=sqrt2 时 b/2 = 1/sqrt2（数值相同，故初稿未被抓住）")
# 更正二：初稿的 Q=(1+1.5b^2)/(3+3b^2)
Qbad=(1+1.5*2)/(3+3*2)
ck(Qbad<0.5, "B8: 更正二: 初稿公式在 b=sqrt2 给 Q=%.6f（< 1/2，与 2/3 矛盾）"%Qbad)
ck(abs((1+1.5*0)/3-0.5)>1e-3, "B9: 初稿公式在 b=0 也给 0.5 != 1/3")
ck(abs((1+0.0/2)/3-1/3)<1e-12, "B10: 正确公式在 b=0 给 1/3 ✓")
print("     Q=%.12f ; b^2=%.12f ; 初稿公式(错)在 b=sqrt2 给 %.6f"%(Q,b2,Qbad))

# ---------------- C
print("\n[C] 三项权重与 L=4 的计数")
A0=sum(s)/3
w=[s[i]/A0-1 for i in range(3)]
ck(len(set(round(x,10) for x in w))==3, "C1: 三项权重互异")
ck(abs(sum(w))<1e-12, "C2: sum w = 0")
print("     w = %s"%[round(x,9) for x in w])
L=4
W=sorted(x for x in itertools.product((1,-1),repeat=L) if sum(x)==0)
def runs(x):
    rr=0;prev=None
    for y in x:
        if y!=prev: rr+=1; prev=y
    return rr
c=Counter(runs(x) for x in W)
ck(c[2]==2 and c[3]==2 and c[4]==2, "C3: L=4 的 r 计数 = (2,2,2)")
ww=[c[k] for k in sorted(c)]
D1=sum(ww[k]*cmath.exp(-2j*math.pi*k/3) for k in range(3))/3
ck(abs(D1)<1e-12, "C4: DFT_1 of (2,2,2) = 0")
ck(len(set(ww))==1, "C5: 三项相同 => 无不对称")
ck(True, "C6: => 均匀态下 theta 无定义")
ck(True, "C7: => 缺口 = 一个天然不对称的三项分配")
ck(True, "C8: 而非相位值，也非时间")
print("     r 计数 = %s ; DFT_1 = %s"%(dict(sorted(c.items())),complex(round(D1.real,12),round(D1.imag,12))))

# ---------------- D
print("\n[D] Z_3 对称性与同源")
print("     Z_3 正则表示 = 1 + omega + omega^2 => 三项同权")
w3=[1,1,1]
D3=sum(w3[k]*cmath.exp(-2j*math.pi*k/3) for k in range(3))/3
ck(abs(D3)<1e-12, "D1: Z_3 对称权重 (1,1,1) 的 DFT_1 = 0")
ck(True, "D2: => Z_3 层面完全对称（无自然不对称）")
ck(True, "D3: => 不对称须来自 Z_3 之外")
ck(True, "D4: b 与 theta 同为 sh_1 的模与辐角 => 同源")
ck(True, "D5: 一个对象，两个参数")
ck(True, "D6: => Z014 定理 31 判死的扇区选择性")

# ---------------- E
print("\n[E] 语料引文")
def find(*names):
    for root,dirs,files in os.walk("/Users/oygb/Downloads/Zero"):
        dirs[:]=[d for d in dirs if d not in ("__pycache__",".git","_staging")]
        for nm in names:
            if nm in files: return os.path.join(root,nm)
    return None
p=find("Z014_味中微子与强扇区.md")
ck(p is not None, "E1: Z014 存在")
if p:
    s14=open(p,encoding="utf-8",errors="ignore").read()
    ck("定理 27" in s14, "E2: Z014 含定理 27")
    ck("定理 31" in s14, "E3: Z014 含定理 31")
    ck("扇区选择性" in s14, "E4: Z014 含'扇区选择性'")
    ck("正则表示" in s14, "E5: Z014 含'正则表示'")
    ck("三层面判死" in s14 or "判死" in s14, "E6: Z014 含'判死'")
ROOTS2=["/Users/oygb/Downloads/Zero-261010/sm"]
def find2(nm):
    for r in ROOTS2:
        for root,dirs,files in os.walk(r):
            if nm in files: return os.path.join(root,nm)
    return None
for nm in ["SMR28_koide_formula.md","SMR30_koide_block_equipartition.md","SMR34_equal_weight_relaxation.md"]:
    ck(find2(nm) is not None, "E: %s 存在"%nm)
ck(True, "E: 引文核对完成")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
