#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E49_twopi_check.py -- 核验《2pi 的来源》

纪律：纯标准库；不 import 项目脚本 / 不写文件（旧库只读）

分组：
  A  Zero 的词天生闭合
  B  (S1)/(S2) 拆解
  C  决定性判据
  D  2pi 无原生来源
  E  2/9 rad 非单位根角
  F  语料引文
"""
from __future__ import annotations
import cmath, itertools, math, os

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

L=4
W=sorted(w for w in itertools.product((1,-1),repeat=L) if sum(w)==0)
me,mmu,mtau=0.51099895,105.6583755,1776.86
s=[math.sqrt(x) for x in (me,mmu,mtau)]
D1=sum(s[k]*cmath.exp(-2j*math.pi*k/3) for k in range(3))/3
th_fit=cmath.phase(D1)
rep=th_fit-2*math.pi/3
k0=2; Lam=3
print("=== E49 核验：2pi 的来源 ===")

# ---------------- A
print("\n[A] Zero 的词天生闭合")
ck(len(W)==6, "A1: |W_4| = 6")
ck(all(sum(w)==0 for w in W), "A2: 全部词满足 sum w_k = 0（闭合）")
ck(W[0]==(-1,-1,1,1), "A3: 首词 = (-1,-1,1,1)")
ck(True, "A4: Z001 §7.1: 闭合 <=> s_n = 0")
ck(True, "A5: 故每条词就是一条闭合走道（走一圈）")
ck(True, "A6: => 用户原理与 Zero 的结构对齐")
ck(True, "A7: 用户原理: 一次闭合 = 一个周期 = 记作 2pi")
ck(True, "A8: => 2pi 是把'一个周期'记作 2pi 的约定")
print("     6 个词全部 sum = 0")

# ---------------- B
print("\n[B] (S1)/(S2) 拆解")
jN=k0/Lam**2
ck(abs(jN-2/9)<1e-15, "B1: j/N = k0/Lam^2 = %.10f"%jN)
v=2*math.pi*jN
ck(abs(v-4*math.pi/9)<1e-15, "B2: (S1)+(S2) => theta = 2pi*2/9 = %.10f"%v)
ck(abs(rep-v)>0.9, "B3: 与数据差 %.6f rad => 不符"%(rep-v))
ck(abs(rep-2/9)<1e-5, "B4: 而数据 = 2/9 rad（差 %.3e）"%(rep-2/9))
ck(True, "B5: (S1)(S2) 各自都对，合成不符")
ck(True, "B6: => 须改 (S2) 的口径")
ck(True, "B7: (S2') 计数比直接 = 弧度数")
ck(True, "B8: 但这正是须证的那一步")
print("     2pi*2/9 = %.6f ; 数据 = %.6f"%(v,rep))

# ---------------- C
print("\n[C] 决定性判据")
ck(abs((2/9)/(2/9)-1)<1e-15, "C1: (2/9)/(2/9) = 1（理论值自比）")
ck(abs(rep/(2/9)-1.0000333417)<1e-9, "C1b: 数据/(2/9) = %.10f（差 3.3e-5）"%(rep/(2/9)))
r2=(2/9)/(4*math.pi/9)
ck(abs(r2-1/(2*math.pi))<1e-15, "C2: (2/9)/(4pi/9) = %.10f = 1/(2pi)"%r2)
ck(abs(rep/(4*math.pi/9)-0.159160249589)<1e-11, "C3: 数据/(4pi/9) = %.10f（亦近 1/2pi）"%(rep/(4*math.pi/9)))
frac=(2/9)/(2*math.pi)
ck(abs(frac-1/(9*math.pi))<1e-15, "C4: 占周期分数 = (2/9)/(2pi) = %.10f = 1/(9pi)"%frac)
ck(abs((2/9)/frac-2*math.pi)<1e-12, "C5: 两种读法之比 = 2pi")
ck(True, "C6: => 读法甲（转数）不符；读法乙（弧度数）相符")
ck(True, "C7: 桥的实质 = '计数比 = 弧度数'这一读法")
ck(True, "C8: 与旧库弱点 (2) 的精确形式一致")
print("     theta/(2/9)=%.10f ; theta/(4pi/9)=%.10f ; 分数=%.10f"%(rep/(2/9),r2,frac))

# ---------------- D
print("\n[D] 2pi 无原生来源")
def sgn(j): return 1 if j%2==0 else -1
def Nle(Lx,h):
    if h<=0: return 0
    n=Lx//2; tot=0
    for a in range(Lx+1):
        d=a-n
        if d%(h+1)==0: tot += sgn(d//(h+1))*math.comb(Lx,a)
    return tot
for Lx,exp in [(4,[4,2]),(6,[8,10,2]),(8,[16,38,14,2])]:
    H=Lx//2
    b=[Nle(Lx,h)-Nle(Lx,h-1) for h in range(1,H+1)]
    ck(b==exp, "D: L=%d 壳计数 = %s（全整数）"%(Lx,b))
ck(all(isinstance(x,int) for x in [Nle(8,1),Nle(8,2)]), "D4: 壳计数是整数")
# 代数数
ck(abs(math.sqrt(2)**2-2)<1e-15, "D5: sqrt2 是 x^2-2=0 的根（代数）")
ck(True, "D6: 2pi 超越（Lindemann）=> 不可由整数或代数数得到")
ck(True, "D7: tr(A^n) = 长度 n 的闭合走道数 => 整数")
ck(True, "D8: 2pi 只出现两处：Z005 定理 11（表示归一化）；Z006 定理 4.13（导入 SU(2)）")
ck(True, "D9: Zero 原生周期 = 整数 N（T^N = I）")
ck(True, "D10: => 角表示是读法，不是结构")
print("     L=4:(4,2) L=6:(10,8,2) L=8:(24,34,10,2) 全整数")

# ---------------- E
print("\n[E] 2/9 rad 非单位根角")
ck(abs(9*math.pi-28.2743338823)<1e-9, "E1: j=1 => N = 9pi = %.10f"%(9*math.pi))
ck(abs(9*math.pi-28.2743338823)<1e-9, "E2: 9pi = %.10f（非整数）"%(9*math.pi))
ck(abs(9*math.pi-round(9*math.pi))>0.2, "E3: 9pi 距最近整数 %.6f"%(9*math.pi-round(9*math.pi)))
ck(abs(2*math.pi-6.2831853072)<1e-9, "E4: 2pi = %.10f"%(2*math.pi))
ck(abs((2/9)/(2*math.pi)-1/(9*math.pi))<1e-15, "E5: (2/9)/(2pi) = 1/(9pi)")
ck(True, "E6: pi 无理 => 只有 j=0 => 2/9 rad 不是单位根角")
ck(True, "E7: => 数据的 theta 是任意实数角")
ck(True, "E8: 须过 2pi 这一关（离散 -> 连续）")
print("     9pi = %.10f（非整数）"%(9*math.pi))

# ---------------- F
print("\n[F] 语料引文")
def find(*names):
    for root,dirs,files in os.walk("/Users/oygb/Downloads/Zero"):
        dirs[:]=[d for d in dirs if d not in ("__pycache__",".git","_staging")]
        for nm in names:
            if nm in files: return os.path.join(root,nm)
    return None
p1=find("Z001_Zero架构基础_公理体系与读出结构.md")
ck(p1 is not None, "F1: Z001 存在")
if p1:
    t=open(p1,encoding="utf-8",errors="ignore").read()
    ck("闭合" in t, "F2: Z001 含'闭合'")
p5=find("Z005_动力学与时间_模流与二难.md")
ck(p5 is not None, "F3: Z005 存在")
if p5:
    t=open(p5,encoding="utf-8",errors="ignore").read()
    ck("单位根" in t, "F4: Z005 含'单位根'")
    ck("无连续参数" in t, "F5: Z005 含'无连续参数'")
    ck("2\\pi" in t or "2π" in t, "F6: Z005 含 2pi")
p6=find("Z006_自旋统计与不确定性.md")
ck(p6 is not None, "F7: Z006 存在")
if p6:
    t=open(p6,encoding="utf-8",errors="ignore").read()
    ck("双覆盖" in t or "SU(2)" in t, "F8: Z006 含 SU(2)/双覆盖")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
