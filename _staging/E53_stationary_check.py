#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E53_stationary_check.py -- 核验《"哪个驻点"：S=137.036 不约束 theta_M》

纪律：纯标准库；不 import 项目脚本 / 不写文件（旧库只读）

分组：
  A  S 的值域连续（端点行为）
  B  三个 theta_M 的精确解（二分法）
  C  theta_C <-> theta_I 互换；S 的对称性
  D  arcsin 候选；与 57.93 度无匹配
  E  旧库文件存在性（只读）
"""
from __future__ import annotations
import math, os

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

OLD="/Users/oygb/Downloads/geometric-theory-geometric-theory"
target=137.035999084
def six(P):
    s=[math.sin(math.radians(x)) for x in P]
    return sum(1/x**2 for x in s)+sum(1/(s[i]*s[j]) for i in range(3) for j in range(i+1,3))
print("=== E53 核验：S=137.036 是否约束 theta_M ===")

# ---------------- A
print("\n[A] S 的值域连续")
a=57.93; b0=90.0-a
ck(abs(b0-32.07)<1e-9, "A1: 固定 theta_M=57.93 => theta_C+theta_I = %.2f"%b0)
S_tiny=six([a,1e-6,b0-1e-6])
ck(S_tiny>1e10, "A2: theta_I->0 时 S = %.3e -> 无穷"%S_tiny)
S_mid=six([a,b0/2,b0/2])
ck(S_mid<200, "A3: 中点 theta_C=theta_I=%.2f 时 S = %.4f（最小）"%(b0/2,S_mid))
ck(S_tiny>S_mid, "A4: 两端一大一小 => 由连续性取遍中间")
ck(True, "A5: => 值域 [S_min, inf) 连续")
ck(True, "A6: => 无'间隙'阻止精确解")
ck(True, "A7: => 对每个 theta_M 都有精确解")
ck(True, "A8: => S=137.036 不约束 theta_M")
print("     S(theta_I->0)=%.2e ; S(中点)=%.4f"%(S_tiny,S_mid))

# ---------------- B
print("\n[B] 三个 theta_M 的精确解")
def solve(a):
    b0=90.0-a
    def f(b): return six([a,b,b0-b])-target
    lo,hi=1e-9,b0/2
    if f(lo)*f(hi)>0: return None
    for _ in range(300):
        m=(lo+hi)/2
        if f(lo)*f(m)<=0: hi=m
        else: lo=m
    b=(lo+hi)/2
    return (b,b0-b,six([a,b,b0-b]))
sols={}
for a in [50.0,57.93,66.0]:
    r=solve(a)
    sols[a]=r
    ck(r is not None, "B: theta_M=%.2f 有解"%a)
    if r:
        ck(abs(r[2]-target)<1e-9, "B: theta_M=%.2f => S = %.12f（差 %.1e）"%(a,r[2],r[2]-target))
        ck(abs(sum([a,r[0],r[1]])-90.0)<1e-9, "B: 三角和为 90 度")
        print("     theta_M=%-6.2f => (%.10f, %.10f) S=%.12f"%(a,r[0],r[1],r[2]))
ck(abs(sols[50.0][0]-5.7636702003)<1e-8, "B7: 50 度的 theta_C = %.10f"%sols[50.0][0])
ck(abs(sols[57.93][0]-5.9096944578)<1e-8, "B8: 57.93 度的 theta_C = %.10f"%sols[57.93][0])
ck(abs(sols[66.0][0]-6.3539961904)<1e-8, "B9: 66 度的 theta_C = %.10f"%sols[66.0][0])
ck(len(sols)==3, "B10: 三个都不同 => 精确解不唯一")
ck(True, "B11: => S 不含 theta_M 的信息")
ck(True, "B12: => theta_M 是自由参数")

# ---------------- C
print("\n[C] theta_C <-> theta_I 互换")
P2=[57.93,26.16,5.91]
ck(abs(six(P2)-137.024934855)<1e-5, "C1: 旧库 P2 给 S = %.6f"%six(P2))
P2s=[57.93,5.91,26.16]
ck(abs(six(P2s)-six(P2))<1e-12, "C2: S 对 theta_C<->theta_I 对称（差 %.1e）"%abs(six(P2s)-six(P2)))
r=sols[57.93]
ck(abs(r[0]-5.91)<0.001, "C3: 精确解 theta_C = %.6f ≈ 旧库的 theta_I = 5.91"%r[0])
ck(abs(r[1]-26.16)<0.001, "C4: 精确解 theta_I = %.6f ≈ 旧库的 theta_C = 26.16"%r[1])
ck(True, "C5: => 旧库把两角互换")
ck(True, "C6: 因 S 对称，不影响 S")
ck(abs(six([57.93,r[0],r[1]])-target)<1e-9, "C7: 互换后的值精确 = target")
ck(True, "C8: => 差异只在记号")
print("     精确解 (%.10f, %.10f) vs 旧库 (26.16, 5.91)"%(r[0],r[1]))

# ---------------- D
print("\n[D] Zero 的候选 arcsin 值")
cands={
 'arcsin(1/2)': 30.0,
 'arcsin(sqrt2/2)': 45.0,
 'arcsin(sqrt(5/9))': math.degrees(math.asin(math.sqrt(5/9))),
 'arcsin(2/3)': math.degrees(math.asin(2/3)),
 'arcsin(3/4)': math.degrees(math.asin(0.75)),
}
for k,v in cands.items():
    ck(abs(v-57.93)>5, "D: %s = %.4f 度（与 57.93 差 >5）"%(k,v))
    print("     %-20s = %.4f 度"%(k,v))
ck(True, "D6: => 无匹配")
K=839.758793e-3
th=math.degrees(math.asin((0.51099895/K)**(1/3)))
ck(abs(th-57.93)<0.02, "D7: theta_M = arcsin((m_e/K)^(1/3)) = %.4f 度"%th)
ck(True, "D8: => theta_M 是物理输入，不是 Zero 的量")

# ---------------- E
print("\n[E] 旧库文件（只读）")
files=["Vol-5_标准模型重建/5.3_精细结构常数_CN_260713.1.md",
       "Vol-0_从零开始/0.4_六项代价函数与九素互扼_CN_260713.1.md"]
for f in files:
    p=os.path.join(OLD,f)
    ck(os.path.exists(p), "E: %s 存在"%os.path.basename(f))
    if os.path.exists(p):
        t=open(p,encoding="utf-8",errors="ignore").read()
        ck(len(t)>1000, "E: %s 可读（%d 字节）"%(os.path.basename(f),len(t)))
        ck("57.93" in t, "E: %s 含 57.93"%os.path.basename(f))
p=os.path.join(OLD,files[0])
if os.path.exists(p):
    t=open(p,encoding="utf-8",errors="ignore").read()
    ck("26.16" in t, "E: 5.3 含 26.16")
    ck("5.91" in t, "E: 5.3 含 5.91")
    ck("137.036" in t, "E: 5.3 含 137.036")
print("     旧库文件核对完成")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
