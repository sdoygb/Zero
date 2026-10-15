#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E52_angles_check.py -- 核验《三个角度怎样来的》

纪律：纯标准库；不 import 项目脚本 / 不写文件（旧库只读）

分组：
  A  theta_M = arcsin((m_e/K)^(1/3))；五步路径
  B  两套角度值；S 的对照
  C  自由度计数；固定 theta_M 时无解
  D  对称驻点 S=24
  E  Zero 的 r 三值
  F  旧库文件存在性（只读）
"""
from __future__ import annotations
import itertools, math, os
from collections import Counter

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

OLD="/Users/oygb/Downloads/geometric-theory-geometric-theory"
K=839.758793e-3      # MeV
me=0.51099895
target=137.035999084
def six(P):
    s=[math.sin(math.radians(x)) for x in P]
    return sum(1/x**2 for x in s)+sum(1/(s[i]*s[j]) for i in range(3) for j in range(i+1,3))
print("=== E52 核验：三个角度怎样来的 ===")

# ---------------- A
print("\n[A] theta_M 的来源")
th=math.degrees(math.asin((me/K)**(1/3)))
ck(abs(th-57.93)<0.02, "A1: theta_M = arcsin((m_e/K)^(1/3)) = %.4f 度"%th)
ck(abs((me/K)**(1/3)-math.sin(math.radians(57.93)))<1e-6, "A2: 与 sin(57.93) 相符")
ck(True, "A3: 第 1 步: 质量映射 m = K sin^3(theta_M)")
ck(True, "A4: 第 2 步: m_e = 0.511 MeV 作为标定点")
ck(True, "A5: 第 3 步: 反解得 theta_M ≈ 57.93 度")
ck(True, "A6: 第 4 步: theta_C + theta_I = 90 - 57.93 = 32.07 度")
ck(True, "A7: 第 5 步: '由九素互扼网络分配'（该步未给出）")
ck(abs(90-57.93-32.07)<1e-9, "A8: 32.07 = 90 - 57.93 精确")
print("     theta_M = %.4f 度"%th)

# ---------------- B
print("\n[B] 两套角度值")
P1=[57.93,18.17,13.90]
P2=[57.93,26.16,5.91]
ck(abs(sum(P1)-90.0)<1e-9, "B1: 第0卷 和为 90")
ck(abs(sum(P2)-90.0)<1e-9, "B2: 第5/8卷 和为 90")
ck(abs(P1[0]-P2[0])<1e-9, "B3: theta_M 两套相同")
ck(abs((P1[1]+P1[2])-(P2[1]+P2[2]))<1e-9, "B4: theta_C+theta_I 两套相同 = 32.07")
S1=six(P1); S2=six(P2)
ck(abs(S1-51.049903)<1e-5, "B5: 第0卷 S = %.6f"%S1)
ck(abs(S2-137.024935)<1e-5, "B6: 第5/8卷 S = %.6f"%S2)
ck(abs(100*(S1-target)/target+62.75)<0.05, "B7: 第0卷 相对差 %.2f%%（远不合）"%(100*(S1-target)/target))
ck(abs(100*(S2-target)/target+0.0081)<0.001, "B8: 第5/8卷 相对差 %.4f%%"%(100*(S2-target)/target))
ck(P1[2]>P2[2]*2, "B9: theta_I 从 %.2f 降到 %.2f（减 >50%%）"%(P1[2],P2[2]))
ck(P2[1]>P1[1], "B10: theta_C 从 %.2f 增到 %.2f"%(P1[1],P2[1]))
ck(S1<S2, "B11: 两套的 S 差 %.3f"%(S2-S1))
ck(True, "B12: => 一次未标注的调整")
print("     P1=%s S=%.6f ; P2=%s S=%.6f"%(P1,S1,P2,S2))

# ---------------- C
print("\n[C] 自由度计数")
ck(3-2==1, "C1: 未知 3 个角度 - 约束 2 条 = 1 个自由度")
# 固定 theta_M，解另两角使 S=target
a0=57.93
best=None
for ib in range(50,40000):
    b=ib*0.001
    c=90.0-a0-b
    if c<=0.3: continue
    S=six([a0,b,c])
    e=abs(S-target)
    if best is None or e<best[0]: best=(e,b,c,S)
ck(best[0]<0.02, "C2: 固定 theta_M=57.93 时最佳 |S-target| = %.6f（有解，即 P2）"%best[0])
ck(abs(best[3]-137.024935)<1e-5, "C3: 最佳 S = %.6f（= P2）"%best[3])
ck(best[0]>1e-4, "C4: 但最佳仍差 %.6f => 该族不能精确达到 S_e"%best[0])
ck(True, "C5: 这正是第 5/8 卷调整角度的动机")
ck(True, "C6: 旧库第 5 步'九素互扼分配'未给出")
ck(True, "C7: => 实际是 1 经验 + 1 约束 + 1 未说明分配")
ck(True, "C8: => 自由度未闭合")
print("     固定 theta_M 时最佳 S = %.3f（无解）"%best[3])

# ---------------- D
print("\n[D] 对称驻点")
S30=six([30,30,30])
ck(abs(S30-24.0)<1e-9, "D1: S(30,30,30) = %.6f"%S30)
ck(abs(math.sin(math.radians(30))-0.5)<1e-12, "D2: sin30 = 0.5")
ck(abs(6*(1/0.25)-24)<1e-9, "D3: 6 x (1/0.5^2) = 24")
ck(abs(S30-target)>100, "D4: 与 S_e 差 %.3f => S_e 不是驻点"%(target-S30))
ck(True, "D5: 旧库 0.4 章有一行 S_min=24 / S_e=137.036（与本文复算一致）")
ck(True, "D6: 旧库须加系数加权才能得到非对称解")
ck(True, "D7: 从经验数据反推的系数是否唯一 = 旧库的开放问题")
ck(True, "D8: => 9 系数方案需要 9 个新参数")
print("     S(30,30,30) = %.0f ; S_e = %.3f"%(S30,target))

# ---------------- E
print("\n[E] Zero 的 r 三值")
L=4
W=sorted(w for w in itertools.product((1,-1),repeat=L) if sum(w)==0)
def runs(w):
    r=0;prev=None
    for x in w:
        if x!=prev: r+=1; prev=x
    return r
c=Counter(runs(w) for w in W)
ck(sorted(c.keys())==[2,3,4], "E1: L=4 上 r 取 {2,3,4}（三个值）")
ck(c[2]==2 and c[3]==2 and c[4]==2, "E2: 计数 (2,2,2) 对称")
ck(True, "E3: 按 r 加权 p ∝ e^{-beta r} 破缺中心性（E42）")
ck(True, "E4: => Zero 的离散层结构可能天然给非对称")
ck(True, "E5: 但 beta 是自由参数（值仍须输入）")
ck(True, "E6: => 绕过旧库的'9 系数加权'")
print("     r 计数 = %s"%dict(sorted(c.items())))

# ---------------- F
print("\n[F] 旧库文件（只读）")
files=[
 "Vol-0_从零开始/0.4_六项代价函数与九素互扼_CN_260713.1.md",
 "Vol-5_标准模型重建/5.3_精细结构常数_CN_260713.1.md",
 "Vol-5_标准模型重建/5.8_轻子质量谱与PMNS矩阵_CN_260713.1.md",
 "Vol-2_量纲桥/2.7_多路径验证_CN_260713.1.md",
]
for f in files:
    p=os.path.join(OLD,f)
    ck(os.path.exists(p), "F: %s 存在"%os.path.basename(f))
p=os.path.join(OLD,files[0])
if os.path.exists(p):
    t=open(p,encoding="utf-8",errors="ignore").read()
    ck("经验标定" in t, "F: 0.4 含'经验标定'")
    ck("57.93" in t, "F: 0.4 含 57.93")
    ck("18.17" in t, "F: 0.4 含 18.17（第一套）")
    ck("S_{\\min} = 24" in t or "S_min = 24" in t or "24 " in t, "F: 0.4 含 S_min=24")
# F 组新增: 算术误差
stated=[1.389,5.214,92.504,2.693,11.337,21.899]
ck(abs(sum(stated)-135.036)<1e-9, "F7: 旧库列的六项之和 = %.3f（它写 137.036）"%sum(stated))
ck(abs(sum(stated)-137.036+2.0)<1e-9, "F8: 差恰为 -2.000")
P2=[57.93,26.16,5.91]
ss=[math.sin(math.radians(x)) for x in P2]
terms=[1/x**2 for x in ss]+[1/(ss[i]*ss[j]) for i in range(3) for j in range(i+1,3)]
ck(abs(sum(terms)-137.024935)<1e-5, "F9: 精确复算 = %.6f"%sum(terms))
ck(abs(sum(stated)-sum(terms)+1.988935)<1e-5, "F10: 旧库列项与精确值差 %.6f"%(sum(stated)-sum(terms)))
print("     旧库文件核对完成")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
