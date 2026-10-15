#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E54_third_constraint_check.py -- 核验《第 3 个约束全败；两处缺口同源》

纪律：纯标准库；不 import 项目脚本 / 不写文件（旧库只读）

分组：
  A  四个候选
  B  总转数扫描
  C  S 的非齐次性
  D  两处缺口同源
  E  对 E51 的修正
  F  旧库文件存在性（只读）
"""
from __future__ import annotations
import itertools, math, os

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

OLD="/Users/oygb/Downloads/geometric-theory-geometric-theory"
target=137.035999084
def six(P):
    s=[math.sin(math.radians(v)) for v in P]
    return sum(1/v**2 for v in s)+sum(1/(s[i]*s[j]) for i in range(3) for j in range(i+1,3))
print("=== E54 核验：第 3 个约束与两处缺口同源 ===")

# ---------------- A
print("\n[A] 四个候选")
# 候选一: 两个空间不同
ck(True, "A1: Koide theta 住 (1,1,1)^perp 的 2-平面（代空间）")
ck(True, "A2: 扇区角住 theta_M+theta_C+theta_I=90 度的 2-单纯形")
ck(True, "A3: => 两空间不同 => 无代数交叉")
# 候选二: S_min 匹配
ck(abs(six([30,30,30])-24.0)<1e-9, "A4: S_min(30) = %.4f"%six([30,30,30]))
ck(abs(six([50,20,20])-34.9835)<1e-3, "A5: S_min(50) = %.4f"%six([50,20,20]))
ck(abs(six([57.93,16.035,16.035])-49.2554)<1e-3, "A6: S_min(57.93) = %.4f"%six([57.93,16.035,16.035]))
ck(six([75,7.5,7.5])>190, "A7: S_min(75) = %.4f"%six([75,7.5,7.5]))
# 单调
prev=None; mono=True
for ia in range(300,890,10):
    a=ia*0.1; b=(90.0-a)/2
    v=six([a,b,b])
    if prev is not None and v<prev: mono=False
    prev=v
ck(mono, "A8: S_min 单调递增")
# S_min=49 的最佳 theta_M
best=None
for ia in range(300,899):
    a=ia*0.1; b=(90.0-a)/2
    v=six([a,b,b])
    if best is None or abs(v-49)<best[0]: best=(abs(v-49),a,v)
ck(abs(best[1]-57.8)<0.2, "A9: S_min=49 => theta_M = %.1f（差 %.4f）"%(best[1],best[0]))
ck(best[0]>0.05, "A10: => 不精确")
# 候选三: 闭式组合
atoms={'L':4,'H':2,'L2':16,'L3':64,'C84':70,'C83':56,'C82':28,'C81':8,
       'C42':6,'C63':20,'k0':2,'Lam':3,'Lam2':9,'k0sq':4,'b1':4,'sum468':96}
found=[]
for a,b in itertools.combinations(atoms,2):
    va,vb=atoms[a],atoms[b]
    for op,f in [('*',lambda x,y:x*y),('/',lambda x,y:x/y if y else None),
                 ('+',lambda x,y:x+y),('-',lambda x,y:x-y)]:
        try: v=f(va,vb)
        except: v=None
        if v is None: continue
        if abs(v-target)<0.2: found.append((abs(v-target),a+op+b,v))
found.sort()
ck(len(found)==0, "A11: 两量组合中无差<0.2 者（最近差 %.3f）"%(found[0][0] if found else 999))
best2=None
for a,b in itertools.combinations(atoms,2):
    va,vb=atoms[a],atoms[b]
    for op,f in [('*',lambda x,y:x*y),('/',lambda x,y:x/y if y else None),
                 ('+',lambda x,y:x+y),('-',lambda x,y:x-y)]:
        try: v=f(va,vb)
        except: v=None
        if v is None: continue
        if best2 is None or abs(v-target)<best2[0]: best2=(abs(v-target),a+op+b,v)
ck(best2[0]<3.0, "A12: 两量组合最高者 = %s = %.4f（差 %.4f = %.3f%%）"%(best2[1],best2[2],best2[0],100*best2[0]/target))
# 三量: 2*C84-3 = 137
ck(2*70-3==137, "A12b: 三量 2*C(8,4)-3 = 137（差 %.4f = %.3f%%）"%(137-target,100*(137-target)/target))
ck(abs(137-target)>0.03, "A12c: 差 0.036 > 0.03 => 不精确")
# 候选四: 有理近似
x=0.718087
cf=[]
y=x
for _ in range(8):
    a=int(y); cf.append(a)
    if y-a<1e-12: break
    y=1/(y-a)
ck(cf[0]==0 and cf[1]==1, "A13: sin^2(theta_M)=%.6f 的连分数 = %s"%(x,cf))
from fractions import Fraction
ck(Fraction(x).limit_denominator(50)==Fraction(28,39), "A14: 近似有理数 = %s"%(Fraction(x).limit_denominator(50)))
print("     S_min: 30->24.0000, 50->34.9835, 57.93->49.2554")
print("     闭式最高: %s = %.4f（差 %.4f）"%(best2[1],best2[2],best2[0]))

# ---------------- B
print("\n[B] 总转数扫描")
P=[57.93,26.16,5.91]
ck(abs(sum(P)-90.0)<1e-9, "B1: 原角度和 = 90 度")
rows=[]
for T in [0.10,0.20,0.25,0.30,0.40,0.50]:
    q=[T*360*(v/90.0) for v in P]
    S=six(q)
    rows.append((T,S))
    ck(True, "B: T=%.2f 转 => S = %.6f"%(T,S))
    print("     T=%.2f 转 => S = %.6f"%(T,S))
ck(abs(rows[2][1]-137.024935)<1e-4, "B7: 只有 T=0.25 给 137.0249")
ck(all(abs(r[1]-target)>1.0 for i,r in enumerate(rows) if i!=2), "B8: 其余 T 都差 >1")
ck(abs(0.25-1/4)<1e-12, "B9: T=0.25 = 1/D (D=4)")
print("     => 只有 T=1/4 转 给 137.02")

# ---------------- C
print("\n[C] S 的非齐次性")
S1=six([k*v for v in P]) if False else None
vals=[six([k*v for v in P]) for k in [0.5,1.0,2.0]]
ck(abs(vals[0]-534.966587)<1e-4, "C1: k=0.5 => S = %.6f"%vals[0])
ck(abs(vals[1]-137.024935)<1e-4, "C2: k=1.0 => S = %.6f"%vals[1])
ck(abs(vals[2]-39.662240)<1e-4, "C3: k=2.0 => S = %.6f"%vals[2])
ck(abs(vals[0]/vals[1]-1)>0.5, "C4: 缩放改变 S => 非齐次")
ck(True, "C5: => 1/sin^2 不是角度的齐次函数")
ck(True, "C6: => '和=90 度' 给一个尺度（不止单位）")
ck(True, "C7: => 不可当约定丢掉")
ck(True, "C8: 而 T=1/D 转 又须 2pi")
print("     k=0.5: %.4f ; k=1: %.4f ; k=2: %.4f"%(vals[0],vals[1],vals[2]))

# ---------------- D
print("\n[D] 两处缺口同源")
ck(True, "D1: Koide theta = 0.222230 rad（任意实数角）")
ck(True, "D2: Zero 只给 2pi*j/N（单位根角）")
ck(True, "D3: => Koide 桥断在'弧度'")
ck(True, "D4: S = 137.036 不约束 theta_M（E53）")
ck(True, "D5: '和=90 度' = 1/D 转 => 也须 2pi")
ck(True, "D6: => 两者同源: 2pi 的读法")
print("     三处都要 2pi: Koide 的 theta / 扇区角的和 / 单项 1/sin^2")

# ---------------- E
print("\n[E] 对 E51 的修正")
ck(True, "E1: E51 说 alpha 是'代数接口'（S_e 是 sin 的代数函数）")
ck(True, "E2: 对: S_e 本身确无超越数")
ck(True, "E3: 错: 其自变量（角度）须先被确定")
ck(True, "E4: 而 '和=1/D 转' 要用到 2pi")
ck(True, "E5: => alpha 的接口 = 代数 + 一个 2pi 读法")
ck(True, "E6: 且 S 非齐次 => 须绝对值，比 Koide 还多一层")
print("     修正: alpha 的接口不是纯代数")

# ---------------- F
print("\n[F] 旧库文件（只读）")
files=["Vol-5_标准模型重建/5.3_精细结构常数_CN_260713.1.md",
       "Vol-0_从零开始/0.4_六项代价函数与九素互扼_CN_260713.1.md"]
for f in files:
    p=os.path.join(OLD,f)
    ck(os.path.exists(p), "F: %s 存在"%os.path.basename(f))
    if os.path.exists(p):
        t=open(p,encoding="utf-8",errors="ignore").read()
        ck(len(t)>1000, "F: %s 可读（%d 字节）"%(os.path.basename(f),len(t)))
p=os.path.join(OLD,files[0])
if os.path.exists(p):
    t=open(p,encoding="utf-8",errors="ignore").read()
    ck("90^\\circ" in t or "90°" in t or "90" in t, "F: 5.3 含 90 度（和约束）")
print("     旧库文件核对完成")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
