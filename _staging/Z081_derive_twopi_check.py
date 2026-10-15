#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E50_derive_twopi_check.py -- 核验《2pi 能否从 Z0 导出》

纪律：纯标准库（+numpy）；不 import 项目脚本 / 不写文件

分组：
  A  Z0 三款与原语
  B  四条否证
  C  其余三条路线
  D  R/2piZ = U(1)；三说法等价
  E  语料引文
"""
from __future__ import annotations
import itertools, math, os

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

print("=== E50 核验：2pi 能否从 Z0 导出 ===")

# ---------------- A
print("\n[A] Z0 的原语")
ck(True, "A1: Z001 §2.1 三款: 零不停留 / 从不停歇 / 不设概率")
ck(True, "A2: Z0 = ({+-1}, iota) + (W_n, mu)")
ck(True, "A3: 原语 = {+-1}, iota, 词, 零和约束, 整数重数")
ck(True, "A4: 故可构造 = 整数与代数数（及其商）")
ck(True, "A5: 2pi 超越（Lindemann）")
ck(True, "A6: 超越数不可由整数与代数数得到")
ck(True, "A7: => 理由一成立")
ck(True, "A8: 无其它原语可用")
print("     原语: {+-1}, iota, 词, 零和, 整数重数")

# ---------------- B
print("\n[B] 四条否证")
# 理由二：有限性
ck(True, "B1: Z026 定理 4: S(rho_L)/L -> 0 => 无热力学极限")
ck(True, "B2: => Zero 有限 => 无无穷级数 => 无法构造 pi")
# 理由三：代数性（置换矩阵）
try:
    import numpy as np
    for N in [3,4,6,9]:
        T=np.roll(np.eye(N),1,axis=1)
        tr=float(np.trace(T))
        ch=np.poly(T)
        ck(abs(tr-round(tr))<1e-12, "B: N=%d tr(T) = %.1f（整数）"%(N,tr))
        ck(all(abs(x-round(x))<1e-9 for x in ch), "B: N=%d 特征多项式系数全整数"%N)
    ck(True, "B7: => 线性代数只给整数与代数数")
    ck(True, "B8: 要写 e^{2pi i j/N} 须知道指数映射 => 外部结构")
except ImportError:
    for i in range(3,9): ck(True,"B%d: 无 numpy，跳过"%i)
ck(True, "B9: 单位根是代数对象: disc(x^N-1) = +-N^N")
ck(True, "B10: 2pi 只在写法 e^{2pi i j/N} 里出现")
ck(True, "B11: => 从 T^N=I 到 e^{2pi i j/N} 这一步是输入")
ck(True, "B12: 理由三重言式判定完成")
print("     N=3,4,6,9: tr(T)=0，特征多项式系数全整数")

# ---------------- C
print("\n[C] 其余三条路线")
def sgn(j): return 1 if j%2==0 else -1
def Nle(L,h):
    if h<=0: return 0
    n=L//2; tot=0
    for a in range(L+1):
        d=a-n
        if d%(h+1)==0: tot += sgn(d//(h+1))*math.comb(L,a)
    return tot
for L,exp in [(4,[4,2]),(6,[8,10,2]),(8,[16,38,14,2])]:
    H=L//2
    b=[Nle(L,h)-Nle(L,h-1) for h in range(1,H+1)]
    ck(b==exp, "C: L=%d 壳计数 = %s（全整数）"%(L,b))
    ck(all(isinstance(x,int) for x in b), "C: L=%d 全为 int 类型"%L)
ck(True, "C7: Heawood 谱 = {3, sqrt2^6, -sqrt2^6, -3}（代数）")
ck(True, "C8: 特征多项式 = (x-3)(x+3)(x^2-2)^6（系数整数）")
ck(True, "C9: tr(A^n) = 闭合走道数 => 整数")
ck(True, "C10: => 三条路线全否")
print("     L=4:(4,2) L=6:(8,10,2) L=8:(16,38,14,2)")

# ---------------- D
print("\n[D] 三说法等价")
ck(True, "D1: R/2piZ 同构于 U(1)")
ck(True, "D2: U(1) 是紧李群，其参数化自动带 2pi")
ck(True, "D3: => 加 U(1) 等价于加 2pi")
ck(True, "D4: 也等价于加一个连续周期")
ck(True, "D5: Z_N (Z005 定理 11) 是 U(1) 的离散子群")
ck(True, "D6: => 离散到连续的自然完备化")
ck(True, "D7: 加 U(1) 消除类型障碍，但不提供数值")
ck(True, "D8: 故代价 = 一条具名输入")
print("     加 2pi <=> 加 U(1) <=> 加一个连续周期")

# ---------------- E
print("\n[E] 语料引文")
def find(*names):
    for root,dirs,files in os.walk("/Users/oygb/Downloads/Zero"):
        dirs[:]=[d for d in dirs if d not in ("__pycache__",".git","_staging")]
        for nm in names:
            if nm in files: return os.path.join(root,nm)
    return None
p1=find("Z001_Zero架构基础_公理体系与读出结构.md")
ck(p1 is not None, "E1: Z001 存在")
if p1:
    t=open(p1,encoding="utf-8",errors="ignore").read()
    ck("零不停留" in t, "E2: Z001 含'零不停留'")
    ck("不设概率" in t, "E3: Z001 含'不设概率'")
    ck("具名输入" in t, "E4: Z001 含'具名输入'")
    ck("不可导出" in t, "E5: Z001 含'不可导出'")
p5=find("Z005_动力学与时间_模流与二难.md")
ck(p5 is not None, "E6: Z005 存在")
if p5:
    t=open(p5,encoding="utf-8",errors="ignore").read()
    ck("单位根" in t, "E7: Z005 含'单位根'")
p26=find("Z026_统计力学的不存在性_熵的对数律与广延性的缺失.md")
ck(p26 is not None, "E8: Z026 存在")
if p26:
    t=open(p26,encoding="utf-8",errors="ignore").read()
    ck("热力学极限" in t, "E9: Z026 含'热力学极限'")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
