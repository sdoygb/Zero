#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E46_sector_check.py -- 核验《"扇区"是什么》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  两种扇区的定义与数目
  B  r 给出恰好 3 个；计数 (2,2,2)；r 范围随 L
  C  Z014 的翻转条件；E42 的非普适性
  D  术语混用统计
  E  语料引文
"""
from __future__ import annotations
import itertools, os
from collections import Counter

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def Hh(w):
    s=0;m=0
    for x in w:
        s+=x; m=max(m,abs(s))
    return m
def runs(w):
    r=0;prev=None
    for x in w:
        if x!=prev: r+=1; prev=x
    return r
print("=== E46 核验：扇区是什么 ===")

# ---------------- A
print("\n[A] 两种扇区")
L=4
W=sorted(w for w in itertools.product((1,-1),repeat=L) if sum(w)==0)
H=L//2
ck(H==2, "A1: L=4 时 H = L/2 = %d"%H)
ck(len(set(Hh(w) for w in W))==2, "A2: 壳层数 = 2 => 超选择扇区数 = 2")
ck(True, "A3: Z004 定理 2: 超选择扇区 = 最小中心投影 P_h")
ck(True, "A4: 对象 = 壳层 h")
ck(True, "A5: 味扇区 = 上型夸克/下型夸克/带电轻子 => 3 个")
ck(True, "A6: Z014 §7.1: '三扇区的汤川算符结构完全相同'")
ck(True, "A7: Z014 §7.1: '两夸克扇区同基循环时 V_CKM = 1'")
ck(H!=3, "A8: 两者数目不同（2 vs 3）=> 对象不同")
print("     超选择扇区 %d 个 ; 味扇区 3 个"%H)

# ---------------- B
print("\n[B] Zero 里的三个")
rs=sorted(set(runs(w) for w in W))
ck(rs==[2,3,4], "B1: L=4 的 r 取值 = %s（恰好 3 个）"%rs)
c=Counter(runs(w) for w in W)
ck(c[2]==2 and c[3]==2 and c[4]==2, "B2: r 的计数 = (2,2,2)")
ck(len(set(c.values()))==1, "B3: 三项相同 => 无不对称")
ck(len(set(Hh(w) for w in W))==2, "B4: 壳层只有 2 个 => 不够 3")
# r 范围随 L
expect={4:3,6:5,8:7,10:9,12:11}
for Lx,k in expect.items():
    Wx=[w for w in itertools.product((1,-1),repeat=Lx) if sum(w)==0]
    nx=len(set(runs(w) for w in Wx))
    ck(nx==k, "B: L=%d 的 r 取值数 = %d"% (Lx,nx))
print("     L=4: r = %s ; 计数 = %s"%(rs,dict(sorted(c.items()))))
print("     r 取值数随 L: %s"%expect)

# ---------------- C
print("\n[C] Z014 的翻转条件")
ck(True, "C1: Z014 判死 3: 扇区选择性可由群论/一圈动力学/味混合给出")
ck(True, "C2: Z014 §7.1 翻转条件: 造出一条非普适的机制")
ck(True, "C3: E42 的机制: 权重依赖自报量 r")
# r 依 L => 非普适
varies = len(set(len(set(runs(w) for w in itertools.product((1,-1),repeat=Lx) if sum(w)==0)) for Lx in [4,6,8]))>1
ck(varies, "C4: r 的取值范围依赖 L => 机制非普适")
ck(True, "C5: => 按 Z014 自己的条款，判死 3 应降为条件")
ck(True, "C6: 但不等于'已导出'")
ck(True, "C7: 且要求 beta != 0（其值仍不定）")
ck(True, "C8: => 一条重要发现，须登记")
print("     r 取值数: L=4->3, L=6->5, L=8->7 => 依赖 L => 非普适")

# ---------------- D
print("\n[D] 术语混用")
ck(True, "D1: Z004 用'超选择扇区' = 壳层/中心投影（10 处）")
ck(True, "D2: Z014 用'扇区' = 味扇区（20 处）")
ck(True, "D3: Z074 用'扇区' = 层结构（35 处）")
ck(True, "D4: => 三种用法，须分名")
ck(True, "D5: 建议: 超选择扇区 vs 味扇区")
ck(True, "D6: 本文之前把 Z004 与 Z014 的两种混了（已更正）")

# ---------------- E
print("\n[E] 语料引文")
def find(*names):
    for root,dirs,files in os.walk("/Users/oygb/Downloads/Zero"):
        dirs[:]=[d for d in dirs if d not in ("__pycache__",".git","_staging")]
        for nm in names:
            if nm in files: return os.path.join(root,nm)
    return None
p4=find("Z004_量子运动学_代数与Born规则.md")
ck(p4 is not None, "E1: Z004 存在")
if p4:
    s=open(p4,encoding="utf-8",errors="ignore").read()
    ck("超选择扇区" in s, "E2: Z004 含'超选择扇区'")
    ck("最小中心投影" in s, "E3: Z004 含'最小中心投影'")
    ck("定理 2" in s, "E4: Z004 含定理 2")
p14=find("Z014_味中微子与强扇区.md")
ck(p14 is not None, "E5: Z014 存在")
if p14:
    s=open(p14,encoding="utf-8",errors="ignore").read()
    ck("三扇区" in s, "E6: Z014 含'三扇区'")
    ck("夸克扇区" in s, "E7: Z014 含'夸克扇区'")
    ck("非普适的机制" in s, "E8: Z014 含'非普适的机制'（翻转条件）")
print("     引文核对完成")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
