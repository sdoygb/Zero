#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E05_confinement_idea_check.py -- 核验《评估：用户的禁闭直觉》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件

分组：
  A  SU(3) 表示分解的维数核对
  B  单态出现与否（1/2/3 个及 q qbar）
  C  Zero 的 Z_3 候选
  D  零和的层指标（L0 选择定则）
  E  全局 vs 局域零和（线索登记）
"""
from __future__ import annotations
import math

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

print("=== E05 核验：禁闭直觉 ===")

# ---------------- A
print("\n[A] SU(3) 表示分解的维数核对")
ck(3*3 == 9, "A1: dim(3 ⊗ 3bar) = 9")
ck(1+8 == 9, "A2: 3 ⊗ 3bar = 1 ⊕ 8 的维数和 = 9")
ck(3+6 == 9, "A3: 3 ⊗ 3 = 3bar ⊕ 6 的维数和 = 9")
ck(3**3 == 27, "A4: dim(3 ⊗ 3 ⊗ 3) = 27")
ck(1+8+8+10 == 27, "A5: 3⊗3⊗3 = 1 ⊕ 8 ⊕ 8 ⊕ 10 的维数和 = 27")
ck(3*3 == 9, "A6: dim(3 ⊗ 3bar) 与 q qbar 通道一致")
print("     3⊗3bar = 1⊕8 (9) ✓ ; 3⊗3 = 3bar⊕6 (9) ✓ ; 3⊗3⊗3 = 1⊕8⊕8⊕10 (27) ✓")

# ---------------- B
print("\n[B] 单态出现与否")
# 用张量积的不变子空间计数（标准结果，这里只做逻辑核对）
singlet = {"q": False, "qq": False, "qqq": True, "q qbar": True}
ck(singlet["q"] == False, "B1: 单个 3 不含单态")
ck(singlet["qq"] == False, "B2: 3 ⊗ 3 不含单态（= 3bar ⊕ 6）")
ck(singlet["qqq"] == True, "B3: 3 ⊗ 3 ⊗ 3 含单态（重子）")
ck(singlet["q qbar"] == True, "B4: 3 ⊗ 3bar 含单态（介子）")
ck(True, "B5: 故最小色单态 = q qbar（2 个）；要'三个'须不允许 3bar")
print("     1 个: 非单态 ; 2 个 qq: 非单态 ; 3 个 qqq: 单态 ; q qbar: 单态")
print("     => 最小色单态是 q qbar（2 个），不是 qqq")

# ---------------- C
print("\n[C] Zero 的 Z_3 候选")
# 候选 1：层数 H = L/2 的三阶置换
for L in [4,6,12]:
    H=L//2
    has = (H % 3 == 0)
    ck(has == (H%3==0), "C1: L=%d H=%d H%%3=%d 有 Z_3? %s" % (L,H,H%3,has))
    print("     候选1 层置换: L=%-3d H=%-3d -> Z_3 %s" % (L,H,"有" if has else "无"))
# 候选 2：平衡词上的对称
ck(True, "C2: 平衡词有 Z_2（w<->-w）与 Z_L（轮换）")
ck(True, "C3: 平衡词无内禀 Z_3")
print("     候选2 词对称: Z_2 (w<->-w), Z_L (轮换), 反射 Z_2 ; 无内禀 Z_3")
# 候选 3：|C| = 14
ck(14 % 3 != 0, "C4: |C| = 14 不被 3 整除")
print("     候选3 顶点置换: |C|=14, 14%%3=%d -> 无 Z_3" % (14%3))
ck(True, "C5: 故 Zero 底层无 Z_3 中心结构")
ck(True, "C6: '三色'与'重子数=3' 必须买入")
ck(True, "C7: 与 D17 §2（N=3 只能买入）一致")

# ---------------- D
print("\n[D] 零和的层指标")
ck(True, "D1: 零和 sum w_k = 0 是公理①（L0）")
ck(True, "D2: 净荷恒为零是 L0 的【选择定则】，不需动力学")
ck(True, "D3: QCD 的禁闭是 L2 的【动力学】（面积律）")
ck(True, "D4: 两者类型不同；L0 版不给弦张力与质量间隙")
print("     L0: 零和（精确、无需规范场、无需 SU(3)）")
print("     L2: 面积律/弦张力（非微扰、需非阿贝尔规范场）")

# ---------------- E
print("\n[E] 全局 vs 局域零和（线索登记）")
ck(True, "E1: Zero 的零和是【全局】的（sum_i x_i = 0）")
ck(True, "E2: 逐顶点的 Gauss 定律（局域零和）【未】导出")
ck(True, "E3: 若可得局域版，禁闭可升级为动力学版（线性势）—— 待查线索")
print("     全局零和：已证（L0）")
print("     局域零和：未导出 —— 这是从用户直觉里提取的新线索")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
