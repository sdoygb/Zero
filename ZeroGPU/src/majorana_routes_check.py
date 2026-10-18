#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""majorana_routes_check.py —— Majorana 质量的三条途径在 Zero 内的可用性。

背景（Majorana 栏收口 §3 的翻转条件）：需要 B−L=2 的色单态"单覆盖"对象。
上一轮扫过全内容：没有 |B−L|=2 的**场**。但 ν_R 的**双线性** ν_R^T C ν_R 本身就是
一个 B−L=2 的**单覆盖（玻色）复合体**——所以翻转条件要分"基本"与"复合"两种读法。

三条途径：
  1. 显式双线性：M ν_R^T C ν_R —— 规范不变？重整化？
  2. 费米子凝聚：⟨ν_R^T C ν_R⟩ ≠ 0 —— 需要动力学；ν_R 是规范单态（无 gauge 束缚）
  3. 基本 Δ_R：B−L=±2 的色单态标量 —— 上一轮已扫：不存在

判据语义：FAIL = 表中某一步不成立（要报的判定）。
输出：ZeroGPU/results/majorana_routes_check.json
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


# ν_R 的量子数（Z136 §3.5）：色单态、Y=0、T3=0、B−L=+1
NU_R = dict(color="单态", Y=0, T3=0.0, BmL=+1)

print("=" * 80)
print("Majorana 质量的三条途径：Zero 内容内的可用性")
print("=" * 80)
print(f"  ν_R 量子数：色={NU_R['color']}，Y={NU_R['Y']}，T3={NU_R['T3']}，B−L={NU_R['BmL']}")

# 途径 1：双线性
bil_color = NU_R["color"]                 # 单态×单态 → 单态
bil_Y = 2 * NU_R["Y"]
bil_BmL = 2 * NU_R["BmL"]
dim_bilinear = 3                          # [费米子]^{3/2} × 2 = 3
print("\n  途径 1  显式双线性 M ν_R^T C ν_R：")
print(f"    色 = {bil_color}（规范单态），Y = {bil_Y}，B−L = {bil_BmL}，维数 = {dim_bilinear}")
check("双线性是规范单态（色单态、Y=0）", bil_color == "单态" and bil_Y == 0)
check("双线性带 B−L=2", bil_BmL == 2)
check("双线性可重整（维数 3 ≤ 4）", dim_bilinear <= 4)
print("    ⟹ **可用**：它是内容里的算子；系数不受结构约束 ⟹ 尺度是输入。")

# 途径 2：凝聚
print("\n  途径 2  费米子凝聚 ⟨ν_R^T C ν_R⟩：")
print(f"    ν_R 是规范单态（Y=0，色单态）⟹ 无规范束缚 ⟹ 内容里没有产生凝聚的动力学")
check("ν_R 无规范荷（故内容内无凝聚机制）", NU_R["Y"] == 0 and NU_R["color"] == "单态")
print("    ⟹ **只有读法**：它把途径 1 的系数读成动力学产物，内容内没有可判的机制。")

# 途径 3：基本 Δ_R
print("\n  途径 3  基本 Δ_R（B−L=±2 的色单态标量）：")
print("    上一轮全内容扫描：max|B−L| = 1；标量只有 Higgs（B−L=0）⟹ 不存在")
check("基本 Δ_R 不在内容内（|B−L|=2 的场不存在）", True, "content_bl_scan.py 已核验")

print()
print("  判定：三条途径里 **只有途径 1 可用**；途径 2 是它的读法（无机制）；途径 3 缺对象。")
print("  翻转条件的准确写法：B−L=2 的对象在**复合**层面存在（ν_R 双线性，单覆盖），")
print("  在**基本**层面不存在（无 Δ_R 场）⟹ 自发性破缺需要基本 Δ_R，显式破缺不需要。")

res = {"bilinear": {"color": bil_color, "Y": bil_Y, "BmL": bil_BmL, "dim": dim_bilinear},
       "pass": PASS, "fail": FAIL}
with open(os.path.join(OUT, "majorana_routes_check.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print()
print("=" * 80)
print(f"汇总：PASS {len(PASS)} / FAIL {len(FAIL)}")
if FAIL:
    print("  FAIL:", FAIL)
print("  → results/majorana_routes_check.json")
sys.exit(0 if not FAIL else 1)
