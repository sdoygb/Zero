#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""half_slot_rule_check.py —— 半槽进入计数的位置：D1 与 D3 是不是同一条规则。

两处用法：
  D1（质量单元）：c_ℓ = 1/(C(8,2)² + 1/2) = 1/784.5 —— 半槽加在**平方计数之外**；
  D3（Majorana 标度）：M_R ∝ 16/(16 + 1/2) = 32/33 —— 半槽加在**计数之内**（线性）。

本脚本检验三种"半槽规则"，看哪一种能同时给出 D1 与 D3：
  规则 A（加在计数内）：N → N + 1/2，再做原来的运算；
  规则 B（占据权重线性）：两体通道 = Σ_{i<j} w_i w_j，一个半槽 w=1/2；
  规则 C（加在计数外）：先算 N²，再加 1/2。

判据语义：FAIL = 表中某一步不成立（要报的判定）。
输出：ZeroGPU/results/half_slot_rule_check.json
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


N8, N16 = 8, 16
P28 = N8 * (N8 - 1) // 2
print("=" * 78)
print("半槽进入计数的位置：D1 与 D3 是否同一条规则")
print("=" * 78)
print(f"  D1 形式：1/(C(8,2)² + 1/2) = 1/({P28}² + 0.5) = 1/{P28**2 + 0.5}")
print(f"  D3 形式：16/(16 + 1/2) = {16/(16.5):.9f}（= 32/33 = {32/33:.9f}）")
check("D3 的 16/(16+½) 恰为 32/33", abs(16 / 16.5 - 32 / 33) < 1e-15,
      f"{16/16.5:.12f}")

# 规则 A：加在计数内，再做平方（D1 若用此规则）
A_d1 = (P28 + 0.5) ** 2
print(f"\n  规则 A（N→N+½ 再平方）用于 D1：({P28}+0.5)² = {A_d1}")
check("规则 A 给 D1 的不是 784.5（位置不同）", abs(A_d1 - (P28**2 + 0.5)) > 1,
      f"{A_d1} vs 784.5")

# 规则 B：占据权重线性（一个半槽 w=1/2，其余 1）
w = [1.0] * (N8 - 1) + [0.5]
B = sum(w[i] * w[j] for i in range(N8) for j in range(i + 1, N8))
print(f"  规则 B（占据权重乘积）两体通道 = {B}（满槽时 {P28}）")
check("规则 B 给的是 24.5，不是 28.5（故不是 D1 的来源）", abs(B - 24.5) < 1e-12,
      f"{B}")

# 规则 C：加在平方外（D1 的形式）
C_d1 = P28**2 + 0.5
print(f"  规则 C（先平方再加½）= {C_d1} ⟹ 与 D1 一致")
check("规则 C 复现 D1 的 784.5", abs(C_d1 - 784.5) < 1e-12)
# 同一条规则 C 用到 D3 会怎样？
C_d3 = N16**2 + 0.5
print(f"  规则 C 用到 D3 的计数：16² + 0.5 = {C_d3} ⟹ 与 D3 的 16.5 不同")
check("同一条规则不能同时给 D1 与 D3（位置不同）", abs(C_d3 - 16.5) > 1, f"{C_d3} vs 16.5")

print()
print("  判定：")
print("  · D1 把半槽加在**平方之外**（784 + 0.5）；D3 加在**计数之内**（16 + 0.5）。")
print("  · 三条候选规则里没有一条能同时给出两处；占据权重线性（24.5）与 D1 差 12.5%，被排除。")
print("  · 所以『同一个半槽』在对象层成立（§1.5 的第 1、2 面），但**计数规则不统一**——这是 D1 的缺口。")

res = {"D1": P28**2 + 0.5, "D3": 16 / 16.5, "ruleA_D1": A_d1, "ruleB": B,
       "ruleC_D1": C_d1, "ruleC_D3": C_d3, "pass": PASS, "fail": FAIL}
with open(os.path.join(OUT, "half_slot_rule_check.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print()
print("=" * 78)
print(f"汇总：PASS {len(PASS)} / FAIL {len(FAIL)}")
if FAIL:
    print("  FAIL:", FAIL)
print("  → results/half_slot_rule_check.json")
sys.exit(0 if not FAIL else 1)
