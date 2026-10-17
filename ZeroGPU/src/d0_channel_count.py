#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""d0_channel_count.py —— D0 的通道数从载体图直接读出（推导稿 §0）。

载体（Z136）：
  Γ8  = K4□K2      （8 槽，第一层）
  Γ16 = K4□K2□K2   （16 槽，第二层）
通道数：
  一体（自配对，Majorana）：|V| 个槽位。
  两体（两个不同槽，Yukawa）：C(|V|,2) 个槽对子；质量单元是振幅²，故取 (C(|V|,2))²。

检查：由 |V| 读出的 N 与链里用的 N 一致；两个 D0 归一化数值正确。

输出：ZeroGPU/results/d0_channel_count.json
用法：python3 ZeroGPU/src/d0_channel_count.py
"""
from __future__ import annotations

import itertools
import json
import math
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


def carrier(dims):
    """K4□K2^dims 的顶点集（只数顶点，用于计数）。"""
    return [(v,) + a for v in range(4) for a in itertools.product([0, 1], repeat=dims)]


V8 = carrier(1)
V16 = carrier(2)
print(f"   |V(Γ8)|={len(V8)}   |V(Γ16)|={len(V16)}")

n8, n16 = len(V8), len(V16)
C8 = math.comb(n8, 2)
C16 = math.comb(n16, 2)
N_mass = C8**2          # 两体 × 两腿
N_MR = n16              # 一体
print(f"   C(8,2)={C8}  →  质量单元 N={N_mass}（=(C(8,2))²）")
print(f"   16 槽  →  M_R 的 N={N_MR}（一体）")

check("质量单元 N=784 由 C(|V(Γ8)|,2)² 读出", N_mass == 784, f"{N_mass}")
check("M_R 的 N=16 由 |V(Γ16)| 读出", N_MR == 16, f"{N_MR}")

d0_mass = N_mass / (N_mass + 0.5)
d0_MR = N_MR / (N_MR + 0.5)
print(f"   D0 质量单元因子 N/(N+1/2) = {d0_mass:.9f}  (1/784.5 的倒数关系)")
print(f"   D0 M_R 因子 = {d0_MR:.9f} = 32/33")
check("质量单元归一化 = 1/(N+1/2)", abs(1 / (N_mass + 0.5) - 1 / 784.5) < 1e-15)
check("M_R 因子 = 32/33", abs(d0_MR - 32 / 33) < 1e-15)

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "V8": n8, "V16": n16, "N_mass": N_mass,
       "N_MR": N_MR, "d0_mass": d0_mass, "d0_MR": d0_MR}
with open(os.path.join(OUT, "d0_channel_count.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/d0_channel_count.json")
sys.exit(0 if not FAIL else 1)
