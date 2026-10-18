#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""amplitude_sublayer_reading.py —— 上/下型振幅能否用三个亚层比表示。

三个亚层比（升序）：r = (r1, r2, r3) = (3/2, 2, 8/3)，等比、公比 4/3（Z129）。
问：
  A. A_up = 16/9 是否 = r3/r1（＝(4/3)²）？
  B. A_up² = (4/3)^4 是否只是它的平方？
  C. A_down² = 7.197684 能否用 r 的乘积/商（小整数幂）表示到 0.1% 以内？

判据语义：FAIL = 表中某一步不成立（要报的判定）。
输出：ZeroGPU/results/amplitude_sublayer_reading.json
"""
from __future__ import annotations

import json
import math
import os
import sys
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


R = (3 / 2, 2.0, 8 / 3)          # 三个亚层比（升序）
A_UP = 16 / 9                     # 候选闭式（A 尺）
A_DOWN2 = 2.682849902237757 ** 2  # 拟合 A_down²
print("=" * 78)
print("上/下型振幅与三个亚层比：r = (3/2, 2, 8/3)，公比 4/3")
print("=" * 78)
print(f"  r3/r1 = {R[2]/R[0]:.9f}   (4/3)² = {(4/3)**2:.9f}   16/9 = {A_UP:.9f}")
check("r3/r1 = (4/3)² = 16/9（精确）",
      abs(R[2] / R[0] - (4 / 3) ** 2) < 1e-12 and abs((4 / 3) ** 2 - A_UP) < 1e-12)
check("A_up = 16/9 与拟合 1.777833 差 <0.01%",
      abs(100 * (A_UP / 1.7778326303083105 - 1)) < 0.01,
      f"{100*(A_UP/1.7778326303083105-1):+.4f}%")
check("A_up² = (4/3)^4 —— 只是 A_up 的平方",
      abs(A_UP ** 2 - (4 / 3) ** 4) < 1e-12, f"{A_UP**2:.9f}")
print(f"  A_up² = {A_UP**2:.9f}；(4/3)^4 = {(4/3)**4:.9f}")

# C：用 r 的小整数幂组合去找 A_down²
best = (1e9, None)
for e in product(range(-3, 4), repeat=3):
    if all(x == 0 for x in e):
        continue
    v = R[0] ** e[0] * R[1] ** e[1] * R[2] ** e[2]
    d = abs(v / A_DOWN2 - 1)
    if d < best[0]:
        best = (d, e, v)
print(f"\n  A_down² = {A_DOWN2:.9f}")
print(f"  最优组合 r1^{best[1][0]}·r2^{best[1][1]}·r3^{best[1][2]} = {best[2]:.9f}"
      f"  相对差 {100*best[0]:+.4f}%")
check("用 r 的小整数幂（|e|≤3）没有 0.1% 以内的 A_down² 表示",
      best[0] > 1e-3, f"最优 {100*best[0]:.4f}%")
print()
print("  判定：")
print("  · A_up = r3/r1 = (4/3)² —— 有结构读法（两个亚层比之比）；(4/3)^4 只是它的平方。")
print("  · A_down² = 36/5 在 r 的小整数幂组合里没有 0.1% 以内的表示 ⟹ 无亚层比读法。")

res = {"r": R, "A_up": A_UP, "A_up_is_r3_r1": abs(R[2] / R[0] - A_UP) < 1e-12,
       "A_down2": A_DOWN2, "best_combo": [best[1], best[2], best[0]],
       "pass": PASS, "fail": FAIL}
with open(os.path.join(OUT, "amplitude_sublayer_reading.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print()
print("=" * 78)
print(f"汇总：PASS {len(PASS)} / FAIL {len(FAIL)}")
if FAIL:
    print("  FAIL:", FAIL)
print("  → results/amplitude_sublayer_reading.json")
sys.exit(0 if not FAIL else 1)
