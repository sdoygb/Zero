#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scale_half_slot_check.py —— 质量单元候选 c_ℓ=1/(C(8,2)²+1/2)（推导稿 §3.3）。

候选（半槽修正）：c_ℓ = 1/(28² + 1/2) = 1/784.5，28 = 8 槽载体的槽对子数，
1/2 = 十六槽里唯一半满轨道的占据份额（与中微子的 m_D 因子同一个 1/2）。

检查：
  A. 三代带电轻子（用相位律的 f_i、v）与测量值的偏差；
  B. 与邻近候选 28²、28²±1/4、28²±3/4、28²+1 的对照（窗口宽度）；
  C. 由三个质量各自反解的 c_i（看形状残差与标度窗口）。

输出：ZeroGPU/results/scale_half_slot_check.json
用法：python3 ZeroGPU/src/scale_half_slot_check.py
"""
from __future__ import annotations

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


V = 246.22
M_DATA = [0.000510999, 0.1056584, 1.77686]
A = math.sqrt(2)
d0 = math.pi / 12 - 2 / 9
de = 3 * math.pi / 4 - d0
f = [1 + A * math.cos(de + 2 * math.pi * k / 3) for k in range(3)]
C = math.comb(8, 2)

c_cand = 1 / (C**2 + 0.5)
pred = [V * c_cand * fk**2 for fk in f]
devs = [100 * (p / d - 1) for p, d in zip(pred, M_DATA)]
print(f"   c=1/(28²+1/2)={c_cand:.9e}")
print(f"   预言（GeV）: {[round(x, 8) for x in pred]}")
print(f"   数据（GeV）: {M_DATA}")
print(f"   偏差%: {[round(x, 4) for x in devs]}")
check("三代偏差全部 <0.01%", max(abs(x) for x in devs) < 0.01,
      f"max={max(abs(x) for x in devs):.4f}%")

# c_i 反解
c_i = [m / (V * fk**2) for m, fk in zip(M_DATA, f)]
print(f"   各自反解 c_i: {[f'{x:.9e}' for x in c_i]}")
print(f"   窗口: {min(c_i):.9e} – {max(c_i):.9e}  相对宽 {100*(max(c_i)/min(c_i)-1):.4f}%")
check("候选落在三个 c_i 的窗口内", min(c_i) <= c_cand <= max(c_i))

print("\n   邻近候选：")
rows = []
for name, N2 in [("28²", C**2), ("28²+1/4", C**2 + 0.25), ("28²+1/2", C**2 + 0.5),
                 ("28²+3/4", C**2 + 0.75), ("28²+1", C**2 + 1)]:
    cc = 1 / N2
    dd = [100 * (V * cc * fk**2 / m - 1) for fk, m in zip(f, M_DATA)]
    rows.append({"form": name, "N2": N2, "max_dev_pct": max(abs(x) for x in dd)})
    print(f"     {name:9s} N²={N2:.3f}  最大偏差 {max(abs(x) for x in dd):.4f}%")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "c_candidate": c_cand, "N2": C**2 + 0.5,
       "pred_GeV": pred, "data_GeV": M_DATA, "dev_pct": devs,
       "c_i": c_i, "alternatives": rows}
with open(os.path.join(OUT, "scale_half_slot_check.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/scale_half_slot_check.json")
sys.exit(0 if not FAIL else 1)
