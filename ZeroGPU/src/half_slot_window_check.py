#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""half_slot_window_check.py —— 半槽候选在自然计数族里的唯一性（推导稿 §3.4）。

族：N² = C(n,k)² + x，n≤16，x∈{0,±1/4,±1/2,±3/4,±1,±3/2}；也查 C(n,k)²+y·C(n,k)。
判据：c=1/N² 是否落在数据反解窗口 [min c_i, max c_i]（三个带电轻子各给一个 c_i）。

输出：ZeroGPU/results/half_slot_window_check.json
用法：python3 ZeroGPU/src/half_slot_window_check.py
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
M = [0.000510999, 0.1056584, 1.77686]
A = math.sqrt(2)
d0 = math.pi / 12 - 2 / 9
de = 3 * math.pi / 4 - d0
f = [1 + A * math.cos(de + 2 * math.pi * k / 3) for k in range(3)]
c_i = [m / (V * fk**2) for m, fk in zip(M, f)]
lo, hi = min(c_i), max(c_i)
print(f"   c_i = {[f'{x:.9e}' for x in c_i]}")
print(f"   窗口 [{lo:.9e}, {hi:.9e}]，宽 {100*(hi/lo-1):.4f}%")

xs = [0, 0.25, -0.25, 0.5, -0.5, 0.75, -0.75, 1.0, -1.0, 1.5, -1.5]
hits = []
for n in range(4, 17):
    for k in range(1, n):
        C = math.comb(n, k)
        for x in xs:
            N2 = C * C + x
            if N2 <= 0:
                continue
            c = 1 / N2
            if lo <= c <= hi:
                hits.append({"n": n, "k": k, "C": C, "x": x, "c": c})
print(f"   C(n,k)²+x 命中窗口：{len(hits)} 个")
for h in hits:
    print(f"     n={h['n']} k={h['k']} C={h['C']} x={h['x']:+.2f}  c={h['c']:.9e}")

hits2 = []
for n in range(4, 17):
    for k in range(1, n):
        C = math.comb(n, k)
        for y in xs:
            N2 = C * C + y * C
            if N2 <= 0:
                continue
            c = 1 / N2
            if lo <= c <= hi:
                hits2.append({"n": n, "k": k, "y": y, "c": c})
print(f"   C(n,k)²+y·C(n,k) 命中窗口：{len(hits2)} 个")

check("唯一命中是 (8,2)²+1/2（与对称等价的 (8,6) 同值）",
      len({(h["C"], h["x"]) for h in hits}) == 1 and hits[0]["C"] == 28 and hits[0]["x"] == 0.5,
      f"{[(h['n'],h['k'],h['x']) for h in hits]}")
check("C²+y·C 形式没有命中", len(hits2) == 0)

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "window": [lo, hi], "hits": hits, "hits2": hits2}
with open(os.path.join(OUT, "half_slot_window_check.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/half_slot_window_check.json")
sys.exit(0 if not FAIL else 1)
