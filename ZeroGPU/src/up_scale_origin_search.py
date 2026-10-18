#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""up_scale_origin_search.py —— 上型标度 u_up ≈ 72 的来历搜索。

目标：
  u_up（拟合）= 71.6020534；
  有效振幅 A_up·√u_up = 15.0437（上型扇区的"实际大小"）。
做法：用**角色原子**（每个在库里有出处）搜单原子幂与三互异原子式，残差 <1%；
并把 72 的几种结构读法（8·9、24·3、27·(8/3)、16·(9/2)）列出来对照。

判据语义：FAIL = 表中某一步不成立（要报的判定）。
输出：ZeroGPU/results/up_scale_origin_search.json
"""
from __future__ import annotations

import json
import math
import os
import sys
from itertools import combinations, product

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


U_FIT = 71.6020533976
EFF = 1.7778326303083105 * math.sqrt(U_FIT)
ATOMS = {
    "2": 2.0, "3": 3.0, "7": 7.0, "8": 8.0, "12": 12.0, "14": 14.0, "16": 16.0,
    "24": 24.0, "27": 27.0, "28": 28.0, "72": 72.0,
    "4/3": 4 / 3, "3/2": 3 / 2, "8/3": 8 / 3, "2/9": 2 / 9, "1/2": 0.5,
    "π": math.pi, "√2": math.sqrt(2), "√3": math.sqrt(3),
}
print("=" * 84)
print(f"上型标度：u_up = {U_FIT:.6f}；有效振幅 A·√u = {EFF:.6f}")
print("=" * 84)


def search(target, tol=1e-2):
    one = []
    for nm, v in ATOMS.items():
        for e in range(-6, 7):
            if e == 0:
                continue
            x = v ** e
            rel = abs(x / target - 1)
            if rel < tol:
                one.append((rel, nm, e, x))
    one.sort()
    tri = []
    names = list(ATOMS)
    for combo in combinations(names, 3):
        for exps in product((-2, -1, 1, 2), repeat=3):
            x = 1.0
            for nm, e in zip(combo, exps):
                x *= ATOMS[nm] ** e
            rel = abs(x / target - 1)
            if rel < tol:
                tri.append((rel, combo, exps, x))
    tri.sort()
    return one, tri


for label, target in (("u_up", U_FIT), ("A·√u", EFF)):
    one, tri = search(target)
    print(f"\n  {label} = {target:.6f}（容差 1%）")
    print(f"    单原子幂：{len(one)} 个 hit" + ("" if not one else f"，最佳 {one[0][1]}^{one[0][2]:+d} = {one[0][3]:.6f}（{100*one[0][0]:+.3f}%）"))
    print(f"    三互异原子：{len(tri)} 个 hit")
    for rel, combo, exps, x in tri[:6]:
        expr = "·".join(f"{n}^{e:+d}" for n, e in zip(combo, exps))
        print(f"      {expr:34s} = {x:.6f}   差 {100*rel:+.4f}%")

print("\n  72 的几种结构读法：")
for nm, v in (("8·9", 8 * 9), ("24·3", 24 * 3), ("27·(8/3)", 27 * (8 / 3)),
              ("16·(9/2)", 16 * 4.5)):
    print(f"    {nm:10s} = {v:.6f}   与 u_up 差 {100*(v/U_FIT-1):+.3f}%")

one_u, tri_u = search(U_FIT, tol=1e-2)
check("u_up 的**单原子**最佳候选是 72（+0.556%）",
      len(one_u) == 1 and one_u[0][1] == "72", f"{one_u[0][1]} {100*one_u[0][0]:+.3f}%")
check("三因子搜索**没有判别力**（1% 内命中 ≥10 个 ⟹ 残差小不构成证据）",
      len(tri_u) >= 10, f"1% 内 {len(tri_u)} 个命中，最佳 {100*tri_u[0][0]:+.3f}%")
check("A·√u 的 1% 内候选：15.0437 与 15 = 16−1 差 0.29%",
      abs(EFF / 15 - 1) < 0.005, f"{100*(EFF/15-1):+.3f}%")
print()
print("  判定：u_up 的最佳简单候选是 72（8·9 或 24·3 或 27·(8/3)），差 0.55%；")
print("        有效振幅 15.0437 接近 15 = 16−1（差 0.29%），但两者都只是读法、无机制。")

res = {"u_fit": U_FIT, "eff": EFF, "one_u": one_u[:5], "tri_u": tri_u[:10],
       "pass": PASS, "fail": FAIL}
with open(os.path.join(OUT, "up_scale_origin_search.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print()
print("=" * 84)
print(f"汇总：PASS {len(PASS)} / FAIL {len(FAIL)}")
if FAIL:
    print("  FAIL:", FAIL)
print("  → results/up_scale_origin_search.json")
sys.exit(0 if not FAIL else 1)
