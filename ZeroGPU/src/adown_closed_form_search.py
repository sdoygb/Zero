#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""adown_closed_form_search.py —— A_down² 的结构表示搜索（换一组结构原子再试）。

目标：A_down² = 7.197683598（拟合值）。它已有候选 36/5 = 7.2（残差 +0.032%），
但 36/5 = 6²/5 里的 5 在库里**没有结构角色**。本脚本用"有角色的原子"重搜：

角色原子（每个都在库里有出处）：
  2（ℤ₂）、3（代/色）、7（满轨道数）、8（轨道数/8 槽载体）、16（槽）、
  27（E6 的 27）、28（C(8,2) 槽对子）、72（8·9）、12（相位格）、14（Heawood）、24（|S4| 相格）、
  4/3（亚层比公比）、3/2、8/3（亚层比）、2/9（δ₀）、1/2（半槽）、π、√2、√3

搜索：三因子式 a^p·b^q·c^r（p,q,r ∈ {−2,−1,1,2}，原子互不相同），残差 <0.1%。
判据语义：FAIL = 表中某一步不成立（要报的判定）。
输出：ZeroGPU/results/adown_closed_form_search.json
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


TARGET = 2.682849902237757 ** 2
ATOMS = {
    "2": 2.0, "3": 3.0, "7": 7.0, "8": 8.0, "12": 12.0, "14": 14.0, "16": 16.0,
    "24": 24.0, "27": 27.0, "28": 28.0, "72": 72.0,
    "4/3": 4 / 3, "3/2": 3 / 2, "8/3": 8 / 3, "2/9": 2 / 9, "1/2": 0.5,
    "π": math.pi, "√2": math.sqrt(2), "√3": math.sqrt(3),
}
print("=" * 80)
print(f"A_down² = {TARGET:.9f}（候选 36/5 = 7.2，残差 +0.032%）的结构表示搜索")
print("=" * 80)
print(f"  角色原子：{', '.join(sorted(ATOMS))}")

def search(target):
    """两套搜索：单原子幂 |e|≤6；三互异原子 |e|≤2。返回 (单原子命中, 三因子命中)。"""
    one = []
    for nm, v in ATOMS.items():
        for e in range(-6, 7):
            if e == 0:
                continue
            x = v ** e
            rel = abs(x / target - 1)
            if rel < 1e-3:
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
            if rel < 1e-3:
                tri.append((rel, combo, exps, x))
    tri.sort()
    return one, tri


ALL = {}
for sec, target in (("up", 1.7778326303083105 ** 2), ("down", TARGET)):
    one, tri = search(target)
    ALL[sec] = {"target": target, "one": one[:5], "tri": tri[:5]}
    print(f"\n  ── {sec}: 目标 {target:.9f} ──")
    print(f"     单原子幂命中 {len(one)} 个：")
    for rel, nm, e, x in one[:5]:
        print(f"       {nm}^{e:+d} = {x:.9f}   残差 {100*rel:+.4f}%")
    print(f"     三互异原子命中 {len(tri)} 个：")
    for rel, combo, exps, x in tri[:5]:
        expr = "·".join(f"{n}^{e:+d}" for n, e in zip(combo, exps))
        print(f"       {expr:34s} = {x:.9f}   残差 {100*rel:+.4f}%")

check("36/5 = 6²/5 需要原子 5 与 6——两者都不在角色原子表里",
      "5" not in ATOMS and "6" not in ATOMS)
print("\n  读数：上型的角色原子命中是 (4/3)^4（残差 −0.006%）；下型的角色原子命中是 16·√2/π（+0.067%），")
print("  而数字候选 36/5 更近（+0.032%）但 5 无角色。⟹ 下型仍是『数字近、读法无』。")

res = {"numbers": {"up": 1.7778326303083105 ** 2, "down": TARGET},
       "role_atom_hits": ALL, "pass": PASS, "fail": FAIL}
with open(os.path.join(OUT, "adown_closed_form_search.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print()
print("=" * 80)
print(f"汇总：PASS {len(PASS)} / FAIL {len(FAIL)}")
if FAIL:
    print("  FAIL:", FAIL)
print("  → results/adown_closed_form_search.json")
sys.exit(0 if not FAIL else 1)
