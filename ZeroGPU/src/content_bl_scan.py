#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""content_bl_scan.py —— 全内容 B−L 谱扫描：有没有 B−L=±2 的 Δ_R 型对象。

内容（Z136 §3.5）：每代 16 槽（费米子，每槽一个 Weyl）＋ 玻色子 {γ, Z, W±, H, g}。
本脚本列出全内容的 B−L 取值，判定：
  · 是否存在 |B−L| = 2 的对象（Δ_R 型的色单态标量）；
  · 若不存在 ⟹ B−L 破缺没有载体 ⟹ 只能走"显式破缺"读法（尺度是输入，Z014 定理 16）。

判据语义：FAIL = 表中某一步不成立（要报的判定）。
输出：ZeroGPU/results/content_bl_scan.json
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


# ── 费米子（16 槽，Z136 §3.5）──
# 名称 -> (色, Y, T3, L)；B = 1/3（夸克）或 0；B−L = B − L
FERMIONS = {
    "ν_L": dict(color="单态", Y=-1, T3=+0.5, L=+1),
    "e_L": dict(color="单态", Y=-1, T3=-0.5, L=+1),
    "e_R": dict(color="单态", Y=-2, T3=0.0, L=+1),
    "ν_R": dict(color="单态", Y=0, T3=0.0, L=-1),
    "u_L^i": dict(color="三重", Y=+1 / 3, T3=+0.5, L=0),
    "d_L^i": dict(color="三重", Y=+1 / 3, T3=-0.5, L=0),
    "u_R^i": dict(color="三重", Y=+4 / 3, T3=0.0, L=0),
    "d_R^i": dict(color="三重", Y=-2 / 3, T3=0.0, L=0),
}
BOSONS = {
    "γ": dict(color="单态", B=0, L=0),
    "Z": dict(color="单态", B=0, L=0),
    "W±": dict(color="单态", B=0, L=0),
    "H": dict(color="单态", B=0, L=0),
    "g": dict(color="八重", B=0, L=0),
}

print("=" * 78)
print("全内容 B−L 谱")
print("=" * 78)
ferm_bl = {}
for nm, q in FERMIONS.items():
    B = 1 / 3 if q["color"] == "三重" else 0
    ferm_bl[nm] = B - q["L"]
    print(f"  {nm:<6} 色={q['color']:<3} Y={q['Y']:+.3f} T3={q['T3']:+.1f} L={q['L']:+.0f}"
          f"  B−L = {ferm_bl[nm]:+.3f}")
bos_bl = {}
for nm, q in BOSONS.items():
    bos_bl[nm] = q["B"] - q["L"]
    print(f"  {nm:<6} 色={q['color']:<3}（玻色子）  B−L = {bos_bl[nm]:+.3f}")

allobj = {**ferm_bl, **bos_bl}
mx = max(abs(v) for v in allobj.values())
print(f"  全内容 |B−L| 的最大值 = {mx:.3f}")
check("全内容里没有 |B−L| = 2 的对象（Δ_R 型不存在）", mx < 2, f"max |B−L| = {mx:.3f}")
check("玻色子（单覆盖）全部 B−L = 0", all(v == 0 for v in bos_bl.values()),
      str(sorted(set(bos_bl.values()))))
print()
print("  ⟹ B−L 破缺在内容里没有载体。")
print("  ⟹ 能用的是『显式破缺』读法：Majorana 项是显式项，系数不固定 ⟹ 尺度是输入（Z014 定理 16）。")
print("  ⟹ 翻转条件：若在 Zero 里造出一个 B−L=±2 的色单态单覆盖对象，第一读法（规范破缺）即打开。")

res = {"fermion_BmL": ferm_bl, "boson_BmL": bos_bl, "max_abs": mx,
       "pass": PASS, "fail": FAIL}
with open(os.path.join(OUT, "content_bl_scan.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print()
print("=" * 78)
print(f"汇总：PASS {len(PASS)} / FAIL {len(FAIL)}")
if FAIL:
    print("  FAIL:", FAIL)
print("  → results/content_bl_scan.json")
sys.exit(0 if not FAIL else 1)
