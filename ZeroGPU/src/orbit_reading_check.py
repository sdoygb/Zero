#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""orbit_reading_check.py —— D0 的读法在"通道级"还是"轨道级"？（链稿 §0.6）

Z136 的计数：8 条 Z2^iso 轨道 = 7 满 + 1 半（占据 1/2）。
两种读法：
  轨道级：满轨道占比 7/(7+1/2) = 14/15 = 0.9333
  通道级：满槽占比   16/(16+1/2) = 32/33 = 0.9697
用中微子标度判别：M_R = √(M_I M_X)·因子，看哪个因子的 m_ν3 与数据（50.10 meV）相容。

输出：ZeroGPU/results/orbit_reading_check.json
用法：python3 ZeroGPU/src/orbit_reading_check.py
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


V, M_Z, M_TAU = 246.22, 91.1876, 1.77696645
M_I = M_Z * math.exp(4 * math.pi)
M_X = 1.009e13
M_R0 = math.sqrt(M_I * M_X)
R = (73 - 28 * math.sqrt(6)) / 25

factors = {
    "轨道级 7/(7+1/2)=14/15": 7 / 7.5,
    "通道级 16/(16+1/2)=32/33": 16 / 16.5,
    "无修正 1": 1.0,
}
rows = []
for name, f in factors.items():
    MR = M_R0 * f
    m3 = (M_TAU / 2) ** 2 / MR * 1e12
    rows.append({"reading": name, "factor": f, "m_nu3": m3, "dev_pct": 100 * (m3 / 50.10 - 1)})
    print(f"   {name:28s} 因子={f:.6f}  m_ν3={m3:.2f} meV  偏差 {100*(m3/50.10-1):+.2f}%")

best = min(rows, key=lambda r: abs(r["dev_pct"]))
print(f"\n   数据最相容的是：{best['reading']}（{best['dev_pct']:+.2f}%）")
check("通道级（32/33）与数据相容（<0.5%）", abs(rows[1]["dev_pct"]) < 0.5,
      f"{rows[1]['dev_pct']:+.2f}%")
check("轨道级（14/15）被数据排除（>2%）", abs(rows[0]["dev_pct"]) > 2,
      f"{rows[0]['dev_pct']:+.2f}%")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "rows": rows, "best": best["reading"]}
with open(os.path.join(OUT, "orbit_reading_check.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/orbit_reading_check.json")
sys.exit(0 if not FAIL else 1)
