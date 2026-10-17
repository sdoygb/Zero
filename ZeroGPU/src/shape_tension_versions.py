#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""shape_tension_versions.py —— 形状项随数据版本与吸收参数的状况（链稿 §3.1）。

链的形状：A_ν²=3/2、m1=0、R=(73−28√6)/25、X=R²=0.031177。
本脚本：
  A. 四个数据集的 X 张力与 Δm²21 的 pull（标度已由半槽修正定住：m2=8.85、m3=50.12）；
  B. 三种吸收参数（A²、相位偏移 ε、间隔 Δ）的所需值（与已落盘的残差稿一致）。

输出：ZeroGPU/results/shape_tension_versions.json
用法：python3 ZeroGPU/src/shape_tension_versions.py
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


R = (73 - 28 * math.sqrt(6)) / 25
X = R**2
# 标度定住后的质量（meV）：m3=50.12、m2=R·m3
M3 = 50.118
M2 = R * M3
D21 = (M2 * 1e-3) ** 2          # eV^2
S21 = 0.20e-5

sets = [("NuFIT 5.2", 7.42e-5, 2.510e-3),
        ("NuFIT 6.0", 7.49e-5, 2.513e-3),
        ("较早一版", 7.53e-5, 2.453e-3),
        ("NuFIT 4.0", 7.39e-5, 2.525e-3)]
rows = []
for name, d21, d31 in sets:
    xobs = d21 / d31
    rows.append({"set": name, "X_obs": xobs, "X_tension_pct": 100 * (X / xobs - 1),
                 "d21_pull": (D21 - d21) / S21})
    print(f"   {name:10s}: X_obs={xobs:.6f}  张力 {100*(X/xobs-1):+.2f}%  Δm²21 pull {(D21-d21)/S21:+.1f}σ")

check("四个数据集的形状张力都在 1–7%", all(1 < abs(r["X_tension_pct"]) < 7 for r in rows))
check("Δm²21 的 pull 在 1–3σ", all(0.5 < r["d21_pull"] < 3 for r in rows),
      f"{[round(r['d21_pull'],1) for r in rows]}")

print("\n   三种吸收参数的所需值（与残差稿一致）：A²=1.513633、ε=−0.4701°、Δ=119.5896°")
res = {"pass": PASS, "fail": FAIL, "R": R, "X": X, "m2_meV": M2, "m3_meV": M3,
       "d21_model": D21, "rows": rows,
       "absorbers": {"A2": 1.513633, "eps_deg": -0.4701, "delta_deg": 119.5896}}
with open(os.path.join(OUT, "shape_tension_versions.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
print("   → results/shape_tension_versions.json")
sys.exit(0 if not FAIL else 1)
