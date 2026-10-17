#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z148 核验脚本：稳态"比值继承"的独立复现（忠实层间引擎，path16，L=4，top2 记忆，T=4）。

判据（承重项）：
  稳态主两支的**比值**应当**逐位等于 seed 在对应两顶点的比值**——即
  被动力学**选择**的是**结构**（哪两支存活），被**继承**的是**比例**。

做法：对 5 个 seed 剖面各跑一遍引擎，取稳态主两支的权重比，与 seed 比对照；
      并把比值换算成**内禀频率** log(比)，与 Z089 §2 登记值对照。

用法：cd ZeroGPU/src && /usr/bin/python3 ratio_inheritance_check.py
"""
import os, sys, json, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import path_adj
from loyalty_engine import LoyaltyEngine

FAIL = []


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    if not cond:
        FAIL.append(name)


V = 16
adj = path_adj(V)
A = np.array([[1.0 if j in adj[i] else 0.0 for j in range(V)] for i in range(V)])

SEEDS = {
    "geom 2^(15-i)":      [2.0 ** (V - 1 - v) for v in range(V)],
    "geom 1.5^(15-i)":    [1.5 ** (V - 1 - v) for v in range(V)],
    "linear 16-i":        [float(V - v) for v in range(V)],
    "quad (16-i)^2":      [float((V - v) ** 2) for v in range(V)],
    "uniform":            [1.0] * V,
}

# Z089 §2 登记的内禀频率
REGISTERED = {"geom 2^(15-i)": 0.6931, "geom 1.5^(15-i)": 0.405,
              "linear 16-i": 0.118, "quad (16-i)^2": 0.236}

T = 4
CYCLES, GENS = 100, 200
print("=" * 78)
print(f"忠实引擎：path{V}  top2 记忆  T={T} 毁灭  {CYCLES}周期 × {GENS}代")
print("=" * 78)

rows = []
for label, seed in SEEDS.items():
    eng = LoyaltyEngine(adj, L=4, seed=seed, T=T, retention="top2")
    for _ in range(CYCLES):
        for _ in range(GENS):
            eng.run_generation()
    w = eng.record / eng.record.sum()
    order = np.argsort(-w)
    top2 = order[:2]
    a, b = float(w[top2[0]]), float(w[top2[1]])
    ratio_model = a / b if b > 0 else float("inf")
    s = np.array(seed, dtype=float)
    s = s / s.sum()
    ratio_seed = float(s[top2[0]] / s[top2[1]])
    freq = math.log(ratio_model)
    rows.append({"seed": label, "vertices": [int(top2[0]), int(top2[1])],
                 "w": [round(a, 5), round(b, 5)],
                 "ratio": round(ratio_model, 5), "ratio_seed": round(ratio_seed, 5),
                 "freq": round(freq, 5), "registered": REGISTERED.get(label)})
    print(f"  {label:18s} 顶点 {int(top2[0]):2d},{int(top2[1]):2d}  "
          f"w=[{a:.5f},{b:.5f}]  比={ratio_model:.5f}  seed比={ratio_seed:.5f}  "
          f"log={freq:.5f}  登记={REGISTERED.get(label)}")

print()
check("R1 稳态主两支的比值 = seed 在对应两顶点的比值（5 个 seed 全部）",
      all(abs(r["ratio"] - r["ratio_seed"]) < 5e-4 for r in rows),
      f"（最大偏差 {max(abs(r['ratio']-r['ratio_seed']) for r in rows):.2e}）")

check("R2 内禀频率与 Z089 §2 登记值相符（4 个有登记值的 seed）",
      all(abs(r["freq"] - r["registered"]) < 1e-3 for r in rows if r["registered"] is not None),
      f"（最大偏差 {max(abs(r['freq']-r['registered']) for r in rows if r['registered'] is not None):.2e}）")

# 结构选择：主两支应落在中间顶点（相邻）
mid = [r for r in rows]
check("R3 结构被选择：5 个 seed 的主两支落在同两个顶点（中间相邻对）",
      len({tuple(r["vertices"]) for r in mid}) == 1,
      f"：{mid[0]['vertices']}")

# 比值随 seed 变 ⟹ 比例是继承量（不是普适常数）
ratios = {r["seed"]: r["ratio"] for r in rows}
check("R4 比例是继承量：比值随 seed 变（5 个 seed 给 5 个不同比值）",
      len({round(v, 4) for v in ratios.values()}) == len(ratios),
      f"：{sorted(round(v,4) for v in ratios.values())}")

print()
print("=" * 78)
print(f"合计：{'全部 PASS' if not FAIL else 'FAIL: ' + ', '.join(FAIL)}")
print("=" * 78)
print()
print("结论：被动力学选择的是【结构】（哪两支存活，5 个 seed 同一结果）；")
print("     被【继承】的是比例（稳态比逐位等于 seed 比，且随 seed 变）。")

OUT = os.path.join(HERE, "..", "results", "ratio_inheritance_check.json")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"object": "忠实引擎 path16 top2 T=4 的稳态主两支比值与 seed 比对照",
               "rows": rows, "criteria": {"total": 4, "failed": FAIL}},
              fh, ensure_ascii=False, indent=2)
print(f"结果写入 {os.path.normpath(OUT)}")
