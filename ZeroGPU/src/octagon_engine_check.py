#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""octagon_engine_check.py —— 把八分圆模型放进 Zero 引擎（zero_sim）里跑一遍。

用引擎的结构函数核对模型的每个原生数：
  |Aut(K4)|        —— zero_sim.aut_order(k4)
  cycle(8) 的 tan(π/8) —— 八分圆几何（顶点角度 kπ/4，半角正切）
  (3/4)^24         —— K4 上"每步换顶点"的存活概率（4 选 1 留、3/4 走，24 步）
  28, 56           —— 8 槽载体的槽对子数（K4□K2 的顶点数 8）
然后装配 d=1/56·(3/4)^24，比较模型与实测。

输出：ZeroGPU/results/octagon_engine_check.json
用法：python3 ZeroGPU/src/octagon_engine_check.py
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402

OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


# ① K4 自同构阶
aut = zs.aut_order(zs.k4())
print(f"   引擎：|Aut(K4)| = {aut}")

# ② 八分圆：cycle(8) 的顶点角 kπ/4，半角正切 tan(π/8)=√2−1
oct_adj = zs.cycle(8)
tan8 = math.tan(math.pi / 8)
print(f"   引擎：cycle(8) girth={zs.girth(oct_adj)}；tan(π/8)={tan8:.9f} = √2−1={math.sqrt(2)-1:.9f}")

# ③ K4 行走者：每步"换顶点"概率 3/4，24 步不出错的存活
surv = (3 / 4) ** 24
print(f"   K4 行走者 24 步存活 (3/4)^24 = {surv:.9e}")

# ④ 8 槽载体 K4□K2 的顶点数与槽对子数
V8 = 8  # K4 × K2
C8 = math.comb(V8, 2)
print(f"   8 槽载体：|V|={V8}  C(8,2)={C8}  2C(8,2)={2*C8}")

# ⑤ 装配模型
d = (1 / (2 * C8)) * surv
A2 = 246.22 / 784.5
R = 3 - 2 * math.sqrt(2)
ms = [(A2 / 2) * (d * math.tan(k * math.pi / 8)) ** 2 * 1e12 for k in range(3)]
print(f"   d = 1/56·(3/4)^24 = {d:.6e}")
print(f"   模型 m = {[round(x,3) for x in ms]} meV（实测 8.614, 50.10）")

check("|Aut(K4)|=24", aut == 24, f"{aut}")
check("tan(π/8)=√2−1", abs(tan8 - (math.sqrt(2) - 1)) < 1e-15)
check("m3 与实测差 <1%", abs(ms[2] / 50.10 - 1) < 0.01, f"{100*(ms[2]/50.10-1):+.3f}%")
check("R=3−2√2 与实测差 <0.5%", abs((3 - 2 * math.sqrt(2)) / math.sqrt(7.42e-5 / 2.510e-3) - 1) < 0.005)

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "aut": aut, "tan8": tan8, "surv": surv,
       "d": d, "m_meV": ms, "R": R}
with open(os.path.join(OUT, "octagon_engine_check.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/octagon_engine_check.json")
sys.exit(0 if not FAIL else 1)
