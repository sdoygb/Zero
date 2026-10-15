#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_delta0_search.py —— 从非中心态(2:1, seed 比)逼 δ₀ = 0.0396 rad
======================================================================================
目标: 电子偏离零点的量 δ₀ = 3π/4 − δ ≈ 0.0396 rad (0.039578 精确反解).

只枚举"有物理意义的自然组合"(非中心态量 + 基础数), 不海量凑数.
诚实报告: 命中(精确) / 近似(几%) / 无(per-界 数据).
"""
from __future__ import annotations
import os, json, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

D0 = 0.039578   # 精确反解: 3π/4 − δ, δ=132.7324°
P_("=" * 96)
P_("从非中心态(2:1, seed 比)逼 δ₀=0.0396 rad · k4_delta0_search.py")
P_("=" * 96)
P_(f"目标 δ₀ = {D0:.6f} rad = {np.degrees(D0):.4f}°")
P_(f"非中心态: 稳态 ρ=(2/3,1/3) [2:1]; seed=[8,4,1.5,1]; 相邻比 [2, 8/3, 3/2]")
P_("")

out = {"target": D0}
cands = []  # (rel_err, label, value)

def add(label, val):
    if val == 0 or not np.isfinite(val):
        return
    rel = abs(val - D0) / D0
    cands.append((rel, label, float(val)))

# ============ 类别 1: 1/(k·π), m/(k·π), 基础分数 ============
P_("【类别 1: π 相关 + 基础分数】")
for k in range(1, 25):
    add(f"1/({k}π)", 1/(k*np.pi))
for k in range(2, 13):
    add(f"1/({k}·2π)", 1/(k*2*np.pi))
for m in [1,2,3,5,7,11,13]:
    for n in range(2, 41):
        add(f"{m}/{n}", m/n)

# ============ 类别 2: log(非中心比) 相关 ============
P_("\n【类别 2: log 结构 (模流频率)】")
logs = {"log2": np.log(2), "log(3/2)": np.log(1.5), "log(8/3)": np.log(8/3),
        "槽差0.288": 0.288, "log3": np.log(3), "log(4/3)": np.log(4/3)}
for name, v in logs.items():
    for k in [2,3,4,6,8,12,16,24,17.5]:
        add(f"{name}/{k}", v/k)
        add(f"{name}/({k}π)", v/(k*np.pi))

# ============ 类别 3: 非中心度量的 arctan/arcsin 与平方 ============
P_("\n【类别 3: 非中心度量 (p−q)/(p+q) 等】")
nc = (2-1)/3   # 1/3 相对偏离
for nm, v in [("arctan(1/3)", np.arctan(nc)), ("arcsin(1/3)", np.arcsin(nc)),
              ("arctan(1/8)", np.arctan(1/8)), ("arcsin(1/8)", np.arcsin(1/8)),
              ("arctan(1/√2)", np.arctan(1/np.sqrt(2)))]:
    add(nm, v)

# ============ 类别 4: seed/相邻比 组合 ============
P_("\n【类别 4: seed 相邻比组合】")
seed = [8, 4, 1.5, 1]
vals = [2, 3, 4, 8, 8/3, 3/2, 2/3, 1/3, np.sqrt(2), 8/4, 4/1.5, 1.5]
for a in vals:
    for b in vals:
        for op in ["*", "/", "-"]:
            try:
                if op == "*": v = a*b
                elif op == "/": v = a/b
                else: v = abs(a-b)
                add(f"{a:.4g}{op}{b:.4g}", v)
            except: pass

# ============ 排序 + 报告 ============
cands.sort()
P_("\n【Top 12 最接近 δ₀ 的自然候选】")
for rel, label, val in cands[:12]:
    flag = "✓精确" if rel < 1e-3 else ("△近似" if rel < 0.01 else "·")
    P_(f"  {flag} {label:22s} = {val:.6f}   (偏差 {rel*100:.2f}%)")
P_("")
out["candidates"] = [{"label": l, "value": v, "rel_err_pct": round(r*100,2)} for r, l, v in cands[:12]]

# 结论
best = cands[0]
P_("=" * 96)
P_("【诚实结论】")
if best[0] < 1e-3:
    P_(f"  ✓ 找到精确自然候选: {best[1]} = {best[2]:.6f}")
elif best[0] < 0.01:
    P_(f"  △ 最接近: {best[1]} = {best[2]:.6f} (偏差 {best[0]*100:.2f}%) —— 近似, 非精确")
else:
    P_(f"  ✗ 无精确自然候选; 最好 {best[1]}={best[2]:.6f} (偏差 {best[0]*100:.2f}%)")
    P_(f"  ⟹ δ₀ 可能是 per-界 数据 (cosmos 判定), 2:1/seed 简单闭式逼不出")
P_(f"  特别检查: 1/(8π) = {1/(8*np.pi):.6f} (seed 8 意义; 偏差 {(abs(1/(8*np.pi)-D0)/D0)*100:.2f}%)")

os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_delta0_search.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_delta0_search.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_delta0_search.json / .txt")
