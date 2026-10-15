#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_delta0_remnant.py —— 从"δ 非满射→剩余"逼 δ₀ (归零下限) 的多模型探索
======================================================================================
机制(最老理论"零之动"): δ 非满射 → 留下剩余 → 归零有下限 δ₀.
本脚本在 Zero 底座上试多种"剩余/迭代/模流"模型 → δ₀=0.0396. 诚实报告.
"""
from __future__ import annotations
import os, json, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

D0 = 0.039578
P_("=" * 96)
P_("从'δ 非满射→剩余'逼 δ₀ · k4_delta0_remnant.py")
P_("=" * 96)
P_(f"目标 δ₀ = {D0:.6f} rad = {np.degrees(D0):.4f}°")
P_("机制: δ 非满射 → 剩余(归不到0的量); 无数迭代 → 非中心态是函数")
P_("")

cands = []
def add(label, val):
    if not np.isfinite(val) or val <= 0:
        return
    cands.append((abs(val-D0)/D0, label, float(val)))

# Zero 基础量
q = 5/9            # 稳态纯度 (2/3)²+(1/3)²
p2, p1 = 2/3, 1/3   # 两支
log2 = np.log(2); l32 = np.log(3/2); l83 = np.log(8/3)
Hmod = np.array([np.log(3/2), np.log(3)])   # 模流谱
tau = 2*np.pi

# ============ A. 剩余率 (δ 非满射) ============
P_("【A. 剩余率: δ 非满射 → 剩余 = 归不到0的量】")
add("1−q (稳态剩余重量)", 1-q)                 # 4/9
add("√(1−q)", np.sqrt(1-q))
add("(1−q)/(1+q)", (1-q)/(1+q))
add("小支 1/3", p1)
add("单支偏离中心 1/6", p2-1/2)
add("(2/3−1/3)/(2/3+1/3)", (p2-p1)/(p2+p1))    # 1/3
add("−log(1−q)", -np.log(1-q))
for f, nm in [(1-q, "4/9"), (1/3, "1/3"), (1/6, "1/6")]:
    for fn, fnn in [(np.arctan,"arctan"), (np.arcsin,"arcsin")]:
        add(f"{fnn}({nm})", fn(f))
        add(f"{fnn}({nm})²", fn(f)**2)
P_("")

# ============ B. 迭代收缩 / Bott 8 步 ============
P_("【B. 迭代收缩 / Bott 周期 8 步】")
for k in range(1, 33):
    add(f"1/{k} (迭代 k 步)", 1/k)
    add(f"log2/{k}", log2/k)
    add(f"1/{k}²", 1/k**2)
add("log2/2π", log2/tau)
add("log2/(8·2π)", log2/(8*tau))
add("1/(8·2π)", 1/(8*tau))
add("1/(2π·√8)", 1/(tau*np.sqrt(8)))
P_("")

# ============ C. 模流衰减 ============
P_("【C. 模流谱衰减】")
for h, nm in [(Hmod[0], "log3/2"), (Hmod[1], "log3"), (log2, "log2"),
              (Hmod[1]-Hmod[0], "log3−log3/2")]:
    add(nm, h)
    add(f"{nm}/(2π)", h/tau)
    add(f"{nm}/8", h/8)
    add(f"{nm}²/(2π)", h**2/tau)
    add(f"e^({nm})", np.exp(h))
add("(e^{log2}−1)/(e^{log2}+1)", (np.exp(log2)-1)/(np.exp(log2)+1))  # 1/3
P_("")

# ============ D. seed 三项组合 ============
P_("【D. seed 相邻比 [2, 8/3, 3/2] 组合】")
s = [2.0, 8/3, 3/2]
for a, b in itertools.combinations(s, 2):
    add(f"|{a:.4g}−{b:.4g}|", abs(a-b))
    add(f"|{a:.4g}−{b:.4g}|/{a+b:.4g}", abs(a-b)/(a+b))
    add(f"min/max", min(a,b)/max(a,b))
add("1/(8·3)", 1/24)
add("(1/3)·(1/8)", 1/24)
add("(1/3)²·(1/8)", (1/3)**2/8)
add("(1/6)²", (1/6)**2)
P_("")

# ============ E. 自洽: δ₀ 与 Koide/质量的关系 ============
P_("【E. 自洽/极值: δ₀ 由某条件定】")
# 电子质量最小 → δ₀ 是"最小偏离". 三代相角 δ+2πi/3, e 最小偏离
# δ₀ 若 = 三代偏离的几何平均? 
me, mmu, mtau = 0.511, 105.66, 1776.86
add("√(me/mτ)", np.sqrt(me/mtau))
add("me/mτ", me/mtau)
add("(1/3)·me 相关", 0)
# 零和 S1+S2+S3=0 (三层自指): 三代质量关系
P_("  (三代质量 Koide 已定 δ; δ₀ 是位置参数, 不与质量直接定 δ₀)")
P_("")

# ============ 排序报告 ============
cands = [c for c in cands if c[1] != "0"]
cands.sort()
P_("【Top 15 最接近 δ₀ 的候选】")
for rel, label, val in cands[:15]:
    flag = "✓" if rel < 1e-3 else ("△" if rel < 0.01 else "·")
    P_(f"  {flag} {label:32s} = {val:.7f}  (偏差 {rel*100:.2f}%)")
out = {"target": D0, "candidates": [{"label": l, "value": v, "rel_err_pct": round(r*100,2)} for r,l,v in cands[:15]]}
best = cands[0]
P_("=" * 96)
P_("【诚实判断】")
if best[0] < 1e-3:
    P_(f"  ✓ 精确命中: {best[1]} = {best[2]:.6f}")
elif best[0] < 0.01:
    P_(f"  △ 最接近: {best[1]} = {best[2]:.6f} (偏差 {best[0]*100:.2f}%)")
else:
    P_(f"  ✗ 无精确; 最好 {best[1]}={best[2]:.6f} ({best[0]*100:.2f}%)")
    P_("  ⟹ δ₀ 非简单'剩余率/迭代/模流'闭式; 需更深机制或保持 per-界 数据")
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_delta0_remnant.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_delta0_remnant.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_delta0_remnant.json / .txt")
