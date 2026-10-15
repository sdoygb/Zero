#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_delta0_floor.py —— δ₀ = 非中心态固有最小偏离 ("归零下限", "0 不重叠")
======================================================================================
用户图像: 非中心态让 0 的位置不重叠 → 归零有下限 δ₀, 怎样归都归不到 0,
         只能归到"宇宙不平衡那个点"(物理最底线) → δ₀ 应由非中心态(2/3,1/3)
         的固有结构唯一决定, 无自由参数.

探索: 非中心态固有量(权重偏离/模流谱/log) → 相位 δ₀=0.0396 rad.
只枚举有物理意义的自然组合, 不海量凑数.
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

D0 = 0.039578   # δ₀ 精确
P_("=" * 96)
P_("δ₀ = 非中心态固有最小偏离 ('归零下限') · k4_delta0_floor.py")
P_("=" * 96)
P_(f"目标 δ₀ = {D0:.6f} rad = {np.degrees(D0):.4f}°")
P_("非中心态: 稳态 ρ=(2/3,1/3), 中心 (1/2,1/2)")
P_("")

cands = []
def add(label, val):
    if not np.isfinite(val) or val <= 0:
        return
    cands.append((abs(val-D0)/D0, label, float(val)))

# ============ 非中心态固有量 (权重偏离) ============
w0 = 2/3 - 1/2   # = 1/6 每支偏离中心
wtot = (2/3-1/3)  # = 1/3 两支差
P_("【权重偏离 (固有)】")
add("1/6 (单支偏离中心)", w0)
add("1/3 (两支总差)", wtot)
add("1/6²", w0**2)
add("1/3²", wtot**2)
add("1/6³", w0**3)
for f, nm in [(w0, "1/6"), (wtot, "1/3"), (w0**2, "1/36"), (wtot**2, "1/9")]:
    for op, fn in [("arctan", np.arctan), ("arcsin", np.arcsin)]:
        add(f"{op}({nm})", fn(f))
P_("")

# ============ 模流/谱 固有量 ============
P_("【模流谱 / log 结构】")
log2 = np.log(2); l32 = np.log(3/2); l83 = np.log(8/3)
for nm, v in [("log2", log2), ("log(3/2)", l32), ("log(8/3)", l83),
              ("log2−log(3/2)", log2-l32), ("log3−log(3/2)", np.log(3)-l32)]:
    add(nm, v)
    add(f"{nm}/(2π)", v/(2*np.pi))
    add(f"{nm}/(8π)", v/(8*np.pi))
    add(f"{nm}/(π/2)", v/(np.pi/2))
    add(f"{nm}²", v**2)
P_("")

# ============ 组合: 权重偏离 × 闭合/seed 相关 ============
P_("【权重偏离 × 闭合(2π)/seed(8) 组合】")
for w, wnm in [(w0, "1/6"), (wtot, "1/3"), (w0**2, "1/36")]:
    for k, knm in [(np.pi, "π"), (2*np.pi, "2π"), (np.pi/2, "π/2"), (8, "seed8"), (4, "seed4"), (8*np.pi, "8π")]:
        add(f"{wnm}/({knm})", w/k)
        add(f"{wnm}·({knm}⁻¹)", w/k)
P_("")

# ============ 排序报告 ============
cands.sort()
P_("【Top 14 最接近 δ₀ 的自然候选】")
for rel, label, val in cands[:14]:
    flag = "✓" if rel < 1e-3 else ("△" if rel < 0.01 else "·")
    P_(f"  {flag} {label:28s} = {val:.7f}  (偏差 {rel*100:.2f}%)")
P_("")
out = {"target": D0, "candidates": [{"label": l, "value": v, "rel_err_pct": round(r*100,2)} for r,l,v in cands[:14]]}

best = cands[0]
P_("=" * 96)
P_("【诚实判断】")
if best[0] < 1e-3:
    P_(f"  ✓ 精确命中: {best[1]} = {best[2]:.6f}")
elif best[0] < 0.01:
    P_(f"  △ 最接近: {best[1]} = {best[2]:.6f} (偏差 {best[0]*100:.2f}%) —— 近似非精确")
else:
    P_(f"  ✗ 无精确; 最好 {best[1]}={best[2]:.6f} ({best[0]*100:.2f}%)")
    P_("  ⟹ 非中心态简单固有量映射不到 δ₀; 需更深的'0 不重叠'机制")
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_delta0_floor.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_delta0_floor.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_delta0_floor.json / .txt")
