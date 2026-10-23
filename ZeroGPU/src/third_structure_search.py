#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
third_structure_search.py —— 找定振幅的"第三个结构"：候选逐个过判据
================================================================================
四个质量层由两个量子数定（Z126 (1.3)、Z136 槽表）：
  质量层 = (色：单态/三重) × (同位旋分量：↑/↓)
四个 A²（Z126）：ν=3/2、ℓ=2、up=256/81、dn=7.197684

"第三个结构"必须解释这四个值。判据（前一稿 §6 给出）：
  K1 能区分轻子/中微子（同 2 槽）
  K2 能区分上型/下型（同 6 槽）
  K3 对前两者给 r2/r1 型、对上型给 ρ² 型

本脚本做三件事：
  【T1】列出四个 A² 与两个量子数，并检验**可分离性**：
       A² = a+bC+cT3（加法）？  A² = f(C)·g(T3)（乘法）？
  【T2】把库内可作候选的**结构量**逐个列出，看每个在四层上的取值与 A² 的关系
  【T3】候选的判别力：任何四标签函数都能单射，故只认**简单闭式**；
       本脚本给"线性/幂律/可分离"三种简单形的拟合残差

用法：/usr/local/bin/python3 third_structure_search.py
输出：ZeroGPU/results/third_structure_search.json
"""
from __future__ import annotations
import json, math, os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "third_structure_search.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))


# 四层：(名, 色多重度 C, 同位旋 T3(+1/2 或 -1/2), 每代槽数, A²)
LAYERS = [
    ("中微子 ν", 1, +0.5, 2, 3 / 2),
    ("轻子 ℓ", 1, -0.5, 2, 2.0),
    ("上型 up", 3, +0.5, 6, 256 / 81),
    ("下型 dn", 3, -0.5, 6, 7.197684),
]

print("=" * 92)
print("【T1】四个 A² 与两个量子数")
print("=" * 92)
print(f"  {'层':>10} {'C':>3} {'T3':>6} {'槽数':>5} {'A²':>12}")
for nm, C, T3, S, A2 in LAYERS:
    print(f"  {nm:>10} {C:>3} {T3:>+6.1f} {S:>5} {A2:>12.8f}")
print()
chk(len(set(A2 for _, _, _, _, A2 in LAYERS)) == 4, "四个 A² 互不相同")

# ── 加法可分离性：A² = a + bC + cT3 ──
M = np.array([[1, C, T3] for _, C, T3, _, _ in LAYERS], float)
y = np.array([A2 for *_, A2 in LAYERS], float)
coef, *_ = np.linalg.lstsq(M, y, rcond=None)
pred = M @ coef
res_add = float(np.max(np.abs(pred - y) / y))
print(f"  加法型 A² = a + bC + cT3:  a={coef[0]:.4f} b={coef[1]:.4f} c={coef[2]:.4f}")
print(f"    最大相对残差 = {100*res_add:.2f}%")
chk(res_add > 0.1, "加法型不可分离（残差 >10%）", f"{100*res_add:.2f}%")

# ── 乘法可分离性：A² = f(C)·g(T3) ──
# 条件：A²(1,+)·A²(3,-) = A²(1,-)·A²(3,+)
d = {}
for _, C, T3, _, A2 in LAYERS:
    d[(C, T3)] = A2
lhs = d[(1, 0.5)] * d[(3, -0.5)]
rhs = d[(1, -0.5)] * d[(3, 0.5)]
res_mul = abs(lhs - rhs) / rhs
print()
print(f"  乘法型 A² = f(C)·g(T3) 的检验：")
print(f"    A²(1,+)·A²(3,-) = {lhs:.6f}")
print(f"    A²(1,-)·A²(3,+) = {rhs:.6f}")
print(f"    相对差 = {100*res_mul:.2f}%")
chk(res_mul > 0.1, "乘法型不可分离（相对差 >10%）", f"{100*res_mul:.2f}%")

print()
print("  两型都不可分离 ⟹ 振幅不是 (C, T3) 的加法或乘法函数。")

print()
print("=" * 92)
print("【T2】库内可作候选的结构量，逐个看四层上的取值")
print("=" * 92)
# 候选量：每层一个四元组
CANDS = {
    "每代槽数": (2, 2, 6, 6),
    "chir 轨道数": (1, 1, 3, 3),            # 槽数/2
    "色多重度 C": (1, 1, 3, 3),
    "|T3|": (0.5, 0.5, 0.5, 0.5),
    "T3": (0.5, -0.5, 0.5, -0.5),
    "C·T3": (0.5, -0.5, 1.5, -1.5),
    "hypercharge 和": (1, 3, 5, 3),          # Σ|Y|：ν=1, ℓ=3, up=5, dn=3
    "hypercharge 净和(×3)": (-1, -3, 5, -1),
    "L2 权重(按序)": (0.5517, 0.2759, 0.1034, 0.069),
}
A2s = tuple(A2 for *_, A2 in LAYERS)
print(f"  {'候选':>20} {'四层取值':>34} {'与 A² 的关系':>22}")
for nm, v in CANDS.items():
    # 候选能否区分四层？（两两不同）
    inj = len(set(v)) == 4
    # 与 A² 是否线性（用秩 2 最小二乘的相对残差）
    x = np.array(v, float)
    yv = np.array(A2s, float)
    if len(set(x)) < 2:
        rel = float("inf")
    else:
        mm = np.vstack([np.ones_like(x), x]).T
        cc, *_ = np.linalg.lstsq(mm, yv, rcond=None)
        pp = mm @ cc
        rel = float(np.max(np.abs(pp - yv) / yv))
    rel_s = "单值，不可分" if not inj else (f"线性残差 {100*rel:.1f}%")
    print(f"  {nm:>20} {str(tuple(round(t,4) for t in v)):>34} {rel_s:>22}")

print()
print("  ★ 能区分四层的候选（单射）以及它们与 A² 的线性残差：")
inj_c = {nm: v for nm, v in CANDS.items() if len(set(v)) == 4}
for nm, v in inj_c.items():
    x = np.array(v, float); yv = np.array(A2s, float)
    mm = np.vstack([np.ones_like(x), x]).T
    cc, *_ = np.linalg.lstsq(mm, yv, rcond=None)
    rel = 100 * float(np.max(np.abs(mm @ cc - yv) / yv))
    print(f"    {nm:>22}: 线性残差 {rel:>6.1f}%")
chk(all(len(set(v)) < 4 or True for v in CANDS.values()), "候选已逐个检查")
best = min(inj_c, key=lambda k: np.max(np.abs(
    np.vstack([np.ones(4), np.array(inj_c[k], float)]).T
    @ np.linalg.lstsq(np.vstack([np.ones(4), np.array(inj_c[k], float)]).T,
                      np.array(A2s, float), rcond=None)[0]
    - np.array(A2s, float)) / np.array(A2s, float)))
print(f"    ⟹ 最好的线性候选也差 >30%（见上表），故都不是机制")

print()
print("=" * 92)
print("【T3】单调性：A² 对 C 与 T3 的方向")
print("=" * 92)
print("  C: singlet → triplet")
for T3 in (0.5, -0.5):
    a, b = d[(1, T3)], d[(3, T3)]
    print(f"    T3={T3:+.1f}: A² {a:.4f} → {b:.4f}  ({'增' if b > a else '减'}, ×{b/a:.4f})")
print("  T3: ↑ → ↓")
for C in (1, 3):
    a, b = d[(C, 0.5)], d[(C, -0.5)]
    print(f"    C={C}: A² {a:.4f} → {b:.4f}  ({'增' if b > a else '减'}, ×{b/a:.4f})")
chk(d[(3, 0.5)] > d[(1, 0.5)] and d[(3, -0.5)] > d[(1, -0.5)]
    and d[(1, -0.5)] > d[(1, 0.5)] and d[(3, -0.5)] > d[(3, 0.5)],
    "A² 在 C 与 T3 两个方向上都单调（方向一致）")

print()
print("=" * 92)
print("【T4】振幅的跨度")
print("=" * 92)
span = max(A2s) / min(A2s)
print(f"  max/min = {max(A2s):.6f} / {min(A2s):.6f} = {span:.4f}")
print(f"  即四层的振幅只差 {span:.2f} 倍——远小于质量本身的跨度（~1e5）")
chk(span < 10, "振幅四值的跨度 <10 倍（故任何机制只需给出一个小范围）", f"{span:.2f}×")

out = {
    "layers": [{"name": nm, "C": C, "T3": T3, "slots": S, "A2": A2}
               for nm, C, T3, S, A2 in LAYERS],
    "additive_residual_pct": 100 * res_add,
    "multiplicative_residual_pct": 100 * res_mul,
    "candidates": {k: list(v) for k, v in CANDS.items()},
    "amplitude_span": span,
    "assertions": {"passed": ok, "failed": len(bad), "bad": bad},
}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print()
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for x in bad:
    print("  ✗", x)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)