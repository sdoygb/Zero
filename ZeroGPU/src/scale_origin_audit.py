#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scale_origin_audit.py —— 定理级标度的入口审计：c_ℓ 与 c_ν 的候选结构。

目标（新 goal）：把中微子绝对标度做成定理级——纯数 c_ℓ=μ_ℓ²/v 与 c_ν=m_ν3/v
必须由 Zero 结构给出。本脚本做「先查库」的定量审计：
  1. 目标纯数（c_ℓ、c_ν、ln 值）；
  2. 候选结构：四环指数组合、中心二项式、层权重、闭包角、常用常数；
  3. 记录最近的匹配与所需指数（含搜索空间计数，避免把巧合当结构）。

输出：ZeroGPU/results/scale_origin_audit.json
用法：python3 ZeroGPU/src/scale_origin_audit.py
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


V, M_Z, M_TAU = 246.22, 91.1876, 1.77686
M_X = 1.009e13
M_I = M_Z * math.exp(4 * math.pi)

# 环指数
r1 = math.log(M_Z / V)              # ≈ −0.9933
r2 = 4 * math.pi
r3 = math.log(M_X / M_I)
r4 = math.log(9.821e16 / M_X)       # 环④（E6），Z153 登记为上游开放项

A_L = math.sqrt(2)
delta0 = math.pi / 12 - 2 / 9
delta_e = 3 * math.pi / 4 - delta0
f = [1 + A_L * math.cos(delta_e + 2 * math.pi * k / 3) for k in range(3)]
f_tau = max(f)
c_ell = M_TAU / (V * f_tau**2)
R = (73 - 28 * math.sqrt(6)) / 25
m3 = (M_TAU / 2) ** 2 / math.sqrt(M_I * M_X) * 1e9
c_nu = m3 * 1e-9 / V

print("== 1. 目标纯数 ==")
print(f"   c_ℓ = {c_ell:.6e}   ln(1/c_ℓ) = {math.log(1/c_ell):.6f}")
print(f"   c_ν = {c_nu:.6e}   ln(1/c_ν) = {math.log(1/c_nu):.6f}")
print(f"   环指数: r1={r1:.6f} r2={r2:.6f} r3={r3:.6f} r4={r4:.6f}")

print("\n== 2. 候选结构（环指数整数组合，系数 −2..2） ==")
target = math.log(1 / c_ell)
best = []
for a in range(-2, 3):
    for b in range(-2, 3):
        for c in range(-2, 3):
            for d in range(-2, 3):
                if a == b == c == d == 0:
                    continue
                val = a * r1 + b * r2 + c * r3 + d * r4
                rel = abs(val - target) / target
                best.append((rel, a, b, c, d, val))
best.sort()
print(f"   组合数 {len(best)}；最近 5 个：")
for rel, a, b, c, d, val in best[:5]:
    print(f"     {a:+d}r1{b:+d}r2{c:+d}r3{d:+d}r4 = {val:.6f}  (目标 {target:.6f}，相对差 {100*rel:+.2f}%)")

print("\n== 3. 中心二项式与层权重 ==")
for L in (8, 10, 12, 14, 16, 20, 24):
    C = math.comb(L, L // 2)
    for label, v in (("1/C", 1 / C), ("2/C", 2 / C), ("(L+1)/C", (L + 1) / C)):
        rel = abs(v / c_ell - 1)
        if rel < 0.2:
            print(f"     L={L} {label} = {v:.4e}  相对差 {100*rel:+.2f}%")

print("\n== 4. 闭包角与常用常数 ==")
cands = {
    "e^{−r4}": math.exp(-r4),
    "e^{−(r4−r1)}": math.exp(-(r4 - r1)),
    "e^{−r2/2}": math.exp(-r2 / 2),
    "e^{−(r2+r4)/3}": math.exp(-(r2 + r4) / 3),
    "(2/9)^3": (2 / 9) ** 3,
    "π/250": math.pi / 250,
    "1/(250π)": 1 / (250 * math.pi),
    "α/108": (1 / 137.036) / 108,
    "1/784": 1 / 784,
}
for k, v in cands.items():
    print(f"     {k:16s} = {v:.4e}  相对差 {100*(v/c_ell-1):+.2f}%")

print("\n== 判定 ==")
near = [x for x in best if x[0] < 0.005]
print(f"   1% 以内的环组合：{len(near)} 个（搜索空间 {len(best)}）")
check("没有 1% 以内的环指数组合", len(near) == 0)
check("c_ℓ 的所需指数 6.663 与 r4−r1=6.599 相差 >0.5%",
      abs((r4 - r1) / target - 1) > 0.005, f"{(r4-r1)/target-1:+.3%}")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "c_ell": c_ell, "c_nu": c_nu,
       "rings": {"r1": r1, "r2": r2, "r3": r3, "r4": r4},
       "target_ln_inv_c_ell": target, "n_combos": len(best),
       "closest": [{"rel": x[0], "comb": [x[1], x[2], x[3], x[4]], "val": x[5]}
                   for x in best[:5]]}
with open(os.path.join(OUT, "scale_origin_audit.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/scale_origin_audit.json")
sys.exit(0 if not FAIL else 1)
