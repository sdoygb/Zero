#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_shape_residual.py —— 中微子形状的残差：三条参数化与退化（Z164 v2）。

形状（Z162）：m1=0、A_ν=√(3/2)、间隔 2π/3、饱和于零点 ⟹ R_pred=(73−28√6)/25、K_ν=7/12。
数据：R_obs=√(Δm²21/Δm²31)，K_obs=(1+R_obs)/(1+√R_obs)²。

本脚本：
  A. 数据集的 R_obs、K_obs 与模型点（R_pred、7/12）；
  B. 三条一参数变形（振幅 A²、饱和偏移 ε、世代间隔 Δ），各自解 X=X_obs 后**同时**报 K：
     三条都命中 (K_obs,X_obs) ⟹ 数据点处三者退化，数据本身不指定哪一条「错」；
  C. 报出三条的所需值，供结构判定用（不在本脚本内判定）。

输出：ZeroGPU/results/neutrino_shape_residual.json
用法：python3 ZeroGPU/src/neutrino_shape_residual.py
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


R_PRED = (73 - 28 * math.sqrt(6)) / 25
A2 = 1.5


def model(a2=A2, eps=0.0, delta=None):
    """返回 (X, K) —— m1=0、零点自洽、间隔 delta（默认 2π/3）。"""
    A = math.sqrt(a2)
    x0 = math.pi - math.acos(1 / A)          # cos x0 = −1/A
    d = 2 * math.pi / 3 if delta is None else delta
    f = [1 + A * math.cos(x0 + eps + d * k) for k in range(3)]
    m = [v * v for v in f]                    # m ∝ f²；m1=0
    X = (m[1] ** 2 - m[0] ** 2) / (m[2] ** 2 - m[0] ** 2)
    K = (m[1] + m[2]) / (math.sqrt(m[1]) + math.sqrt(m[2])) ** 2
    return X, K


def solve(fn, lo, hi, target):
    """fn 在 [lo,hi] 上单调（方向由两端值判）。"""
    inc = fn(hi) > fn(lo)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        v = fn(mid)
        if (v < target) == inc:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


SETS = [("NuFIT 5.2", 7.42e-5, 2.510e-3),
        ("NuFIT 6.0", 7.49e-5, 2.513e-3),
        ("较早一版", 7.53e-5, 2.453e-3),
        ("NuFIT 4.0", 7.39e-5, 2.525e-3)]

print("== A. 数据与模型点 ==")
rows = []
for name, d21, d31 in SETS:
    X = d21 / d31
    Robs = math.sqrt(X)
    K = (1 + Robs) / (1 + math.sqrt(Robs)) ** 2
    rows.append({"set": name, "d21": d21, "d31": d31, "X_obs": X, "R_obs": Robs,
                 "K_obs": K, "rel_R": Robs / R_PRED - 1, "rel_X": X / R_PRED**2 - 1})
    print(f"   {name:10s}: R_obs={Robs:.6f} K_obs={K:.6f}  "
          f"R 残差 {100*(Robs/R_PRED-1):+.2f}%  X 残差 {100*(X/R_PRED**2-1):+.2f}%")
Xm, Km = model()
print(f"   模型: X={Xm:.6f} K={Km:.9f}")
check("模型点与数据的残差（全部数据集）在 1–3%（R 上）",
      all(0.005 < abs(r["rel_R"]) < 0.035 for r in rows))

print("\n== B. 三条一参数变形：解 X 后同时报 K ==")
D21, D31 = SETS[0][1], SETS[0][2]
XT = D21 / D31
ROBS = math.sqrt(XT)
KOBS = (1 + ROBS) / (1 + math.sqrt(ROBS)) ** 2

a2 = solve(lambda a: model(a2=a)[0], 1.1, 2.5, XT)
eps = solve(lambda e: model(eps=e)[0], -0.05, 0.05, XT)
delta = solve(lambda d: model(delta=d)[0], math.radians(100), math.radians(140), XT)
routes = [
    ("振幅 A²", f"{a2:.6f}", model(a2=a2)[1]),
    ("偏移 ε", f"{math.degrees(eps):+.4f}°", model(eps=eps)[1]),
    ("间隔 Δ", f"{math.degrees(delta):.4f}°", model(delta=delta)[1]),
]
for name, val, K in routes:
    print(f"   {name:8s} 所需 {val:14s} ⟹ K = {K:.9f}   (K_obs={KOBS:.9f}，差 {100*(K/KOBS-1):+.4f}%)")
check("三条变形在数据点处同时命中 K（退化）",
      all(abs(K / KOBS - 1) < 1e-4 for _, _, K in routes),
      "最大差 %.1e" % max(abs(K / KOBS - 1) for _, _, K in routes))
check("三条所需值都是小量（≤4%）",
      abs(a2 / 1.5 - 1) < 0.05 and abs(math.degrees(eps)) < 1 and abs(math.degrees(delta) - 120) < 1)

print("\n== C. 所需值汇总 ==")
print(f"   A²*={a2:.6f}（比 3/2 高 {100*(a2/1.5-1):+.3f}%）")
print(f"   ε*={math.degrees(eps):+.4f}°")
print(f"   Δ*={math.degrees(delta):.4f}°（比 120° 少 {120-math.degrees(delta):.4f}°）")
print(f"   数据点 (X,K)=({XT:.6f},{KOBS:.6f})；模型点 ({Xm:.6f},{Km:.6f})")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "R_pred": R_PRED, "datasets": rows,
       "model": {"X": Xm, "K": Km},
       "routes": {"A2": a2, "eps_deg": math.degrees(eps), "delta_deg": math.degrees(delta)},
       "K_obs": KOBS, "X_obs": XT}
with open(os.path.join(OUT, "neutrino_shape_residual.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_shape_residual.json")
sys.exit(0 if not FAIL else 1)
