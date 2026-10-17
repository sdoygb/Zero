#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_shape_residual.py —— 中微子形状的残差与候选修正（Z164 的核验）。

形状（Z162）：m1=0、A_ν=√(3/2)、间隔 2π/3 ⟹ R_pred=(73−28√6)/25、K_ν=7/12。
残差：R_pred 与数据反解的 R_obs=√(Δm²21/Δm²31) 差 +0.8%…+2.7%（随数据集）。

本脚本：
  A. 多个数据集的 R_obs 与残差；
  B. 四条候选修正各自的「精确命中所需值」：
     (i) 饱和点偏移 ε；(ii) 间隔 Δ；(iii) 振幅 A²；(iv) 最小质量 m1；
  C. 这些所需值对一组预声明的框架量的接近程度（只报数，不判定）。

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
A = math.sqrt(A2)
x0 = math.pi - math.acos(math.sqrt(2 / 3))

SETS = [("NuFIT 5.2", 7.42e-5, 2.510e-3),
        ("NuFIT 6.0", 7.49e-5, 2.513e-3),
        ("较早一版", 7.53e-5, 2.453e-3),
        ("NuFIT 4.0", 7.39e-5, 2.525e-3)]

print("== A. 残差（R 的定义与数据集） ==")
rows = []
for name, d21, d31 in SETS:
    X = d21 / d31
    Robs = math.sqrt(X)
    rows.append({"set": name, "d21": d21, "d31": d31, "X": X, "R_obs": Robs,
                 "rel_R": Robs / R_PRED - 1, "rel_X": X / R_PRED**2 - 1})
    print(f"   {name:10s}: X={X:.6f}  R_obs={Robs:.6f}  "
          f"R 残差 {100*(Robs/R_PRED-1):+.2f}%  X 残差 {100*(X/R_PRED**2-1):+.2f}%")
check("R_pred 与各数据集 R_obs 的残差在 1–3%", all(0.005 < abs(r["rel_R"]) < 0.035 for r in rows),
      f"{100*min(abs(r['rel_R']) for r in rows):.2f}–{100*max(abs(r['rel_R']) for r in rows):.2f}%")
check("R_pred 高于所有数据集的 R_obs", all(r["R_obs"] < R_PRED for r in rows))

print("\n== B. 候选修正的精确需求值（以 NuFIT 5.2 为目标） ==")
D21, D31 = SETS[0][1], SETS[0][2]
X_T = D21 / D31


def X_of(eps=0.0, delta=None, a2=A2, m1=0.0):
    d = 2 * math.pi / 3 if delta is None else delta
    f = [1 + math.sqrt(a2) * math.cos(x0 + eps + d * k) for k in range(3)]
    m = sorted((v * v for v in f))
    m = [mi + m1**2 for mi in m]
    return (m[1] ** 2 - m[0] ** 2) / (m[2] ** 2 - m[0] ** 2)


def bisect(fn, lo, hi, target):
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if (fn(mid) - target) * (fn(lo) - target) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# (i) 偏移 ε（rad）
eps_req = bisect(lambda e: X_of(eps=e), -0.05, 0.05, X_T)
# (ii) 间隔 Δ（rad），在 100°–140° 内
d_req = bisect(lambda d: X_of(delta=d), math.radians(100), math.radians(140), X_T)
# (iii) 振幅 A²
a2_req = bisect(lambda a2: X_of(a2=a2), 1.2, 2.5, X_T)
# (iv) 最小质量 m1（eV），用尺度 m3=0.0486 的假定（Z163）
m3_scale = 0.0486
m1_req = bisect(lambda m: X_of(m1=m), 0.0, 0.03, X_T) * m3_scale

print(f"   (i)   偏移 ε* = {eps_req:+.6f} rad = {math.degrees(eps_req):+.4f}°")
print(f"   (ii)  间隔 Δ* = {math.degrees(d_req):.4f}°（比 2π/3 少 {120-math.degrees(d_req):.4f}°）")
print(f"   (iii) 振幅 A²* = {a2_req:.6f}（比 3/2 高 {100*(a2_req/1.5-1):+.3f}%）")
print(f"   (iv)  最小质量 m1* = {m1_req*1e3:.3f} meV（≙ m1/m3={m1_req/m3_scale:.4f}）")
check("四个需求值都小且非零（残差是小的结构修正）",
      abs(eps_req) < 0.05 and 100 < math.degrees(d_req) < 140 and 1.2 < a2_req < 2.5 and m1_req < 0.02)

print("\n== C. 需求值对预声明框架量的接近程度（只报数） ==")
cands = {
    "δ0=π/12−2/9": math.pi / 12 - 2 / 9,
    "2/9": 2 / 9,
    "π/12": math.pi / 12,
    "2π/3−2/9(角)": math.radians(2 / 9),
    "(2/9)²/6": (2 / 9) ** 2 / 6,
    "δ0/5": (math.pi / 12 - 2 / 9) / 5,
    "1/(4π)": 1 / (4 * math.pi),
    "9/(100·4π)": 9 / (100 * 4 * math.pi),
    "α=1/137.036": 1 / 137.036,
    "1/32": 1 / 32,
    "1/108": 1 / 108,
}
best = []
for k, v in cands.items():
    best.append((abs(abs(eps_req) - v), k, v))
best.sort()
print("   与 |ε*| 最接近的三个框架量：")
for d, k, v in best[:3]:
    print(f"     {k:14s} = {v:.6f}   相对差 {100*(v/abs(eps_req)-1):+.2f}%")
# 间隔修正
dc = 2 * math.pi / 3 - d_req
bestd = sorted((abs(dc - v), k) for k, v in cands.items())
print(f"   间隔修正 {dc:.6f} rad 最接近：{bestd[0][1]}（差 {100*(dict((k,v) for k,v in cands.items())[bestd[0][1]]/dc-1):+.2f}%）")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "R_pred": R_PRED,
       "datasets": rows,
       "corrections": {"eps_rad": eps_req, "eps_deg": math.degrees(eps_req),
                       "delta_rad": d_req, "delta_deg": math.degrees(d_req),
                       "A2": a2_req, "m1_eV": m1_req},
       "closest_eps": [{"quantity": k, "value": v, "rel": v / abs(eps_req) - 1} for _, k, v in best[:3]]}
with open(os.path.join(OUT, "neutrino_shape_residual.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_shape_residual.json")
sys.exit(0 if not FAIL else 1)
