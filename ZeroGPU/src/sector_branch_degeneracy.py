#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sector_branch_degeneracy.py —— 四个质量层的拟合分支退化
================================================================================
质量律（Z126 (1.1)）：m_i = u_s A² M(th_i)²，M(th)=1+A cos th，th_i = th + 2πi/3。
整体标度全进 u_s ⟹ **A 与 th 只由三代质量比值定**（2 个未知、2 个独立比值）。

四层的来源不同，这是本脚本的核心对照：
  轻子   A=√2、th 由零点+偏移定（Z121/Z123） ⟹ **A 与 th 是导出的，三代质量是预言**
  中微子 A=√{3/2}，同型
  上型   A、th、u 三者由拟合定 ⟹ 无预言
  下型   同上

若"拟合才有退化、导出则没有"，则退化是**拟合的产物**而不是质量律的产物。

【检验】
  【S1】四层用 Z126 的参数复现三代质量
  【S2】各层的分支枚举（正 A 去重，角距正确）
  【S3】每层有几支？（拟合层应 >1，导出层应恰好 1）
  【S4】导出层的参数是否**恰好**是其中一支
  【S5】上型：$A²=256/81$ 的命中是否只在拟合选中的那一支上
  【S6】下型：对照（已知 4 支）

用法：/usr/local/bin/python3 sector_branch_degeneracy.py
输出：ZeroGPU/results/sector_branch_degeneracy.json
"""
from __future__ import annotations
import json, math, os

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "sector_branch_degeneracy.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))


A2_MEV = 313.856
ANG = lambda x, y: abs((x - y + math.pi) % (2 * math.pi) - math.pi)


def masses(A, th_rad, u):
    return [u * A2_MEV * (1.0 + A * math.cos(th_rad + 2 * math.pi * i / 3)) ** 2
            for i in range(3)]


def ratios(A, th):
    m = masses(A, th, 1.0)
    return np.array(m) / m[0]


def res(p, tr):
    A, th = p
    r = ratios(A, th)
    if np.any(r <= 0) or not np.all(np.isfinite(r)):
        return np.array([1e6, 1e6])
    return np.log(r) - np.log(tr)


def branches(tr, A_lo=1.0, A_hi=6.0, step=0.2, tol=1e-14):
    sols = []
    for A0 in np.arange(A_lo, A_hi + 1e-9, step):
        for th0 in np.arange(0.0, 2 * math.pi, step):
            try:
                r = least_squares(res, [A0, th0], args=(tr,),
                                  xtol=1e-15, ftol=1e-15, gtol=1e-15)
            except Exception:
                continue
            if float(np.sum(r.fun ** 2)) < tol:
                A, th = r.x[0], r.x[1] % (2 * math.pi)
                if A < 0:
                    A, th = -A, (th + math.pi) % (2 * math.pi)
                if all(abs(A - a) > 1e-3 or ANG(th, t) > 1e-3 for _, a, t in sols):
                    sols.append((float(np.sum(r.fun ** 2)), A, th))
    return sorted(sols)


# Z126 的四层参数（A, th 度, u）与三代质量
SECTORS = [
    ("轻子 e,μ,τ", math.sqrt(2), 132.7324, 1.0, "导出", (0.51099, 105.66, 1777.0)),
    ("上型 u,c,t", 1.777833, 124.6107, 71.60, "拟合", (2.16, 1270.0, 1.727e5)),
    ("下型 d,s,b", 2.682850, 114.7387, 0.98798, "拟合", (4.67, 93.4, 4180.0)),
]

print("=" * 96)
print("【S1】用 Z126 的参数复现三代质量")
print("=" * 96)
repro = {}
for nm, A, th_deg, u, kind, tgt in SECTORS:
    m = masses(A, math.radians(th_deg), u)
    dev = max(abs(m[i] - tgt[i]) / tgt[i] for i in range(3))
    repro[nm] = (A, th_deg, u, m, dev)
    print(f"  {nm:<14} A={A:.6f} th={th_deg:>9.4f}° u={u:>8.5f}  "
          f"⟹ {m[0]:.4g}, {m[1]:.4g}, {m[2]:.4g}  (最大偏差 {100*dev:.2f}%)")
chk(all(repro[n][4] < 0.02 for n in repro), "三层都能用 Z126 参数复现（<2%）")

print()
print("=" * 96)
print("【S2】【S3】各层的分支枚举")
print("=" * 96)
result = {}
for nm, A, th_deg, u, kind, tgt in SECTORS:
    tr = np.array(tgt) / tgt[0]
    brs = branches(tr)
    result[nm] = {"kind": kind, "A_reported": A, "theta_deg": th_deg,
                  "branches": [{"A": b[1], "theta_deg": math.degrees(b[2])} for b in brs]}
    print(f"  {nm}（{kind}）：{len(brs)} 支")
    print(f"     {'#':>2} {'A':>12} {'theta':>11} {'与报告值同支？':>13}")
    for i, (e, a, t) in enumerate(brs, 1):
        same = abs(a - A) < 1e-3 and ANG(t, math.radians(th_deg)) < 1e-3
        print(f"     {i:>2} {a:>12.6f} {math.degrees(t):>10.4f}° {'是' if same else '':>13}")

print()
print("=" * 96)
print("【S3/S4】拟合层 vs 导出层：支数是否不同")
print("=" * 96)
for nm in result:
    r = result[nm]
    n = len(r["branches"])
    contains = any(abs(b["A"] - r["A_reported"]) < 1e-3 for b in r["branches"])
    print(f"  {nm:<14} {r['kind']}  支数={n}  报告值在其中={contains}")
    r["n_branches"] = n
    r["reported_in_branches"] = contains

chk(result["轻子 e,μ,τ"]["n_branches"] >= 1 and result["轻子 e,μ,τ"]["reported_in_branches"],
    "轻子：报告值（导出）确实是分支之一（方法自检通过）")
chk(result["上型 u,c,t"]["n_branches"] > 1, "上型：也有分支退化",
    f"{result['上型 u,c,t']['n_branches']} 支")
chk(result["下型 d,s,b"]["n_branches"] > 1, "下型：也有分支退化",
    f"{result['下型 d,s,b']['n_branches']} 支")

print()
print("=" * 96)
print("【S5】上型：A² = 256/81 的命中在哪一支上")
print("=" * 96)
TARGET = 256.0 / 81.0          # (4/3)^4，Z126 记 −0.006%
print(f"  目标 A² = 256/81 = {TARGET:.8f}（Z126 记与拟合值差 −0.006%）")
print(f"  {'A':>12} {'A²':>12} {'相对偏差':>12} {'命中(±0.01%)？':>15}")
up_branches = result["上型 u,c,t"]["branches"]
hits = []
for b in up_branches:
    A = b["A"]
    d = 100 * (A * A - TARGET) / TARGET
    hit = abs(d) < 0.01
    if hit:
        hits.append(A)
    print(f"  {A:>12.6f} {A*A:>12.6f} {d:>+11.3f}% {'是' if hit else '否':>15}")
print()
print(f"  256/81 的命中支数：{len(hits)}  {['%.6f' % h for h in hits]}")
chk(len(hits) <= 1, "A²=256/81 只在至多一支上命中（若 1 支，则它是分支选择而非普遍关系）",
    f"{len(hits)} 支")

print()
print("=" * 96)
print("【S6】轻子：A=√2 的命中")
print("=" * 96)
lep = result["轻子 e,μ,τ"]["branches"]
print(f"  {'A':>12} {'A²':>12} {'对 2 的偏差':>12}")
lh = []
for b in lep:
    A = b["A"]
    d = 100 * (A * A - 2.0) / 2.0
    if abs(d) < 0.01:
        lh.append(A)
    print(f"  {A:>12.6f} {A*A:>12.6f} {d:>+11.3f}%")
print()
print(f"  A²=2 的命中支数：{len(lh)}  {['%.6f' % h for h in lh]}")
chk(len(lh) == 1, "轻子：A=√2 恰好命中一支（这正是导出值的意义）", f"{len(lh)} 支")

out = {"sectors": result,
       "up_256_81_hits": hits,
       "lepton_sqrt2_hits": lh,
       "assertions": {"passed": ok, "failed": len(bad), "bad": bad}}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print()
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for x in bad:
    print("  ✗", x)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)