#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
adown_phase.py —— （甲）θ 侧：四支的相位偏移有没有结构
================================================================================
背景：`四分支重跑前三重判据_结论翻转.md` —— 振幅侧四支都无结构命中。
本篇查 θ 侧：零点 θ* = π - arccos(1/A) 由 A 定；最轻代相位 θ 是拟合量；
偏移 δ = θ - θ* 因此是**被定的**，不是自由参数。

轻子侧的参照（Z123）：δ0 = π/12 - 2/9 rad = 2.2676°，且被称"完全内生"。

问：（1）四支的 δ 是多少？
    （2）有没有一支的 δ 等于 δ0（若 Z123 的机制普遍适用）？
    （3）δ 或 θ* 是否落在任何结构角上（带偶然期望）？

【检验】
  【P1】四支的 δ（度与弧度）
  【P2】δ 对 δ0 的偏离
  【P3】δ 是否在任何"结构角"上（小族 + 偶然期望）
  【P4】θ*/2π 与 θ*/π 的有理近似
  【P5】符号：轻子是负偏移，哪些支是负的

用法：/usr/local/bin/python3 adown_phase.py
输出：ZeroGPU/results/adown_phase.json
"""
from __future__ import annotations
import json, math, os
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "adown_phase.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))


# 四支（A, θ 度）——来自 adown_uncertainty.py 的分支枚举
BRANCHES = [
    ("1", 2.420560, 111.4853),
    ("2", 2.682850, 114.7387),
    ("3", 1.545499, 126.3126),
    ("4", 1.711737, 129.5455),
]

# 轻子的内生偏移（Z123）：δ0 = π/12 - 2/9 弧度
DELTA0_RAD = math.pi / 12 - 2 / 9
DELTA0_DEG = math.degrees(DELTA0_RAD)

print("=" * 92)
print("【P1】四支的零点与偏移")
print("=" * 92)
print(f"  {'#':>2} {'A':>10} {'θ (最轻代)':>12} {'θ*=180°-arccos(1/A)':>21} "
      f"{'δ (度)':>10} {'δ (弧)':>10}")
rows = []
for tag, A, th_deg in BRANCHES:
    tp = 180.0 - math.degrees(math.acos(1.0 / A))
    d_deg = th_deg - tp
    d_rad = math.radians(d_deg)
    rows.append({"tag": tag, "A": A, "theta_deg": th_deg, "thetastar_deg": tp,
                 "delta_deg": d_deg, "delta_rad": d_rad})
    print(f"  {tag:>2} {A:>10.6f} {th_deg:>11.4f}° {tp:>20.4f}° "
          f"{d_deg:>+9.4f}° {d_rad:>+10.6f}")
chk(all(abs(r["thetastar_deg"]) < 180 for r in rows), "四支的零点都在 (0,180°] 内")

print()
print("=" * 92)
print("【P2】δ 对 δ0（轻子内生偏移）的偏离")
print("=" * 92)
print(f"  δ0 = π/12 - 2/9 = {DELTA0_DEG:.4f}° = {DELTA0_RAD:.6f} rad")
print()
print(f"  {'#':>2} {'δ (度)':>10} {'|δ| - δ0':>11} {'同号？':>7} {'相对 δ0 的偏离':>15}")
for r in rows:
    same = (r["delta_deg"] > 0) == (DELTA0_DEG > 0)
    rel = 100 * (abs(r["delta_deg"]) - abs(DELTA0_DEG)) / abs(DELTA0_DEG)
    print(f"  {r['tag']:>2} {r['delta_deg']:>+9.4f}° {abs(r['delta_deg']) - DELTA0_DEG:>+10.4f}° "
          f"{'是' if same else '否':>7} {rel:>+14.1f}%")
best = min(rows, key=lambda r: abs(abs(r["delta_deg"]) - DELTA0_DEG))
chk(abs(abs(best["delta_deg"]) - DELTA0_DEG) / DELTA0_DEG > 0.1,
    "没有一支的 |δ| 接近 δ0（最接近也差 >10%）",
    f"支{best['tag']}: {abs(best['delta_deg']):.4f}° vs {DELTA0_DEG:.4f}° "
    f"({100*(abs(best['delta_deg'])-DELTA0_DEG)/DELTA0_DEG:+.1f}%)")

print()
print("=" * 92)
print("【P3】δ 是否落在结构角上")
print("=" * 92)
# 小族：框架里出现过的角（弧度表示，单位 2π 的分数或 π 的分数）
cands = {
    "π/12 - 2/9 (轻子 δ0)": DELTA0_RAD,
    "π/45 (=4°)": math.pi / 45,
    "π/60 (=3°)": math.pi / 60,
    "π/90 (=2°)": math.pi / 90,
    "2/9 (rad)": 2 / 9,
    "π/12 (=15°)": math.pi / 12,
    "π/180 (=1°)": math.pi / 180,
}
print(f"  {'候选角':>24} {'值(度)':>10} " + " ".join(f"{'支' + r['tag']:>11}" for r in rows))
for nm, v in cands.items():
    cells = []
    for r in rows:
        dev = 100 * (abs(r["delta_rad"]) - v) / v
        cells.append(f"{dev:>+10.2f}%")
    print(f"  {nm:>24} {math.degrees(v):>9.4f}° " + " ".join(cells))
print()
w = 0.001        # 容差 0.1%
exp = len(cands) * len(rows) * 2 * w
print(f"  偶然期望（{len(cands)} 候选 × {len(rows)} 支 × 窗口 ±0.1%） = {exp:.2f} 个")
hits = []
for r in rows:
    for nm, v in cands.items():
        if abs(abs(r["delta_rad"]) - v) / v < w:
            hits.append((r["tag"], nm, v))
print(f"  窗口内命中：{len(hits)} 个  {[(h[0], h[1]) for h in hits]}")
chk(len(hits) <= 1, "±0.1% 窗口内的命中不超过 1 个（与偶然同量级）",
    f"{len(hits)} 个，偶然 {exp:.2f}")

print()
print("=" * 92)
print("【P4】θ* 作为 2π 与 π 的分数（有理近似，分母 ≤ 24）")
print("=" * 92)
print(f"  {'#':>2} {'θ* (度)':>10} {'θ*/2π':>10} {'最简分数':>12} {'相对偏差':>10}")
for r in rows:
    x = r["thetastar_deg"] / 360.0
    fr = Fraction(x).limit_denominator(24)
    print(f"  {r['tag']:>2} {r['thetastar_deg']:>9.4f}° {x:>10.6f} "
          f"{fr.numerator}/{fr.denominator:>3} {100*(x-float(fr))/float(fr):>+9.3f}%")

print()
print("=" * 92)
print("【P5】符号与轻子的一致性")
print("=" * 92)
neg = [r["tag"] for r in rows if r["delta_deg"] < 0]
pos = [r["tag"] for r in rows if r["delta_deg"] > 0]
print(f"  轻子：δ0 < 0（负偏移）")
print(f"  四支：负偏移 = {neg}；正偏移 = {pos}")
chk(len(neg) == 2 and len(pos) == 2, "正负各 2 支（符号不唯一，故符号也不能选支）",
    f"负 {neg}／正 {pos}")

out = {
    "delta0_deg": DELTA0_DEG,
    "branches": rows,
    "candidate_angles": {k: v for k, v in cands.items()},
    "hits": [{"branch": h[0], "candidate": h[1]} for h in hits],
    "chance_expectation": exp,
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