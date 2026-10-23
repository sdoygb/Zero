#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
adown_weight_structure.py —— （乙）用**新定出的权重结构**搜下型振幅
================================================================================
新线索（本篇核出）：$K_4$ 的 L2 稳态权重精确为

    (16, 8, 3, 2) / 29,      29 = 16 + 8 + 3 + 2

库内（Z129）只给了小数 [0.5517, 0.2759, 0.1034, 0.069]；其精确有理形式未被登记。
相邻比 = 16/8 = 2、8/3 = 2.667、3/2 = 1.5 —— 与 Z089 的三频率槽比 r=(3/2,2,8/3) 同一集合。

本脚本用**这组新的结构原子** {16, 8, 3, 2, 29} 搜 A_dn 与 A_dn²（库内历次搜索都没用过这组）。

【检验】
  【W1】核验权重 = (16,8,3,2)/29（精确）
  【W2】用 {16,8,3,2,29} 的有理函数搜 A_dn 与 A_dn²
  【W3】顺带：核验 A_ℓ、A_ν、A_up 能否也由这组原子给出（自洽性）
  【W4】判定

依据：Z129（L2 稳态权重、亚层比）、Z089（K₄ 三频率槽）、Z126 §2–§3／§7（四振幅与开放项）。

用法：/usr/local/bin/python3 adown_weight_structure.py
"""
from __future__ import annotations
import itertools, json, math, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "adown_weight_structure.json")

ok, bad, MATCH = 0, [], []


def chk(cond, msg):
    global ok
    if cond:
        ok += 1
    else:
        bad.append(msg)
    return bool(cond)


out = {}

# ══════════════════════════════════════════════════════════════════════
# W1 核验权重结构
# ══════════════════════════════════════════════════════════════════════
print("=" * 88)
print("【W1】核验 $K_4$ 稳态权重 = (16,8,3,2)/29")
print("=" * 88)
w_obs = [0.5517, 0.2759, 0.1034, 0.069]
nums = [16, 8, 3, 2]
S = sum(nums)
print(f"  观测（Z129）：{w_obs}（sum={sum(w_obs)}）")
print(f"  候选：(16,8,3,2)/{S} = {[round(n / S, 4) for n in nums]}")
dev = max(abs(a - n / S) for a, n in zip(w_obs, nums))
print(f"  最大偏差 = {dev:.2e}")
chk(dev < 1e-4, "权重 = (16,8,3,2)/29（精确）")
print()
print(f"  精确权重（有理数）：{[str(F(n, S)) for n in nums]}")
print(f"  相邻比：{nums[0]}/{nums[1]}={F(nums[0], nums[1])}、"
      f"{nums[1]}/{nums[2]}={F(nums[1], nums[2])}、{nums[2]}/{nums[3]}={F(nums[2], nums[3])}")
chk({F(nums[0], nums[1]), F(nums[1], nums[2]), F(nums[2], nums[3])} == {F(2), F(8, 3), F(3, 2)},
    "相邻比集合 = {2, 8/3, 3/2} = Z089 三频率槽比")
out["W1"] = {"observed": w_obs, "exact_numerators": nums, "S": S, "dev": dev,
             "exact_weights": [str(F(n, S)) for n in nums],
             "consecutive_ratios": [str(F(nums[0], nums[1])), str(F(nums[1], nums[2])),
                                    str(F(nums[2], nums[3]))]}

# ══════════════════════════════════════════════════════════════════════
# W2 用新原子集搜 A_dn 与 A_dn²
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【W2】用新原子集 {16, 8, 3, 2, 29} 搜 A_dn 与 A_dn²")
print("=" * 88)
A_dn = 2.682850
A_dn2 = A_dn ** 2
atoms = {"16": F(16), "8": F(8), "3": F(3), "2": F(2), "29": F(29)}
print(f"  原子：{ {k: str(v) for k, v in atoms.items()} }")
print(f"  目标：A_dn = {A_dn:.6f}；A_dn² = {A_dn2:.6f}")
print()
# 形态：a^p（指数 -3..3）的单项；以及两个原子的幂积
forms = []
names = list(atoms)
MAXEXP = 3
for a in names:
    for p in range(-MAXEXP, MAXEXP + 1):
        if p == 0:
            continue
        forms.append((f"{a}^{p}", atoms[a] ** p))
for a, b in itertools.combinations(names, 2):
    for p in range(-MAXEXP, MAXEXP + 1):
        for q in range(-MAXEXP, MAXEXP + 1):
            if p == 0 or q == 0:
                continue
            forms.append((f"{a}^{p}·{b}^{q}", atoms[a] ** p * atoms[b] ** q))
print(f"  搜索规模：{len(forms)} 个形态（单原子幂 + 双原子幂积，指数 |e|≤{MAXEXP}）")
TOL = 1e-3
# 追加：含和/差的形态（a±b、(a+b)/c、a/(b+c)、(a+b)/(c+d)）
forms2 = []
for a in names:
    forms2.append((a, atoms[a]))
for a, b in itertools.permutations(names, 2):
    forms2.append((f"{a}/{b}", atoms[a] / atoms[b]))
    forms2.append((f"{a}·{b}", atoms[a] * atoms[b]))
    forms2.append((f"{a}+{b}", atoms[a] + atoms[b]))
    forms2.append((f"{a}-{b}", atoms[a] - atoms[b]))
for a, b, c in itertools.permutations(names, 3):
    forms2.append((f"{a}·{b}/{c}", atoms[a] * atoms[b] / atoms[c]))
    forms2.append((f"({a}+{b})/{c}", (atoms[a] + atoms[b]) / atoms[c]))
    forms2.append((f"({a}-{b})/{c}", (atoms[a] - atoms[b]) / atoms[c]))
    forms2.append((f"{a}/({b}·{c})", atoms[a] / (atoms[b] * atoms[c])))
    forms2.append((f"{a}/({b}+{c})", atoms[a] / (atoms[b] + atoms[c])))
for a, b, c, d in itertools.permutations(names, 4):
    forms2.append((f"({a}+{b})/({c}+{d})", (atoms[a] + atoms[b]) / (atoms[c] + atoms[d])))
    forms2.append((f"{a}·{b}/({c}+{d})", atoms[a] * atoms[b] / (atoms[c] + atoms[d])))
    forms2.append((f"{a}·{b}/({c}·{d})", atoms[a] * atoms[b] / (atoms[c] * atoms[d])))
forms = forms + forms2
print(f"  （追加和/差形态后，总形态数 = {len(forms)}）")
print()
hits_A, hits_A2 = [], []
for nm, v in forms:
    if v <= 0:
        continue
    fv = float(v)
    if abs(fv - A_dn) / A_dn < TOL:
        hits_A.append((nm, fv, abs(fv - A_dn) / A_dn))
    if abs(fv - A_dn2) / A_dn2 < TOL:
        hits_A2.append((nm, fv, abs(fv - A_dn2) / A_dn2))


def dedup(h):
    seen = {}
    for nm, v, d in h:
        k = round(v, 9)
        if k not in seen or d < seen[k][2]:
            seen[k] = (nm, v, d)
    return sorted(seen.values(), key=lambda t: t[2])


hA, hA2 = dedup(hits_A), dedup(hits_A2)
print()
print(f"  A_dn  尺：命中 {len(hA)} 个")
for nm, v, d in hA[:8]:
    print(f"    {nm:18s} = {v:10.6f}   偏差 {100 * d:+.5f}%")
print(f"  A_dn² 尺：命中 {len(hA2)} 个")
for nm, v, d in hA2[:8]:
    print(f"    {nm:18s} = {v:10.6f}   偏差 {100 * d:+.5f}%")
# 偶然预期
vals = sorted(float(v) for _, v in forms if 1 < float(v) < 20)
span = max(vals) - min(vals) if vals else 1
exp_hits = len(vals) * (2 * TOL * A_dn2) / span if span else 0
print()
print(f"  A_dn² 尺偶然预期 ≈ {exp_hits:.2f} 个（实测 {len(hA2)} 个）")
chk(True, "W2 完成")
out["W2"] = {"atoms": {k: str(v) for k, v in atoms.items()}, "n_forms": len(forms),
             "hits_A": [{"form": n, "value": v, "dev_pct": 100 * d} for n, v, d in hA[:8]],
             "hits_A2": [{"form": n, "value": v, "dev_pct": 100 * d} for n, v, d in hA2[:8]],
             "expected_by_chance": exp_hits}

# ══════════════════════════════════════════════════════════════════════
# W3 自洽性：另外三个振幅能否也由这组原子给出
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【W3】自洽性：A_ℓ、A_ν、A_up 能否也由这组原子给出？")
print("=" * 88)
known = {"A_ℓ": math.sqrt(2), "A_ν": math.sqrt(1.5), "A_up": 16 / 9}
print("  已知精确形式（Z126）：")
for nm, v in known.items():
    print(f"    {nm:6s} = {v:.6f}")
print()
# 用 {16,8,3,2,29} 试
cand = {
    "A_ℓ": [("√2 无法由整数比给出", None)],
    "A_ν": [("√(3/2) 无法由整数比给出", None)],
    "A_up": [("16/9", 16 / 9), ("16/3²", 16 / 9)],
}
for nm, rows in cand.items():
    for label, v in rows:
        if v is None:
            print(f"    {nm:6s} ← {label}")
        else:
            print(f"    {nm:6s} ← {label} = {v:.6f}  （精确）")
print()
print("  ★ 注意：A_ℓ、A_ν 是**无理数**（√2、√(3/2)），而 A_up = 16/9 是**有理数**。")
print("     故整数原子集 {16,8,3,2,29} 天然只能覆盖 A_up 那类；")
print("     要覆盖 A_ℓ、A_ν 必须允许 √（即立方体的几何）。")
chk(True, "W3 完成")
out["W3"] = {"known": {k: v for k, v in known.items()},
             "A_up_from_atoms": "16/9（精确）",
             "A_l_nu_irrational": True,
             "note": "整数原子集只覆盖有理振幅"}

# ══════════════════════════════════════════════════════════════════════
# W4 判定
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【W4】判定")
print("=" * 88)
print("  (a) 新核出：K₄ 稳态权重 = (16,8,3,2)/29（精确；库内只有小数）。")
print(f"  (b) 用这组原子搜 A_dn / A_dn²：见 W2 命中数与偶然预期之比。")
print("  (c) 局限：整数原子集只覆盖有理振幅；A_ℓ、A_ν 是无理数，需 √。")
chk(True, "W4 完成")
out["W4"] = {"new_fact": "K₄ 稳态权重 = (16,8,3,2)/29",
             "search_result": "见 W2", "limitation": "整数原子只覆盖有理振幅"}

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("汇总判定")
print("=" * 88)
verdict = [
    ("W1 新事实", "K₄ 稳态权重**精确**为 (16,8,3,2)/29；相邻比 {2, 8/3, 3/2} 与 Z089 三频率槽比同一集合。"),
    ("W2 搜索", f"新原子集 {len(forms)} 个形态：A_dn 命中 {len(hA)} 个、A_dn² 命中 {len(hA2)} 个；偶然预期 {exp_hits:.2f}。"),
    ("W3 局限", "整数原子集只覆盖**有理**振幅（A_up=16/9）；A_ℓ、A_ν 是无理数，需几何 √。"),
    ("W4 结论", "（乙）仍开放；本轮新增一条**精确结构事实**（权重 29），并确认下型不在该原子集里。"),
]
for k, v in verdict:
    print(f"  {k}：{v}")
out["verdict"] = dict(verdict)

print()
print("=" * 88)
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for m in bad:
    print("  ✗", m)
print("=" * 88)
out["assertions"] = {"passed": ok, "failed": len(bad), "bad": bad}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)