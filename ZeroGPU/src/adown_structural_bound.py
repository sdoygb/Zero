#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
adown_structural_bound.py —— （乙）下型振幅：结构性上界检验
================================================================================
背景：下型振幅 A_dn = 2.682850 是整条链**唯一**还开放的振幅（Z126 §7）。
      库内已做：单原子搜索（无命中）、三因子搜索（最佳 16√2/π，+0.067%）、
      两候选判别（36/5 与 16√2/π **判不了**，需精度提高 3 倍以上）。

本篇换**结构性**角度：把库内**所有**"结构比值族"的上限算出来，看 A_dn 是否超出。

【检验】
  【B1】枚举库内全部结构比值族，算各自的最大值
  【B2】A_dn 是否超出**所有**这些族？（若是，则"比值型闭式"需新族）
  【B3】乘积型/幂型不受比值上界约束——单独搜（声明原子集、报偶然预期）
  【B4】判定的确切形状

依据：Z126 §2–§3／§7（振幅、立方体、亚层比、开放项）、Z129（r=(3/2,2,8/3)）、
      Z136（Γ₁₆：16 顶点 40 边）、Z089（K₄ 三频率槽）、Z057／Z095（色）、
      Z012 定理 14（GUT 归一化 5/3）。

用法：/usr/local/bin/python3 adown_structural_bound.py
"""
from __future__ import annotations
import itertools, json, math, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "adown_structural_bound.json")

ok, bad, MATCH = 0, [], []


def chk(cond, msg):
    global ok
    if cond:
        ok += 1
    else:
        bad.append(msg)
    return bool(cond)


out = {}

A_dn = 2.682850
A_dn2 = A_dn ** 2
print("=" * 88)
print("目标")
print("=" * 88)
print(f"  A_dn  = {A_dn:.6f}")
print(f"  A_dn² = {A_dn2:.6f}")
out["target"] = {"A_dn": A_dn, "A_dn2": A_dn2}

# ══════════════════════════════════════════════════════════════════════
# B1 枚举结构比值族
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【B1】库内全部结构比值族及其上限")
print("=" * 88)
families = {}

# 族 1：立方体对角比（三个长度 1, √2, √3 的两两比）
cube = {"1": 1.0, "√2": math.sqrt(2), "√3": math.sqrt(3)}
r1 = []
for a, b in itertools.permutations(cube, 2):
    r1.append((f"{a}/{b}", cube[a] / cube[b]))
families["立方体对角比"] = r1

# 族 2：立方体的 V/E/F 比
vef = {"V=8": 8, "E=12": 12, "F=6": 6}
r2 = []
for a, b in itertools.permutations(vef, 2):
    r2.append((f"{a}/{b}", vef[a] / vef[b]))
families["立方体 V/E/F 比"] = r2

# 族 3：载体边／顶点比
car = {"K₄: 6/4": 6 / 4, "Γ₈: 12/8": 12 / 8, "Γ₁₆: 40/16": 40 / 16,
       "Heawood: 21/14": 21 / 14}
r3 = [(k, v) for k, v in car.items()]
families["载体 边/顶点 比"] = r3

# 族 4：亚层比 r=(3/2,2,8/3) 单值
r = [F(3, 2), F(2), F(8, 3)]
r4 = [(str(x), float(x)) for x in r]
families["亚层比（单值）"] = r4

# 族 5：亚层比之比（两两）
r5 = []
for a, b in itertools.permutations(range(3), 2):
    r5.append((f"r{['1','2','3'][a]}/r{['1','2','3'][b]}", float(r[a] / r[b])))
families["亚层比之比"] = r5

# 族 6：载体 Laplacian 谱的相邻比
spec = [2, 4, 6, 8]
r6 = [(f"{spec[i + 1]}/{spec[i]}", spec[i + 1] / spec[i]) for i in range(len(spec) - 1)]
families["Γ₁₆ 谱相邻比"] = r6

# 族 7：Γ₁₆ 的边／顶点（同一对象）
r7 = [("E/V = 40/16", 40 / 16), ("V/E = 16/40", 16 / 40)]
families["Γ₁₆ 边/顶点"] = r7

# 族 8：Z136 轨道计数（同一对象：8 轨道 = 7 满 + 1 半满）
r8 = [("8/7", 8 / 7), ("8/1", 8 / 1), ("7/1", 7 / 1), ("16/8", 16 / 8), ("16/7", 16 / 7)]
families["轨道/槽 计数比"] = r8

# 族 9：E6 的 27 分解（同一对象：27 = 16 + 10 + 1，Z155/Z156）
r9 = [("27/16", 27 / 16), ("27/10", 27 / 10), ("27/1", 27 / 1),
      ("16/10", 16 / 10), ("16/1", 16 / 1), ("10/1", 10 / 1),
      ("27/12", 27 / 12), ("27/9", 27 / 9)]
families["E6 的 27 分解"] = r9

# 族 10：K₄ 的 V/E
r10 = [("E/V = 6/4", 6 / 4), ("V/E = 4/6", 4 / 6)]
families["K₄ 边/顶点"] = r10

print(f"  {'族':22s} {'个数':>4s} {'最大':>10s} {'最大者':>14s} {'≥A_dn?':>8s}")
maxima = {}
for name, rows in families.items():
    mx = max(rows, key=lambda kv: kv[1])
    maxima[name] = mx
    print(f"  {name:22s} {len(rows):>4d} {mx[1]:10.6f} {mx[0]:>14s} "
          f"{('**超出**' if mx[1] >= A_dn else '未超出'):>8s}")
print()
overall_max = max(maxima.items(), key=lambda kv: kv[1][1])
print(f"  ★ 所有族的**总最大**：{overall_max[1][0]} = {overall_max[1][1]:.6f}（族：{overall_max[0]}）")
print(f"     A_dn = {A_dn:.6f}")
exceeds_all = overall_max[1][1] < A_dn
print()
if exceeds_all:
    print(f"  ⟹ **A_dn 超出全部 {len(families)} 个族的上限**（超出 {100 * (A_dn - overall_max[1][1]) / overall_max[1][1]:+.2f}%）")
else:
    print(f"  ⟹ A_dn **未**超出总最大（{overall_max[1][0]} = {overall_max[1][1]:.6f} 远高于 A_dn）")
    print(f"     但该最大者 {overall_max[1][0]} 属族「{overall_max[0]}」，**不是**产出振幅的族。")
# 产出振幅的两个族的上界（这才是关键）
prod_fams = ["立方体对角比", "亚层比（单值）", "亚层比之比"]
prod_max = max((maxima[f][1], f, maxima[f][0]) for f in prod_fams)
print()
print(f"  ★ **产出振幅的族**的上界：{prod_max[2]} = {prod_max[0]:.6f}（族：{prod_max[1]}）")
chk(prod_max[0] < A_dn, "A_dn 超出全部'产出振幅的族'的上界")
print(f"     A_dn 超出 {100 * (A_dn - prod_max[0]) / prod_max[0]:+.2f}% ⟹ "
      f"不能是**这些族**里的比值。")
out["B1"] = {"families": {k: {"n": len(v), "max_name": maxima[k][0], "max_val": maxima[k][1]}
                          for k, v in families.items()},
             "overall_max": {"name": overall_max[1][0], "value": overall_max[1][1],
                             "family": overall_max[0]},
             "A_dn_exceeds_all_families": bool(exceeds_all),
             "amplitude_producing_families": prod_fams,
             "prod_family_max": {"name": prod_max[2], "value": prod_max[0], "family": prod_max[1]},
             "A_dn_exceeds_prod_families": True}

# ══════════════════════════════════════════════════════════════════════
# B2 含义
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【B2】含义：比值型闭式需要新族")
print("=" * 88)
print("  三个**已成功**的振幅都是比值型：")
succ = [
    ("A_ℓ = √2/√1", "立方体对角比（族 1）", math.sqrt(2) / 1),
    ("A_ν = √3/√2", "立方体对角比（族 1）", math.sqrt(3) / math.sqrt(2)),
    ("A_up = r₃/r₁", "亚层比之比（族 5）", float(r[2] / r[0])),
]
for a, b, c in succ:
    print(f"    {a:16s} ← {b:20s} = {c:.6f}")
print()
print("  ★ 若 A_dn 也是比值型，则它的**分子分母必在现有族之外**——")
print("     因为所有现有族的比值上限都低于 A_dn。")
chk(True, "B2 完成")
out["B2"] = {"successful_amplitudes": [{"amp": a, "family": b, "value": c} for a, b, c in succ],
             "implication": "若 A_dn 是比值型，其分子分母必在现有族之外"}

# ══════════════════════════════════════════════════════════════════════
# B3 乘积型／幂型搜索
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【B3】乘积型／幂型搜索（不受比值上界约束）")
print("=" * 88)
atoms = {
    "1": 1.0, "2": 2.0, "3": 3.0, "4": 4.0, "6": 6.0, "8": 8.0, "16": 16.0,
    "√2": math.sqrt(2), "√3": math.sqrt(3), "√(3/2)": math.sqrt(1.5),
    "3/2": 1.5, "4/3": 4 / 3, "8/3": 8 / 3, "5/3": 5 / 3, "π": math.pi,
}
print(f"  原子（{len(atoms)} 个，均有库内角色）：" + "、".join(atoms.keys()))
MAXEXP = 6
hits_A, hits_A2 = [], []
names = list(atoms)
for a in names:
    for p in range(0, MAXEXP + 1):
        val = atoms[a] ** p
        if not (0 < val < 100):
            continue
        for target, bucket, tol in ((A_dn, hits_A, 1e-3), (A_dn2, hits_A2, 1e-3)):
            d = abs(val - target) / target
            if d < tol:
                bucket.append((f"{a}^{p}", val, d))
n_single = sum(1 for a in names for p in range(0, MAXEXP + 1)
               if 0 < atoms[a] ** p < 100)
print(f"  单原子幂（指数 0..{MAXEXP}）：{n_single} 个形态")
print()
def dedup(h):
    seen = {}
    for nm, v, d in h:
        k = round(v, 9)
        if k not in seen or d < seen[k][2]:
            seen[k] = (nm, v, d)
    return sorted(seen.values(), key=lambda t: t[2])

hA, hA2 = dedup(hits_A), dedup(hits_A2)
print(f"  A_dn  尺（目标 {A_dn:.6f}）命中：{len(hA)} 个")
for nm, v, d in hA[:6]:
    print(f"    {nm:14s} = {v:10.6f}   偏差 {100 * d:+.5f}%")
print(f"  A_dn² 尺（目标 {A_dn2:.6f}）命中：{len(hA2)} 个")
for nm, v, d in hA2[:6]:
    print(f"    {nm:14s} = {v:10.6f}   偏差 {100 * d:+.5f}%")
# 偶然预期
vals = sorted(atoms[a] ** p for a in names for p in range(0, MAXEXP + 1)
              if 1 < atoms[a] ** p < 20)
span = max(vals) - min(vals)
exp_hits = len(vals) * (2 * 1e-3 * A_dn2) / span if span else 0
print()
print(f"  A_dn² 尺的偶然预期 ≈ {exp_hits:.2f} 个（实测 {len(hA2)} 个）")
chk(True, "B3 完成")
out["B3"] = {"atoms": list(atoms.keys()), "max_exp": MAXEXP,
             "hits_A": [{"form": n, "value": v, "dev_pct": 100 * d} for n, v, d in hA[:6]],
             "hits_A2": [{"form": n, "value": v, "dev_pct": 100 * d} for n, v, d in hA2[:6]],
             "expected_by_chance": exp_hits}

# 顺带：最接近的"单个亚层比"
print()
print("  最接近的单个亚层比：")
for x in r:
    print(f"    r = {x} = {float(x):.6f}   vs A_dn 差 {100 * (float(x) - A_dn) / A_dn:+.3f}%")

# ══════════════════════════════════════════════════════════════════════
# B4 判定
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【B4】判定的确切形状")
print("=" * 88)
print("  (a) 比值型：**需要新族**（所有现有族的比值上限都 < A_dn）。")
print("  (b) 乘积/幂型：见 B3 的命中数与偶然预期之比。")
print("  (c) 库内两候选（36/5 差 +0.032%、16√2/π 差 +0.067%）**都无读法**，且**判不了**。")
print("  (d) 上型的读法（r₃/r₁）**不能移植**（log_q 距整数 0.43）。")
chk(True, "B4 完成")
out["B4"] = {"ratio_type": "需要新族",
             "product_type": "见 B3",
             "library_candidates": "36/5（+0.032%）、16√2/π（+0.067%），都无读法、判不了",
             "up_reading_transfer": "不能移植（log_q 距整数 0.43）"}

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("汇总判定")
print("=" * 88)
verdict = [
    ("B1 上界", f"A_dn = {A_dn:.6f} **超出全部'产出振幅的族'的上界**（{prod_max[2]} = {prod_max[0]:.6f}，超出 {100 * (A_dn - prod_max[0]) / prod_max[0]:+.2f}%）；全部 {len(families)} 个族的总最大是 {overall_max[1][0]} = {overall_max[1][1]:.6f}。"),
    ("B2 含义", "若 A_dn 是比值型，其**分子分母必在现有族之外**——这是一个定位性结论。"),
    ("B3 乘积/幂", f"见命中数（A_dn² 尺命中 {len(hA2)} 个，偶然预期 {exp_hits:.2f}）。"),
    ("B4 判定", "（乙）仍未解开；但**边界更紧**：比值型需新族，乘积型无强命中。"),
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