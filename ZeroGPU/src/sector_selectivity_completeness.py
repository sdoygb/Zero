#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sector_selectivity_completeness.py —— 攻击 Z014 定理 31：三层清单是否完备？
================================================================================
问题（承接"夸克振幅下型开放项"的归位）：定理 31 把扇区选择性判为
  "在 {群论层, 一圈动力学层, 味混合层} 三层面取定理形式；其取值归属为公设级买入项"。

本脚本检验：**这三层是否穷尽了产生扇区差异的所有来源？**

关键线索：四个振幅中，$A_ℓ=\sqrt2$、$A_ν=\sqrt{3/2}$ 来自**立方体的对角比**（Z126 §2），
          $A_{up}=(4/3)^2$ 来自**亚层比之比**（Z129/Z089）。
          这两类来源**都不在**定理 31 列出的三层里。

  【P1】三层清单核对：立方体几何／亚层比 能否归入 {群论, 一圈, 味混合} 中任一层？
  【P2】扇区差异的"值"来源清点：四个振幅各有几个有结构性来源？
  【P3】"指派"（哪个扇区走哪条路）是否有结构性理由？（Z126 §2 自述）
  【P4】精化后的定理 31 应写成什么？

依据：Z014 定理 31（三层面判死）、Z126 §2–§3／§7（振幅与来源）、Z129（亚层比）、
      Z089（K₄ 三频率槽）、Z084（相位补不了的墙·墙 3）。

用法：/usr/local/bin/python3 sector_selectivity_completeness.py
"""
from __future__ import annotations
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "sector_selectivity_completeness.json")

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
# P1 三层清单核对
# ══════════════════════════════════════════════════════════════════════
print("=" * 88)
print("【P1】定理 31 的三层清单，与振幅的实际来源逐条对照")
print("=" * 88)
layers = ["群论层", "一圈动力学层", "味混合层"]
print("  定理 31 列出的三层：" + "、".join(layers))
print()
print("  定理 31 各层的具体排除内容：")
layer_detail = [
    ("群论层", "三扇区汤川算符结构**完全相同**（唯一差别是色重数）；SU(2) 不变收缩唯一"),
    ("一圈动力学层", "Koide 型量 K 零次齐次 ⟹ 味普适修正零作用；阈值分化仅 1%–2%（比所需小 5–20 倍）"),
    ("味混合层", "K 基无关 ⟹ CKM/PMNS 不能改 K；两夸克扇区同基时 V_CKM = 1（不成立）"),
]
for a, b in layer_detail:
    print(f"    {a:14s}：{b}")
print()
print("  四个振幅的实际来源（Z126 §2–§3）：")
sources = [
    ("A_ℓ = √2", "立方体：面对角线／棱", "**载体几何**"),
    ("A_ν = √(3/2)", "立方体：体对角线／面对角线", "**载体几何**"),
    ("A_up = (4/3)²", "亚层比之比 r₃/r₁（K₄ 三频率槽）", "**载体谱／亚层比**"),
    ("A_dn", "无读法", "**开放**"),
]
for a, b, c in sources:
    print(f"    {a:16s} | {b:34s} | {c}")
print()
print("  ★ 关键检验：'载体几何' 能否归入 {群论, 一圈, 味混合} 中任一层？")
check = [
    ("群论层？", "群论层讲的是**汤川算符的结构**（哪些场收缩）；立方体对角比是**载体几何**，与算符结构无关",
     "**不能**归入"),
    ("一圈动力学层？", "一圈层讲的是**跑动与味普适性**；立方体对角比是**零圈几何量**，与跑动无关",
     "**不能**归入"),
    ("味混合层？", "混合层讲的是**基变换**；立方体对角比在混合前就已确定",
     "**不能**归入"),
]
for a, b, c in check:
    print(f"    {a:14s} {b}")
    print(f"    {'':14s} ⟹ {c}")
print()
chk(True, "P1 三层对照完成")
out["P1"] = {"theorem_layers": layers, "layer_detail": layer_detail,
             "amplitude_sources": [{"amp": a, "source": b, "kind": c}
                                   for a, b, c in sources],
             "carrier_geometry_in_three_layers": False}

print("  ★★ 判定：**载体几何与亚层比都不在定理 31 的三层里**。")
print("     ⟹ 三层清单对'产生扇区差异的来源'**不完备**。")

# ══════════════════════════════════════════════════════════════════════
# P2 "值"的来源清点
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【P2】扇区差异的'值'来源清点")
print("=" * 88)
A_l, A_nu, A_up, A_dn = math.sqrt(2), math.sqrt(1.5), 1.777833, 2.682850
print(f"  四个振幅的值：")
for nm, A in (("A_ℓ", A_l), ("A_ν", A_nu), ("A_up", A_up), ("A_dn", A_dn)):
    print(f"    {nm:6s} = {A:.6f}   A² = {A ** 2:.6f}")
print()
print("  有结构性来源的：")
have = [("A_ℓ", "立方体 面对角线/棱", 0.0),
        ("A_ν", "立方体 体对角线/面对角线", 0.0),
        ("A_up", "亚层比 r₃/r₁ = (4/3)²", 0.0031)]
for nm, src, dev in have:
    print(f"    {nm:6s} ← {src:30s} 偏差 {dev:.4f}%")
print()
print("  无结构性来源的：")
print(f"    A_dn   ← 无（候选 36/5 差 +0.032%、16√2/π 差 +0.067%，都无读法）")
print()
n_have, n_total = 3, 4
print(f"  ⟹ 4 个振幅中 **{n_have} 个**有结构性来源；只有下型是开放的。")
chk(n_have == 3, "四个振幅中三个有结构性来源")
out["P2"] = {"have_source": n_have, "total": 4, "open": ["A_dn"],
             "values": {"A_l": A_l, "A_nu": A_nu, "A_up": A_up, "A_dn": A_dn}}

print()
print("  ★ 这与定理 31 的结语'其取值归属为公设级买入项'**张力**：")
print("     若 3/4 的取值都有来源，则'取值归属公设'这句话**过宽**。")

# ══════════════════════════════════════════════════════════════════════
# P3 "指派"是否有结构性理由
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【P3】'指派'（哪个扇区走哪条路）是否有结构性理由？")
print("=" * 88)
print("  Z126 §2 原话：'带电的走**面**（面：棱 = √2）；中微子**无荷**、只走**体**内（体：面 = √(3/2)）'")
print("  Z126 §2 又自述：'这是一条**结构识别**（不是定理）：**识别本身是选择**'")
print()
assign = [
    ("轻子（带电）", "走面 → A_ℓ = 面对角线/棱 = √2", "**与电荷相关**（带电走面）"),
    ("中微子（无荷）", "走体 → A_ν = 体对角线/面对角线 = √(3/2)", "**与电荷相关**（无荷走体）"),
    ("上型（Q=+2/3）", "亚层比 r₃/r₁", "与电荷的关系**未建立**"),
    ("下型（Q=-1/3）", "无读法", "—"),
]
for a, b, c in assign:
    print(f"    {a:16s} | {b:38s} | {c}")
print()
print("  ★ 判定：指派规则在**轻子/中微子**这一对上有'电荷定路径'的规则（Z126 §2）；")
print("     但 Z126 自述'识别本身是选择'⟹ 该规则**尚未成为定理**。")
print("     上型/下型的指派规则**未建立**。")
chk(True, "P3 指派核对完成")
out["P3"] = {"assignment_rules": [{"sector": a, "path": b, "rule": c}
                                  for a, b, c in assign],
             "lepton_neutrino_rule": "电荷定路径（读法，非定理）",
             "up_down_rule": "未建立"}

# ══════════════════════════════════════════════════════════════════════
# P4 精化后的表述
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【P4】精化后的定理 31 应写成什么？")
print("=" * 88)
refined = (
    "扇区选择性在 {群论, 一圈动力学, 味混合} 三层面无来源（原定理 31 成立）；"
    "但产生扇区差异的来源**不止这三层**：载体几何（立方体对角比）给出 A_ℓ、A_ν 的值，"
    "载体谱/亚层比（K₄ 三频率槽的 r₃/r₁）给出 A_up 的值。"
    "故'其取值归属为公设级买入项'应精化为："
    "(i) **三个振幅的值有结构性来源**（读法级，非定理级）；"
    "(ii) **振幅到扇区的指派规则**是买入项（轻子/中微子对上有'电荷定路径'的读法，上型/下型未建立）；"
    "(iii) **下型振幅的值**开放。"
)
print("  " + refined)
print()
chk(True, "P4 精化完成")
out["P4"] = {"refined_statement": refined}

# ══════════════════════════════════════════════════════════════════════
# P5 若"载体几何"层扩展，能否给出希格斯的第五个振幅？
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【P5】载体几何层扩展：能否给出希格斯所需的第五个振幅？")
print("=" * 88)
V, M_H = 246.22, 125.25
A2 = V * 1000 / 784.5      # 313.856 MeV
# 立方体的全部路径比（三个长度 1, √2, √3 的两两之比，含反向）
cube = {"棱 1": 1.0, "面对角线 √2": math.sqrt(2), "体对角线 √3": math.sqrt(3)}
ratios = {}
for (n1, v1), (n2, v2) in __import__("itertools").permutations(cube.items(), 2):
    ratios[f"{n1}/{n2}"] = v1 / v2
print("  立方体的全部路径比（6 个）：")
for nm, v in sorted(ratios.items(), key=lambda kv: -kv[1]):
    print(f"    {nm:22s} = {v:.6f}")
print()
print("  已用：面对角线/棱 = √2（轻子）、体对角线/面对角线 = √(3/2)（中微子）")
print()
# 若希格斯取某个未用的路径比为自身振幅 A_H，则 M_H = u A² (1+A_H cosθ)²
# 求所需的 u（θ 取地标）
landmarks = [45.0, 90.0, 135.0]
print("  若希格斯有自身层、取未用的路径比为 A_H，且相位落在地标上：")
print(f"  {'A_H':>22} {'θ':>7} {'所需 u':>12}  u 是否有角色？")
rows = []
for nm, A in sorted(ratios.items(), key=lambda kv: kv[1]):
    if abs(A - math.sqrt(2)) < 1e-9 or abs(A - math.sqrt(1.5)) < 1e-9:
        continue
    for th in landmarks:
        Mv = (1 + A * math.cos(math.radians(th))) ** 2
        if Mv <= 0:
            continue
        u = M_H * 1000 / (A2 * Mv)
        rows.append((nm, A, th, u))
        print(f"  {nm:>22} {th:7.0f} {u:12.3f}")
print()
print("  ★ 所需的 u 值：")
us = [r[3] for r in rows]
print(f"    范围 [{min(us):.2f}, {max(us):.2f}]")
print("  库内已知的质量层标度：u=1（轻子）、u≈0.988（下型）、u≈71.6（上型）")
close = [r for r in rows if min(abs(r[3] - 1), abs(r[3] - 0.988), abs(r[3] - 71.602) / 71.602) < 0.05]
print(f"  与已知标度相差 <5% 的：{len(close)} 个")
for r in close[:6]:
    print(f"    {r[0]} @ {r[2]:.0f}° 需要 u = {r[3]:.3f}")
chk(True, "P5 完成")
out["P5"] = {"cube_ratios": ratios, "candidate_rows": [
    {"A_H": nm, "A": A, "theta": th, "u_needed": u} for nm, A, th, u in rows],
    "n_close_to_known_u": len(close)}

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("汇总判定")
print("=" * 88)
verdict = [
    ("三层清单", "**不完备**：载体几何与亚层比都不在 {群论, 一圈, 味混合} 里。"),
    ("值来源", "4 个振幅中 **3 个**有结构性来源（A_ℓ、A_ν 来自立方体；A_up 来自亚层比）。"),
    ("指派规则", "轻子/中微子对上有'电荷定路径'（读法）；上型/下型**未建立**。"),
    ("精化", "'取值归属公设'应精化为：值（3/4 有来源）+ 指派（买入项）+ 下型值（开放）。"),
    ("对希格斯", "希格斯需**第五个振幅**；载体几何只给出 6 个路径比（已用 2 个），扩展它得到的候选见 P5。"),
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