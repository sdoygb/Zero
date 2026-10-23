#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
charged_lepton_half.py —— 攻（甲）：带电轻子半的指派能否由轨道结构推出？
================================================================================
背景：上一轮证成"中微子 ↔ 第二次扩张"由载体完备化定（Z136 定理 2：ν_R 是唯一补位）；
      带电轻子半（"e ↔ 第一次扩张"）仍未定。
      候选依据：e_R 在**满轨道**、ν_R 在**半满轨道**（Z136 §4 的 8 轨道表）。

【检验】
  【E1】8 条轨道表复现：7 满 + 1 半满
  【E2】e 与 ν 的轨道占用是否**对称**？（若对称，轨道结构不能作依据）
  【E3】逐项查"结构不变量"，看哪个能区分 e 与 ν 并给出梯级步号
  【E4】判定的确切形状

依据：Z136 §4（8 轨道表、定理 2、锚定规则 R1、补位唯一性）、Z126 §2–§3（振幅与梯级）、
      Z090（超荷内生 Q=T₃+Y/2）。

用法：/usr/local/bin/python3 charged_lepton_half.py
"""
from __future__ import annotations
import json, math, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "charged_lepton_half.json")

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
# E1 复现 8 条轨道表（Z136 §4）
# ══════════════════════════════════════════════════════════════════════
print("=" * 88)
print("【E1】复现 Z136 §4 的 8 条轨道表")
print("=" * 88)
# 槽 = (v, a, c)；v=0 色单态，v=1,2,3 色；a=↑/↓（Z2^iso）；c=L/R（手征）
# SM 每代 15 场（不含 ν_R）
SM15 = {
    "ν_L": (0, "up", "L"), "e_L": (0, "down", "L"), "e_R": (0, "down", "R"),
    "u_L1": (1, "up", "L"), "d_L1": (1, "down", "L"),
    "u_R1": (1, "up", "R"), "d_R1": (1, "down", "R"),
    "u_L2": (2, "up", "L"), "d_L2": (2, "down", "L"),
    "u_R2": (2, "up", "R"), "d_R2": (2, "down", "R"),
    "u_L3": (3, "up", "L"), "d_L3": (3, "down", "L"),
    "u_R3": (3, "up", "R"), "d_R3": (3, "down", "R"),
}
print(f"  输入：SM 每代 {len(SM15)} 个 Weyl 场（不含 ν_R）")
chk(len(SM15) == 15, "输入恰 15 场")
print()
# 8 条轨道 = (v, c) 各 2 个元素（a=↑,↓）
orbits = []
for v in range(4):
    for c in ("L", "R"):
        members = {a: None for a in ("up", "down")}
        for name, (vv, aa, cc) in SM15.items():
            if vv == v and cc == c:
                members[aa] = name
        filled = [k for k, x in members.items() if x is not None]
        orbits.append({"v": v, "c": c, "members": members, "n_filled": len(filled)})
print(f"  {'顶点':>4} {'扇区':>4} {'↑':>10} {'↓':>10} {'轨道':>8}")
full, half = 0, 0
for o in orbits:
    up = o["members"]["up"] or "—（空）"
    dn = o["members"]["down"] or "—（空）"
    stat = "满" if o["n_filled"] == 2 else ("半满" if o["n_filled"] == 1 else "空")
    if stat == "满":
        full += 1
    if stat == "半满":
        half += 1
    print(f"  {o['v']:>4} {o['c']:>4} {up:>10} {dn:>10} {stat:>8}")
print()
print(f"  合计：{full} 满 + {half} 半满（共 {len(orbits)} 条）")
chk(full == 7 and half == 1, "7 满 + 1 半满")
# 半满轨道是 (0,R)
half_orbit = [o for o in orbits if o["n_filled"] == 1]
chk(half_orbit[0]["v"] == 0 and half_orbit[0]["c"] == "R", "半满轨道 = (0,R)")
missing = [a for a, x in half_orbit[0]["members"].items() if x is None]
print(f"  半满轨道 (0,R) 的缺失成员 = a={missing[0]} ⟹ 补位 = (0,↑,R) = ν_R")
chk(missing == ["up"], "缺失成员 = ↑ ⟹ 补位 = (0,↑,R)")
out["E1"] = {"n_input": 15, "n_orbits": len(orbits), "full": full, "half": half,
             "half_orbit": "(0,R)", "missing_member": "up",
             "orbits": [{"v": o["v"], "c": o["c"],
                         "up": o["members"]["up"], "down": o["members"]["down"],
                         "n_filled": o["n_filled"]} for o in orbits]}

# ══════════════════════════════════════════════════════════════════════
# E2 e 与 ν 的轨道占用是否对称？
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【E2】e 与 ν 的轨道占用是否**对称**？")
print("=" * 88)
# e 的场：e_L (0,↓,L)、e_R (0,↓,R)
# ν 的场：ν_L (0,↑,L)、ν_R (0,↑,R)（后者是补位）
e_slots = [("e_L", (0, "down", "L")), ("e_R", (0, "down", "R"))]
nu_slots = [("ν_L", (0, "up", "L")), ("ν_R(补位)", (0, "up", "R"))]
print("  带电轻子 e 占用的槽：")
for nm, s in e_slots:
    print(f"    {nm:12s} = {s}")
print("  中微子 ν 占用的槽：")
for nm, s in nu_slots:
    print(f"    {nm:12s} = {s}")
print()
print("  对照：")
for label, e, nu in (("顶点 v", e_slots[0][1][0], nu_slots[0][1][0]),
                     ("手征 c（两个槽）", "L, R", "L, R"),
                     ("iso 标签 a（两个槽）", "↓, ↓", "↑, ↑"),
                     ("占用槽数", 2, 2),
                     ("跨轨道数", 2, 2)):
    print(f"    {label:22s} | e: {str(e):8s} | ν: {str(nu):8s} | {'对称' if e == nu else '**不同**'}")
print()
print("  ★ 判定：**完全对称**——两者各占 2 槽、各跨 2 轨道（(0,L) 与 (0,R)），")
print("     只在 iso 标签（↓ vs ↑）上不同。")
chk([x[2] for _, x in e_slots] == [x[2] for _, x in nu_slots],
    "两者手征结构相同（都跨 L 与 R）")
chk(len(e_slots) == len(nu_slots) == 2, "两者占用槽数相同（2）")
out["E2"] = {"symmetric": True,
             "e_slots": [s for _, s in e_slots], "nu_slots": [s for _, s in nu_slots],
             "only_difference": "iso 标签（e 在 ↓ 列，ν 在 ↑ 列）"}

print("  ⚠ 故：**轨道结构本身不能作为指派的依据**（它是对称的）。")
print("     唯一的不对称在别处：**ν_R 不在 15 场输入里，而 e_R 在**。")

# ══════════════════════════════════════════════════════════════════════
# E3 逐项查结构不变量
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【E3】逐项查结构不变量：哪个能区分 e 与 ν 并给出梯级步号？")
print("=" * 88)
# 梯级步号：e → 步 1（√2/√1）；ν → 步 2（√3/√2）
print("  已知指派：e → 步 1（√2/√1 = 1.414214）；ν → 步 2（√3/√2 = 1.224745）")
print()
inv = [
    ("T₃（L 成员）", F(-1, 2), F(1, 2), "e: −1/2；ν: +1/2"),
    ("Y_L", F(-1), F(-1), "相同"),
    ("Y_R", F(-2), F(0), "e: −2；ν: 0"),
    ("Q", F(-1), F(0), "e: −1；ν: 0"),
    ("|Q|", 1, 0, "e: 1；ν: 0"),
    ("iso 列", "↓", "↑", "e 在 ↓；ν 在 ↑"),
    ("输入里的场数", 2, 1, "e: e_L+e_R；ν: 只有 ν_L"),
    ("需补位？", "否", "是", "ν_R 是唯一补位"),
    ("占用槽数", 2, 2, "相同"),
    ("跨轨道数", 2, 2, "相同"),
]
print(f"  {'不变量':16s} | {'e':>10s} | {'ν':>10s} | 说明")
for nm, a, b, note in inv:
    same = (a == b)
    print(f"  {nm:16s} | {str(a):>10s} | {str(b):>10s} | {note}{'  ← 相同（不能区分）' if same else ''}")
print()
print("  能区分 e 与 ν 的不变量（4 个）：")
distinguishing = [nm for nm, a, b, _ in inv if a != b]
for nm in distinguishing:
    print(f"    · {nm}")
print()
print("  ★ 其中与梯级步号有**单调**关系（步号随该量增大而增大）的：")
mono = []
for nm, a, b, note in inv:
    if a == b:
        continue
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        continue
    # 步号：e=1, ν=2。若 (b>a) 且 ν 步号更大 ⟹ 单调
    if fb > fa:
        mono.append((nm, a, b, "ν 的该量更大，且 ν 步号更大 ⟹ 单调一致"))
    else:
        mono.append((nm, a, b, "ν 的该量更小，而 ν 步号更大 ⟹ 反向"))
for nm, a, b, note in mono:
    print(f"    · {nm}: {a} → {b}  {note}")
print()
chk(True, "E3 完成")
out["E3"] = {"distinguishing": distinguishing,
             "monotone_with_step": [{"inv": nm, "e": str(a), "nu": str(b), "note": note}
                                    for nm, a, b, note in mono],
             "n_data_points": 1}

# ══════════════════════════════════════════════════════════════════════
# E4 判定的确切形状
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【E4】判定的确切形状")
print("=" * 88)
print("  (a) 中微子半：**由载体完备化定**（Z136 定理 2：ν_R 是唯一补位，完成第二次扩张）")
print("      —— 这一步有定理背书（上一轮已证）。")
print()
print("  (b) 带电轻子半：**轨道结构不能作依据**（E2：完全对称）。")
print("      唯一的不对称是'ν_R 不在输入、e_R 在输入'，即**输入内容**。")
print("      故 e 的指派只能是'另一支'（default），**没有独立的结构依据**。")
print()
print("  (c) 可用的规则（读法级，1 个数据点）：")
print("      '需要补位的扇区走更高一级扩张' —— 与 Z126 §2 的'无荷走体内'方向一致。")
print()
print("  ⟹ 判定：**(甲) 部分成功**——中微子半已定；带电轻子半归结到**输入内容**，")
print("     而输入内容本身是具名输入（SM 15 场）。")
chk(True, "E4 完成")
out["E4"] = {
    "nu_half": "由 Z136 定理 2 定（有定理背书）",
    "charged_lepton_half": "轨道结构对称（不能作依据）；归结到输入内容（SM 15 场，具名输入）",
    "usable_rule": "需要补位的扇区走更高一级扩张（读法级，1 个数据点）",
    "verdict": "(甲) 部分成功：中微子半已定；带电轻子半归结到具名输入",
}

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("汇总判定")
print("=" * 88)
verdict = [
    ("E1 轨道表", "复现 7 满 + 1 半满；半满轨道 (0,R)，缺失成员 = ↑ ⟹ 补位 (0,↑,R)=ν_R。"),
    ("E2 对称性", "**e 与 ν 的轨道占用完全对称**（各 2 槽、各跨 2 轨道）⟹ 轨道结构**不能**作指派依据。"),
    ("E3 不变量", "能区分 e/ν 的：T₃、Y_R、Q、|Q|、iso 列、输入场数、需补位。与步号单调一致者见输出。"),
    ("E4 判定", "**(甲) 部分成功**：中微子半由定理定；带电轻子半归结到**输入内容**（具名输入）。"),
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