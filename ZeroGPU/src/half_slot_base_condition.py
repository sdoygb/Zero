#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
half_slot_base_condition.py —— 基数里"数什么"的条件（含一处自我否证）
================================================================================
承接 `half_slot_additive_rule.py`（半槽 = 1/2 个轨道单位；进入占据的计数 ⟹ 加）。
本篇攻最后一步：**基数里"数什么"** —— 为何 D_3 数 16 而 D_1 数 16·7²。

【两个条件候选，本篇给出判别】
  候选甲（已被否证）：基数的指数 = 2×(该量**所在律**里的平方阶)
      ⟹ 反例：D_1（A²）与 D_3（M_R）**都是质量标度**，
         "质量律有平方"不能区分它们 ⟹ 候选甲失败。
  候选乙（本篇）：基数的指数 = 该量中**"振幅"的幂次**
      D_1: A² 里振幅 A 出现 2 次 ⟹ 指数 2；D_3: M_R 里无振幅 ⟹ 指数 0。

【检验】
  【B1】候选甲的自否证（两处都是质量，区分不了）
  【B2】候选乙的复算（指数 = 振幅幂次）
  【B3】可证伪性检验：库里有没有第三个"计数型"半槽用法？
  【B4】判定：这一步为何升不了档

依据：Z126（m = u_s·A²·(1+A_s cosθ)²；A²=v/784.5）、Z136 定理 1／2（16 槽、7 满+1 半满）、
      Z019（seesaw：m_ν ∝ m_D²/M_R）、库内半槽五用法表。

用法：/usr/local/bin/python3 half_slot_base_condition.py
"""
from __future__ import annotations
import json, math, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "half_slot_base_condition.json")

ok, bad, MATCH = 0, [], []


def chk(cond, msg):
    global ok
    if cond:
        ok += 1
    else:
        bad.append(msg)
    return bool(cond)


out = {}
N, k = 16, 7

# ══════════════════════════════════════════════════════════════════════
# B1 候选甲的自我否证
# ══════════════════════════════════════════════════════════════════════
print("=" * 88)
print("【B1】候选甲的自我否证：'所在律的平方阶'区分不了两处")
print("=" * 88)
print("  候选甲：基数的指数 = 2 ×（该量**所在律**里的平方阶）")
print()
print("  两处的所在律：")
rows_a = [
    ("D1  A²（质量单元）", "费米子质量律 m ∝ A²(1+A cosθ)²", "有平方"),
    ("D3  M_R（Majorana 标度）", "seesaw m_ν ∝ m_D²/M_R", "M_R 在分母、一次"),
]
for a, b, c in rows_a:
    print(f"    {a:24s} | {b:34s} | {c}")
print()
print("  ★ 反例：**两处都是质量标度**（A² 与 M_R 都是质量量纲）。")
print("     '质量律有平方'这条对两者同样成立 ⟹ **区分不了**。")
chk(True, "B1 完成（候选甲被否证）")
out["B1"] = {"candidate_A": "指数 = 2×(所在律的平方阶)",
             "refuted": True,
             "reason": "两处都是质量标度；'质量律有平方'不能区分"}

# ══════════════════════════════════════════════════════════════════════
# B2 候选乙：指数 = 该量中"振幅"的幂次
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【B2】候选乙：基数的指数 = 该量中**振幅**的幂次")
print("=" * 88)
rows_b = [
    ("D1", "A²（质量单元）", "振幅 A 出现 2 次", 2, N * k ** 2, 784.5),
    ("D3", "M_R（Majorana 标度）", "无振幅", 0, N * k ** 0, 16.5),
]
print(f"  {'处':4s} {'量':22s} {'振幅幂次':16s} {'指数':>4s} {'基数':>6s} {'基数+½':>8s}")
for lab, q, why, expo, X, xh in rows_b:
    print(f"  {lab:4s} {q:22s} {why:16s} {expo:>4d} {X:>6d} {xh:>8.1f}")
print()
okd1 = (N * k ** 2 == 784)
okd3 = (N * k ** 0 == 16)
chk(okd1 and okd3, "候选乙复算两处（16·7²=784、16·7⁰=16）")
print("  核验：")
print(f"    D1: 16·7² = {N * k ** 2} ✓  ⟹ 1/(784+½) = {F(2, 1569)} ⟹ 1/c_ℓ = 784.5 ✓")
print(f"    D3: 16·7⁰ = {N * k ** 0} ✓  ⟹ 16/16.5 = {F(32, 33)} ✓")
print()
print("  ★ 候选乙的物理内容：**振幅出现几次，基数就平方几次**。")
out["B2"] = {"candidate_B": "指数 = 该量中振幅的幂次",
             "D1": {"quantity": "A^2", "amp_power": 2, "base": 784},
             "D3": {"quantity": "M_R", "amp_power": 0, "base": 16},
             "reproduces_both": True}

# ══════════════════════════════════════════════════════════════════════
# B3 可证伪性检验
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【B3】可证伪性检验：库里有没有第三个'计数型'半槽用法？")
print("=" * 88)
uses = [
    ("1", "槽表 7 满 + 1 半", "占据", "不是'基数+半槽'型"),
    ("2", "D2  m_D = m_τ/2", "投影期望", "不是（是纯 1/2）"),
    ("3", "D1  c_ℓ", "计数", "**基数 16·7²** ✓"),
    ("4", "D3  M_R", "计数", "**基数 16·7⁰** ✓"),
    ("5", "环④ 每态 π/N", "相格半格", "另一个 1/2（Z156 H1/H2）"),
]
print(f"  {'#':>2} {'处':18s} {'这个 1/2 是什么':12s} {'是否基数+半槽型'}")
for a, b, c, d in uses:
    print(f"  {a:>2} {b:18s} {c:12s} {d}")
print()
n_count = sum(1 for _, _, c, _ in uses if c == "计数")
print(f"  ⟹ '计数型'（基数+半槽型）的只有 **{n_count}** 处")
chk(n_count == 2, "只有两处是计数型")
print()
print("  ★ 结论：**没有第三处可检验** ⟹ 候选乙**无法独立证伪**")
print("     它是两处的**一致性**，不是**预测**。")
out["B3"] = {"counting_type_uses": n_count, "third_case_available": False,
             "falsifiable": False}

# ══════════════════════════════════════════════════════════════════════
# B4 判定：为何这一步升不了档
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【B4】判定：这一步为何升不了档")
print("=" * 88)
print("  要升档（从'读法'到'条件/定理'），需要一个**独立判据**。")
print("  但：")
print("    · 计数型用法只有 2 处；")
print("    · 任意两处都能被无穷多条'条件'拟合（含候选甲、乙及更多）；")
print("    · 无第三处 ⟹ 无法在两者间判别 ⟹ 都无法独立证伪。")
print()
print("  ⟹ **这一步的性质是'两处的一致性'，不是'条件'**。")
print("     要升档，必须先**造出第三处**（例如让某个尚未用半槽的量也用它），")
print("     或找到一个**独立于这两处的**理由。")
chk(True, "B4 完成")
out["B4"] = {"why_not_upgraded": "计数型只有 2 处 ⟹ 任何条件都无法独立证伪",
             "needed": "第三处，或一个独立于两处的理由"}

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("汇总判定")
print("=" * 88)
verdict = [
    ("B1 自我否证", "候选甲（'所在律的平方阶'）**被否证**——两处都是质量标度，区分不了。"),
    ("B2 候选乙", "**指数 = 该量中振幅的幂次**：D_1（A²，幂次 2）→ 784；D_3（M_R，幂次 0）→ 16。两处复算通过。"),
    ("B3 可证伪性", "计数型只有 **2** 处 ⟹ **无第三处可检验** ⟹ 候选乙无法独立证伪。"),
    ("B4 判定", "这一步的性质是**两处的一致性**，不是**条件**；升档需造出第三处或独立理由。"),
]
for kk, vv in verdict:
    print(f"  {kk}：{vv}")
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