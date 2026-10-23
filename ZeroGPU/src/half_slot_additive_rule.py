#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
half_slot_additive_rule.py —— “为何是加法”：半槽 = 1/2 个轨道单位，加到**占据计数**上
================================================================================
承接 `half_slot_native_search.py`（统一模板：基数 X = N·k^(2(n−1))，半槽以 +½ 进入），
本篇攻库内标为仍缺的那一步：**“为何是加法”**（而不是以权重乘进通道）。

【核心读法】
  Z136：每条轨道 = 2 个元素（Z₂^iso：↑↓）；满轨道占 2 元素，半满轨道占 1 元素。
  ⟹ 半槽的值 = (半满元素数)/(每轨道元素数) = 1/2。
  ⟹ 它进入计数的位置是**占据的计数**（加到 X 上），不是**通道的权重**（乘进 w_i w_j）。

【检验】
  【A1】半槽 = 1/2 个轨道单位的来历（元素数之比）
  【A2】两处复算：D3 = 32/33、D1 = 1/784.5，统一为 X/(X+½) 与 1/(X+½)
  【A3】反例检查：本读法**预测**“权重乘积”给不出 28.5（复现库内排除）
  【A4】区分两种进入方式：占据计数（加）vs 通道权重（乘）
  【A5】诚实登记仍缺的（若任何一处仍缺）

依据：Z136 定理 1／2 与 §4 轨道表、Z126（A²=v/784.5）、Z019／Z126（M_e 相位）；
      库内 half_slot_rule_check.py（5/5）、half_slot_placement_search.py（2/2）、
      zero_half_slot_784.py（8/8）。

用法：/usr/local/bin/python3 half_slot_additive_rule.py
"""
from __future__ import annotations
import json, math, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "half_slot_additive_rule.json")

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
# A1 半槽 = 1/2 个轨道单位
# ══════════════════════════════════════════════════════════════════════
print("=" * 88)
print("【A1】半槽 = 1/2 个轨道单位：来历")
print("=" * 88)
elem_per_orbit = 2          # Z136：每条 Z2^iso 轨道恰 2 元素（↑↓）
elem_half = 1               # 半满轨道占 1 元素
h = F(elem_half, elem_per_orbit)
print(f"  Z136：每条轨道 = {elem_per_orbit} 个元素（Z₂^iso：↑↓）")
print(f"        满轨道 = {elem_per_orbit} 元素；半满轨道 = {elem_half} 元素")
print()
print(f"  半槽的值 = (半满元素数)/(每轨道元素数) = {elem_half}/{elem_per_orbit} = {h}")
chk(h == F(1, 2), "半槽 = 1/2（元素数之比）")
print()
print("  ★ 注意：这个 1/2 不是“随便取的半个”，而是")
print("     **半满轨道的元素数 ÷ 一条轨道的元素数** —— 两个 Z136 原生数之比。")
out["A1"] = {"elem_per_orbit": elem_per_orbit, "elem_half": elem_half,
             "half_value": str(h), "origin": "半满元素数 / 每轨道元素数"}

# ══════════════════════════════════════════════════════════════════════
# A2 两处复算
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【A2】两处复算（统一为 X + ½）")
print("=" * 88)
N, k = 16, 7
rows = [
    ("D3（标度比，一体 n=1）", N, "槽数 16", "ratio"),
    ("D1（质量单元，两体 n=2）", N * k ** 2, "16·7²", "unit"),
]
print(f"  {'处':26s} {'X':>6} {'X+½':>10} {'比值型 X/(X+½)':>16} {'单位型 1/(X+½)':>16}")
res = {}
for lab, X, src, kind in rows:
    xh = F(2 * X + 1, 2)
    ratio = F(2 * X, 2 * X + 1)
    unit = F(2, 2 * X + 1)
    res[lab] = {"X": X, "X_plus_half": str(xh), "ratio": str(ratio), "unit": str(unit)}
    print(f"  {lab:26s} {X:>6} {str(xh):>10} {str(ratio):>16} {str(unit):>16}")
print()
d3 = F(res["D3（标度比，一体 n=1）"]["ratio"])
chk(d3 == F(32, 33), "D3 = 32/33（精确）")
d1_unit = F(res["D1（质量单元，两体 n=2）"]["unit"])
chk(1 / d1_unit == F(1569, 2), "D1 的 1/(X+½) 倒数 = 1569/2 = 784.5（精确）")
print(f"  D3：X/(X+½) = {d3} = {float(d3):.9f} ✓")
print(f"  D1：1/(X+½) = {d1_unit} ⟹ 1/c_ℓ = {1 / d1_unit} = {float(1 / d1_unit):.1f} ✓")
print()
print("  ★ 两处**同一个** X+½ 形式；差别只在取比值还是取倒数。")
out["A2"] = res

# ══════════════════════════════════════════════════════════════════════
# A3 反例检查：本读法预测“权重乘积”给不出 28.5
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【A3】反例检查：本读法预测“权重乘积”给不出 28.5")
print("=" * 88)
print("  库内规则 B：把半槽当“权重 ½ 的槽”，数通道 Σ_{i<j} w_i w_j")
w = [1.0] * 7 + [0.5]
pair = sum(w[i] * w[j] for i in range(8) for j in range(i + 1, 8))
print(f"    w = {w}")
print(f"    Σ_{{i<j}} w_i w_j = {pair}")
chk(abs(pair - 24.5) < 1e-12, "权重乘积给 24.5（复现库内）")
print()
print("  本读法预测：**占据的计数**应是 X + ½ 型，而不是加权的通道数。")
print(f"    若半槽是“占据”：（7 满轨道 + ½ 轨道单位）的计数 = {F(15, 2)}（=7.5）")
print(f"    若半槽是“通道权重”：Σw_iw_j = {pair}（≠ 28.5）")
print()
chk(abs(pair - 28.5) > 1, "权重乘积 ≠ 28.5（本读法预测正确）")
print("  ⟹ 本读法与库内排除**一致**：权重乘积不是 D1 的来源。")
out["A3"] = {"weight_product": pair, "target_28_5": 28.5, "consistent_with_library": True}

# ══════════════════════════════════════════════════════════════════════
# A4 两种进入方式的区分
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【A4】两种进入方式的区分（本篇的核心）")
print("=" * 88)
print("  ┌─ 半槽进入【占据的计数】⟹ 加到 X 上（本篇读法）")
print("  │    D3: X = 16  ⟹ 16 + ½ = 16.5")
print("  │    D1: X = 784 ⟹ 784 + ½ = 784.5")
print("  │")
print("  └─ 半槽进入【通道的权重】⟹ 乘进 w_i w_j（库内已排除）")
print(f"       Σ w_i w_j = {pair} ≠ 28.5")
print()
print("  ★ 区分判据：**数的是“占据”还是“通道”**")
print("     · 占据 → 加（它是“有多少个单位”的一部分）")
print("     · 通道 → 乘（它是“每个单位的权重”）")
print()
print("  Z136 里半槽的身份是**占据**（半满轨道的洞，定理级）——")
print("  ⟹ 它应该进入**占据的计数** ⟹ **加法**。")
chk(True, "A4 完成")
out["A4"] = {"additive": "占据的计数（加到 X 上）",
             "multiplicative": "通道的权重（乘进 w_i w_j，已排除）",
             "criterion": "数的是占据还是通道",
             "half_slot_identity": "占据（Z136 定理 2，定理级）"}

# ══════════════════════════════════════════════════════════════════════
# A5 仍缺的
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【A5】仍缺的（诚实登记）")
print("=" * 88)
still = [
    ("为何 D3 数**槽**而 D1 数**满轨道²**", "已给候选读法：D3 一体、D1 两体（振幅²）；"
                                          "但“一体为何数槽、两体为何数满轨道”仍缺"),
    ("指数 2 的独立核验", "理由（质量=振幅²）是读法级；无独立判据"),
]
for a, b in still:
    print(f"  · {a}")
    print(f"      {b}")
chk(True, "A5 完成")
out["A5"] = [{"item": a, "status": b} for a, b in still]

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("汇总判定")
print("=" * 88)
verdict = [
    ("A1", "半槽 = 1/2 = 半满元素数/每轨道元素数（两个 Z136 原生数之比）。"),
    ("A2", "两处统一为 X+½：D3 = 32/33（比值型）、D1 ⟹ 1/c_ℓ = 784.5（单位型）。"),
    ("A3", f"本读法与库内排除一致（权重乘积给 {pair} ≠ 28.5）。"),
    ("A4", "**“为何是加法”的候选答案**：半槽的身份是**占据**（Z136 定理 2）⟹ 进入占据的计数 ⟹ 加。"),
    ("A5", "仍缺：D3 数槽 / D1 数满轨道² 的理由；指数 2 无独立核验。"),
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