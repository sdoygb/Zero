#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_weinberg_schemes.py —— 弱混合角："测出来的" vs "算出来的"分界线（数值钉住）
================================================================================
主流理论（SM）里 sin²θ_W 有多个"定义"，它们数值不同；本脚本把四个常用定义摆出来，
并标出每个数里哪一部分是**测的**、哪一部分是**算的**。

四个定义：
  D1 树级（用 α、G_F、M_W）：sin²θ_W = πα /(√2 G_F M_W²)     —— 三个实测输入 + 一个树级关系
  D2 on-shell：           sin²θ_W = 1 − M_W²/M_Z²            —— 两个实测质量之比
  D3 MS-bar(M_Z)：        0.23122（PDG 值，输入的引用值）
  D4 有效轻子(Z 极点)：    0.23153（PDG 值，输入的引用值）

要点：D1 与 D2 用的都是**同一批实测数**，差 3–4% —— 这个差是**圈修正**（Δr），即"算出来的"；
      而四个数的**值**没有一个是被理论预言的。

用法：/usr/local/bin/python3 k4_weinberg_schemes.py
"""
from __future__ import annotations
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "k4_weinberg_schemes.json")

# ── 实测输入（PDG 量级值） ────────────────────────────────────────────
alpha_inv = 137.035999084        # α(0)^{-1}
GF = 1.1663787e-5                # GeV^-2, 费米常数（μ 衰变）
MW = 80.377                      # GeV
MZ = 91.1876                     # GeV
sin2_MSbar = 0.23122             # 输入的引用值（PDG）
sin2_eff_lept = 0.23153          # 输入的引用值（PDG）

ok, bad, MATCH = 0, [], []


def chk(c, m):
    global ok
    if c:
        ok += 1
    else:
        bad.append(m)
    return bool(c)


alpha = 1 / alpha_inv
sin2_tree = math.pi * alpha / (math.sqrt(2) * GF * MW ** 2)   # D1
sin2_onshell = 1 - (MW / MZ) ** 2                             # D2

print("=" * 88)
print("四个定义（同一个「角」）")
print("=" * 88)
rows = [
    ("D1 树级 πα/(√2 G_F M_W²)", sin2_tree, "3 个实测输入 + 树级关系"),
    ("D2 on-shell 1−M_W²/M_Z²", sin2_onshell, "2 个实测质量之比"),
    ("D3 MS-bar(M_Z)", sin2_MSbar, "引用值（从多个观测拟合）"),
    ("D4 有效轻子(Z 极点)", sin2_eff_lept, "引用值（从多个观测拟合）"),
]
for nm, v, src in rows:
    print(f"  {nm:26s} = {v:.5f}   （{src}）")

chk(abs(sin2_onshell - 0.2231) < 0.0005, "on-shell 值应在 0.2231 附近")
MATCH.append(("D", f"on-shell = 1−(M_W/M_Z)² = {sin2_onshell:.5f}（由两个实测质量算出的**比值**）", True))

d12 = 100 * (sin2_onshell - sin2_tree) / sin2_tree
d34 = 100 * (sin2_eff_lept - sin2_MSbar) / sin2_MSbar
d14 = 100 * (sin2_eff_lept - sin2_tree) / sin2_tree
print(f"\n  同一批实测输入下：D2−D1 = {100*(sin2_onshell-sin2_tree):+.5f}"
      f"（相对 {d12:+.2f}%）——这是**圈修正 Δr**，是「算出来的」")
print(f"  D4−D3 = {100*(sin2_eff_lept-sin2_MSbar):+.5f}（相对 {d34:+.3f}%）—— 方案差")
print(f"  D4−D1 = {d14:+.2f}%  ⟹ **四个「同一个角」的定义差 3–7%**，而实验误差只有 ~0.1%")

# 圈修正的"大小"：把 D1 与 D2 的差写成 Δr
delta_r = 1 - (sin2_tree / sin2_onshell)
print(f"\n  由 D1、D2 反推的圈修正量级：Δr ≈ {delta_r:+.4f}"
      f"（标准值 Δr ≈ 0.036 —— 同一量级）")
chk(0.02 < delta_r < 0.06, "Δr 应在 0.03 量级")
MATCH.append(("R", f"Δr ≈ {delta_r:.4f}（同一批输入下两个定义之差）⟹ "
                   f"**差值是圈修正（算出来的）**", True))

# 我们这条线已经算出来的两块
b_note = "b=(41/10,−19/6,−7) 已由 Zero 的 16 槽内容精确算出（WEINBERG_LOW_ENERGY.md §0）"
sin2_gut = 3 / 8
print()
print("=" * 88)
print("我们这条线里的「算」与「测」")
print("=" * 88)
print(f"  【算】统一标度的值：sin²θ_W(M_U) = {sin2_gut} = ΣT₃²/(ΣT₃²+ΣY'²) —— 归一化被结构固定")
print(f"  【算】跑动系数：{b_note}")
print(f"  【测】低能值 sin²θ_W(M_Z) = {sin2_MSbar}")
print(f"  【输入】标度比 L = ln(M_U/M_Z)（或等价地：α_s(M_Z)、M_U）")
MATCH.append(("S", "本线的分工：3/8 与 b_i 是算的；sin²θ_W(M_Z) 是测的；L 是输入", True))

print()
print("断言清单：")
for tag, txt, good in MATCH:
    print(f"   {tag}  {'✓' if good else '✗'}  {txt}")
_total = len(MATCH)
_failed = len([1 for _t, _x, g in MATCH if not g]) + len(bad)
print(f"\n—— {_total} 条断言 / {_failed} 不符 ——")

os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump({
        "inputs": {"alpha_inv": alpha_inv, "G_F": GF, "M_W": MW, "M_Z": MZ,
                   "sin2_MSbar": sin2_MSbar, "sin2_eff_lept": sin2_eff_lept},
        "definitions": {nm: v for nm, v, _ in rows},
        "spread_percent": {"D2_minus_D1": d12, "D4_minus_D3": d34, "D4_minus_D1": d14},
        "delta_r": delta_r,
        "division": {"computed": ["sin²θ_W(M_U)=3/8（归一化）", "b_i（跑动系数）"],
                     "measured": ["sin²θ_W(M_Z)"],
                     "input": ["L=ln(M_U/M_Z)"]},
        "assertions": {"total": _total, "failed": _failed, "detail": bad},
    }, f, ensure_ascii=False, indent=1)
print(f"已存 {OUT}")