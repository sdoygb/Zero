#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routeA_27_slot_criterion.py —— 路线 A 第 0 步：把"中间谱"写成**可判定的槽位判据**
================================================================================
问题（承接 _staging/推导/L的下一步_环③中间谱内生.md 路线 A）：
  Z155 的中间谱 = SM + 3 代色单态 SU(2) 双重态，是**具名输入**。
  它能否从 Zero 的 16 槽载体**再一次扩张**中内生？

本脚本做三件事（都是**可复算的**，不含新假设）：

  【一】计数核验：E6 的 27 = 16 + 10 + 1；Z136 已把 16 落槽；剩下 11 个态待落。
  【二】β 位移核验：中间谱对 (b_L, b_Y) 的贡献必须是 (+2, +2)（Z155 (2.2)）。
       由此**反定**中间谱的态数 —— 这是一个把"具名输入"变成"可判定"的步骤。
  【三】把判据写成一条显式命题，供后续在 16 槽载体上验证。

纪律：本脚本发现的单位混淆（b_Y vs b_1）在 §3 登记，防复用。

用法：/usr/local/bin/python3 routeA_27_slot_criterion.py
"""
from __future__ import annotations
import json, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "routeA_27_slot_criterion.json")

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
# §1  计数：27 = 16 + 10 + 1
# ══════════════════════════════════════════════════════════════════════
print("=" * 88)
print("§1  计数核验：E6 的 27 ⊃ SO(10) 的 16，且 16 已由 Z136 落槽")
print("=" * 88)
print(f"  E6 的 27 = 16 + 10 + 1  ⟹  计数 16+10+1 = {16 + 10 + 1}")
chk(16 + 10 + 1 == 27, "27 = 16 + 10 + 1 计数相符")
print(f"  Z136：每代 16 槽 = SO(10) 旋量 5̄ + 10 + 1  ⟹  计数 5+10+1 = {5 + 10 + 1}")
chk(5 + 10 + 1 == 16, "16 槽 = 5̄+10+1 计数相符")
print("  ⟹ 27 分解为 (Z136 已落槽的 16) + (待落的 10+1 = 11 个态)")
print("     中间谱的候选对象在剩下的 11 个态里")
out["counting"] = {"27_decomp": "16+10+1", "Z136_decomp": "5bar+10+1",
                    "remaining_states": 11}

# ══════════════════════════════════════════════════════════════════════
# §2  β 位移反定态数
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("§2  β 位移反定中间谱的态数（Z155 (2.2): b_L→b_L+2, b_Y→b_Y+2）")
print("=" * 88)
# 一个复的色单态 SU(2) 双重态（N 或 L'，各含 2 个复分量）
#   b_L 贡献：(2/3)·T(2) = (2/3)(1/2) = 1/3
#   b_1 贡献：(2/3)·(3/5)·Σ(Y/2)² = (2/3)(3/5)(1/2) = 1/5
# 注意 b_Y 与 b_1 的关系：b_Y = (5/3)·b_1（Z155 (2.1)）
bL_per = F(2, 3) * F(1, 2)
b1_per = F(2, 3) * F(3, 5) * F(1, 2)
print(f"  单个复双重态：b_L 贡献 = (2/3)·T(2) = {bL_per}")
print(f"                b_1 贡献 = (2/3)·(3/5)·Σ(Y/2)² = {b1_per}")
print(f"  换算：b_Y = (5/3)·b_1  ⟹  b_Y 位移 = (5/3)·(b_1 位移)")

# 要求：b_L 位移 = +2，b_Y 位移 = +2
n_from_bL = F(2) / bL_per
n_from_bY = F(2) / (F(5, 3) * b1_per)
print(f"\n  由 b_L 位移 = +2 反定：需要 {n_from_bL} 个复双重态")
print(f"  由 b_Y 位移 = +2 反定：需要 {n_from_bY} 个复双重态")
chk(n_from_bL == n_from_bY, "两个位移给出一致的态数")
print(f"  ⟹ 两者一致：需要 {n_from_bL} 个复的色单态 SU(2) 双重态")
n_per_gen = n_from_bL / 3
print(f"  ⟹ 每代 {n_per_gen} 个（= N 与 L'，各 1 个）")
chk(n_per_gen == 2, "每代 2 个色单态双重态")
out["beta_shift"] = {"bL_per_doublet": str(bL_per), "b1_per_doublet": str(b1_per),
                     "n_doublets_total": str(n_from_bL), "n_per_generation": str(n_per_gen)}

# 复现 Z155 (2.1)-(2.3) 的数值
print("\n  复现 Z155 (2.1)–(2.3)：")
bL_sm, bY_sm = F(-19, 6), F(5, 3) * F(41, 10)
bR_sm = F(3, 4) * (bY_sm - bL_sm / 3)
bR_int = F(3, 4) * ((bY_sm + 2) - (bL_sm + 2) / 3)
print(f"    SM 段  b = (-7, {bL_sm}={float(bL_sm):.4f}, {bR_sm}={float(bR_sm):.4f})"
      f"  ← Z155 报 (-7, -3.167, 5.917)")
print(f"    中间段 b = (-7, {bL_sm + 2}={float(bL_sm + 2):.4f}, {bR_int}={float(bR_int):.4f})"
      f"  ← Z155 报 (-7, -1.167, 6.917)")
chk(abs(float(bL_sm) + 3.167) < 0.001 and abs(float(bR_sm) - 5.917) < 0.001,
    "SM 段 b 复现 Z155 (2.3)")
chk(abs(float(bL_sm + 2) + 1.167) < 0.001 and abs(float(bR_int) - 6.917) < 0.001,
    "中间段 b 复现 Z155 (2.3)")
out["Z155_reproduced"] = {"b_SM": [-7, float(bL_sm), float(bR_sm)],
                          "b_int": [-7, float(bL_sm + 2), float(bR_int)]}

# ══════════════════════════════════════════════════════════════════════
# §3  判据（供后续在 16 槽载体上验证）
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("§3  路线 A 的判据（显式命题）")
print("=" * 88)
criterion = (
    "在 K_4 的 16 槽载体（Z136：V(K_4)×Z_2^iso×Z_2^chir）之外，"
    "存在一次新的 Z_2 扩张，其新增槽位恰好给出："
    "每代 2 个色单态 SU(2) 双重态（Y/2 = ±1/2，手征配对），共 3 代 6 个；"
    "且这 6 个态对 (b_L, b_1) 的贡献恰为 (+2, +6/5)（等价于 b_Y 位移 +2）。"
)
print("  " + criterion)
print()
print("  ★ 该命题的可判定性：β 位移是**精确有理数**，故判据是离散的、非拟合的。")
print("  ★ 若成立 ⟹ 环③ 从'具名输入'变'算的'；若否 ⟹ 环③ 是一条内容层 no-go。")
out["criterion"] = criterion

# ══════════════════════════════════════════════════════════════════════
# §4  登记的判据错误（防复用）
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("§4  本轮出现并已修正的判据错误（登记，防复用）")
print("=" * 88)
print("  错误：把 b_Y 与 b_1 混为一谈。Z155 (2.1) 用 b_Y = (5/3)·b_1，")
print("        故 'b_Y 抬 +2' 等价于 'b_1 抬 +6/5'，不是 'b_1 抬 +2'。")
print("        直接按 b_1 比会得'需要 10 个双重态'的假结论（正确为 6 个）。")
out["errors_fixed"] = ["b_Y 与 b_1 混淆：b_Y = (5/3)b_1，故 +2 的 b_Y 位移 = +6/5 的 b_1 位移"]

# ══════════════════════════════════════════════════════════════════════
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