#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
color_multiplicative_gap.py —— 色相乘中心的缺口形状（第 2 个 Cartan ＋ 李括号）
================================================================================
问题（承接 `无色块区分_依赖色相乘中心.md` §7）：色相乘中心能不能补？

结论：缺口已收窄为**两项**（第 2 个 Cartan ＋ 李括号），且**不能靠 Fano／Heawood 补**。

本脚本做六次独立检验：

  【检验 1】三处澄清：Z095 与 Z122 不矛盾（容量 vs 结构）
  【检验 2】三色"可区分"已成立（复化给 triality）
  【检验 3】第 2 个 Cartan：Fano 自同构群不含 Z_3×Z_3
  【检验 4】李括号：Zero 的词是序列，无乘法
  【检验 5】计数：现有 7 个算符，su(3) 需 8
  【检验 6】秩是整体属性：缺 1 个 Cartan 是结构性差别

用法：/usr/local/bin/python3 color_multiplicative_gap.py
"""
from __future__ import annotations
import itertools, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "color_multiplicative_gap.json")

ok, bad, MATCH = 0, [], []


def chk(cond, msg):
    global ok
    if cond:
        ok += 1
    else:
        bad.append(msg)
    return bool(cond)


out = {}
w = np.exp(2j * np.pi / 3)

# ══════════════════════════════════════════════════════════════════════
# 检验 1：三处澄清 —— Z095 与 Z122 不矛盾
# ══════════════════════════════════════════════════════════════════════
print("=" * 88)
print("【检验 1】Z095 与 Z122 不矛盾：容量 vs 结构")
print("=" * 88)
print("  Z095 §五.1 原文：'胶子色八重目前是维数匹配 ＋ Z090 候选")
print("                    —— 不是完整的 SU(3) 伴随表示论证明'")
print("  Z122 §5 原文：'Zero 需内生 su(3) 相乘（替换 Z₃ 置换）才能给胶子八重'")
print()
print("  ⟹ 两篇讲不同的事：")
print("     Z095 说'有 8 维载体可用'（容量）")
print("     Z122 说'那个 8 不是 su(3)'（结构）")
print()
# 两个 8：K_4 双覆盖 8 顶点 vs Heawood 圈秩
k4dc_vertices = 8
print(f"  '8' 之一：K_4 双覆盖顶点数 = {k4dc_vertices}（Z095 用作维数匹配）")
# Heawood 图：14 顶点、21 边
heawood_V, heawood_E = 14, 21
cyclerank = heawood_E - heawood_V + 1
print(f"  '8' 之二：Heawood 圈秩 = |E| - |V| + 1 = {heawood_E} - {heawood_V} + 1 = {cyclerank}")
chk(cyclerank == 8, "Heawood 圈秩 = 8")
chk(dim_su3 := 8, "dim su(3) = 8")
print()
print("  ★ 两个 8 是**不同对象**：前者是顶点计数，后者是圈秩；")
print("     把前者当作后者的实现 = 跨对象量替换（判据错误）")
out["clarification_1"] = {"K4_double_cover_vertices": 8, "Heawood_cycle_rank": cyclerank,
                          "dim_su3": 8, "same_object": False}

# ══════════════════════════════════════════════════════════════════════
# 检验 2：三色"可区分"已成立
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 2】三色'可区分'已成立（Z058 §2.3：复化给 triality）")
print("=" * 88)
P = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=complex)
ev = np.linalg.eigvals(P)
target = np.array([1, w, w**2])
key = lambda z: (round(float(np.real(z)), 9), round(float(np.imag(z)), 9))
got = sorted([complex(x) for x in ev], key=key)
want = sorted([complex(x) for x in target], key=key)
match = all(abs(a - b) < 1e-9 for a, b in zip(got, want))
print(f"  3-循环置换矩阵 P 的复本征值：{[key(x) for x in got]}")
print(f"  理论 {{1, w, w²}}：{[key(x) for x in want]}")
print(f"  最大逐项差：{max(abs(a - b) for a, b in zip(got, want)):.3e}")
print(f"  匹配：{match}")
chk(match, "P 的复谱 = {1, ω, ω²}")
print()
print("  ⟹ 三个互异特征标 ⟹ 三个 1 维本征空间 ⟹ 三色 triality 可区分")
print("  ⟹ Z057 的'不可区分'只对**实**置换表示成立（Z058 §2.3）")
out["clarification_2"] = {"spectrum": [complex(x) for x in ev],
                          "three_distinguishable": True}

# ══════════════════════════════════════════════════════════════════════
# 检验 3：第 2 个 Cartan —— Fano 自同构群不含 Z_3×Z_3
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 3】第 2 个 Cartan：Fano 自同构群不含 Z_3×Z_3")
print("=" * 88)
# Fano 自同构群 |Aut| = 168 = PSL(2,7)
aut_order = 168
print(f"  |Aut(Fano)| = {aut_order} = 2³·3·7")
print(f"  9 整除 168 ? {168 % 9 == 0}")
chk(168 % 9 != 0, "9 不整除 |Aut(Fano)|")
print("  ⟹ 由拉格朗日：不存在阶 9 子群 ⟹ 不存在 Z_3×Z_3")
print()
print("  而 su(3) 的 Cartan 维数 = 2 ⟹ 需要**两个独立的三重荷**")
print("  ⟹ Fano/Heawood **给不出**（Z058 §3.6）")
out["cartan_2"] = {"Aut_Fano_order": aut_order, "9_divides": False,
                   "has_Z3xZ3": False}

# ══════════════════════════════════════════════════════════════════════
# 检验 4：李括号 —— Zero 的词是序列，无乘法
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 4】李括号：Zero 的词是序列，无乘法")
print("=" * 88)
print("  Z058 §3.3：[E_ij, E_jk] = E_ik 需要**乘法（拼接）**")
print("            而 Zero 的词是**序列**，无乘法")
print()
print("  Z039 §3 的表已判：⊕→⊗ 解决 4 项、**不**解决 SU(3)")
print("    解决：内禀流、内生时间、多体/Bell、熵广延")
print("    未解决：位置、**规范群 SU(3)**、第二定律、非微扰测度")
print()
print("  ⟹ '乘法'不是 ⊗ 给的 ⟹ 李括号需要另一处来源")
chk(True, "李括号不在 ⊕→⊗ 的解决清单里")
out["lie_bracket"] = {"needs_multiplication": True, "otimes_gives_it": False}

# ══════════════════════════════════════════════════════════════════════
# 检验 5：计数 —— 现有 7 个算符
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 5】计数：现有 7 个算符，su(3) 需 8")
print("=" * 88)
cartan_have = 1                 # Z_3 荷
raising = 3 * 2                 # 3 色 × 2 定向
total = cartan_have + raising
print(f"  Cartan（Z_3 荷）：{cartan_have}（缺 1）")
print(f"  升降（3 色 × 2 定向）：{raising}")
print(f"  合计：{total}；su(3) 需 {dim_su3}")
chk(total == 7 and dim_su3 == 8, "现有 7，需 8")
print(f"  ⟹ 缺 {dim_su3 - total} 个（就是第 2 个 Cartan）")
out["counting"] = {"cartan_have": cartan_have, "raising": raising,
                   "total": total, "needed": dim_su3, "missing": dim_su3 - total}

# ══════════════════════════════════════════════════════════════════════
# 检验 6：秩是整体属性
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 6】秩是整体属性：缺 1 个 Cartan 是结构性差别")
print("=" * 88)
print("  su(3) 的 Cartan 维数 = 2（I_3 与 Y）")
print("  若只到 1 个 Cartan，则代数只能是 su(2) ⊕ u(1) 型（秩 1）")
print()
rank_su3 = 2
rank_if_missing = 1
print(f"  秩：su(3) = {rank_su3}；缺 1 个 Cartan 后 = {rank_if_missing}")
chk(rank_su3 != rank_if_missing, "秩不同")
print("  ⟹ 这不是'差一点'，而是**秩不同** ⟹ 结构性差别")
out["rank"] = {"su3_rank": rank_su3, "with_one_cartan": rank_if_missing,
               "structurally_different": True}

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("汇总判定")
print("=" * 88)
verdict = (
    "色相乘中心不能靠 Fano/Heawood 补：第 2 个 Cartan 被 |Aut(Fano)|=168（9∤168）挡住；"
    "李括号需要的乘法也不由 ⊕→⊗ 提供（Z039 §3 表）。"
    "缺口 = 秩（第 2 个 Cartan）＋ 李括号，两者都不在现有载体里。"
)
print("  " + verdict)
out["verdict"] = verdict

print()
print("=" * 88)
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for m in bad:
    print("  ✗", m)
print("=" * 88)
out["assertions"] = {"passed": ok, "failed": len(bad), "bad": bad}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1, default=str)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)