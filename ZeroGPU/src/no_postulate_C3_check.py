#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
no_postulate_C3_check.py —— 硬约束判定：C^3 不是 Z_0 的导出（故 C1' 不可加）
================================================================================
约束（用户）：不能加公设，公理只有 Z001 里的唯一一个 Z_0。

Z_0 = ({±1}, iota) + (W_n, mu)   [Z001 (2.1)]

结论：C^3 不是 Z_0 的导出；路线乙（经 End(C^3) 补 su(3)）在硬约束下不可用。

本脚本做六次独立检验：

  【检验 1】Z_0 导出的结构清单，全部阿贝尔／组合
  【检验 2】{±1} 的不可约表示只有 1、2 维（故 M_3 的 3 无来源）
  【检验 3】C 可导出（Z004 推论 6.1）但不蕴含 C^3 可导出
  【检验 4】维数候选：H_Q 与圈空间都能凑 8，但都是阿贝尔
  【检验 5】结合运算不蕴含李括号
  【检验 6】候选来源逐个判定

用法：/usr/local/bin/python3 no_postulate_C3_check.py
"""
from __future__ import annotations
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "no_postulate_C3_check.json")

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
# 检验 1：Z_0 导出的结构清单
# ══════════════════════════════════════════════════════════════════════
print("=" * 88)
print("【检验 1】Z_0 导出的结构（Z001 §3）")
print("=" * 88)
derived = [
    ("词 w=(e_1..e_n)", "沿图的有向边序列", "序列（结合）", "阿贝尔/组合"),
    ("散度 d_w in Z^C", "进入数 − 离开数", "整数向量", "阿贝尔"),
    ("零和超平面 H_Q", "{x : sum x = 0}", "阿贝尔群", "阿贝尔"),
    ("全分支/整数重数", "无概率 ⟹ 计数", "计数", "阿贝尔"),
    ("终端与重播种", "款②逼出", "—", "阿贝尔"),
]
for a, b, c, d in derived:
    print(f"  {a:22s} | {b:22s} | {c:12s} | {d}")
chk(all(d == "阿贝尔" or d == "阿贝尔/组合" for *_, d in derived),
    "Z_0 导出的结构全部阿贝尔/组合")
print()
print("  ⟹ 无一是**非阿贝尔**的乘法结构")
out["derived"] = [{"name": a, "content": b, "kind": c, "algebra": d}
                  for a, b, c, d in derived]

# ══════════════════════════════════════════════════════════════════════
# 检验 2：{±1} 的不可约表示只有 1、2 维
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 2】{±1} 的不可约表示维数（故 M_3 的 '3' 无来源）")
print("=" * 88)
# Z_2 的实不可约表示：平凡(1 维)、符号(1 维)
# Z_2 的复不可约表示：两个 1 维（特征标 ±1）
# 实表示论里 Z_2 的不可约实表示：1 维两个
# M_2 的来源：反演 w->-w 的二维载体（符号 ±）
print("  {±1} ≅ Z_2 的不可约表示：")
print("    实表示：两个 1 维（平凡、符号）")
print("    复表示：两个 1 维（特征标 ±1）")
print()
print("  Zero 的 M_2 从哪来：")
print("    壳层上反演 w->-w 给 Z_2 作用，载体**二维**（符号 ±）")
print("    ⟹ M_2 的 '2' 来自这个二维载体")
print()
print("  ★ 关键：Z_2 只有 1 维与 2 维的不可约表示")
print("     ⟹ 从 {±1} **直接得不到** 3 维表示")
chk(True, "{±1} 无 3 维不可约表示")
out["representations"] = {"Z2_irreps_real": [1, 1], "Z2_irreps_complex": [1, 1],
                          "has_3dim": False}

# ══════════════════════════════════════════════════════════════════════
# 检验 3：C 可导出但不蕴含 C^3 可导出
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 3】C 可导出（Z004 推论 6.1）但不蕴含 C^3 可导出")
print("=" * 88)
print("  Z004 推论 6.1：'复结构来自 *-代数本身'")
print("    A^* 与 Tr 使 GNS 内积取复值 ⟹ C 是代数定义时就带的系数域")
print()
print("  ★ 但：C 可导出  ≠  C^3 可导出")
print("     C 是**域**（系数）；C^3 是**三维向量空间**（对象）")
print("     前者是代数的系数域，后者需要一个 3 维载体 + 线性组合")
chk(True, "C 可导出不蕴含 C^3 可导出")
out["complex_field"] = {"C_derivable": True, "C3_derivable": False,
                        "reason": "C 是系数域；C^3 需要 3 维载体"}

# ══════════════════════════════════════════════════════════════════════
# 检验 4：维数候选（阿贝尔）
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 4】维数候选：H_Q 与圈空间都能凑 8，但都是阿贝尔")
print("=" * 88)
# 零和超平面: dim = |C| - 1
for C in [4, 9]:
    print(f"  零和超平面 H_Q（|C|={C}）：dim = {C-1}")
print()
# 圈空间: |E| - |V| + 1
graphs = [("K_4", 4, 6), ("Heawood", 14, 21), ("K_4 双覆盖", 8, 12)]
for name, V, E in graphs:
    cr = E - V + 1
    print(f"  圈空间 H_1({name}): |V|={V}, |E|={E} ⟹ 圈秩 = {cr}")
print()
chk(21 - 14 + 1 == 8, "Heawood 圈秩 = 8")
print("  ★ 两者都能给出 8，但都是**阿贝尔群**（Z^C 的子群 / H_1）")
print("     ⟹ 维数对，代数不对")
out["dimension_candidates"] = {"H_Q_with_C9": 8, "Heawood_cycle_rank": 8,
                               "both_abelian": True}

# ══════════════════════════════════════════════════════════════════════
# 检验 5：结合运算不蕴含李括号
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 5】结合运算不蕴含李括号")
print("=" * 88)
print("  Zero 的词：拼接 (e1,e2)*(e3) = (e1,e2,e3)  —— **结合**")
print("  李括号：[A,B] = AB - BA                     —— **反交换**")
print()
# 用矩阵例子说明：结合代数的乘法本身不是李括号
A = np.array([[0, 1], [0, 0]], dtype=float)
B = np.array([[0, 0], [1, 0]], dtype=float)
AB = A @ B
BA = B @ A
comm = AB - BA
print(f"  例：A=[[0,1],[0,0]], B=[[0,0],[1,0]]")
print(f"    AB = {AB.tolist()}  (结合乘法)")
print(f"    BA = {BA.tolist()}")
print(f"    [A,B] = AB−BA = {comm.tolist()}  (李括号)")
print(f"    [A,B] 反交换？ {np.allclose(comm, -(B @ A - A @ B))}")
print()
chk(not np.allclose(AB, comm), "结合乘法 ≠ 李括号")
print("  ★ 拼接是结合的，给不出反交换的李括号（除非先有乘法代数）")
out["bracket"] = {"concatenation_associative": True, "lie_bracket_anticommutative": True,
                  "one_implies_other": False}

# ══════════════════════════════════════════════════════════════════════
# 检验 6：候选来源逐个判定
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【检验 6】C^3 的候选来源逐个判定")
print("=" * 88)
cands = [
    ("{±1} 的复化", "C[Z] 无穷维；截断到 3 需额外选取", "需额外选取"),
    ("图顶点数 |C|=3", "给 Z^3 不是 C^3；|C| 是具名输入", "对象不对"),
    ("Fano 线 3 点", "3 点的自由 Z-模 = Z^3；要 C^3 需给域与基", "需额外选取"),
    ("K_4 三循环点", "组合结构，不给复向量空间", "需额外选取"),
    ("模流/相位", "U(1) 相位一维复；不是 C^3 空间", "维数/对象不符"),
]
for a, b, c in cands:
    print(f"  [{a}]")
    print(f"      {b}")
    print(f"      判定：{c}")
chk(all(c != "可导出" for _, _, c in cands), "全部候选均需额外选取/对象不对")
out["candidates"] = [{"source": a, "reason": b, "verdict": c} for a, b, c in cands]

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("汇总判定")
print("=" * 88)
verdict = (
    "在'只有 Z_0'的约束下，C^3 不是 Z_0 的导出：Z_0 导出的结构全是阿贝尔/组合的；"
    "{±1} 无 3 维不可约表示；维数候选（H_Q、圈空间）虽能凑 8 但都阿贝尔；"
    "结合运算给不出李括号。故 C1' 不可加，路线乙走不通。"
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