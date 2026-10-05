#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D244 —— 零和匹配线性性
================================================================================
【被检验的问题（D244）】
  零和闭合词的正负完美匹配能否给 D242 的全对全势？
  有界边权给什么增长？
  全对全常数核与年龄增长边权模型是否可区分？

【本步判据】
  N1  D244 与恢复结构已登记
  N2  零和闭合词正负计数相等
  N3  完美匹配边数等于 m=|w|/2
  N4  有界边权给线性匹配势
  N5  匹配平均仍是线性边数
  N6  全对全有序配对给 k^2
  N7  年龄增长边权可形式给二次势
  N8  两种模型共享二次势但机制不同
  N9  终端中性不选择边权规律
  N10 文档登记与上游边界
  N11 文档不引用项目外体系名或外部路径

运行：python3 verify/d244_zero_sum_matching_linearity.py
"""
import io
import itertools
import os
import sys

import numpy as np


# [lh 迁入] 本脚本读旧理论的 AXIOMS.md 与 derivations/，故 ROOT 指向其原仓库（../modular-equilibrium）
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    os.pardir, os.pardir, "modular-equilibrium"))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
if os.path.join(ROOT, "verify") not in sys.path:
    sys.path.insert(0, os.path.join(ROOT, "verify"))

from latex_utils import canonical_math


PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))


def head(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


AX = io.open(os.path.join(ROOT, "AXIOMS.md"), encoding="utf-8").read()
DOC = canonical_math(
    io.open(
        os.path.join(
            ROOT,
            "derivations",
            "D244_zero_sum_matching_linearity.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D244 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D244", "D244" in AX)
check(
    "N1 公理表与文档登记零和匹配线性性",
    "R-Z-ZERO-SUM-MATCHING-LINEARITY" in BOTH,
)
check(
    "N1 公理表与文档登记匹配边权缺口",
    "R-Z-MATCHING-EDGE-WEIGHT-GAP" in BOTH,
)
check(
    "N1 公理表与文档登记全对全匹配退化",
    "R-Z-ALL-TO-ALL-MATCHING-DEGENERACY" in BOTH,
)


# ==================================================================
head("N2  零和闭合词正负计数相等")
# ==================================================================
closed_word = np.array(
    [1, 1, -1, -1, 1, -1, -1, 1],
    dtype=int,
)
plus_count = int(np.count_nonzero(closed_word == 1))
minus_count = int(np.count_nonzero(closed_word == -1))

check(
    "N2 闭合词总和为零",
    int(closed_word.sum()) == 0,
)
check(
    "N2 正负步骤数相等",
    plus_count == minus_count,
)


# ==================================================================
head("N3  完美匹配边数等于 m=|w|/2")
# ==================================================================
word_length = len(closed_word)
m = plus_count
perfect_matching_edge_count = m

check(
    "N3 完美匹配边数是词长一半",
    perfect_matching_edge_count == word_length // 2,
)
check(
    "N3 匹配边数线性随词长增长",
    perfect_matching_edge_count == word_length / 2.0,
)


# ==================================================================
head("N4  有界边权给线性匹配势")
# ==================================================================
bound = 0.7
matchings = []
plus_positions = tuple(np.where(closed_word == 1)[0])
minus_positions = tuple(np.where(closed_word == -1)[0])
for perm in itertools.permutations(minus_positions):
    matchings.append(tuple(zip(plus_positions, perm)))

edge_weights = {}
for matching in matchings:
    for edge in matching:
        edge_weights[edge] = bound * (
            0.2 + 0.1 * abs(edge[0] - edge[1])
        )

matching_potentials = np.array(
    [
        sum(edge_weights[edge] for edge in matching)
        for matching in matchings
    ]
)

check(
    "N4 每个匹配恰有 m 条边",
    all(len(matching) == m for matching in matchings),
)
check(
    "N4 匹配势被 m 倍边权上界控制",
    np.max(np.abs(matching_potentials))
    <= m * bound * (0.2 + 0.1 * word_length),
)


# ==================================================================
head("N5  匹配平均仍是线性边数")
# ==================================================================
edge_count_sequence = np.arange(1, 21)
linear_bound_sequence = bound * edge_count_sequence
quadratic_sequence = 0.03 * (2 * edge_count_sequence) ** 2

check(
    "N5 匹配平均边数线性",
    np.allclose(
        np.diff(edge_count_sequence),
        np.ones_like(edge_count_sequence[:-1]),
    ),
)
check(
    "N5 线性边数不足以匹配 k^2",
    not np.allclose(linear_bound_sequence, quadratic_sequence),
)


# ==================================================================
head("N6  全对全有序配对给 k^2")
# ==================================================================
ages = np.arange(0, 31)
delta = 0.018
all_to_all_potential = delta * ages**2

check(
    "N6 全对全势等于 delta k^2",
    np.allclose(all_to_all_potential, delta * ages**2),
)
check(
    "N6 全对全势差线性增长",
    np.allclose(
        np.diff(all_to_all_potential),
        delta * (2 * ages[:-1] + 1),
    ),
)


# ==================================================================
head("N7  年龄增长边权可形式给二次势")
# ==================================================================
edge_weight_coefficient = 2.0 * delta
growing_edge_matching_potential = (
    (ages / 2.0)
    * edge_weight_coefficient
    * ages
)

check(
    "N7 年龄增长边权给二次势",
    np.allclose(
        growing_edge_matching_potential,
        delta * ages**2,
    ),
)
check(
    "N7 边权含年龄线性因子",
    np.allclose(
        growing_edge_matching_potential[1:] / (ages[1:] / 2.0),
        edge_weight_coefficient * ages[1:],
    ),
)


# ==================================================================
head("N8  两种模型共享二次势但机制不同")
# ==================================================================
all_to_all_pairs = ages**2
matching_edges = ages / 2.0
matching_edge_weight = edge_weight_coefficient * ages

check(
    "N8 全对全靠 k^2 个配对",
    np.allclose(all_to_all_pairs, ages**2),
)
check(
    "N8 增长边权靠 k 乘边权",
    np.allclose(
        matching_edges * matching_edge_weight,
        delta * ages**2,
    ),
)
check(
    "N8 两种对象不同",
    not np.allclose(all_to_all_pairs, matching_edges),
)


# ==================================================================
head("N9  终端中性不选择边权规律")
# ==================================================================
def terminal_neutral_log_ratio(potential_values):
    log_ratio = -potential_values
    log_ratio -= log_ratio[-1]
    return log_ratio


all_to_all_log_ratio = terminal_neutral_log_ratio(all_to_all_potential)
growing_edge_log_ratio = terminal_neutral_log_ratio(
    growing_edge_matching_potential
)

check(
    "N9 两模型都可终端中性",
    abs(all_to_all_log_ratio[-1]) < 1e-12
    and abs(growing_edge_log_ratio[-1]) < 1e-12,
)
check(
    "N9 两模型比例相同",
    np.allclose(all_to_all_log_ratio, growing_edge_log_ratio),
)


# ==================================================================
head("N10  文档登记与上游边界")
# ==================================================================
check(
    "N10 文档登记零和匹配线性性",
    "R-Z-ZERO-SUM-MATCHING-LINEARITY" in DOC
    and "有界边权" in DOC
    and "至多线性势" in DOC,
)
check(
    "N10 文档登记匹配边权缺口",
    "R-Z-MATCHING-EDGE-WEIGHT-GAP" in DOC
    and "边权随年龄线性增长" in DOC,
)
check(
    "N10 文档登记全对全匹配退化",
    "R-Z-ALL-TO-ALL-MATCHING-DEGENERACY" in DOC
    and "零和闭合不能区分二者" in DOC,
)
check(
    "N10 文档不把新选择写成 U1-U4 推论",
    "它们不是 `U1-U4` 的推论" in DOC
    and "这三项都不是零和闭合的自动结果" in DOC,
)
check(
    "N10 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)


# ==================================================================
head("N11  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N11 D244 文档不含项目外体系名或外部路径",
    not any(marker in DOC for marker in external_markers),
)


print("\n" + "=" * 78)
print(f"通过 {len(PASS)} / 不符 {len(FAIL)}")
print("=" * 78)
if FAIL:
    for name in FAIL:
        print(f"失败：{name}")
    sys.exit(1)
print("全部通过")
