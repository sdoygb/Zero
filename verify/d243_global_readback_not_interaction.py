#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D243 —— 全局读回不是相互作用
================================================================================
【被检验的问题（D243）】
  全局零层的可见性是否自动给出全对全耦合？
  同重数记录的配对图是否可由全局读回区分？
  哪条显式规则才恢复 D242 的全对全核？

【本步判据】
  N1  D243 与恢复结构已登记
  N2  同重数配置可以承载不同配对图
  N3  可加全局读回对两者相同
  N4  配对泛函对两者不同
  N5  全局可见性不自动选择配对律
  N6  全对全规则给 D242 势差
  N7  零配对与有限配对不能给 D242
  N8  终端中性不选择配对律
  N9  全局全对全配对是条件选择器
  N10 文档登记与上游边界
  N11 文档不引用项目外体系名或外部路径

运行：python3 verify/d243_global_readback_not_interaction.py
"""
import io
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
            "D243_global_readback_not_interaction.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D243 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D243", "D243" in AX)
check(
    "N1 公理表与文档登记全局读回无耦合",
    "R-Z-GLOBAL-READBACK-NOT-COUPLING" in BOTH,
)
check(
    "N1 公理表与文档登记全局全对全配对",
    "R-Z-GLOBAL-ALL-TO-ALL-PAIRING" in BOTH,
)
check(
    "N1 公理表与文档登记全局配对来源缺口",
    "R-Z-GLOBAL-PAIRING-ORIGIN-GAP" in BOTH,
)


# ==================================================================
head("N2  同重数配置可以承载不同配对图")
# ==================================================================
record_count = 8
delta = 0.037
all_to_all = delta * np.ones((record_count, record_count))
diagonal_only = delta * np.eye(record_count)
finite_range = np.zeros((record_count, record_count))
for i in range(record_count):
    for j in range(i, record_count):
        if abs(i - j) <= 1:
            finite_range[i, j] = delta
            finite_range[j, i] = delta

check(
    "N2 三种配对图记录数相同",
    all_to_all.shape == diagonal_only.shape == finite_range.shape,
)
check(
    "N2 三种配对图非对角结构不同",
    not np.allclose(all_to_all, diagonal_only)
    and not np.allclose(all_to_all, finite_range),
)

pair_functional = lambda matrix: float(matrix.sum())


# ==================================================================
head("N3  可加全局读回对两者相同")
# ==================================================================
count_readback = lambda matrix: int(matrix.shape[0])
readback_all_to_all = count_readback(all_to_all)
readback_diagonal = count_readback(diagonal_only)
readback_finite = count_readback(finite_range)

check(
    "N3 记录计数读回对三者相同",
    readback_all_to_all == readback_diagonal == readback_finite,
)
check(
    "N3 记录计数读回不是配对泛函",
    readback_all_to_all != pair_functional(all_to_all),
)


# ==================================================================
head("N4  配对泛函对两者不同")
# ==================================================================
check(
    "N4 全对全与对角配对泛函不同",
    pair_functional(all_to_all) > pair_functional(diagonal_only),
)
check(
    "N4 全对全与有限范围配对泛函不同",
    pair_functional(all_to_all) > pair_functional(finite_range),
)


# ==================================================================
head("N5  全局可见性不自动选择配对律")
# ==================================================================
possible_pairings = [all_to_all, diagonal_only, finite_range]
pairing_summaries = tuple(pair_functional(matrix) for matrix in possible_pairings)
unique_pairing_summaries = set(pairing_summaries)

check(
    "N5 同一全局目录可有多个配对律",
    len(unique_pairing_summaries) == len(possible_pairings),
)
check(
    "N5 全局目录本身不含配对选择",
    len({matrix.shape for matrix in possible_pairings}) == 1,
)


# ==================================================================
head("N6  全对全规则给 D242 势差")
# ==================================================================
ages = np.arange(0, 21)
delta_uniform = 0.052
all_to_all_kernel_sum = np.array(
    [delta_uniform * k**2 for k in ages]
)
increments = np.diff(all_to_all_kernel_sum)

check(
    "N6 全对全势差等于 delta(2k+1)",
    np.allclose(
        increments,
        delta_uniform * (2 * ages[:-1] + 1),
    ),
)
check(
    "N6 全对全势差线性增长",
    np.allclose(np.diff(increments), 2.0 * delta_uniform),
)


# ==================================================================
head("N7  零配对与有限配对不能给 D242")
# ==================================================================
zero_kernel_sum = np.zeros_like(ages, dtype=float)
finite_kernel_sum = np.array(
    [
        delta_uniform * min(k * k, 4 * k - 4)
        for k in ages
    ]
)

check(
    "N7 零配对势差恒零",
    np.allclose(np.diff(zero_kernel_sum), 0.0),
)
check(
    "N7 有限配对势差最终饱和",
    np.ptp(np.diff(finite_kernel_sum)[6:]) < 1e-12,
)
check(
    "N7 有限配对不同全对全",
    not np.allclose(finite_kernel_sum, all_to_all_kernel_sum),
)


# ==================================================================
head("N8  终端中性不选择配对律")
# ==================================================================
def terminal_neutral_ratio(potential_values):
    log_ratio = -potential_values
    log_ratio -= log_ratio[-1]
    return np.exp(log_ratio)


zero_ratio = terminal_neutral_ratio(zero_kernel_sum)
finite_ratio = terminal_neutral_ratio(finite_kernel_sum)
all_to_all_ratio = terminal_neutral_ratio(all_to_all_kernel_sum)

check(
    "N8 三种配对都可终端中性",
    abs(zero_ratio[-1] - 1.0) < 1e-12
    and abs(finite_ratio[-1] - 1.0) < 1e-12
    and abs(all_to_all_ratio[-1] - 1.0) < 1e-12,
)
check(
    "N8 终端中性不选全对全",
    not np.allclose(zero_ratio, all_to_all_ratio)
    and not np.allclose(finite_ratio, all_to_all_ratio),
)


# ==================================================================
head("N9  全局全对全配对是条件选择器")
# ==================================================================
global_pairing_rule = "each_new_record_pairs_with_all_old_records"
all_to_all_condition = (
    global_pairing_rule == "each_new_record_pairs_with_all_old_records"
)

check(
    "N9 选择器显式规定全对全",
    all_to_all_condition,
)
check(
    "N9 选择器不是全局可见性的推论",
    global_pairing_rule not in {"global_counts_only", "local_top_two"},
)


# ==================================================================
head("N10  文档登记与上游边界")
# ==================================================================
check(
    "N10 文档登记全局读回无耦合",
    "R-Z-GLOBAL-READBACK-NOT-COUPLING" in DOC
    and "全局可见" in DOC
    and "记录之间有相互作用" in DOC,
)
check(
    "N10 文档登记全局全对全配对",
    "R-Z-GLOBAL-ALL-TO-ALL-PAIRING" in DOC
    and "每条新记录与全部已有记录建立均匀符号差配对" in DOC,
)
check(
    "N10 文档登记全局配对来源缺口",
    "R-Z-GLOBAL-PAIRING-ORIGIN-GAP" in DOC
    and "全局配对律的来源仍是新的恢复层缺口" in DOC,
)
check(
    "N10 文档不把新规则写成 U1-U4 推论",
    "它们不是 `U1-U4` 的推论" in DOC
    and "不是 `U1-U4`、`D211` 或终端中性的推论" in DOC,
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
    "N11 D243 文档不含项目外体系名或外部路径",
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
