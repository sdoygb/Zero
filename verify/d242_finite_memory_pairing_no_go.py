#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D242 —— 有限记忆配对无解
================================================================================
【被检验的问题（D242）】
  有限范围配对核能否生成 D238 的线性危险率？
  D222 的两层历史保留属于哪一类记忆？
  什么耦合范围才能给二次生成势？

【本步判据】
  N1  D242 与恢复结构已登记
  N2  平移不变配对核给前缀增量公式
  N3  有限范围核给饱和增量
  N4  饱和危险率给指数型比例
  N5  D222 两层保留是有限范围记忆
  N6  全对全常数核给 delta(2k+1)
  N7  全对全核给 D238 抛物型比例
  N8  幂律核给对数增长
  N9  终端中性不选择耦合范围
  N10 有限范围与全对全核给出不同年龄比例
  N11 文档登记与上游边界
  N12 文档不引用项目外体系名或外部路径

运行：python3 verify/d242_finite_memory_pairing_no_go.py
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
            "D242_finite_memory_pairing_no_go.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D242 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D242", "D242" in AX)
check(
    "N1 公理表与文档登记有限范围无解",
    "R-Z-FINITE-RANGE-PAIR-NOGO" in BOTH,
)
check(
    "N1 公理表与文档登记全对全耦合",
    "R-Z-ALL-TO-ALL-AGE-COUPLING" in BOTH,
)
check(
    "N1 公理表与文档登记长程记忆缺口",
    "R-Z-LONG-RANGE-MEMORY-GAP" in BOTH,
)


# ==================================================================
head("N2  平移不变配对核给前缀增量公式")
# ==================================================================
size = 80
delta = 0.021
ages = np.arange(0, size + 1)


def potential_from_kernel(kernel, k):
    matrix = np.empty((k, k), dtype=float)
    for i in range(k):
        for j in range(k):
            matrix[i, j] = kernel(abs(i - j))
    return float(matrix.sum())


def increments_from_kernel(kernel, max_k):
    potentials = np.array(
        [potential_from_kernel(kernel, k) for k in range(max_k + 1)]
    )
    return np.diff(potentials)


def finite_kernel(range_cutoff):
    return lambda r: delta if r <= range_cutoff else 0.0


finite_increments = increments_from_kernel(finite_kernel(2), 8)
finite_formula = np.array(
    [
        delta
        + 2 * sum(delta for r in range(1, k + 1) if r <= 2)
        for k in range(8)
    ]
)

check(
    "N2 有限核增量符合距离求和公式",
    np.allclose(finite_increments, finite_formula),
)
check(
    "N2 有限核早期增量非零",
    finite_increments[0] > 0.0,
)


# ==================================================================
head("N3  有限范围核给饱和增量")
# ==================================================================
finite_long_increments = increments_from_kernel(finite_kernel(2), 40)

check(
    "N3 有限范围后增量饱和",
    np.allclose(finite_long_increments[5:], finite_long_increments[5]),
)
check(
    "N3 饱和增量不等于线性增长",
    np.ptp(finite_long_increments[5:]) < 1e-12,
)


# ==================================================================
head("N4  饱和危险率给指数型比例")
# ==================================================================
saturated_hazard = finite_long_increments[-1]
finite_log_ratio = -np.cumsum(
    np.concatenate([[0.0], finite_long_increments])
)
finite_log_ratio -= finite_log_ratio[-1]
finite_ratio = np.exp(finite_log_ratio)
quadratic_target = np.exp(delta * (size**2 - ages**2))

check(
    "N4 有限核比例终端归零",
    abs(finite_log_ratio[-1]) < 1e-12,
)
check(
    "N4 饱和危险率非线性增长",
    abs(saturated_hazard - delta * (2 * size + 1)) > 1e-2,
)
check(
    "N4 有限核比例不同于 D238",
    np.max(
        np.abs(finite_ratio - quadratic_target[: len(finite_ratio)])
    )
    > 1e-2,
)


# ==================================================================
head("N5  D222 两层保留是有限范围记忆")
# ==================================================================
top_two_memory_range = 2
top_two_kernel = finite_kernel(top_two_memory_range)
top_two_increments = increments_from_kernel(top_two_kernel, 12)

check(
    "N5 两层保留只给有限核",
    top_two_kernel(top_two_memory_range + 1) == 0.0,
)
check(
    "N5 两层保留增量饱和",
    np.allclose(top_two_increments[4:], top_two_increments[4]),
)


# ==================================================================
head("N6  全对全常数核给 delta(2k+1)")
# ==================================================================
all_to_all_kernel = lambda r: delta
all_to_all_increments = increments_from_kernel(
    all_to_all_kernel,
    size,
)
target_increments = delta * (2 * np.arange(size) + 1)

check(
    "N6 全对全核增量等于 delta(2k+1)",
    np.allclose(all_to_all_increments, target_increments),
)
check(
    "N6 全对全核增量线性增长",
    np.allclose(np.diff(all_to_all_increments), 2.0 * delta),
)


# ==================================================================
head("N7  全对全核给 D238 抛物型比例")
# ==================================================================
all_to_all_potential = np.array(
    [potential_from_kernel(all_to_all_kernel, k) for k in ages]
)
all_to_all_log_ratio = -all_to_all_potential
all_to_all_log_ratio -= all_to_all_log_ratio[-1]
all_to_all_ratio = np.exp(all_to_all_log_ratio)

check(
    "N7 全对全核恢复 D238 比例",
    np.allclose(all_to_all_ratio, quadratic_target),
)
check(
    "N7 全对全核终端中性",
    abs(all_to_all_ratio[-1] - 1.0) < 1e-12,
)


# ==================================================================
head("N8  幂律核给对数增长")
# ==================================================================
power_law_kernel = lambda r: delta / r if r > 0 else delta
power_law_increments = increments_from_kernel(
    power_law_kernel,
    80,
)
power_law_linear_fit = np.polyfit(
    np.arange(10, 80),
    power_law_increments[10:80],
    1,
)
all_to_all_linear_fit = np.polyfit(
    np.arange(10, 80),
    all_to_all_increments[10:80],
    1,
)

check(
    "N8 幂律核增量仍增长",
    power_law_increments[-1] > power_law_increments[10],
)
check(
    "N8 幂律核不是精确线性增长",
    abs(power_law_linear_fit[0]) < abs(all_to_all_linear_fit[0]) / 2.0,
)


# ==================================================================
head("N9  终端中性不选择耦合范围")
# ==================================================================
finite_log_ratio_compensated = -finite_log_ratio + finite_log_ratio[-1]
all_to_all_log_ratio_compensated = (
    -all_to_all_log_ratio + all_to_all_log_ratio[-1]
)

check(
    "N9 有限核可终端补偿",
    abs(finite_log_ratio_compensated[-1]) < 1e-12,
)
check(
    "N9 全对全核可终端补偿",
    abs(all_to_all_log_ratio_compensated[-1]) < 1e-12,
)
check(
    "N9 终端补偿不选择核范围",
    not np.allclose(
        finite_log_ratio_compensated,
        all_to_all_log_ratio_compensated[
            : len(finite_log_ratio_compensated)
        ],
    ),
)


# ==================================================================
head("N10  有限范围与全对全核给出不同年龄比例")
# ==================================================================
finite_ratio_long = np.exp(
    finite_log_ratio - finite_log_ratio[-1]
)
all_to_all_ratio_long = np.exp(
    all_to_all_log_ratio - all_to_all_log_ratio[-1]
)

check(
    "N10 两核比例显著不同",
    np.max(
        np.abs(
            finite_ratio_long
            - all_to_all_ratio_long[: len(finite_ratio_long)]
        )
    )
    > 1e-2,
)
check(
    "N10 全对全核给二次对数曲率",
    np.ptp(np.diff(all_to_all_log_ratio, 2)) < 1e-12,
)


# ==================================================================
head("N11  文档登记与上游边界")
# ==================================================================
check(
    "N11 文档登记有限记忆无解",
    "R-Z-FINITE-RANGE-PAIR-NOGO" in DOC
    and "有限记忆" in DOC
    and "饱和危险率" in DOC,
)
check(
    "N11 文档登记全对全耦合",
    "R-Z-ALL-TO-ALL-AGE-COUPLING" in DOC
    and "全对全常数配对核" in DOC,
)
check(
    "N11 文档登记长程记忆缺口",
    "R-Z-LONG-RANGE-MEMORY-GAP" in DOC
    and "长程全对全耦合来源仍缺" in DOC,
)
check(
    "N11 文档不把新耦合写成 U1-U4 推论",
    "它们不是 `U1-U4` 的推论" in DOC
    and "D238 需要长程记忆或全对全耦合" in DOC,
)
check(
    "N11 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)


# ==================================================================
head("N12  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N12 D242 文档不含项目外体系名或外部路径",
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
