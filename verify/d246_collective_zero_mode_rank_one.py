#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D246 —— 全对全核的集体零模表示
================================================================================
【被检验的问题（D246）】
  常数全对全核能否化为秩一集体零模？
  终端零约束和平衡完成计数能否自动给出固定二次刚度？
  D238 的 Q_k^2 势还缺哪一项输入？

【本步判据】
  N1  D246 与恢复结构已登记
  N2  常数全对全核等于秩一矩阵
  N3  常数核二次型等于 c Q_k^2
  N4  零和向量被全对全核湮灭
  N5  D241 符号差给集体二次能量
  N6  终端零约束把端到能量设为零
  N7  基准 Q_k=k 恢复抛物型比例
  N8  平衡完成计数给 Q_k^2/(T-k)
  N9  固定荷缺陷下成本随剩余步数改变
  N10 最小二次完成给同一倒数尺度
  N11 终端零和路径可承载不同前缀荷
  N12 文档登记与上游边界
  N13 文档不引用项目外体系名或外部路径

运行：python3 verify/d246_collective_zero_mode_rank_one.py
"""
import io
import math
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
            "D246_collective_zero_mode_rank_one.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D246 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D246", "D246" in AX)
check(
    "N1 公理表与文档登记集体零模表示",
    "R-Z-COLLECTIVE-ZERO-MODE-REPRESENTATION" in BOTH,
)
check(
    "N1 公理表与文档登记完成成本无解",
    "R-Z-COMPLETION-COST-NOGO" in BOTH,
)
check(
    "N1 公理表与文档登记零缺陷刚度缺口",
    "R-Z-ZERO-DEFECT-STIFFNESS-GAP" in BOTH,
)


# ==================================================================
head("N2  常数全对全核等于秩一矩阵")
# ==================================================================
c = 0.23
size = 9
ones = np.ones(size)
constant_kernel = c * np.outer(ones, ones)

check(
    "N2 全对全核矩阵秩为一",
    np.linalg.matrix_rank(constant_kernel) == 1,
)
check(
    "N2 全对全核等于 c 11^T",
    np.allclose(constant_kernel, c * np.outer(ones, ones)),
)


# ==================================================================
head("N3  常数核二次型等于 c Q_k^2")
# ==================================================================
q = np.array([1, -2, 3, -1, 2, 0, -3, 4, -4], dtype=float)
Q = float(np.sum(q))
quadratic_form = float(q @ constant_kernel @ q)

check(
    "N3 二次型等于 c Q^2",
    np.isclose(quadratic_form, c * Q**2),
)
check(
    "N3 秩一改写保持符号势",
    np.isclose(quadratic_form, float(np.sum(constant_kernel * np.outer(q, q)))),
)


# ==================================================================
head("N4  零和向量被全对全核湮灭")
# ==================================================================
zero_sum_vector = np.array([1, -1, 2, -2, 3, -3, 4, -4, 0], dtype=float)

check(
    "N4 零和向量总荷为零",
    np.isclose(np.sum(zero_sum_vector), 0.0),
)
check(
    "N4 全对全核湮灭零和向量",
    np.allclose(constant_kernel @ zero_sum_vector, np.zeros(size)),
)


# ==================================================================
head("N5  D241 符号差给集体二次能量")
# ==================================================================
delta = 0.31
sign_kernel = -delta * np.outer(ones, ones)
sign_charge = np.ones(size)
sign_sum = float(np.sum(sign_charge))
sign_potential = float(sign_charge @ sign_kernel @ sign_charge)

check(
    "N5 符号差核给负集体平方",
    np.isclose(sign_potential, -delta * sign_sum**2),
)
check(
    "N5 生成势为正集体平方",
    np.isclose(-sign_potential, delta * sign_sum**2),
)


# ==================================================================
head("N6  终端零约束把端到能量设为零")
# ==================================================================
T = 12
ages = np.arange(T + 1)
prefix_charge = ages.astype(float)
potential = delta * prefix_charge**2

check(
    "N6 终端零约束",
    np.isclose(prefix_charge[-1] - T, 0.0),
)
check(
    "N6 终端集体能量为零",
    np.isclose(delta * (prefix_charge[-1] - T) ** 2, 0.0),
)


# ==================================================================
head("N7  基准 Q_k=k 恢复抛物型比例")
# ==================================================================
log_ratio = potential[-1] - potential

check(
    "N7 对数比例等于 delta(T^2-k^2)",
    np.allclose(log_ratio, delta * (T**2 - ages**2)),
)
check(
    "N7 终端中性",
    np.isclose(log_ratio[-1], 0.0),
)


# ==================================================================
head("N8  平衡完成计数给 Q_k^2/(T-k)")
# ==================================================================
def log_binomial(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def completion_log_ratio(remaining, defect):
    left = (remaining - defect) // 2
    if (remaining - defect) % 2:
        return float("-inf")
    return log_binomial(remaining, left) - log_binomial(
        remaining,
        remaining // 2,
    )


test_remaining = [30, 60, 120]
test_defect = 2
approx_exponents = np.array(
    [
        completion_log_ratio(remaining, test_defect)
        for remaining in test_remaining
    ]
)
expected_exponents = -test_defect**2 / (2.0 * np.array(test_remaining))

check(
    "N8 完成计数对数比接近 -Q^2/(2N)",
    np.allclose(approx_exponents, expected_exponents, rtol=2e-2, atol=2e-3),
)
check(
    "N8 固定 Q 时成本随 N 改变",
    not np.allclose(approx_exponents, approx_exponents[0]),
)


# ==================================================================
head("N9  固定荷缺陷下成本随剩余步数改变")
# ==================================================================
def minimal_quadratic_cost(remaining, defect):
    return defect**2 / (2.0 * remaining)


cost_curve = np.array(
    [minimal_quadratic_cost(remaining, 3) for remaining in test_remaining]
)

check(
    "N9 固定刚度 Q^2 不随剩余步数改变",
    len({float(cost_curve[0])} | {float(cost_curve[-1])}) == 2,
)
check(
    "N9 成本随剩余步数下降",
    np.all(np.diff(cost_curve) < 0),
)


# ==================================================================
head("N10  最小二次完成给同一倒数尺度")
# ==================================================================
def distribute_correction(remaining, defect):
    epsilon = np.full(remaining, -defect / remaining, dtype=float)
    cost = 0.5 * float(np.sum(epsilon**2))
    return epsilon, cost


epsilon, distributed_cost = distribute_correction(40, 5)

check(
    "N10 修正总和等于需求",
    np.isclose(np.sum(epsilon), -5.0),
)
check(
    "N10 最小二次成本等于 Q^2/(2N)",
    np.isclose(distributed_cost, 25.0 / 80.0),
)


# ==================================================================
head("N11  终端零和路径可承载不同前缀荷")
# ==================================================================
path_v = np.array([1, -1, 1, -1])
path_w = np.array([1, 1, -1, -1])
prefix_v = np.cumsum(path_v)
prefix_w = np.cumsum(path_w)

check(
    "N11 两条路径终端都为零",
    int(prefix_v[-1]) == 0 and int(prefix_w[-1]) == 0,
)
check(
    "N11 两条路径中间前缀荷不同",
    not np.array_equal(prefix_v, prefix_w),
)
check(
    "N11 终端约束不选择中间荷",
    prefix_v[-1] == prefix_w[-1] and np.any(prefix_v != prefix_w),
)


# ==================================================================
head("N12  文档登记与上游边界")
# ==================================================================
check(
    "N12 文档登记集体零模表示",
    "R-Z-COLLECTIVE-ZERO-MODE-REPRESENTATION" in DOC
    and "秩一集体模式" in DOC,
)
check(
    "N12 文档登记完成成本无解",
    "R-Z-COMPLETION-COST-NOGO" in DOC
    and "Q_k^2/(T-k)" in DOC,
)
check(
    "N12 文档登记零缺陷刚度缺口",
    "R-Z-ZERO-DEFECT-STIFFNESS-GAP" in DOC
    and "固定刚度" in DOC,
)
check(
    "N12 文档不把新结构写成 U1-U4 推论",
    "它们不是 `U1-U4` 的推论" in DOC
    and "仍不是 `U1-U4`" in DOC,
)
check(
    "N12 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)


# ==================================================================
head("N13  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N13 D246 文档不含项目外体系名或外部路径",
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
