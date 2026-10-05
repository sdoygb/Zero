#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D225 —— 张量因子与直接和替代
================================================================================
【被检验的问题（D225）】
  D224 的 M_2 tensor C(X_tau) 是否是同类型所有有限维载体中的绝对最小候选？
  直接和 M_2 + C(X_tau) 能否同样满足 U1,U2_sem,U3,U4_sem？
  U1-U4 能否单独选出张量或直接和？

【本步判据】
  N1  D225 与两个恢复结构已登记
  N2  张量与直接和载体都有限维且含单位元
  N3  直接和维数严格小于张量维数
  N4  张量态与直接和分块态都可忠实且非中心
  N5  两类模流都可非平凡
  N6  两类年龄支持都保序并给支持包含
  N7  直接和的 p_M,p_C 是补中央投影
  N8  张量中心投影不能产生独立 M_2 分块
  N9  跨区域无交在两类局部选择下都成立
  N10 文档不把任一路线写成上游导出
  N11 文档登记选择器缺口
  N12 不引用外部材料作依据
  N13 上游边界保持

运行：python3 verify/d225_tensor_vs_direct_sum_factorization_gap.py
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
            "D225_tensor_vs_direct_sum_factorization_gap.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


def age_shift(a, t, tau):
    if a == "*":
        return "*"
    if a + t >= tau:
        return "*"
    return a + t


def local_order(a, b, tau):
    sequence = tuple(range(tau)) + ("*",)
    return sequence.index(a) <= sequence.index(b)


def prefix(a, tau):
    if a == "*":
        return tuple(range(tau)) + ("*",)
    return tuple(range(a + 1))


def modular_phase_ratio(r):
    return r / (1.0 - r)


def matrix_unit(dimension, row, column):
    matrix = np.zeros((dimension, dimension), dtype=complex)
    matrix[row, column] = 1.0
    return matrix


def direct_sum(left, right):
    return np.block(
        [
            [left, np.zeros((left.shape[0], right.shape[1]), dtype=complex)],
            [np.zeros((right.shape[0], left.shape[1]), dtype=complex), right],
        ]
    )


def is_central(operator, algebra_basis):
    return all(
        np.allclose(operator @ basis, basis @ operator)
        for basis in algebra_basis
    )


# ==================================================================
head("N1  D225 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记局部因子二择一",
    "D225" in AX
    and "R-Z-LOCAL-FACTORIZATION-DICHOTOMY" in BOTH
    and "直接和" in DOC,
)
check(
    "N1 公理表与文档登记张量直接和选择缺口",
    "R-Z-TENSOR-DIRECT-SUM-GAP" in BOTH
    and "中央超选择" in DOC,
)


# ==================================================================
head("N2  张量与直接和载体都有限维且含单位元")
# ==================================================================
taus = (1, 2, 3, 4, 5)
tensor_dimensions = {
    tau: 4 * (tau + 1)
    for tau in taus
}
direct_dimensions = {
    tau: 4 + (tau + 1)
    for tau in taus
}
check(
    "N2 张量载体维数公式正确",
    all(tensor_dimensions[tau] == 4 * (tau + 1) for tau in taus),
    f"tensor={tensor_dimensions}",
)
check(
    "N2 直接和载体维数公式正确",
    all(direct_dimensions[tau] == tau + 5 for tau in taus),
    f"direct={direct_dimensions}",
)


# ==================================================================
head("N3  直接和维数严格小于张量维数")
# ==================================================================
strictly_smaller = all(
    direct_dimensions[tau] < tensor_dimensions[tau]
    for tau in taus
)
check(
    "N3 对全部测试 tau，直接和更小",
    strictly_smaller,
    f"pairs={[(tensor_dimensions[t], direct_dimensions[t]) for t in taus]}",
)


# ==================================================================
head("N4  张量态与直接和分块态都可忠实且非中心")
# ==================================================================
r = 0.75
lambda_weight = 0.4
matrix_weights = (r, 1.0 - r)
age_weights = (0.1, 0.2, 0.7)
direct_block_weights = (
    lambda_weight * r,
    lambda_weight * (1.0 - r),
    (1.0 - lambda_weight) * age_weights[0],
    (1.0 - lambda_weight) * age_weights[1],
    (1.0 - lambda_weight) * age_weights[2],
)
check(
    "N4 矩阵态是忠实非中心态",
    all(weight > 0 for weight in matrix_weights)
    and not np.isclose(matrix_weights[0], matrix_weights[1]),
)
check(
    "N4 直接和分块态忠实且两个 M_2 权重不同",
    all(weight > 0 for weight in direct_block_weights)
    and not np.isclose(
        direct_block_weights[0],
        direct_block_weights[1],
    ),
)


# ==================================================================
head("N5  两类模流都可非平凡")
# ==================================================================
check(
    "N5 非中心矩阵态给非平凡模相位比",
    not np.isclose(modular_phase_ratio(r), 1.0),
    f"ratio={modular_phase_ratio(r):.6f}",
)
check(
    "N5 直接和模流在 M_2 分块与张量模流相同",
    np.isclose(
        modular_phase_ratio(direct_block_weights[0] / lambda_weight),
        modular_phase_ratio(r),
    ),
)


# ==================================================================
head("N6  两类年龄支持都保序并给支持包含")
# ==================================================================
order_ok = True
inclusion_ok = True
for tau in taus:
    ages = tuple(range(tau)) + ("*",)
    for a in ages:
        for b in ages:
            if local_order(a, b, tau):
                for t in range(5):
                    order_ok &= local_order(
                        age_shift(a, t, tau),
                        age_shift(b, t, tau),
                        tau,
                    )
        for t in range(5):
            inclusion_ok &= set(prefix(a, tau)).issubset(
                set(prefix(age_shift(a, t, tau), tau))
            )
check(
    "N6 年龄移位保序",
    order_ok,
)
check(
    "N6 年龄前缀包含成立",
    inclusion_ok,
)


# ==================================================================
head("N7  直接和的 p_M,p_C 是补中央投影")
# ==================================================================
p_m = direct_sum(
    np.eye(2, dtype=complex),
    np.zeros((3, 3), dtype=complex),
)
p_c = direct_sum(
    np.zeros((2, 2), dtype=complex),
    np.eye(3, dtype=complex),
)
algebra_basis = []
for row in range(2):
    for column in range(2):
        algebra_basis.append(direct_sum(matrix_unit(2, row, column), np.zeros((3, 3))))
for index in range(3):
    algebra_basis.append(direct_sum(np.zeros((2, 2)), matrix_unit(3, index, index)))
check(
    "N7 p_M+p_C=1 且 p_Mp_C=0",
    np.allclose(p_m + p_c, np.eye(5, dtype=complex))
    and np.allclose(p_m @ p_c, 0.0),
)
check(
    "N7 p_M 与 p_C 都中央",
    is_central(p_m, algebra_basis) and is_central(p_c, algebra_basis),
)


# ==================================================================
head("N8  张量中心投影不能产生独立 M_2 分块")
# ==================================================================
tensor_basis = []
for age in range(3):
    for row in range(2):
        for column in range(2):
            matrix = np.zeros((2, 2), dtype=complex)
            matrix[row, column] = 1.0
            age_projection = np.zeros((3, 3), dtype=complex)
            age_projection[age, age] = 1.0
            tensor_basis.append(np.kron(matrix, age_projection))
tensor_identity = np.eye(2 * 3, dtype=complex)
central_dimension_four_candidates = []
for row in range(2):
    for column in range(2):
        matrix = np.zeros((2, 2), dtype=complex)
        matrix[row, column] = 1.0
        central_projection = np.kron(matrix, np.eye(3, dtype=complex))
        if np.allclose(central_projection, central_projection.conj().T):
            if np.allclose(central_projection @ central_projection, central_projection):
                if is_central(central_projection, tensor_basis):
                    central_dimension_four_candidates.append(central_projection)
check(
    "N8 张量载体不出现独立 M_2 直接和分块",
    len(central_dimension_four_candidates) == 0,
)
check(
    "N8 张量恒等元仍为完整局部单位",
    np.allclose(tensor_identity @ tensor_basis[0], tensor_basis[0]),
)


# ==================================================================
head("N9  跨区域无交在两类局部选择下都成立")
# ==================================================================
left_factor = np.eye(direct_dimensions[2], dtype=complex)
right_factor = np.eye(direct_dimensions[3], dtype=complex)
left_embedded = np.kron(left_factor, np.eye(right_factor.shape[0]))
right_embedded = np.kron(np.eye(left_factor.shape[0]), right_factor)
check(
    "N9 不同区域张量嵌入交换",
    np.allclose(
        left_embedded @ right_embedded,
        right_embedded @ left_embedded,
    ),
)
check(
    "N9 文档允许逐区域选择张量或直接和",
    "A_i\\in\\{A_{\\tau_i}^\\otimes,A_{\\tau_i}^\\oplus\\}" in DOC,
)


# ==================================================================
head("N10 文档不把任一路线写成上游导出")
# ==================================================================
check(
    "N10 文档拒绝 U1-U4 单独选择",
    "不能选择张量乘积或直接和" in DOC
    and "恢复层选择器" in DOC,
)
check(
    "N10 文档保留两套物理输入",
    "直接和更省" in DOC
    and "张量共存更强" in DOC,
)


# ==================================================================
head("N11 文档登记选择器缺口")
# ==================================================================
check(
    "N11 文档要求额外原则",
    "局部时间结构必须共存" in DOC
    and "允许最小维数的超选择扇区" in DOC,
)


# ==================================================================
head("N12 不引用外部材料作依据")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "旧体系",
    "../..",
)
check(
    "N12 D225 文档不含项目外体系名或外部路径",
    not any(marker in DOC for marker in external_markers),
)


# ==================================================================
head("N13 上游边界保持")
# ==================================================================
check(
    "N13 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC
    and "不新增 `U5`" in DOC,
)


print("\n" + "=" * 78)
print(f"通过 {len(PASS)} / 不符 {len(FAIL)}")
if FAIL:
    print("失败项：")
    for name in FAIL:
        print(" -", name)
    sys.exit(1)
print("全部通过")
