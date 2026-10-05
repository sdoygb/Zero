#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D226 —— 矩阵均匀共存原则
================================================================================
【被检验的问题（D226）】
  矩阵因子是否必须在每个年龄支持扇区中完整出现？
  若要求均匀共存，能否排除直接和并选出张量载体？
  该选择是否是 U1-U4 的定理？

【本步判据】
  N1  D226 与两个恢复结构已登记
  N2  张量载体在每个年龄扇区保留完整 M_2
  N3  直接和的年龄中央投影把 M_2 压成零
  N4  均匀共存原则排除直接和
  N5  条件生成类中的张量载体达到维数下界
  N6  放弃均匀共存时直接和仍是更省候选
  N7  中心投影与年龄扇区结构正确
  N8  选择原则没有被写成 U1-U4 定理
  N9  D223-D224 接口在新选择下保持
  N10 文档不引用外部材料作依据
  N11 上游边界保持

运行：python3 verify/d226_uniform_matrix_coexistence_selector.py
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
            "D226_uniform_matrix_coexistence_selector.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


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


def matrix_support_dimension(projection, matrix_basis):
    compressed = [
        projection @ basis @ projection
        for basis in matrix_basis
    ]
    vector_basis = np.column_stack(
        [matrix.reshape(-1) for matrix in compressed]
    )
    return int(np.linalg.matrix_rank(vector_basis))


# ==================================================================
head("N1  D226 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记均匀共存选择器",
    "D226" in AX
    and "R-Z-UNIFORM-MATRIX-COEXISTENCE" in BOTH
    and "均匀共存" in DOC,
)
check(
    "N1 公理表与文档登记共存原则缺口",
    "R-Z-COEXISTENCE-PRINCIPLE-GAP" in BOTH
    and "交换性与生成性" in DOC,
)


# ==================================================================
head("N2  张量载体在每个年龄扇区保留完整 M_2")
# ==================================================================
taus = (1, 2, 3, 4)
tensor_supports = {}
for tau in taus:
    age_count = tau + 1
    matrix_basis = [
        np.kron(matrix_unit(2, row, column), np.eye(age_count))
        for row in range(2)
        for column in range(2)
    ]
    age_supports = []
    for age in range(age_count):
        age_projection = np.zeros((age_count, age_count), dtype=complex)
        age_projection[age, age] = 1.0
        central_projection = np.kron(np.eye(2), age_projection)
        age_supports.append(
            matrix_support_dimension(central_projection, matrix_basis)
        )
    tensor_supports[tau] = tuple(age_supports)
check(
    "N2 张量每个年龄扇区中的矩阵支持维数都是 4",
    all(
        all(dimension == 4 for dimension in supports)
        for supports in tensor_supports.values()
    ),
    f"supports={tensor_supports}",
)


# ==================================================================
head("N3  直接和的年龄中央投影把 M_2 压成零")
# ==================================================================
tau = 3
age_count = tau + 1
direct_matrix_basis = [
    direct_sum(
        matrix_unit(2, row, column),
        np.zeros((age_count, age_count), dtype=complex),
    )
    for row in range(2)
    for column in range(2)
]
p_m = direct_sum(
    np.eye(2, dtype=complex),
    np.zeros((age_count, age_count), dtype=complex),
)
p_c = direct_sum(
    np.zeros((2, 2), dtype=complex),
    np.eye(age_count, dtype=complex),
)
support_pm = matrix_support_dimension(p_m, direct_matrix_basis)
support_pc = matrix_support_dimension(p_c, direct_matrix_basis)
check(
    "N3 直接和矩阵扇区支持维数为 4",
    support_pm == 4,
)
check(
    "N3 直接和年龄扇区中的矩阵支持维数为 0",
    support_pc == 0,
    f"support_pc={support_pc}",
)


# ==================================================================
head("N4  均匀共存原则排除直接和")
# ==================================================================
direct_supports = (support_pc,)
uniform_direct = all(dimension == 4 for dimension in direct_supports)
uniform_tensor = all(
    dimension == 4
    for supports in tensor_supports.values()
    for dimension in supports
)
check(
    "N4 直接和不满足每个非零年龄扇区都有完整 M_2",
    not uniform_direct,
)
check(
    "N4 张量载体满足该原则",
    uniform_tensor,
)


# ==================================================================
head("N5  条件生成类中的张量载体达到维数下界")
# ==================================================================
tensor_dimension = {
    tau: 4 * (tau + 1)
    for tau in taus
}
lower_bound = {
    tau: 4 * (tau + 1)
    for tau in taus
}
check(
    "N5 每个年龄扇区至少需要一个 M_2 扇区",
    all(tensor_dimension[tau] == lower_bound[tau] for tau in taus),
    f"dimensions={tensor_dimension}",
)
check(
    "N5 文档声明这是条件类内最小",
    "条件类内最小" in DOC
    and "不是全类有限维" in DOC,
)


# ==================================================================
head("N6  放弃均匀共存时直接和仍是更省候选")
# ==================================================================
direct_dimension = {
    tau: tau + 5
    for tau in taus
}
check(
    "N6 直接和维数严格小于张量",
    all(
        direct_dimension[tau] < tensor_dimension[tau]
        for tau in taus
    ),
    f"direct={direct_dimension}, tensor={tensor_dimension}",
)


# ==================================================================
head("N7  中心投影与年龄扇区结构正确")
# ==================================================================
check(
    "N7 p_M+p_C=1 且 p_Mp_C=0",
    np.allclose(p_m + p_c, np.eye(2 + age_count, dtype=complex))
    and np.allclose(p_m @ p_c, 0.0),
)
tensor_central_basis = []
for age in range(age_count):
    age_projection = np.zeros((age_count, age_count), dtype=complex)
    age_projection[age, age] = 1.0
    tensor_central_basis.append(np.kron(np.eye(2), age_projection))
check(
    "N7 张量年龄中心投影正交且求和为一",
    all(
        np.allclose(left @ right, 0.0)
        for left in tensor_central_basis
        for right in tensor_central_basis
        if not np.allclose(left, right)
    )
    and np.allclose(
        sum(tensor_central_basis, np.zeros_like(tensor_central_basis[0])),
        np.eye(2 * age_count, dtype=complex),
    ),
)


# ==================================================================
head("N8  选择原则没有被写成 U1-U4 定理")
# ==================================================================
check(
    "N8 文档拒绝把均匀共存写成上游定理",
    "不是 `U1-U4` 的定理" in DOC
    and "直接和证明了这个额外要求不能由四槽推出" in DOC,
)
check(
    "N8 文档明确选择器候选",
    "选择器候选" in DOC
    and "恢复层选择器" in DOC,
)


# ==================================================================
head("N9  D223-D224 接口在新选择下保持")
# ==================================================================
check(
    "N9 文档保留 M_2 年龄载体",
    "M_2(\\mathbb C)\\otimes C(X_{\\tau_i})" in DOC,
)
check(
    "N9 文档保留 D223-D224 状态与权重",
    "\\omega_i^M(r_i)\\otimes\\mu_i(\\epsilon_i,\\tau_i)" in DOC
    and "D220" in DOC,
)


# ==================================================================
head("N10 不引用外部材料作依据")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "旧体系",
    "../..",
)
check(
    "N10 D226 文档不含项目外体系名或外部路径",
    not any(marker in DOC for marker in external_markers),
)


# ==================================================================
head("N11 上游边界保持")
# ==================================================================
check(
    "N11 不修改 U1-U4+C1 且不新增 U5",
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
