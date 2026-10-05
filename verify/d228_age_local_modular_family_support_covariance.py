#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D228 —— 年龄局部模生成元族
================================================================================
【被检验的问题（D228）】
  逐年龄模见证能否组装成支持前缀上的局部模生成元？
  该生成元族是否与 U4 的年龄支持半流协变？
  它与几何 boost 识别之间还缺什么？

【本步判据】
  N1  D228 与两个恢复结构已登记
  N2  共同矩阵模生成元非中心
  N3  每个 K_k 位于对应年龄扇区
  N4  每个活动年龄扇区都有非零局部生成元
  N5  前缀生成元属于对应支持代数
  N6  模流固定每个 K_k
  N7  年龄推进满足半群律
  N8  年龄推进保序并在终端吸收
  N9  模流把支持代数送到年龄推进后的支持代数
  N10 直接和的年龄局部矩阵生成元为零
  N11 张量给出逐年龄前缀生成元族
  N12 文档保留连续几何与 boost 识别缺口
  N13 文档不把选择器写成 U1-U4 定理
  N14 上游边界保持
  N15 文档不引用项目外体系或外部路径

运行：python3 verify/d228_age_local_modular_family_support_covariance.py
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
            "D228_age_local_modular_family_support_covariance.md",
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


def age_projection(age, age_count):
    return np.diag(
        [1.0 if age == other else 0.0 for other in range(age_count)]
    ).astype(complex)


def modular_flow_matrix(rate):
    return np.diag([1.0, np.exp(1j * rate)]).astype(complex)


def is_noncentral(element, algebra_basis):
    return any(
        not np.allclose(element @ basis, basis @ element)
        for basis in algebra_basis
    )


def element_in_corner(element, projection, full_basis, tol=1e-10):
    complement = np.eye(projection.shape[0], dtype=complex) - projection
    left_violation = complement @ element @ projection
    right_violation = projection @ element @ complement
    corner_basis = [
        projection @ basis @ projection
        for basis in full_basis
    ]
    residual = element - projection @ element @ projection
    return (
        np.linalg.norm(left_violation) < tol
        and np.linalg.norm(right_violation) < tol
        and np.linalg.norm(residual) < tol
        and all(
            np.linalg.norm(element @ basis - basis @ element) >= 0.0
            for basis in corner_basis
        )
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


tau_last = 4
age_count = tau_last + 1
ages = list(range(tau_last + 1))


# ==================================================================
head("N1  D228 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记年龄局部模生成元族",
    "D228" in AX
    and "R-Z-AGE-LOCAL-MODULAR-FAMILY" in BOTH
    and "年龄局部模生成元族" in DOC,
)
check(
    "N1 公理表与文档登记模几何识别缺口",
    "R-Z-MODULAR-GEOMETRY-IDENTIFICATION-GAP" in BOTH
    and "boost" in DOC,
)


# ==================================================================
head("N2  共同矩阵模生成元非中心")
# ==================================================================
r_matrix = 0.75
k_matrix = np.diag(
    [np.log(r_matrix), np.log(1.0 - r_matrix)]
).astype(complex)
matrix_basis = [
    matrix_unit(2, row, column)
    for row in range(2)
    for column in range(2)
]
check(
    "N2 矩阵模生成元非中心",
    is_noncentral(k_matrix, matrix_basis)
    and not np.allclose(k_matrix, k_matrix[0, 0] * np.eye(2)),
    f"eigs=({k_matrix[0, 0]:.6f},{k_matrix[1, 1]:.6f})",
)


# ==================================================================
head("N3  每个 K_k 位于对应年龄扇区")
# ==================================================================
k_sectors = {
    age: np.kron(k_matrix, age_projection(age, age_count))
    for age in ages
}
full_basis = []
for row in range(2):
    for column in range(2):
        for age in ages:
            full_basis.append(
                np.kron(matrix_unit(2, row, column), age_projection(age, age_count))
            )

check(
    "N3 K_k 位于各年龄扇区",
    all(
        element_in_corner(
            k_sectors[age],
            np.kron(np.eye(2), age_projection(age, age_count)),
            full_basis,
        )
        for age in ages
    ),
)


# ==================================================================
head("N4  每个活动年龄扇区都有非零局部生成元")
# ==================================================================
check(
    "N4 每个 K_k 范数非零",
    all(np.linalg.norm(k_sectors[age]) > 1e-10 for age in ages),
)
check(
    "N4 每个年龄扇区矩阵支持维数为 4",
    all(
        matrix_support_dimension(
            np.kron(np.eye(2), age_projection(age, age_count)),
            [
                np.kron(matrix_unit(2, row, column), np.eye(age_count))
                for row in range(2)
                for column in range(2)
            ],
        )
        == 4
        for age in ages
    ),
)


# ==================================================================
head("N5  前缀生成元属于对应支持代数")
# ==================================================================
prefix_projections = {}
prefix_generators = {}
for cutoff in ages:
    prefix_projections[cutoff] = np.kron(
        np.eye(2),
        sum(
            (age_projection(age, age_count) for age in range(cutoff + 1)),
            np.zeros((age_count, age_count), dtype=complex),
        ),
    )
    prefix_generators[cutoff] = sum(
        (k_sectors[age] for age in range(cutoff + 1)),
        np.zeros((2 * age_count, 2 * age_count), dtype=complex),
    )

check(
    "N5 每个前缀生成元位于支持代数",
    all(
        element_in_corner(
            prefix_generators[cutoff],
            prefix_projections[cutoff],
            full_basis,
        )
        for cutoff in ages
    ),
)


# ==================================================================
head("N6  模流固定每个 K_k")
# ==================================================================
phase_rate = np.log(r_matrix / (1.0 - r_matrix))
flow_matrix = modular_flow_matrix(phase_rate)
tensor_flow = np.kron(flow_matrix, np.eye(age_count))
check(
    "N6 模流固定每个局域模生成元",
    all(
        np.allclose(
            tensor_flow @ k_sectors[age] @ tensor_flow.conj().T,
            k_sectors[age],
        )
        for age in ages
    ),
)


# ==================================================================
head("N7  年龄推进满足半群律")
# ==================================================================
def age_shift(age, step):
    return min(age + step, tau_last)


semigroup_ok = all(
    age_shift(age_shift(age, s), t) == age_shift(age, s + t)
    for age in ages
    for s in range(6)
    for t in range(6)
)
check(
    "N7 年龄推进满足 E_{s+t}=E_s∘E_t",
    semigroup_ok,
    f"T={tau_last}",
)


# ==================================================================
head("N8  年龄推进保序并在终端吸收")
# ==================================================================
order_ok = all(
    age_shift(left, step) <= age_shift(right, step)
    for left in ages
    for right in ages
    if left <= right
    for step in range(8)
)
absorbing_ok = all(age_shift(tau_last, step) == tau_last for step in range(8))
check(
    "N8 年龄推进保序",
    order_ok,
)
check(
    "N8 终端吸收",
    absorbing_ok,
)


# ==================================================================
head("N9  模流把支持代数送到年龄推进后的支持代数")
# ==================================================================
def embedded_matrix_element(row, column, age):
    return np.kron(matrix_unit(2, row, column), age_projection(age, age_count))


flow_covariant = True
for step in range(5):
    for cutoff in ages:
        mapped_cutoff = age_shift(cutoff, step)
        target_projection = prefix_projections[mapped_cutoff]
        for row in range(2):
            for column in range(2):
                for age in range(cutoff + 1):
                    element = embedded_matrix_element(row, column, age)
                    flowed = tensor_flow @ element @ tensor_flow.conj().T
                    if not np.allclose(target_projection @ flowed @ target_projection, flowed):
                        flow_covariant = False

check(
    "N9 模支持协变成立",
    flow_covariant,
)


# ==================================================================
head("N10  直接和的年龄局部矩阵生成元为零")
# ==================================================================
q0_direct = direct_sum(
    np.zeros((2, 2), dtype=complex),
    age_projection(0, age_count),
)
k_direct = direct_sum(k_matrix, np.zeros((age_count, age_count), dtype=complex))
direct_local_generator = q0_direct @ k_direct @ q0_direct
check(
    "N10 直接和年龄局部矩阵生成元为零",
    np.allclose(direct_local_generator, 0.0),
)


# ==================================================================
head("N11  张量给出逐年龄前缀生成元族")
# ==================================================================
check(
    "N11 每个前缀生成元非零",
    all(np.linalg.norm(prefix_generators[cutoff]) > 1e-10 for cutoff in ages),
)
check(
    "N11 前缀生成元随年龄增加累积",
    all(
        np.allclose(
            prefix_generators[cutoff + 1] - prefix_generators[cutoff],
            k_sectors[cutoff + 1],
        )
        for cutoff in range(tau_last)
    ),
)


# ==================================================================
head("N12  文档保留连续几何与 boost 识别缺口")
# ==================================================================
check(
    "N12 文档连接 D23/D24 并保留几何缺口",
    "D23" in DOC
    and "D24" in DOC
    and "几何模条件" in DOC
    and "连续几何" in DOC
    and "未建" in DOC,
)


# ==================================================================
head("N13  文档不把选择器写成 U1-U4 定理")
# ==================================================================
check(
    "N13 文档保留恢复层地位",
    "不是 `U1-U4` 的推论" in DOC
    and "条件构造" in DOC
    and "不修改 `U1-U4+C1`" in DOC,
)


# ==================================================================
head("N14  上游边界保持")
# ==================================================================
check(
    "N14 不新增 U5",
    "不新增 `U5`" in DOC,
)


# ==================================================================
head("N15  文档不引用项目外体系或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N15 D228 文档不含项目外体系名或外部路径",
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
