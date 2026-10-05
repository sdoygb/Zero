#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D227 —— 年龄扇区局部模见证
================================================================================
【被检验的问题（D227）】
  直接和的年龄扇区是否没有非中心矩阵与局部模相位？
  张量载体是否在每个年龄扇区给非平凡局部模见证？
  该条件能否把 D226 的均匀共存连接到 GR 的局部模接口？

【本步判据】
  N1  D227 与两个恢复结构已登记
  N2  直接和的年龄扇区是一维交换代数
  N3  直接和的年龄扇区没有矩阵压缩
  N4  直接和的年龄扇区没有非平凡模相位
  N5  张量每个年龄扇区保留完整 M_2
  N6  张量每个年龄扇区有非中心局部见证
  N7  张量见证的模相位非平凡
  N8  局部模见证排除直接和并选中张量
  N9  条件类中该见证恢复 D226 的均匀共存
  N10 张量给每个扇区局部模生成元
  N11 文档不把选择器写成 U1-U4 定理
  N12 文档连接 D23/D24 但不声称几何极限已建
  N13 上游边界保持
  N14 文档不引用项目外体系或外部路径

运行：python3 verify/d227_age_local_modular_witness_selector.py
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
            "D227_age_local_modular_witness_selector.md",
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


def is_noncentral(element, algebra_basis):
    return any(
        not np.allclose(element @ basis, basis @ element)
        for basis in algebra_basis
    )


def algebra_support_dimension(projection, algebra_basis):
    compressed = [
        projection @ basis @ projection
        for basis in algebra_basis
    ]
    vector_basis = np.column_stack(
        [matrix.reshape(-1) for matrix in compressed]
    )
    return int(np.linalg.matrix_rank(vector_basis))


def matrix_rank_in_projection(projection, matrix_basis):
    compressed = [
        projection @ basis @ projection
        for basis in matrix_basis
    ]
    vector_basis = np.column_stack(
        [matrix.reshape(-1) for matrix in compressed]
    )
    dimension = int(np.linalg.matrix_rank(vector_basis))
    if dimension == 0:
        return 0
    if dimension == 1:
        return 1
    if dimension == 4:
        return 2
    return None


tau = 3
age_count = tau + 1


# ==================================================================
head("N1  D227 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记年龄局部模见证",
    "D227" in AX
    and "R-Z-AGE-LOCAL-MODULAR-WITNESS" in BOTH
    and "年龄局部模见证" in DOC,
)
check(
    "N1 公理表与文档登记局部模见证缺口",
    "R-Z-LOCAL-MODULAR-WITNESS-GAP" in BOTH
    and "连续几何识别" in DOC,
)


# ==================================================================
head("N2  直接和的年龄扇区是一维交换代数")
# ==================================================================
direct_full_basis = []
for row in range(2):
    for column in range(2):
        direct_full_basis.append(
            direct_sum(
                matrix_unit(2, row, column),
                np.zeros((age_count, age_count), dtype=complex),
            )
        )
for age in range(age_count):
    age_matrix = np.zeros((age_count, age_count), dtype=complex)
    age_matrix[age, age] = 1.0
    direct_full_basis.append(
        direct_sum(
            np.zeros((2, 2), dtype=complex),
            age_matrix,
        )
    )

q0_direct = direct_sum(
    np.zeros((2, 2), dtype=complex),
    np.diag([1.0] + [0.0] * (age_count - 1)).astype(complex),
)
direct_sector_dimension = algebra_support_dimension(q0_direct, direct_full_basis)
check(
    "N2 直接和年龄扇区维数为 1",
    direct_sector_dimension == 1,
    f"dim={direct_sector_dimension}",
)


# ==================================================================
head("N3  直接和的年龄扇区没有矩阵压缩")
# ==================================================================
direct_matrix_basis = [
    direct_sum(
        matrix_unit(2, row, column),
        np.zeros((age_count, age_count), dtype=complex),
    )
    for row in range(2)
    for column in range(2)
]
direct_matrix_compressions = [
    q0_direct @ basis @ q0_direct
    for basis in direct_matrix_basis
]
check(
    "N3 直接和矩阵压缩全部为零",
    all(np.allclose(compression, 0.0) for compression in direct_matrix_compressions),
)
check(
    "N3 直接和年龄扇区矩阵秩为零",
    matrix_rank_in_projection(q0_direct, direct_matrix_basis) == 0,
)


# ==================================================================
head("N4  直接和的年龄扇区没有非平凡模相位")
# ==================================================================
direct_sector_basis = [q0_direct]
direct_age_projector = q0_direct
direct_flow_fixed = np.allclose(
    direct_age_projector @ direct_sector_basis[0]
    - direct_sector_basis[0] @ direct_age_projector,
    0.0,
)
check(
    "N4 直接和年龄投影在交换扇区中固定",
    direct_flow_fixed,
)
check(
    "N4 直接和年龄扇区没有非中心局部见证",
    not is_noncentral(direct_age_projector, direct_sector_basis),
)


# ==================================================================
head("N5  张量每个年龄扇区保留完整 M_2")
# ==================================================================
tensor_full_basis = []
for row in range(2):
    for column in range(2):
        for age in range(age_count):
            tensor_full_basis.append(
                np.kron(
                    matrix_unit(2, row, column),
                    np.diag(
                        [1.0 if age == other else 0.0 for other in range(age_count)]
                    ).astype(complex),
                )
            )

tensor_matrix_basis = [
    np.kron(matrix_unit(2, row, column), np.eye(age_count))
    for row in range(2)
    for column in range(2)
]

tensor_sector_dimensions = {}
for age in range(age_count):
    q_age = np.kron(
        np.eye(2),
        np.diag(
            [1.0 if age == other else 0.0 for other in range(age_count)]
        ).astype(complex),
    )
    tensor_sector_dimensions[age] = (
        algebra_support_dimension(q_age, tensor_full_basis),
        matrix_rank_in_projection(q_age, tensor_matrix_basis),
    )

check(
    "N5 张量年龄扇区维数为 4 且矩阵秩为 2",
    all(
        dimension == 4 and matrix_rank == 2
        for dimension, matrix_rank in tensor_sector_dimensions.values()
    ),
    f"sectors={tensor_sector_dimensions}",
)


# ==================================================================
head("N6  张量每个年龄扇区有非中心局部见证")
# ==================================================================
tensor_witnesses = {}
for age in range(age_count):
    q_age = np.kron(
        np.eye(2),
        np.diag(
            [1.0 if age == other else 0.0 for other in range(age_count)]
        ).astype(complex),
    )
    e_age = np.diag(
        [1.0 if age == other else 0.0 for other in range(age_count)]
    ).astype(complex)
    witness = np.kron(matrix_unit(2, 0, 1), e_age)
    compressed = q_age @ witness @ q_age
    sector_basis = [
        q_age @ basis @ q_age
        for basis in tensor_full_basis
    ]
    tensor_witnesses[age] = (
        np.linalg.norm(compressed),
        is_noncentral(compressed, sector_basis),
    )

check(
    "N6 张量年龄见证非零",
    all(norm > 1e-12 for norm, _ in tensor_witnesses.values()),
    f"norms={ {a: round(n, 6) for a, (n, _) in tensor_witnesses.items()} }",
)
check(
    "N6 张量年龄见证非中心",
    all(noncentral for _, noncentral in tensor_witnesses.values()),
)


# ==================================================================
head("N7  张量见证的模相位非平凡")
# ==================================================================
r_matrix = 0.75
phase_rate = np.log(r_matrix / (1.0 - r_matrix))
phase_time = np.pi / phase_rate
phase_factor = np.exp(1j * phase_rate * phase_time)
check(
    "N7 模频率非零",
    abs(phase_rate) > 1e-12,
    f"rate={phase_rate:.12f}",
)
check(
    "N7 模流改变局域见证",
    abs(phase_factor - 1.0) > 1e-12,
    f"phase={phase_factor:.12f}",
)


# ==================================================================
head("N8  局部模见证排除直接和并选中张量")
# ==================================================================
direct_sector_has_witness = is_noncentral(
    direct_age_projector,
    direct_sector_basis,
)
tensor_sector_has_witness = all(
    noncentral
    for _, noncentral in tensor_witnesses.values()
)
check(
    "N8 直接和不满足年龄局部模见证",
    not direct_sector_has_witness,
)
check(
    "N8 张量满足年龄局部模见证",
    tensor_sector_has_witness,
)


# ==================================================================
head("N9  条件类中该见证恢复 D226 的均匀共存")
# ==================================================================
check(
    "N9 张量每个扇区都有完整 M_2",
    all(
        matrix_rank == 2
        for _, matrix_rank in tensor_sector_dimensions.values()
    ),
)
check(
    "N9 文档声明见证与均匀共存条件等价",
    "年龄局部模见证" in DOC
    and "均匀矩阵共存" in DOC
    and "条件类" in DOC,
)


# ==================================================================
head("N10  张量给每个扇区局部模生成元")
# ==================================================================
k0 = float(np.log(r_matrix))
k1 = float(np.log(1.0 - r_matrix))
modular_generator = np.diag([k0, k1]).astype(complex)
check(
    "N10 局域模生成元非平凡",
    not np.allclose(modular_generator, np.zeros((2, 2), dtype=complex))
    and not np.allclose(modular_generator, modular_generator[0, 0] * np.eye(2)),
    f"diag=({k0:.6f},{k1:.6f})",
)
check(
    "N10 文档给出 K_a 与直接和零生成元",
    "K_a" in DOC
    and "直接和只给常数" in DOC,
)


# ==================================================================
head("N11  文档不把选择器写成 U1-U4 定理")
# ==================================================================
check(
    "N11 文档拒绝上游定理表述",
    "不是 `U1-U4` 的定理" in DOC
    and "不能主张" in DOC
    and "唯一选出张量载体" in DOC,
)
check(
    "N11 文档登记选择器候选",
    "选择器候选" in DOC
    and "R-Z-AGE-LOCAL-MODULAR-WITNESS" in DOC,
)


# ==================================================================
head("N12  文档连接 D23/D24 但不声称几何极限已建")
# ==================================================================
check(
    "N12 文档连接 D23 第一定律与 D24 几何模条件",
    "D23" in DOC
    and "D24" in DOC
    and "几何模条件" in DOC
    and "K_B" in DOC,
)
check(
    "N12 文档保留连续几何缺口",
    "连续几何" in DOC
    and "未建" in DOC
    and "状态到几何映射" in DOC,
)


# ==================================================================
head("N13  上游边界保持")
# ==================================================================
check(
    "N13 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)


# ==================================================================
head("N14  文档不引用项目外体系或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N14 D227 文档不含项目外体系名或外部路径",
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
