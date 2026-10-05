#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D229 —— 连续年龄支持的原子障碍
================================================================================
【被检验的问题（D229）】
  有限年龄扇区能否直接取连续极限？
  为什么 C([0,T]) 不能保留非平凡年龄原子？
  最小替代路线是什么？

【本步判据】
  N1  D229 与三个恢复结构已登记
  N2  有限年龄代数的中央投影完备正交
  N3  连续年龄代数没有非平凡幂等函数
  N4  非平凡连续 0/1 转换产生 q^2-q 缺陷
  N5  连续年龄推进满足半群律
  N6  直接连续化不能分解出至少两个非零年龄投影
  N7  有限细化块非空且粗粒化保序
  N8  粗块投影是细块投影之和
  N9  每层每个粗块仍保留张量矩阵扇区
  N10 文档登记连续年龄原子障碍
  N11 文档登记有限细化粗粒化路线
  N12 文档不把连续极限写成已建
  N13 上游边界保持
  N14 文档不引用项目外体系或外部路径

运行：python3 verify/d229_continuous_age_support_obstruction.py
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
            "D229_continuous_age_support_obstruction.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


def matrix_unit(dimension, row, column):
    matrix = np.zeros((dimension, dimension), dtype=complex)
    matrix[row, column] = 1.0
    return matrix


def age_projection(age, age_count):
    return np.diag(
        [1.0 if age == other else 0.0 for other in range(age_count)]
    ).astype(complex)


def block_projection(start, stop, fine_count):
    return sum(
        (age_projection(age, fine_count) for age in range(start, stop)),
        np.zeros((fine_count, fine_count), dtype=complex),
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
head("N1  D229 与三个恢复结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记连续年龄原子障碍",
    "D229" in AX
    and "R-Z-CONTINUOUS-AGE-ATOM-OBSTRUCTION" in BOTH
    and "没有非平凡投影" in DOC,
)
check(
    "N1 公理表与文档登记细化粗粒化路线与极限缺口",
    "R-Z-AGE-REFINEMENT-COARSEGRAINING-ROUTE" in BOTH
    and "R-Z-CONTINUOUS-AGE-LIMIT-GAP" in BOTH
    and "有限细化塔" in DOC,
)


# ==================================================================
head("N2  有限年龄代数的中央投影完备正交")
# ==================================================================
finite_count = 5
finite_projections = [
    age_projection(age, finite_count)
    for age in range(finite_count)
]
check(
    "N2 有限年龄投影正交",
    all(
        np.allclose(left @ right, 0.0)
        for left in finite_projections
        for right in finite_projections
        if not np.allclose(left, right)
    ),
)
check(
    "N2 有限年龄投影求和为一",
    np.allclose(
        sum(finite_projections, np.zeros_like(finite_projections[0])),
        np.eye(finite_count, dtype=complex),
    ),
)


# ==================================================================
head("N3  连续年龄代数没有非平凡幂等函数")
# ==================================================================
grid = np.linspace(0.0, 1.0, 2001)
constant_zero = np.zeros_like(grid)
constant_one = np.ones_like(grid)
nontrivial_candidates = {
    "linear": grid,
    "sine_squared": np.sin(np.pi * grid) ** 2,
}

check(
    "N3 常数函数是幂等函数",
    np.allclose(constant_zero**2 - constant_zero, 0.0)
    and np.allclose(constant_one**2 - constant_one, 0.0),
)
check(
    "N3 非平凡连续候选不幂等",
    all(
        np.max(np.abs(candidate**2 - candidate)) > 1e-3
        for candidate in nontrivial_candidates.values()
    ),
)


# ==================================================================
head("N4  非平凡连续 0/1 转换产生 q^2-q 缺陷")
# ==================================================================
linear_q = grid
midpoint_index = int(np.argmin(np.abs(linear_q - 0.5)))
midpoint_value = linear_q[midpoint_index]
midpoint_defect = midpoint_value**2 - midpoint_value
check(
    "N4 连续 0 到 1 必经过中间值",
    abs(midpoint_value - 0.5) < 1e-12,
    f"q_mid={midpoint_value:.6f}",
)
check(
    "N4 中间值处缺陷为 -1/4",
    abs(midpoint_defect + 0.25) < 1e-12,
    f"defect={midpoint_defect:.6f}",
)


# ==================================================================
head("N5  连续年龄推进满足半群律")
# ==================================================================
T = 4.0


def age_shift(age, step):
    return min(age + step, T)


sample_ages = np.linspace(0.0, T, 101)
semigroup_ok = all(
    abs(age_shift(age_shift(age, s), t) - age_shift(age, s + t)) < 1e-12
    for age in sample_ages
    for s in (0.0, 0.2, 0.5, 1.0, 2.0, 6.0)
    for t in (0.0, 0.3, 0.7, 1.5, 5.0)
)
order_ok = all(
    age_shift(left, step) <= age_shift(right, step) + 1e-12
    for left in sample_ages
    for right in sample_ages
    if left <= right
    for step in (0.0, 0.4, 1.0, 3.0, 7.0)
)
check(
    "N5 连续年龄推进满足半群律",
    semigroup_ok,
)
check(
    "N5 连续年龄推进保序",
    order_ok,
)


# ==================================================================
head("N6  直接连续化不能分解出至少两个非零年龄投影")
# ==================================================================
possible_projection_values = []
for value in np.linspace(0.0, 1.0, 10001):
    if abs(value**2 - value) < 1e-12:
        possible_projection_values.append(value)
check(
    "N6 连续投影值只有 0 和 1",
    np.allclose(np.unique(np.round(possible_projection_values, 12)), [0.0, 1.0]),
    f"values={np.unique(np.round(possible_projection_values, 12))}",
)
check(
    "N6 文档拒绝直接连续年龄扇区",
    "不能直接连续化" in DOC
    and "无法给每个年龄点一个非零扇区" in DOC,
)


# ==================================================================
head("N7  有限细化块非空且粗粒化保序")
# ==================================================================
coarse_count = 4
fine_count = 12
if fine_count % coarse_count != 0:
    raise SystemExit("核验配置错误：fine_count 必须能被 coarse_count 整除")
block_size = fine_count // coarse_count


def coarse_of_fine(age):
    return age // block_size


blocks = []
for coarse in range(coarse_count):
    blocks.append(range(coarse * block_size, (coarse + 1) * block_size))

check(
    "N7 每个粗年龄块非空",
    all(len(list(block)) > 0 for block in blocks),
)
check(
    "N7 粗粒化保序",
    all(
        (left <= right) <= (coarse_of_fine(left) <= coarse_of_fine(right))
        for left in range(fine_count)
        for right in range(fine_count)
    ),
)


# ==================================================================
head("N8  粗块投影是细块投影之和")
# ==================================================================
block_projections = [
    block_projection(block.start, block.stop, fine_count)
    for block in blocks
]
check(
    "N8 每个粗块投影是细投影之和",
    all(
        np.allclose(
            block_projections[coarse],
            sum(
                (
                    age_projection(age, fine_count)
                    for age in blocks[coarse]
                ),
                np.zeros((fine_count, fine_count), dtype=complex),
            ),
        )
        for coarse in range(coarse_count)
    ),
)
check(
    "N8 粗块投影完备",
    np.allclose(
        sum(block_projections, np.zeros_like(block_projections[0])),
        np.eye(fine_count, dtype=complex),
    ),
)


# ==================================================================
head("N9  每层每个粗块仍保留张量矩阵扇区")
# ==================================================================
matrix_basis = [
    np.kron(matrix_unit(2, row, column), np.eye(fine_count))
    for row in range(2)
    for column in range(2)
]
tensor_block_supports = []
for block in block_projections:
    central_projection = np.kron(np.eye(2), block)
    tensor_block_supports.append(
        matrix_support_dimension(central_projection, matrix_basis)
    )
check(
    "N9 每个粗块保留完整矩阵支持",
    all(dimension == 4 for dimension in tensor_block_supports),
    f"supports={tensor_block_supports}",
)


# ==================================================================
head("N10  文档登记连续年龄原子障碍")
# ==================================================================
check(
    "N10 文档登记障碍",
    "R-Z-CONTINUOUS-AGE-ATOM-OBSTRUCTION" in DOC
    and "有限年龄结果不能通过直接把" in DOC,
)


# ==================================================================
head("N11  文档登记有限细化粗粒化路线")
# ==================================================================
check(
    "N11 文档登记细化路线",
    "R-Z-AGE-REFINEMENT-COARSEGRAINING-ROUTE" in DOC
    and "粗粒化映射" in DOC
    and "在第 $N$ 层继续存在" in DOC,
)


# ==================================================================
head("N12  文档不把连续极限写成已建")
# ==================================================================
check(
    "N12 文档保留连续极限缺口",
    "R-Z-CONTINUOUS-AGE-LIMIT-GAP" in DOC
    and "连续年龄极限本身仍未建" in DOC
    and "未解选择器" in DOC,
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
    "N14 D229 文档不含项目外体系名或外部路径",
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
