#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D230 —— L^infinity 年龄支持与投影压缩补单位
================================================================================
【被检验的问题（D230）】
  连续年龄支持能否用 L^infinity 的区间投影构造？
  A(I)=q_I A q_I + C(1-q_I) 是否是含单位的局部子代数？
  嵌套、无交、矩阵因子、模生成元与支持协变是否成立？

【本步判据】
  N1  D230 与三个恢复结构已登记
  N2  年龄投影幂等且正交
  N3  A(I) 含单位
  N4  A(I) 对乘法和共轭闭合
  N5  嵌套支持给出 A(I) 包含于 A(J)
  N6  无交支持交换
  N7  每个正测度支持保留完整矩阵因子
  N8  区间局部模生成元非零且属于 A(I)
  N9  模流固定局部生成元
  N10 年龄前缀推进满足半群律
  N11 模支持协变成立
  N12 有限细化在弱极限中给区间投影加细
  N13 文档登记 L^infinity 构造、补单位规则与测度缺口
  N14 文档不把连续支持写成几何导出
  N15 上游边界保持
  N16 文档不引用项目外体系或外部路径

运行：python3 verify/d230_linfinity_age_support_and_projection_complement.py
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
            "D230_linfinity_age_support_and_projection_complement.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


def matrix_unit(dimension, row, column):
    matrix = np.zeros((dimension, dimension), dtype=complex)
    matrix[row, column] = 1.0
    return matrix


def interval_projection(interval, site_count):
    return np.diag(
        [1.0 if site in interval else 0.0 for site in range(site_count)]
    ).astype(complex)


def kron_matrix_age(matrix, age_matrix):
    return np.kron(matrix, age_matrix)


def in_span(matrix, basis, tol=1e-9):
    columns = np.column_stack([entry.reshape(-1) for entry in basis])
    solution, _, _, _ = np.linalg.lstsq(
        columns,
        matrix.reshape(-1),
        rcond=None,
    )
    residual = columns @ solution - matrix.reshape(-1)
    return np.linalg.norm(residual) < tol


def support_algebra_basis(interval, site_count):
    q_interval = interval_projection(interval, site_count)
    basis = []
    for site in sorted(interval):
        site_projection = interval_projection({site}, site_count)
        for row in range(2):
            for column in range(2):
                basis.append(
                    kron_matrix_age(
                        matrix_unit(2, row, column),
                        site_projection,
                    )
                )
    basis.append(
        kron_matrix_age(np.eye(2), np.eye(site_count) - q_interval)
    )
    return q_interval, basis


def matrix_support_dimension(projection, matrix_basis):
    compressed = [
        projection @ basis @ projection
        for basis in matrix_basis
    ]
    vector_basis = np.column_stack(
        [matrix.reshape(-1) for matrix in compressed]
    )
    return int(np.linalg.matrix_rank(vector_basis))


site_count = 8
sites = set(range(site_count))
I = {0, 1, 2}
J = {0, 1, 2, 3, 4, 5}
K_disjoint = {5, 6, 7}

q_I, basis_I = support_algebra_basis(I, site_count)
q_J, basis_J = support_algebra_basis(J, site_count)
q_K, basis_K = support_algebra_basis(K_disjoint, site_count)


# ==================================================================
head("N1  D230 与三个恢复结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记 L^infinity 支持与补单位构造",
    "D230" in AX
    and "R-Z-LINFINITY-AGE-SUPPORT-ALGEBRA" in BOTH
    and "R-Z-PROJECTION-COMPRESSION-COMPLEMENT" in BOTH
    and "投影压缩加补单位" in DOC,
)
check(
    "N1 公理表与文档登记连续测度缺口",
    "R-Z-CONTINUOUS-AGE-MEASURE-GAP" in BOTH
    and "测度" in DOC,
)


# ==================================================================
head("N2  年龄投影幂等且正交")
# ==================================================================
q_empty = interval_projection(set(), site_count)
q_full = interval_projection(sites, site_count)
check(
    "N2 年龄投影幂等",
    np.allclose(q_I @ q_I, q_I)
    and np.allclose(q_J @ q_J, q_J)
    and np.allclose(q_empty, 0.0)
    and np.allclose(q_full, np.eye(site_count, dtype=complex)),
)
check(
    "N2 无交年龄投影正交",
    np.allclose(q_I @ q_K, 0.0) and np.allclose(q_K @ q_I, 0.0),
)


# ==================================================================
head("N3  A(I) 含单位")
# ==================================================================
global_identity = kron_matrix_age(np.eye(2), np.eye(site_count))
check(
    "N3 A(I) 含全局单位",
    in_span(global_identity, basis_I),
)


# ==================================================================
head("N4  A(I) 对乘法和共轭闭合")
# ==================================================================
product_closed = all(
    in_span(left @ right, basis_I)
    for left in basis_I
    for right in basis_I
)
star_closed = all(
    in_span(left.conj().T, basis_I)
    for left in basis_I
)
check(
    "N4 A(I) 对乘法闭合",
    product_closed,
)
check(
    "N4 A(I) 对共轭闭合",
    star_closed,
)


# ==================================================================
head("N5  嵌套支持给出 A(I) 包含于 A(J)")
# ==================================================================
check(
    "N5 I 包含于 J 且支持代数单调",
    I.issubset(J)
    and all(in_span(entry, basis_J) for entry in basis_I),
)


# ==================================================================
head("N6  无交支持交换")
# ==================================================================
check(
    "N6 无交支持代数交换",
    all(
        np.allclose(left @ right, right @ left)
        for left in basis_I
        for right in basis_K
    ),
)


# ==================================================================
head("N7  每个正测度支持保留完整矩阵因子")
# ==================================================================
matrix_basis = [
    kron_matrix_age(matrix_unit(2, row, column), np.eye(site_count))
    for row in range(2)
    for column in range(2)
]
local_central_projection = kron_matrix_age(np.eye(2), q_I)
matrix_support = matrix_support_dimension(
    local_central_projection,
    matrix_basis,
)
check(
    "N7 区间支持中的矩阵支持维数为 4",
    matrix_support == 4,
    f"support={matrix_support}",
)


# ==================================================================
head("N8  区间局部模生成元非零且属于 A(I)")
# ==================================================================
r_matrix = 0.75
k_matrix = np.diag(
    [np.log(r_matrix), np.log(1.0 - r_matrix)]
).astype(complex)
local_generator = kron_matrix_age(k_matrix, q_I)
check(
    "N8 区间局部模生成元非零",
    np.linalg.norm(local_generator) > 1e-10,
)
check(
    "N8 区间局部模生成元属于 A(I)",
    in_span(local_generator, basis_I),
)


# ==================================================================
head("N9  模流固定局部生成元")
# ==================================================================
phase_rate = np.log(r_matrix / (1.0 - r_matrix))
flow_matrix = np.diag([1.0, np.exp(1j * phase_rate)]).astype(complex)
global_flow = kron_matrix_age(flow_matrix, np.eye(site_count))
check(
    "N9 模流固定局部模生成元",
    np.allclose(
        global_flow @ local_generator @ global_flow.conj().T,
        local_generator,
    ),
)


# ==================================================================
head("N10  年龄前缀推进满足半群律")
# ==================================================================
T = float(site_count - 1)


def prefix_endpoint(a):
    return min(a, T)


def age_shift(a, s):
    return min(a + s, T)


semigroup_ok = all(
    abs(age_shift(age_shift(a, s), t) - age_shift(a, s + t)) < 1e-12
    for a in np.linspace(0.0, T, 41)
    for s in (0.0, 0.25, 0.5, 1.0, 3.0, 10.0)
    for t in (0.0, 0.5, 1.5, 4.0)
)
check(
    "N10 年龄前缀推进满足半群律",
    semigroup_ok,
)


# ==================================================================
head("N11  模支持协变成立")
# ==================================================================
flow_covariant = True
for step in (0, 1, 2, 3):
    expanded_interval = {
        site for site in sites
        if site <= prefix_endpoint(max(I) + step)
    }
    _, expanded_basis = support_algebra_basis(expanded_interval, site_count)
    for entry in basis_I:
        flowed = global_flow @ entry @ global_flow.conj().T
        if not in_span(flowed, expanded_basis):
            flow_covariant = False
check(
    "N11 模流把支持代数送到扩张支持代数",
    flow_covariant,
)


# ==================================================================
head("N12  有限细化在弱极限中给区间投影加细")
# ==================================================================
coarse_block = {0, 1}
fine_children = {0} | {1}
coarse_projection = interval_projection(coarse_block, site_count)
child_sum = sum(
    (
        interval_projection({child}, site_count)
        for child in fine_children
    ),
    np.zeros((site_count, site_count), dtype=complex),
)
check(
    "N12 粗块投影等于细子块投影之和",
    np.allclose(coarse_projection, child_sum),
)
check(
    "N12 文档说明弱极限落在 L^infinity",
    "弱*稠密" in DOC
    and "自然极限是 $L^\\infty$" in DOC,
)


# ==================================================================
head("N13  文档登记 L^infinity 构造、补单位规则与测度缺口")
# ==================================================================
check(
    "N13 文档登记 L^infinity 与补单位规则",
    "L^\\infty" in DOC
    and "q_IAq_I" in DOC
    and "1-q_I" in DOC,
)
check(
    "N13 文档登记测度缺口",
    "R-Z-CONTINUOUS-AGE-MEASURE-GAP" in DOC
    and "测度 $\\mu$ 的来源" in DOC,
)


# ==================================================================
head("N14  文档不把连续支持写成几何导出")
# ==================================================================
check(
    "N14 文档保留几何识别缺口",
    "几何 boost" in DOC
    and "未导出" in DOC
    and "区域到四维几何的嵌入" in DOC,
)


# ==================================================================
head("N15  上游边界保持")
# ==================================================================
check(
    "N15 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)


# ==================================================================
head("N16  文档不引用项目外体系或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N16 D230 文档不含项目外体系名或外部路径",
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
