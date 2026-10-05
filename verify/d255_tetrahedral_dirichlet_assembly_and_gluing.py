#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D255 —— 交换权重的四面体 Dirichlet 组装与面胶合
================================================================================
【被检验的问题（D255）】
  给定四面体顶点嵌入与交换权重，能否唯一装配局部 Dirichlet 张量？
  装配出的张量能否条件反解空间度规？
  相邻单元面胶合成功的充分条件与失败后果是什么？
  细化与单元几何还缺什么？

【本步判据】
  N1  D255 与恢复结构已登记
  N2  四面体体积与边向量
  N3  离散交换能量与 A 矩阵
  N4  装配公式 Q=A/omega
  N5  Q 正定并反解 h
  N6  面胶合成功
  N7  面胶合失败
  N8  细化与单元几何缺口
  N9  文档登记组装结果
  N10 文档不把单元几何写成上游推论
  N11 上游边界保持
  N12 文档不引用项目外体系名或外部路径

运行：python3 verify/d255_tetrahedral_dirichlet_assembly_and_gluing.py
"""
import io
import itertools
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
            "D255_tetrahedral_dirichlet_assembly_and_gluing.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


vertices = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
)
edges = list(itertools.combinations(range(4), 2))
edge_vectors = {
    edge: vertices[edge[1]] - vertices[edge[0]]
    for edge in edges
}
edge_weights = {
    (0, 1): 1.2,
    (0, 2): 0.9,
    (0, 3): 1.1,
    (1, 2): 0.7,
    (1, 3): 0.8,
    (2, 3): 1.0,
}
coordinate_volume = abs(
    np.linalg.det(
        np.column_stack(
            [
                vertices[1] - vertices[0],
                vertices[2] - vertices[0],
                vertices[3] - vertices[0],
            ]
        )
    )
) / 6.0
gradient_matrix = sum(
    edge_weights[edge] * np.outer(edge_vectors[edge], edge_vectors[edge])
    for edge in edges
)
Q_three = gradient_matrix / coordinate_volume
h_recovered = np.linalg.det(Q_three) * np.linalg.inv(Q_three)


def discrete_energy(gradient):
    return 0.5 * sum(
        edge_weights[edge] * float(np.dot(gradient, edge_vectors[edge])) ** 2
        for edge in edges
    )


def continuous_energy(gradient):
    return 0.5 * coordinate_volume * float(gradient @ Q_three @ gradient)


# ==================================================================
head("N1  D255 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D255", "D255" in AX)
check(
    "N1 公理表与文档登记四面体组装",
    "R-Z-TETRAHEDRAL-DIRICHLET-ASSEMBLY" in BOTH,
)
check(
    "N1 公理表与文档登记面胶合",
    "R-Z-DIRICHLET-ASSEMBLY-GLUING" in BOTH,
)
check(
    "N1 公理表与文档登记单元几何缺口",
    "R-Z-CELL-COMPLEX-GEOMETRY-GAP" in BOTH,
)
check(
    "N1 公理表与文档登记细化一致性缺口",
    "R-Z-REFINEMENT-TENSOR-CONSISTENCY-GAP" in BOTH,
)


# ==================================================================
head("N2  四面体体积与边向量")
# ==================================================================
check(
    "N2 单元体积为 1/6",
    np.isclose(coordinate_volume, 1.0 / 6.0),
    f"omega={coordinate_volume:.12f}",
)
check(
    "N2 边数等于四面体边数",
    len(edges) == 6,
)
check(
    "N2 至少三条边向量线性独立",
    np.linalg.matrix_rank(np.array([edge_vectors[edge] for edge in edges])) == 3,
)


# ==================================================================
head("N3  离散交换能量与 A 矩阵")
# ==================================================================
test_gradients = [
    np.array([1.0, 0.0, 0.0]),
    np.array([0.0, 1.0, 0.0]),
    np.array([0.3, -0.2, 0.7]),
]
energy_errors_a = [
    abs(discrete_energy(g_vec) - 0.5 * float(g_vec @ gradient_matrix @ g_vec))
    for g_vec in test_gradients
]

check(
    "N3 离散能量等于 g^T A g/2",
    max(energy_errors_a) < 1.0e-14,
    f"max_error={max(energy_errors_a):.3e}",
)
check(
    "N3 A 由边权重与边向量外积组成",
    np.allclose(
        gradient_matrix,
        sum(
            edge_weights[edge]
            * np.outer(edge_vectors[edge], edge_vectors[edge])
            for edge in edges
        ),
    ),
)


# ==================================================================
head("N4  装配公式 Q=A/omega")
# ==================================================================
energy_errors_q = [
    abs(discrete_energy(g_vec) - continuous_energy(g_vec))
    for g_vec in test_gradients
]

check(
    "N4 Q=A/omega 对多个梯度匹配能量",
    max(energy_errors_q) < 1.0e-14,
    f"max_error={max(energy_errors_q):.3e}",
)
check(
    "N4 文档登记显式装配公式",
    "Q_\\sigma" in DOC and "e_{ab}\\otimes e_{ab}" in DOC,
)


# ==================================================================
head("N5  Q 正定并反解 h")
# ==================================================================
Q_eigenvalues = np.linalg.eigvalsh(Q_three)
h_eigenvalues = np.linalg.eigvalsh(h_recovered)
Q_identity = math.sqrt(float(np.linalg.det(h_recovered))) * np.linalg.inv(h_recovered)

check(
    "N5 Q 为正定矩阵",
    np.all(Q_eigenvalues > 0.0),
    f"eig={np.round(Q_eigenvalues, 8)}",
)
check(
    "N5 h=(det Q)Q^{-1} 为正定矩阵",
    np.all(h_eigenvalues > 0.0),
    f"eig={np.round(h_eigenvalues, 8)}",
)
check(
    "N5 反解满足 Q=sqrt(det h)h^{-1}",
    np.allclose(Q_identity, Q_three),
)


# ==================================================================
head("N6  面胶合成功")
# ==================================================================
shared_face_vectors = np.array(
    [
        edge_vectors[(0, 1)],
        edge_vectors[(0, 2)],
        edge_vectors[(1, 2)],
    ]
)
same_Q = Q_three.copy()
same_h = np.linalg.det(same_Q) * np.linalg.inv(same_Q)
shared_energy_left = sum(
    float(vector @ Q_three @ vector)
    for vector in shared_face_vectors
)
shared_energy_right = sum(
    float(vector @ same_Q @ vector)
    for vector in shared_face_vectors
)

check(
    "N6 相同 Q 在共享面给相同能量配对",
    np.isclose(shared_energy_left, shared_energy_right),
)
check(
    "N6 相同 Q 反解相同 h",
    np.allclose(same_h, h_recovered),
)


# ==================================================================
head("N7  面胶合失败")
# ==================================================================
different_Q = 2.0 * np.eye(3)
different_h = np.linalg.det(different_Q) * np.linalg.inv(different_Q)
shared_energy_mismatch = sum(
    float(vector @ different_Q @ vector) - float(vector @ Q_three @ vector)
    for vector in shared_face_vectors
)

check(
    "N7 不同 Q 在共享面给不同能量配对",
    not np.isclose(shared_energy_mismatch, 0.0),
    f"mismatch={shared_energy_mismatch:.6f}",
)
check(
    "N7 不同 Q 反解不同 h",
    not np.allclose(different_h, h_recovered),
)
check(
    "N7 文档登记胶合失败 no-go",
    "没有单一叶层度规" in DOC
    or "胶合失败" in DOC,
)


# ==================================================================
head("N8  细化与单元几何缺口")
# ==================================================================
check(
    "N8 文档登记单元几何缺口",
    "R-Z-CELL-COMPLEX-GEOMETRY-GAP" in DOC
    and "顶点嵌入" in DOC,
)
check(
    "N8 文档登记细化一致性缺口",
    "R-Z-REFINEMENT-TENSOR-CONSISTENCY-GAP" in DOC
    and "唯一连续极限" in DOC,
)
check(
    "N8 文档不把面胶合写成连续唯一性",
    "面胶合" in DOC
    and "不是" in DOC
    and "唯一连续" in DOC,
)


# ==================================================================
head("N9  文档登记组装结果")
# ==================================================================
check(
    "N9 文档登记局部度规",
    "\\mathrm g_\\sigma" in DOC and "dt^2-h_\\sigma" in DOC,
)
check(
    "N9 文档登记条件定理",
    "条件定理" in DOC and "交换权重" in DOC,
)
check(
    "N9 文档登记输入预算",
    "## §3 输入预算" in DOC,
)


# ==================================================================
head("N10  文档不把单元几何写成上游推论")
# ==================================================================
check(
    "N10 文档明确条件性",
    "不是 `U1-U4` 的推论" in DOC
    and "条件闭合" in DOC,
)
check(
    "N10 文档不把单元几何写成已导出",
    "尚未由上游选择" in DOC
    and "未导出" in DOC,
)


# ==================================================================
head("N11  上游边界保持")
# ==================================================================
check(
    "N11 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)
check(
    "N11 上游边界声明出现",
    "不修改 `U1-U4+C1`" in DOC,
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
    "N12 D255 文档不含项目外体系名或外部路径",
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
