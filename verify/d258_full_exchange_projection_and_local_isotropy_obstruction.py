#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D258 —— 全边交换投影、K4 普遍固定点与局部各向异性障碍
================================================================================
【被检验的问题（D258）】
  全边交换刚度是否是有效电阻嵌入的范围投影？
  任意正权重四顶点完全图是否都能经尺度匹配 P1？
  单位权重四顶点图是否回收 c=576？
  局部四面体的 P1 相容是否只是“选择合适尺度”？
  加入外部边后局部各向同性是否可能被破坏？

【本步判据】
  N1  D258 与恢复结构已登记
  N2  全边交换刚度是范围投影
  N3  任意正权重 K4 的唯一固定点尺度
  N4  单位权重 K4 回收 c=576
  N5  局部 P1 相容的标量条件
  N6  外部边算例破坏局部各向同性
  N7  文档登记局部化缺口
  N8  文档不把局部几何写成自动结论
  N9  上游边界保持
  N10 文档不引用项目外体系名或外部路径

运行：python3 verify/d258_full_exchange_projection_and_local_isotropy_obstruction.py
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


def weighted_laplacian(vertex_count, edges, weights):
    laplacian = np.zeros((vertex_count, vertex_count), dtype=float)
    for edge, weight in zip(edges, weights):
        i, j = edge
        laplacian[i, i] += weight
        laplacian[j, j] += weight
        laplacian[i, j] -= weight
        laplacian[j, i] -= weight
    return laplacian


def resistance_embedding(laplacian):
    """Return c=1 coordinates X with |X_i-X_j|^2=R_ij."""
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian)
    positive = eigenvalues > 1.0e-10
    return eigenvectors[:, positive] / np.sqrt(eigenvalues[positive])[None, :]


def local_coordinates(points, subset):
    """Orthonormal 3D coordinates for the affine span of a local subset."""
    base = points[subset[0]]
    span_vectors = np.column_stack(
        [points[vertex] - base for vertex in subset[1:]]
    )
    left, singular_values, _ = np.linalg.svd(span_vectors, full_matrices=False)
    if int(np.sum(singular_values > 1.0e-10)) != 3:
        raise ValueError("local edge vectors do not span three dimensions")
    basis = left[:, :3]
    return points @ basis, float(abs(np.linalg.det(span_vectors.T @ basis)))


def local_exchange_stiffness(points, subset, edges, weights):
    stiffness = np.zeros((3, 3), dtype=float)
    for edge, weight in zip(edges, weights):
        i, j = edge
        if i in subset and j in subset:
            edge_vector = points[j] - points[i]
            stiffness += weight * np.outer(edge_vector, edge_vector)
    return stiffness


AX = io.open(os.path.join(ROOT, "AXIOMS.md"), encoding="utf-8").read()
DOC = canonical_math(
    io.open(
        os.path.join(
            ROOT,
            "derivations",
            "D258_full_exchange_projection_and_local_isotropy_obstruction.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D258 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D258", "D258" in AX)
check(
    "N1 公理表与文档登记全边交换投影",
    "R-Z-EFFECTIVE-EXCHANGE-PROJECTION" in BOTH,
)
check(
    "N1 公理表与文档登记 K4 普遍固定点",
    "R-Z-K4-UNIVERSAL-FIXED-POINT" in BOTH,
)
check(
    "N1 公理表与文档登记局部各向同性条件",
    "R-Z-LOCAL-SIMPLEX-ISOTROPY-CONDITION" in BOTH,
)
check(
    "N1 公理表与文档登记局部反作用缺口",
    "R-Z-LOCAL-EXCHANGE-BACKREACTION-GAP" in BOTH,
)


# ==================================================================
head("N2  全边交换刚度是范围投影")
# ==================================================================
k4_edges = list(itertools.combinations(range(4), 2))
k4_weights = np.array([0.7, 1.3, 0.9, 1.1, 0.8, 1.2], dtype=float)
k4_laplacian = weighted_laplacian(4, k4_edges, k4_weights)
k4_points = resistance_embedding(k4_laplacian)
k4_stiffness = np.zeros((3, 3), dtype=float)
for edge, weight in zip(k4_edges, k4_weights):
    i, j = edge
    edge_vector = k4_points[j] - k4_points[i]
    k4_stiffness += weight * np.outer(edge_vector, edge_vector)

check(
    "N2 非均匀正权重 K4 的嵌入维数为三",
    k4_points.shape == (4, 3),
    f"shape={k4_points.shape}",
)
check(
    "N2 全边交换刚度等于 I3",
    np.allclose(k4_stiffness, np.eye(3)),
    f"max_error={np.max(np.abs(k4_stiffness - np.eye(3))):.3e}",
)
check(
    "N2 全边刚度不依赖权重是否均匀",
    not np.allclose(k4_weights, k4_weights[0]),
)


# ==================================================================
head("N3  任意正权重 K4 的唯一固定点尺度")
# ==================================================================
k4_volume_unit = abs(
    np.linalg.det(
        np.column_stack(
            [
                k4_points[1] - k4_points[0],
                k4_points[2] - k4_points[0],
                k4_points[3] - k4_points[0],
            ]
        )
    )
) / 6.0
k4_scale = 1.0 / k4_volume_unit ** 2
k4_scaled_points = math.sqrt(k4_scale) * k4_points
k4_scaled_stiffness = k4_scale * k4_stiffness
k4_p1_scaled = (k4_scale ** 1.5) * k4_volume_unit * np.eye(3)

check(
    "N3 c=1 体积为正",
    k4_volume_unit > 0.0,
    f"V1={k4_volume_unit:.12f}",
)
check(
    "N3 尺度 c=V1^{-2} 使全边交换刚度等于连续 P1",
    np.allclose(k4_scaled_stiffness, k4_p1_scaled),
    f"c={k4_scale:.12f}",
)
check(
    "N3 尺度重标度保持嵌入形状",
    np.allclose(
        np.diag(k4_scaled_points @ k4_scaled_points.T),
        k4_scale * np.diag(k4_points @ k4_points.T),
    ),
)


# ==================================================================
head("N4  单位权重 K4 回收 c=576")
# ==================================================================
uniform_edges = list(itertools.combinations(range(4), 2))
uniform_weights = np.ones(6, dtype=float)
uniform_laplacian = weighted_laplacian(4, uniform_edges, uniform_weights)
uniform_points = resistance_embedding(uniform_laplacian)
uniform_volume = abs(
    np.linalg.det(
        np.column_stack(
            [
                uniform_points[1] - uniform_points[0],
                uniform_points[2] - uniform_points[0],
                uniform_points[3] - uniform_points[0],
            ]
        )
    )
) / 6.0
uniform_scale = 1.0 / uniform_volume ** 2

check(
    "N4 单位权重 K4 的 c=1 体积为 1/24",
    np.isclose(uniform_volume, 1.0 / 24.0),
    f"V1={uniform_volume:.12f}",
)
check(
    "N4 单位权重 K4 的固定点尺度为 576",
    np.isclose(uniform_scale, 576.0),
    f"c={uniform_scale:.12f}",
)
check(
    "N4 文档登记 c=576",
    "c=576" in DOC and "V_1" in DOC,
)


# ==================================================================
head("N5  局部 P1 相容的标量条件")
# ==================================================================
scalar_local = np.eye(3)
non_scalar_local = np.diag([0.8, 1.0, 1.0])
scalar_eigenvalues = np.linalg.eigvalsh(scalar_local)
non_scalar_eigenvalues = np.linalg.eigvalsh(non_scalar_local)

check(
    "N5 标量局部刚度满足各向同性条件",
    np.allclose(scalar_eigenvalues, scalar_eigenvalues[0]),
)
check(
    "N5 非标量局部刚度不满足各向同性条件",
    not np.allclose(non_scalar_eigenvalues, non_scalar_eigenvalues[0]),
    f"eig={non_scalar_eigenvalues}",
)
check(
    "N5 文档登记局部相容必要条件",
    "A_S^0" in DOC and "标量矩阵" in DOC,
)


# ==================================================================
head("N6  外部边算例破坏局部各向同性")
# ==================================================================
external_edges = k4_edges + [(0, 4), (1, 4)]
external_weights = np.ones(len(external_edges), dtype=float)
external_laplacian = weighted_laplacian(5, external_edges, external_weights)
external_points = resistance_embedding(external_laplacian)
external_local_points, external_local_volume = local_coordinates(
    external_points,
    (0, 1, 2, 3),
)
external_local_stiffness = local_exchange_stiffness(
    external_local_points,
    (0, 1, 2, 3),
    external_edges,
    external_weights,
)
external_local_eigenvalues = np.linalg.eigvalsh(external_local_stiffness)
external_scale_condition = np.allclose(
    external_local_eigenvalues,
    external_local_eigenvalues[0],
    atol=1.0e-8,
)

check(
    "N6 外部边算例局部体积为正",
    external_local_volume > 0.0,
    f"V1={external_local_volume:.12f}",
)
check(
    "N6 外部边算例局部特征值为 4/5,1,1",
    np.allclose(
        external_local_eigenvalues,
        np.array([0.8, 1.0, 1.0]),
        atol=1.0e-8,
    ),
    f"eig={np.round(external_local_eigenvalues, 10)}",
)
check(
    "N6 外部边算例不存在统一各向同性尺度",
    not external_scale_condition,
)
check(
    "N6 文档登记外部边反作用",
    "外部边" in DOC and "4/5" in DOC and "1,1" in DOC,
)


# ==================================================================
head("N7  文档登记局部化缺口")
# ==================================================================
check(
    "N7 文档登记局部化规则是输入",
    "局部化规则" in DOC and "未导出" in DOC,
)
check(
    "N7 文档登记 K4 正例不选择权重",
    "P1 固定点固定的是尺度，不是交换权重本身" in DOC,
)
check(
    "N7 文档登记输入预算",
    "## §3 输入预算" in DOC,
)


# ==================================================================
head("N8  文档不把局部几何写成自动结论")
# ==================================================================
check(
    "N8 文档明确条件性",
    "不是 `U1-U4` 的推论" in DOC and "条件" in DOC,
)
check(
    "N8 文档否认任意大图自动给局部几何",
    "任意大于四顶点的图自动给局部 P1 几何" in DOC,
)
check(
    "N8 文档否认外部边不影响局部度量",
    "外部边只是局部背景" in DOC,
)


# ==================================================================
head("N9  上游边界保持")
# ==================================================================
check(
    "N9 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)
check(
    "N9 上游边界声明出现",
    "不修改 `U1-U4+C1`" in DOC,
)


# ==================================================================
head("N10  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N10 D258 文档不含项目外体系名或外部路径",
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
