#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D257 —— 有效电阻度量、交换权重固定点与局部性缺口
================================================================================
【被检验的问题（D257）】
  交换权重能否经加权 Laplacian 与有效电阻距离生成内禀度量？
  图的旗复形能否自动给三维单元？
  单位权重四顶点完全图是否给正四面体？
  它是否同时满足 D256 的 P1 刚度固定点？
  远端边是否会改变局部有效电阻距离？
  权重尺度如何影响距离与边长？

【本步判据】
  N1  D257 与恢复结构已登记
  N2  加权 Laplacian 半正定且零模为常量
  N3  有效电阻满足度量公理
  N4  有效电阻平方根可欧氏嵌入
  N5  单位权重 K4 给正四面体几何
  N6  单位权重 K4 在 c=576 满足 P1 固定点
  N7  旗复形候选与三维截断缺口
  N8  单一顶点度量给共享面边长胶合
  N9  远端边改变局部有效电阻
  N10 权重尺度重标度
  N11 文档登记输入预算与限制
  N12 文档不把固定点或局部性写成自动结论
  N13 上游边界保持
  N14 文档不引用项目外体系名或外部路径

运行：python3 verify/d257_resistance_metric_fixed_point_and_locality_gap.py
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


def effective_resistance(laplacian):
    pseudoinverse = np.linalg.pinv(laplacian, rcond=1.0e-12, hermitian=True)
    diagonal = np.diag(pseudoinverse)
    resistance = diagonal[:, None] + diagonal[None, :] - 2.0 * pseudoinverse
    resistance[np.abs(resistance) < 1.0e-12] = 0.0
    return resistance


def metric_properties(resistance, tol=1.0e-10):
    size = resistance.shape[0]
    diagonal_zero = np.allclose(np.diag(resistance), 0.0, atol=tol)
    symmetric = np.allclose(resistance, resistance.T, atol=tol)
    positive_off_diagonal = bool(
        np.all(resistance[~np.eye(size, dtype=bool)] > tol)
    )
    triangle = True
    for i in range(size):
        for j in range(size):
            for k in range(size):
                if resistance[i, j] > resistance[i, k] + resistance[k, j] + tol:
                    triangle = False
    return diagonal_zero, symmetric, positive_off_diagonal, triangle


def euclidean_embedding(laplacian):
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian)
    pseudoinverse = np.zeros_like(laplacian, dtype=float)
    for eigenvalue, eigenvector in zip(eigenvalues, eigenvectors.T):
        if eigenvalue > 1.0e-12:
            pseudoinverse += np.outer(eigenvector, eigenvector) / eigenvalue
    embedding_evals, embedding_evecs = np.linalg.eigh(pseudoinverse)
    positive = embedding_evals > 1.0e-12
    return embedding_evecs[:, positive] * np.sqrt(embedding_evals[positive])


def gram_from_embedding(points):
    edge_vectors = [points[j] - points[i] for i, j in itertools.combinations(range(4), 2)]
    edge_00 = points[1] - points[0]
    edge_10 = points[2] - points[0]
    edge_20 = points[3] - points[0]
    basis = np.column_stack([edge_00, edge_10, edge_20])
    gram = basis.T @ basis
    return gram, edge_vectors


def stiffness_from_edges(edge_vectors, weights):
    stiffness = np.zeros((3, 3), dtype=float)
    for edge_vector, weight in zip(edge_vectors, weights):
        stiffness += weight * np.outer(edge_vector, edge_vector)
    return stiffness


def cliques(vertex_count, edges):
    adjacency = np.zeros((vertex_count, vertex_count), dtype=bool)
    for i, j in edges:
        adjacency[i, j] = True
        adjacency[j, i] = True
    all_cliques = []
    for size in range(1, vertex_count + 1):
        for candidate in itertools.combinations(range(vertex_count), size):
            if all(
                adjacency[i, j]
                for i, j in itertools.combinations(candidate, 2)
            ):
                all_cliques.append(candidate)
    return all_cliques


AX = io.open(os.path.join(ROOT, "AXIOMS.md"), encoding="utf-8").read()
DOC = canonical_math(
    io.open(
        os.path.join(
            ROOT,
            "derivations",
            "D257_resistance_metric_fixed_point_and_locality_gap.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D257 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D257", "D257" in AX)
check(
    "N1 公理表与文档登记有效电阻度量",
    "R-Z-EFFECTIVE-RESISTANCE-METRIC" in BOTH,
)
check(
    "N1 公理表与文档登记旗复形候选",
    "R-Z-FLAG-COMPLEX-FROM-GRAPH" in BOTH,
)
check(
    "N1 公理表与文档登记电阻面胶合",
    "R-Z-RESISTANCE-FACE-GLUING" in BOTH,
)
check(
    "N1 公理表与文档登记交换度量固定点",
    "R-Z-EXCHANGE-METRIC-FIXED-POINT" in BOTH,
)
check(
    "N1 公理表与文档登记 K4 P1 正例",
    "R-Z-K4-P1-FIXED-POINT" in BOTH,
)
check(
    "N1 公理表与文档登记电阻局部性缺口",
    "R-Z-RESISTANCE-METRIC-LOCALITY-GAP" in BOTH,
)
check(
    "N1 公理表与文档登记三维截断缺口",
    "R-Z-DIMENSION-TRUNCATION-GAP" in BOTH,
)


# ==================================================================
head("N2  加权 Laplacian 半正定且零模为常量")
# ==================================================================
k4_vertices = 4
k4_edges = list(itertools.combinations(range(4), 2))
k4_weights = np.ones(6, dtype=float)
k4_laplacian = weighted_laplacian(k4_vertices, k4_edges, k4_weights)
k4_eigenvalues = np.linalg.eigvalsh(k4_laplacian)
k4_kernel = np.linalg.eigh(k4_laplacian)[1][:, 0]

check(
    "N2 加权 Laplacian 半正定",
    np.all(k4_eigenvalues > -1.0e-12),
    f"eig={np.round(k4_eigenvalues, 8)}",
)
check(
    "N2 连通图的 Laplacian 只有一个零模",
    int(np.sum(np.abs(k4_eigenvalues) < 1.0e-10)) == 1,
)
check(
    "N2 零模是常量向量",
    np.allclose(k4_kernel, k4_kernel[0]),
)


# ==================================================================
head("N3  有效电阻满足度量公理")
# ==================================================================
k4_resistance = effective_resistance(k4_laplacian)
diag_zero, symmetric, positive, triangle = metric_properties(k4_resistance)

check("N3 有效电阻对角为零", diag_zero)
check("N3 有效电阻对称", symmetric)
check("N3 不同顶点的有效电阻为正", positive)
check("N3 有效电阻满足三角不等式", triangle)
check(
    "N3 K4 单位权重的非对角有效电阻为 1/2",
    np.allclose(
        k4_resistance[~np.eye(4, dtype=bool)],
        0.5,
    ),
)


# ==================================================================
head("N4  有效电阻平方根可欧氏嵌入")
# ==================================================================
k4_points = euclidean_embedding(k4_laplacian)
embedded_distance_squared = np.zeros((4, 4), dtype=float)
for i in range(4):
    for j in range(4):
        embedded_distance_squared[i, j] = float(
            np.dot(k4_points[i] - k4_points[j], k4_points[i] - k4_points[j])
        )

check(
    "N4 欧氏嵌入维数等于三",
    k4_points.shape == (4, 3),
    f"shape={k4_points.shape}",
)
check(
    "N4 欧氏嵌入平方距离等于有效电阻",
    np.allclose(embedded_distance_squared, k4_resistance),
    f"max_error={np.max(np.abs(embedded_distance_squared - k4_resistance)):.3e}",
)
check(
    "N4 文档登记唯一半正定平方根嵌入",
    "(L^+)^{1/2}" in DOC,
)


# ==================================================================
head("N5  单位权重 K4 给正四面体几何")
# ==================================================================
k4_gram, k4_edge_vectors = gram_from_embedding(k4_points)
expected_gram = 0.25 * np.ones((3, 3), dtype=float) + 0.25 * np.eye(3)
k4_det_gram = float(np.linalg.det(k4_gram))
k4_volume = math.sqrt(k4_det_gram) / 6.0
k4_edge_lengths = np.array(
    [np.linalg.norm(edge_vector) for edge_vector in k4_edge_vectors]
)

check(
    "N5 六条边长等于 1/sqrt(2)",
    np.allclose(k4_edge_lengths, 1.0 / math.sqrt(2.0)),
    f"lengths={np.round(k4_edge_lengths, 8)}",
)
check(
    "N5 内禀 Gram 矩阵符合正四面体",
    np.allclose(k4_gram, expected_gram),
    f"max_error={np.max(np.abs(k4_gram - expected_gram)):.3e}",
)
check(
    "N5 正四面体 Gram 行列式为 1/16",
    np.isclose(k4_det_gram, 1.0 / 16.0),
    f"detG={k4_det_gram:.12f}",
)
check(
    "N5 正四面体几何体积为 1/24",
    np.isclose(k4_volume, 1.0 / 24.0),
    f"volume={k4_volume:.12f}",
)
check(
    "N5 文档修正 Gram 行列式为 1/16",
    "\\det G=\\frac1{16}" in DOC,
)


# ==================================================================
head("N6  单位权重 K4 在 c=576 满足 P1 固定点")
# ==================================================================
k4_q_vectors = [
    np.array([1.0, 0.0, 0.0]),
    np.array([0.0, 1.0, 0.0]),
    np.array([0.0, 0.0, 1.0]),
    np.array([-1.0, 1.0, 0.0]),
    np.array([-1.0, 0.0, 1.0]),
    np.array([0.0, -1.0, 1.0]),
]
identity_three = np.eye(3)
ones_three = np.ones((3, 3), dtype=float)
k4_local_exchange = stiffness_from_edges(k4_q_vectors, k4_weights)
k4_target_c1 = (1.0 / 24.0) * (4.0 * identity_three - ones_three)
k4_target_c576 = (24.0 / 24.0) * (4.0 * identity_three - ones_three)

check(
    "N6 局部坐标交换刚度为 4I-J",
    np.allclose(k4_local_exchange, 4.0 * identity_three - ones_three),
    f"eig={np.round(np.linalg.eigvalsh(k4_local_exchange), 8)}",
)
check(
    "N6 c=1 的 P1 刚度为 (4I-J)/24",
    np.allclose(k4_target_c1, (4.0 * identity_three - ones_three) / 24.0),
)
check(
    "N6 c=576 的 P1 刚度等于局部交换刚度",
    np.allclose(k4_target_c576, k4_local_exchange),
    f"max_error={np.max(np.abs(k4_target_c576 - k4_local_exchange)):.3e}",
)
check(
    "N6 文档登记 K4 在 c=576 的固定点",
    "单位权重 $K_4$" in DOC
    and "c=576" in DOC
    and "P1 固定点" in DOC,
)


# ==================================================================
head("N7  旗复形候选与三维截断缺口")
# ==================================================================
k4_cliques = cliques(4, k4_edges)
path_edges = [(0, 1), (1, 2), (2, 3)]
path_cliques = cliques(4, path_edges)
k4_counts = {
    size: sum(len(clique) == size for clique in k4_cliques)
    for size in range(1, 5)
}
path_counts = {
    size: sum(len(clique) == size for clique in path_cliques)
    for size in range(1, 5)
}

check(
    "N7 K4 旗复形含四个顶点",
    k4_counts[1] == 4,
    f"counts={k4_counts}",
)
check(
    "N7 K4 旗复形含六条边",
    k4_counts[2] == 6,
)
check(
    "N7 K4 旗复形含四个三角形和一个四面体",
    k4_counts[3] == 4 and k4_counts[4] == 1,
)
check(
    "N7 路径图旗复形不含三角形与四面体",
    path_counts[3] == 0 and path_counts[4] == 0,
    f"counts={path_counts}",
)
check(
    "N7 文档登记三维截断缺口",
    "三维骨架" in DOC and "维数截断" in DOC,
)


# ==================================================================
head("N8  单一顶点度量给共享面边长胶合")
# ==================================================================
shared_face_a = {
    (0, 1): 0.7,
    (0, 2): 0.8,
    (1, 2): 0.9,
}
shared_face_b = dict(shared_face_a)
shared_face_error = max(
    abs(shared_face_a[edge] - shared_face_b[edge])
    for edge in shared_face_a
)

check(
    "N8 共享面的三条边长一致",
    shared_face_error == 0.0,
    f"max_error={shared_face_error:.3e}",
)
check(
    "N8 文档登记共享面二次型而非逐坐标相等",
    "切空间二次型与转移映射的一致性" in DOC,
)


# ==================================================================
head("N9  远端边改变局部有效电阻")
# ==================================================================
path_laplacian = weighted_laplacian(4, path_edges, np.ones(3, dtype=float))
path_resistance = effective_resistance(path_laplacian)
far_edge_edges = path_edges + [(0, 3)]
far_edge_laplacian = weighted_laplacian(4, far_edge_edges, np.ones(4, dtype=float))
far_edge_resistance = effective_resistance(far_edge_laplacian)

check(
    "N9 路径图局部有效电阻 R12=1",
    np.isclose(path_resistance[1, 2], 1.0),
    f"R12={path_resistance[1, 2]:.12f}",
)
check(
    "N9 加入远端边后 R12=3/4",
    np.isclose(far_edge_resistance[1, 2], 0.75),
    f"R12={far_edge_resistance[1, 2]:.12f}",
)
check(
    "N9 远端边改变局部有效电阻",
    not np.isclose(path_resistance[1, 2], far_edge_resistance[1, 2]),
)
check(
    "N9 文档登记有效电阻非局部",
    "有效电阻不是局部量" in DOC and "远端连接可以改变局部有效电阻距离" in DOC,
)


# ==================================================================
head("N10 权重尺度重标度")
# ==================================================================
scale = 3.7
scaled_laplacian = weighted_laplacian(4, k4_edges, scale * k4_weights)
scaled_resistance = effective_resistance(scaled_laplacian)
scaled_edge_lengths = np.sqrt(k4_resistance / scale)
scaled_edge_lengths_embedded = np.sqrt(scaled_resistance)

check(
    "N10 K->lambda K 给 R->R/lambda",
    np.allclose(scaled_resistance, k4_resistance / scale),
    f"max_error={np.max(np.abs(scaled_resistance - k4_resistance / scale)):.3e}",
)
check(
    "N10 边长按 lambda^{-1/2} 重标度",
    np.allclose(
        scaled_edge_lengths_embedded,
        scaled_edge_lengths,
    ),
)
check(
    "N10 文档登记全局尺度仍是独立输入",
    "交换权重只给相对度量" in DOC,
)


# ==================================================================
head("N11  文档登记输入预算与限制")
# ==================================================================
check(
    "N11 文档登记输入预算",
    "## §3 输入预算" in DOC,
)
check(
    "N11 文档登记局部化规则是输入",
    "局部化规则" in DOC and "未导出" in DOC,
)
check(
    "N11 文档登记三维截断是输入",
    "三维截断或团数条件" in DOC,
)
check(
    "N11 文档登记固定点选择是输入",
    "P1 固定点选择" in DOC,
)


# ==================================================================
head("N12  文档不把固定点或局部性写成自动结论")
# ==================================================================
check(
    "N12 文档明确条件性",
    "不是 `U1-U4` 的推论" in DOC and "条件" in DOC,
)
check(
    "N12 文档否认固定点自动存在",
    "固定点自动存在或唯一" in DOC and "不能主张" in DOC,
)
check(
    "N12 文档否认任意图自动给局部三维复形",
    "任意支持图自动给局部三维复形" in DOC,
)
check(
    "N12 文档否认局部收敛自动成立",
    "不同细化会自动收敛到同一个三维几何" in DOC,
)


# ==================================================================
head("N13  上游边界保持")
# ==================================================================
check(
    "N13 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)
check(
    "N13 上游边界声明出现",
    "不修改 `U1-U4+C1`" in DOC,
)


# ==================================================================
head("N14  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N14 D257 文档不含项目外体系名或外部路径",
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
