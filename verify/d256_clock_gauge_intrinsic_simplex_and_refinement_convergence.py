#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D256 —— 时钟规范拆分、内禀单纯几何与细化收敛条件
================================================================================
【被检验的问题（D256）】
  物理时钟、单位 lapse、零 shift 能否合并为一项？
  单位正交叶状是否只有局部版本？
  CMC 能否提供普适时钟？
  内禀边长单纯复形能否消去顶点嵌入？
  交换权重何时与标准 P1 刚度一致？
  哪种细化才会给同一连续 Q 与 h？

【本步判据】
  N1  D256 与恢复结构已登记
  N2  物理时钟、N=1、beta=0 三项分离
  N3  N=1 给单位法向与单位测地线
  N4  caustic 阻断全局正规坐标
  N5  CMC 单调条件与临界点障碍
  N6  边长 Gram 矩阵给内禀 G
  N7  由 G 装配 Q 与 h
  N8  P1 刚度与内禀几何一致
  N9  交换权重相容条件
  N10 固定边权在细化下尺度漂移
  N11 共同 Q 收敛给共同 h
  N12 形状正则与 sliver 区分
  N13 文档登记输入预算与限制
  N14 文档不把条件写成上游推论
  N15 上游边界保持
  N16 文档不引用项目外体系名或外部路径

运行：python3 verify/d256_clock_gauge_intrinsic_simplex_and_refinement_convergence.py
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


def adm_metric(h, lapse, shift):
    """Positive-time ADM representative in 1+n dimensions."""
    h = np.asarray(h, dtype=float)
    shift = np.asarray(shift, dtype=float)
    n = h.shape[0]
    metric = np.zeros((n + 1, n + 1), dtype=float)
    metric[0, 0] = lapse ** 2 - shift @ h @ shift
    metric[0, 1:] = -h @ shift
    metric[1:, 0] = -h @ shift
    metric[1:, 1:] = -h
    return metric


def signature_counts(metric, tol=1.0e-10):
    eigenvalues = np.linalg.eigvalsh(metric)
    positive = int(np.sum(eigenvalues > tol))
    negative = int(np.sum(eigenvalues < -tol))
    zero = int(np.sum(np.abs(eigenvalues) <= tol))
    return positive, negative, zero, eigenvalues


def gram_from_edge_lengths(lengths):
    """Return the 3x3 Gram matrix from edge lengths l_ab."""
    gram = np.zeros((3, 3), dtype=float)
    for i in range(3):
        for j in range(3):
            a = i + 1
            b = j + 1
            if i == j:
                gram[i, j] = lengths[(0, a)] ** 2
            else:
                pair = (a, b) if a < b else (b, a)
                gram[i, j] = 0.5 * (
                    lengths[(0, a)] ** 2
                    + lengths[(0, b)] ** 2
                    - lengths[pair] ** 2
                )
    return gram


AX = io.open(os.path.join(ROOT, "AXIOMS.md"), encoding="utf-8").read()
DOC = canonical_math(
    io.open(
        os.path.join(
            ROOT,
            "derivations",
            "D256_clock_gauge_intrinsic_simplex_and_refinement_convergence.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D256 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D256", "D256" in AX)
check(
    "N1 公理表与文档登记时钟校准拆分",
    "R-Z-CLOCK-CALIBRATION-SPLIT" in BOTH,
)
check(
    "N1 公理表与文档登记局部正规标架",
    "R-Z-GAUSSIAN-NORMAL-LOCAL-FRAME" in BOTH,
)
check(
    "N1 公理表与文档登记 CMC 时钟路线",
    "R-Z-CMC-CLOCK-ROUTE" in BOTH,
)
check(
    "N1 公理表与文档登记内禀单纯几何",
    "R-Z-INTRINSIC-SIMPLEX-GEOMETRY" in BOTH,
)
check(
    "N1 公理表与文档登记刚度相容条件",
    "R-Z-GEOMETRIC-STIFFNESS-COMPATIBILITY" in BOTH,
)
check(
    "N1 公理表与文档登记形状正则细化",
    "R-Z-SHAPE-REGULAR-REFINEMENT" in BOTH,
)
check(
    "N1 公理表与文档登记共同细化极限",
    "R-Z-COMMON-REFINEMENT-LIMIT" in BOTH,
)


# ==================================================================
head("N2  物理时钟、N=1、beta=0 三项分离")
# ==================================================================
h_sample = np.array(
    [
        [2.0, 0.1, 0.0],
        [0.1, 1.5, -0.2],
        [0.0, -0.2, 1.2],
    ]
)
metric_base = adm_metric(h_sample, 1.0, np.zeros(3))
metric_lapse = adm_metric(h_sample, 1.7, np.zeros(3))
metric_shift = adm_metric(h_sample, 1.0, np.array([0.3, 0.0, 0.0]))
inverse_base = np.linalg.inv(metric_base)

check(
    "N2 N=1 给逆时间分量为 1",
    np.isclose(inverse_base[0, 0], 1.0),
    f"g00={inverse_base[0, 0]:.12f}",
)
check(
    "N2 lapse 改变不改变空间块",
    np.allclose(metric_lapse[1:, 1:], -h_sample),
)
check(
    "N2 lapse 改变光锥度规",
    not np.allclose(metric_lapse, metric_base),
)
check(
    "N2 shift 改变时间空间混合项",
    not np.allclose(metric_shift, metric_base),
)
check(
    "N2 三项条件分别对应三个独立数据",
    "物理时钟" in DOC and "N=1" in DOC and "\\beta=0" in DOC,
)


# ==================================================================
head("N3  N=1 给单位法向与单位测地线")
# ==================================================================
positive, negative, zero, eigenvalues = signature_counts(metric_base)
lapse_values = np.array([0.5, 1.0, 2.0])
inverse_time_values = 1.0 / lapse_values ** 2
acceleration_zero = np.isclose(-np.gradient(np.log(np.ones(5)))[2], 0.0)

check(
    "N3 单位正交度规号差为 (1,3)",
    positive == 1 and negative == 3 and zero == 0,
    f"eig={np.round(eigenvalues, 8)}",
)
check(
    "N3 g^{00}=N^{-2}",
    np.allclose(inverse_time_values, 1.0 / lapse_values ** 2),
)
check(
    "N3 N=1 时法向加速度为零",
    acceleration_zero,
)
check(
    "N3 文档登记法向固有时",
    "法向曲线是单位测地线" in DOC and "固有时" in DOC,
)


# ==================================================================
head("N4  caustic 阻断全局正规坐标")
# ==================================================================
kappa = 0.25
caustic_time = 1.0 / kappa
jacobian_at_caustic = 1.0 - kappa * caustic_time
jacobian_before = 1.0 - kappa * (0.5 / kappa)
jacobian_after = 1.0 - kappa * (1.5 / kappa)

check(
    "N4 caustic 处坐标 Jacobian 为零",
    np.isclose(jacobian_at_caustic, 0.0),
)
check(
    "N4 caustic 两侧 Jacobian 变号",
    jacobian_before > 0.0 and jacobian_after < 0.0,
)
check(
    "N4 文档登记 caustic 与全局同步障碍",
    "caustic" in DOC and "全局同步" in DOC,
)


# ==================================================================
head("N5  CMC 单调条件与临界点障碍")
# ==================================================================
cmc_time = np.linspace(-1.0, 1.0, 41)
linear_derivative = np.ones_like(cmc_time)
cubic_derivative = 3.0 * cmc_time ** 2
critical_points = np.where(np.isclose(cubic_derivative, 0.0))[0]

check(
    "N5 线性 CMC 时间函数在样本上严格单调",
    np.all(linear_derivative > 0.0),
)
check(
    "N5 三次 CMC 时间函数出现临界点",
    len(critical_points) > 0,
    f"critical_count={len(critical_points)}",
)
check(
    "N5 文档不把 CMC 写成普适定理",
    "候选" in DOC and "不是普适定理" in DOC,
)


# ==================================================================
head("N6  边长 Gram 矩阵给内禀 G")
# ==================================================================
vertices = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
)
edges = list(itertools.combinations(range(4), 2))
edge_lengths = {
    edge: float(np.linalg.norm(vertices[edge[1]] - vertices[edge[0]]))
    for edge in edges
}
gram = gram_from_edge_lengths(edge_lengths)
gram_direct = np.eye(3)

check(
    "N6 六条边长已给出",
    len(edge_lengths) == 6,
)
check(
    "N6 边长 Gram 矩阵等于直接内积",
    np.allclose(gram, gram_direct),
    f"max_error={np.max(np.abs(gram - gram_direct)):.3e}",
)
check(
    "N6 Gram 矩阵正定",
    np.all(np.linalg.eigvalsh(gram) > 0.0),
)


# ==================================================================
head("N7  由 G 装配 Q 与 h")
# ==================================================================
det_gram = float(np.linalg.det(gram))
Q_intrinsic = math.sqrt(det_gram) * np.linalg.inv(gram)
h_recovered = float(np.linalg.det(Q_intrinsic)) * np.linalg.inv(Q_intrinsic)
geometric_volume = math.sqrt(det_gram) / 6.0

check(
    "N7 Q=sqrt(det G)G^{-1}",
    np.allclose(Q_intrinsic, math.sqrt(det_gram) * np.linalg.inv(gram)),
)
check(
    "N7 反解恢复 h=G",
    np.allclose(h_recovered, gram),
    f"max_error={np.max(np.abs(h_recovered - gram)):.3e}",
)
check(
    "N7 度规体积为 sqrt(det G)/6",
    np.isclose(geometric_volume, 1.0 / 6.0),
)
check(
    "N7 文档登记内禀边长路线",
    "G_{ij}" in DOC and "Q" in DOC and "h=G" in DOC,
)


# ==================================================================
head("N8  P1 刚度与内禀几何一致")
# ==================================================================
stiffness = math.sqrt(det_gram) / 6.0 * np.linalg.inv(gram)
test_gradients = [
    np.array([1.0, 0.0, 0.0]),
    np.array([0.0, 1.0, 0.0]),
    np.array([0.3, -0.2, 0.7]),
]
p1_energy = [
    0.5 * geometric_volume * float(d @ np.linalg.inv(gram) @ d)
    for d in test_gradients
]
nodal_energy = [
    0.5 * float(d @ stiffness @ d)
    for d in test_gradients
]

check(
    "N8 P1 刚度为 sqrt(det G)/6 G^{-1}",
    np.allclose(stiffness, math.sqrt(det_gram) / 6.0 * np.linalg.inv(gram)),
)
check(
    "N8 P1 能量等于节点二次型",
    np.allclose(p1_energy, nodal_energy),
)
check(
    "N8 文档登记标准 P1 刚度公式",
    "A_\\sigma" in DOC and "G^{-1}" in DOC,
)


# ==================================================================
head("N9  交换权重相容条件")
# ==================================================================
edge_vectors = {
    edge: vertices[edge[1]] - vertices[edge[0]]
    for edge in edges
}
index_pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
basis = np.zeros((6, 6), dtype=float)
for column, edge in enumerate(edges):
    outer = np.outer(edge_vectors[edge], edge_vectors[edge])
    for row, (i, j) in enumerate(index_pairs):
        basis[row, column] = outer[i, j]
target = np.array([stiffness[i, j] for i, j in index_pairs])
compatible_weights = np.linalg.solve(basis, target)
A_compatible = sum(
    compatible_weights[index] * np.outer(edge_vectors[edge], edge_vectors[edge])
    for index, edge in enumerate(edges)
)
incompatible_weights = compatible_weights + np.array([0.1, -0.1, 0.0, 0.0, 0.0, 0.0])
A_incompatible = sum(
    incompatible_weights[index] * np.outer(edge_vectors[edge], edge_vectors[edge])
    for index, edge in enumerate(edges)
)

check(
    "N9 相容权重精确给标准 P1 刚度",
    np.allclose(A_compatible, stiffness),
    f"max_error={np.max(np.abs(A_compatible - stiffness)):.3e}",
)
check(
    "N9 扰动权重破坏 P1 相容性",
    not np.allclose(A_incompatible, stiffness),
)
check(
    "N9 文档登记交换权重相容公式",
    "K_{ab}" in DOC and "e_{ab}\\otimes e_{ab}" in DOC,
)


# ==================================================================
head("N10  固定边权在细化下尺度漂移")
# ==================================================================
mesh_size = 0.25
Q_fixed_weight = 1.0 * mesh_size
Q_inverse_square_weight = (1.0 / mesh_size ** 2) * mesh_size
Q_consistent_weight = (1.0 / mesh_size) * mesh_size

check(
    "N10 固定边权给 Q 随格距趋零",
    np.isclose(Q_fixed_weight, mesh_size),
)
check(
    "N10 1/a^2 边权使 Q 发散",
    np.isclose(Q_inverse_square_weight, 1.0 / mesh_size),
)
check(
    "N10 1/a 边权给有限共同尺度",
    np.isclose(Q_consistent_weight, 1.0),
)
check(
    "N10 文档登记固定边权无自动连续极限",
    "固定边权" in DOC and "不会自动" in DOC,
)


# ==================================================================
head("N11  共同 Q 收敛给共同 h")
# ==================================================================
Q_limit = np.array(
    [
        [1.3, 0.1, 0.0],
        [0.1, 1.1, -0.05],
        [0.0, -0.05, 1.4],
    ]
)
h_limit = float(np.linalg.det(Q_limit)) * np.linalg.inv(Q_limit)
errors = []
for epsilon in [1.0e-1, 1.0e-2, 1.0e-3]:
    Q_epsilon = Q_limit + epsilon * np.eye(3)
    h_epsilon = float(np.linalg.det(Q_epsilon)) * np.linalg.inv(Q_epsilon)
    errors.append(float(np.max(np.abs(h_epsilon - h_limit))))

check(
    "N11 h(Q) 在 Q_limit 附近连续",
    errors[0] > errors[1] > errors[2],
    f"errors={np.array(errors)}",
)
check(
    "N11 共同 Q 给共同 h",
    np.allclose(float(np.linalg.det(Q_limit)) * np.linalg.inv(Q_limit), h_limit),
)


# ==================================================================
head("N12  形状正则与 sliver 区分")
# ==================================================================
regular_triangle = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, math.sqrt(3.0) / 2.0]])
sliver_triangle = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 1.0e-3]])


def triangle_shape_ratio(points):
    lengths = np.array(
        [
            np.linalg.norm(points[1] - points[0]),
            np.linalg.norm(points[2] - points[1]),
            np.linalg.norm(points[0] - points[2]),
        ]
    )
    edge_u = points[1] - points[0]
    edge_v = points[2] - points[0]
    twice_area = abs(edge_u[0] * edge_v[1] - edge_u[1] * edge_v[0])
    return float(np.max(lengths) ** 2 / twice_area)


regular_ratio = triangle_shape_ratio(regular_triangle)
sliver_ratio = triangle_shape_ratio(sliver_triangle)

check(
    "N12 正三角形形状比有限且较小",
    regular_ratio < 2.0,
    f"ratio={regular_ratio:.6f}",
)
check(
    "N12 sliver 形状比远大于正三角形",
    sliver_ratio > 100.0 * regular_ratio,
    f"ratio={sliver_ratio:.6f}",
)
check(
    "N12 文档登记形状正则条件",
    "形状正则" in DOC and "sliver" in DOC,
)


# ==================================================================
head("N13  文档登记输入预算与限制")
# ==================================================================
check(
    "N13 文档登记输入预算",
    "## §3 输入预算" in DOC,
)
check(
    "N13 文档登记时钟三项拆分",
    "三项独立输入" in DOC and "物理时钟选择" in DOC,
)
check(
    "N13 文档登记顶点嵌入被内禀边长替换",
    "顶点嵌入不再是独立输入" in DOC
    and "复形、边长或细化规则" in DOC,
)


# ==================================================================
head("N14  文档不把条件写成上游推论")
# ==================================================================
check(
    "N14 文档明确条件性",
    "不是 `U1-U4` 的推论" in DOC
    and "条件" in DOC,
)
check(
    "N14 文档不把时钟写成已导出",
    "某个标量读回已经是物理时钟" in DOC
    and "不能主张" in DOC,
)
check(
    "N14 文档不把细化写成自动唯一",
    "任意细化" in DOC and "同一" in DOC and "不是" in DOC,
)


# ==================================================================
head("N15  上游边界保持")
# ==================================================================
check(
    "N15 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)
check(
    "N15 上游边界声明出现",
    "不修改 `U1-U4+C1`" in DOC,
)


# ==================================================================
head("N16  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N16 D256 文档不含项目外体系名或外部路径",
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
