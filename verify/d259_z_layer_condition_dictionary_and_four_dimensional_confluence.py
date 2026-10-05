#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D259 —— Z 层条件字典与四维时空的汇流构造
================================================================================
【被检验的问题（D259）】
  Z 结果类是否可以和历史实现分离？
  五通道零和约束为什么给四个局域方向？
  四通道是否只给三个方向？
  完整交换型是否给各向同性正定型？
  正定型加一条无向时间线是否给 (1,3) 号差？
  跨区域不同历史何时能汇流为同一度规？
  零和置换对称是否能自然选出时间线？

【本步判据】
  N1  D259 与恢复结构已登记
  N2  零和通道秩证书
  N3  完整交换的各向同性投影
  N4  时间线提升给四维洛伦兹号差
  N5  体积尺度的条件归一
  N6  跨区域汇流与转移映射
  N7  零和置换对称不选择唯一时间线
  N8  文档登记条件字典输入与缺口
  N9  文档不把 Z 无时间性写成四维自动导出
  N10 上游边界保持
  N11 文档不引用项目外体系名或外部路径

运行：python3 verify/d259_z_layer_condition_dictionary_and_four_dimensional_confluence.py
"""
import io
import itertools
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


def zero_sum_rank(channel_count):
    constraint = np.ones((1, channel_count), dtype=float)
    singular_values = np.linalg.svd(constraint, compute_uv=False)
    return channel_count - int(np.sum(singular_values > 1.0e-12))


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
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian)
    positive = eigenvalues > 1.0e-10
    return eigenvectors[:, positive] / np.sqrt(eigenvalues[positive])[None, :]


def full_exchange_stiffness(points, edges, weights):
    dimension = points.shape[1]
    stiffness = np.zeros((dimension, dimension), dtype=float)
    for edge, weight in zip(edges, weights):
        i, j = edge
        edge_vector = points[j] - points[i]
        stiffness += weight * np.outer(edge_vector, edge_vector)
    return stiffness


def signature_counts(matrix):
    eigenvalues = np.linalg.eigvalsh(matrix)
    positive = int(np.sum(eigenvalues > 1.0e-10))
    negative = int(np.sum(eigenvalues < -1.0e-10))
    zero = int(np.sum(np.abs(eigenvalues) <= 1.0e-10))
    return positive, negative, zero


AX = io.open(os.path.join(ROOT, "AXIOMS.md"), encoding="utf-8").read()
DOC = canonical_math(
    io.open(
        os.path.join(
            ROOT,
            "derivations",
            "D259_z_layer_condition_dictionary_and_four_dimensional_confluence.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D259 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D259", "D259" in AX)
check(
    "N1 公理表与文档登记结果历史分离",
    "R-Z-INVARIANT-HISTORY-SPLIT" in BOTH,
)
check(
    "N1 公理表与文档登记五通道秩证书",
    "R-Z-FIVE-CHANNEL-RANK-CERTIFICATE" in BOTH,
)
check(
    "N1 公理表与文档登记时间线证书",
    "R-Z-CLOSURE-TIME-LINE-CERTIFICATE" in BOTH,
)
check(
    "N1 公理表与文档登记跨区域汇流",
    "R-Z-CROSS-REGION-CONFLUENCE" in BOTH,
)
check(
    "N1 公理表与文档登记条件四维图册",
    "R-Z-CONDITIONAL-4D-ATLAS" in BOTH,
)
check(
    "N1 公理表与文档登记字典来源缺口",
    "R-Z-DICTIONARY-SOURCE-GAP" in BOTH,
)


# ==================================================================
head("N2  零和通道秩证书")
# ==================================================================
rank_two = zero_sum_rank(2)
rank_three = zero_sum_rank(3)
rank_four = zero_sum_rank(4)
rank_five = zero_sum_rank(5)
rank_six = zero_sum_rank(6)

check("N2 两通道零和秩为一", rank_two == 1, f"rank={rank_two}")
check("N2 三通道零和秩为二", rank_three == 2, f"rank={rank_three}")
check("N2 四通道零和秩为三", rank_four == 3, f"rank={rank_four}")
check("N2 五通道零和秩为四", rank_five == 4, f"rank={rank_five}")
check("N2 六通道零和秩为五", rank_six == 5, f"rank={rank_six}")
check(
    "N2 四通道不足四维而五通道达到四维",
    rank_four == 3 and rank_five == 4,
)
check(
    "N2 最小四维零和证书为五通道",
    min(m for m in range(2, 8) if zero_sum_rank(m) >= 4) == 5,
)


# ==================================================================
head("N3  完整交换的各向同性投影")
# ==================================================================
rng = np.random.default_rng(259)
edges_five = list(itertools.combinations(range(5), 2))
weights_five = rng.uniform(0.5, 1.5, size=len(edges_five))
laplacian_five = weighted_laplacian(5, edges_five, weights_five)
points_five = resistance_embedding(laplacian_five)
stiffness_five = full_exchange_stiffness(points_five, edges_five, weights_five)

check(
    "N3 五通道完整交换嵌入维数为四",
    points_five.shape == (5, 4),
    f"shape={points_five.shape}",
)
check(
    "N3 五通道完整交换刚度等于 I4",
    np.allclose(stiffness_five, np.eye(4)),
    f"max_error={np.max(np.abs(stiffness_five - np.eye(4))):.3e}",
)

sample = rng.normal(size=(8, 5))
sample -= np.mean(sample, axis=1, keepdims=True)
exchange_quadratic = np.array(
    [
        sum((row[j] - row[i]) ** 2 for i in range(5) for j in range(i + 1, 5))
        for row in sample
    ]
)
direct_quadratic = 5.0 * np.sum(sample**2, axis=1)
check(
    "N3 等权完整交换恒等式 Q=5 sum p_i^2",
    np.allclose(exchange_quadratic, direct_quadratic),
    f"max_error={np.max(np.abs(exchange_quadratic - direct_quadratic)):.3e}",
)
check(
    "N3 各向同性矩阵的四个特征值相等",
    np.allclose(np.linalg.eigvalsh(stiffness_five), np.ones(4)),
)


# ==================================================================
head("N4  时间线提升给四维洛伦兹号差")
# ==================================================================
g_positive = np.eye(4)
time_like = np.array([1.0, 0.0, 0.0, 0.0])
g_lorentz = 2.0 * np.outer(time_like, time_like) - g_positive
positive, negative, zero = signature_counts(g_lorentz)

check(
    "N4 时间线提升给对角 (1,-1,-1,-1)",
    np.allclose(g_lorentz, np.diag([1.0, -1.0, -1.0, -1.0])),
    f"g={np.diag(g_lorentz)}",
)
check(
    "N4 号差为一正三负",
    (positive, negative, zero) == (1, 3, 0),
    f"signature=({positive},{negative},{zero})",
)
check("N4 四维洛伦兹体积行列式为负", np.linalg.det(g_lorentz) < 0.0)
check(
    "N4 时间方向保持类时正号",
    np.isclose(time_like @ g_lorentz @ time_like, 1.0),
)


# ==================================================================
head("N5  体积尺度的条件归一")
# ==================================================================
g_reference = np.diag([1.0, 1.0, 1.0, 1.0])
volume_reference = np.sqrt(abs(np.linalg.det(g_reference)))
target_volume = 3.7
conformal_factor = (target_volume / volume_reference) ** 0.25
g_scaled = (conformal_factor**2) * g_reference
volume_scaled = np.sqrt(abs(np.linalg.det(g_scaled)))

check(
    "N5 四维体积按 Omega^4 缩放",
    np.isclose(volume_scaled, target_volume),
    f"volume={volume_scaled:.12f}",
)
check(
    "N5 体积归一后共形类未变",
    np.allclose(g_scaled / g_scaled[0, 0], g_reference),
)
check(
    "N5 时间线归一后提升仍给洛伦兹号差",
    signature_counts(
        2.0 * np.outer(time_like, time_like) - g_scaled
    )
    == (1, 3, 0),
)


# ==================================================================
head("N6  跨区域汇流与转移映射")
# ==================================================================
rotation = np.array(
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -1.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
)
good_metric = rotation.T @ g_lorentz @ rotation
bad_metric = np.diag([1.7, -1.0, -1.0, -1.0])

check(
    "N6 保时线和体积的空间转动保持度规",
    np.allclose(good_metric, g_lorentz),
    f"max_error={np.max(np.abs(good_metric - g_lorentz)):.3e}",
)
check(
    "N6 保时线的空间转动保持时间线",
    np.allclose(rotation @ time_like, time_like),
)
check(
    "N6 错误局部时间标度破坏重叠度规",
    not np.allclose(bad_metric, g_lorentz),
)
check(
    "N6 错误标度改变时间方向 norm",
    not np.isclose(time_like @ bad_metric @ time_like, 1.0),
)


# ==================================================================
head("N7  零和置换对称不选择唯一时间线")
# ==================================================================
permutation = np.array(
    [
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 1.0],
    ]
)
transpositions = []
for i, j in itertools.combinations(range(5), 2):
    matrix = np.eye(5)
    matrix[[i, j], :] = matrix[[j, i], :]
    transpositions.append(matrix)

average_projector = np.mean(
    [permutation_matrix for permutation_matrix in transpositions], axis=0
)
eigenvalues, eigenvectors = np.linalg.eigh(average_projector)
fixed_space = eigenvectors[:, np.abs(eigenvalues - 1.0) < 1.0e-10]
fixed_space_in_zero_sum = fixed_space - np.mean(fixed_space, axis=0, keepdims=True)
candidate_channel_line = np.array([1.0, 0.0, 0.0, 0.0, 0.0])

check(
    "N7 交换通道标签改变局域时间候选人",
    not np.allclose(permutation @ candidate_channel_line, candidate_channel_line),
)
check(
    "N7 全置换固定空间只含均匀方向",
    np.allclose(
        fixed_space @ fixed_space.T,
        np.outer(np.ones(5), np.ones(5)) / 5.0,
    ),
)
check(
    "N7 零和超平面内没有全置换固定时间线",
    np.linalg.norm(fixed_space_in_zero_sum) < 1.0e-10,
)
check(
    "N7 文档登记时间线仍是独立证书",
    "无向时间线是独立闭合证书" in DOC,
)


# ==================================================================
head("N8  文档登记条件字典输入与缺口")
# ==================================================================
check(
    "N8 文档登记输入预算",
    "## §3 输入预算" in DOC,
)
check(
    "N8 文档登记五通道证书未导出",
    "五通道证书本身仍是恢复层输入" in DOC,
)
check(
    "N8 文档登记汇流映射未导出",
    "区域覆盖与转移映射" in DOC and "未导出" in DOC,
)
check(
    "N8 文档登记外边局部化缺口",
    "外部边局部隔离" in DOC and "R-Z-LOCAL-EXCHANGE-BACKREACTION-GAP" in DOC,
)
check(
    "N8 文档登记结果与历史分离",
    "Z 保存闭合结果类，P/E 保存不同历史实现" in DOC,
)


# ==================================================================
head("N9  文档不把 Z 无时间性写成四维自动导出")
# ==================================================================
check(
    "N9 文档明确 Z 无时间性不自动给字典",
    "Z 的无时间性" in DOC and "五通道证书、时间线、尺度或跨区域汇流自动存在" in DOC,
)
check(
    "N9 文档明确 Z 分离于完整几何",
    "Z 层保存的是四维时空的“结果类”，不是构造四维时空的完整字典" in DOC,
)
check(
    "N9 文档明确条件性",
    "它们不是 `U1-U4` 的推论" in DOC and "条件" in DOC,
)
check(
    "N9 文档不声称无条件四维导出",
    "本文不声称四维时空已经无条件导出" in DOC,
)


# ==================================================================
head("N10  上游边界保持")
# ==================================================================
check(
    "N10 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)
check(
    "N10 上游边界声明出现",
    "不修改 `U1-U4+C1`" in DOC,
)


# ==================================================================
head("N11  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N11 D259 文档不含项目外体系名或外部路径",
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
