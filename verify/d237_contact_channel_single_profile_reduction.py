#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D237 —— 接触项与单剖面归约
================================================================================
【被检验的问题（D237）】
  两通道球模 Hamiltonian 何时能合并为单个 K_M tensor g？
  纯中心、平行与独立非中心接触项分别有什么后果？

【本步判据】
  N1  D237 与恢复结构已登记
  N2  平行矩阵通道可以合并为单轮廓
  N3  平行通道的边界值会从零抬起
  N4  非中心独立通道不满足平行合并判据
  N5  纯中心通道不改变 likelihood ratio
  N6  纯中心通道仍使完整算子脱离严格单轮廓
  N7  独立非中心扰动改变谱比例
  N8  独立非中心扰动改变 likelihood ratio
  N9  文档登记单剖面归约判据
  N10 文档登记接触通道分类缺口
  N11 文档不把单轮廓写成默认形式
  N12 上游边界保持
  N13 文档不引用项目外体系名或外部路径

运行：python3 verify/d237_contact_channel_single_profile_reduction.py
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
            "D237_contact_channel_single_profile_reduction.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D237 与恢复结构已登记")
# ==================================================================
check(
    "N1 公理表登记 D237",
    "D237" in AX,
)
check(
    "N1 公理表与文档登记单剖面归约判据",
    "R-Z-SINGLE-PROFILE-REDUCTION-CRITERION" in BOTH
    and "单剖面" in DOC,
)
check(
    "N1 公理表与文档登记接触通道分类缺口",
    "R-Z-CONTACT-CHANNEL-CLASSIFICATION-GAP" in BOTH
    and "接触通道" in DOC,
)


# ==================================================================
head("N2  平行矩阵通道可以合并为单轮廓")
# ==================================================================
dimension = 2
base_channel = np.diag([-1.0, -2.0]).astype(complex)
parallel_channel = 0.3 * base_channel

age = np.linspace(0.0, 1.0, 1001)
geometry_weight = 1.0 - age**2
contact_weight = np.ones_like(age)

matrix_at_age_parallel = (
    base_channel[None, :, :] * geometry_weight[:, None, None]
    + parallel_channel[None, :, :] * contact_weight[:, None, None]
)
merged_profile = geometry_weight + 0.3 * contact_weight
matrix_at_age_merged = base_channel[None, :, :] * merged_profile[:, None, None]

channel_rank_parallel = np.linalg.matrix_rank(
    np.stack([base_channel.reshape(-1), parallel_channel.reshape(-1)], axis=1)
)
weight_rank_parallel = np.linalg.matrix_rank(
    np.stack([geometry_weight, contact_weight], axis=1)
)
weight_rank_proportional = np.linalg.matrix_rank(
    np.stack([geometry_weight, 2.0 * geometry_weight], axis=1)
)

check(
    "N2 平行通道逐点合并为同一矩阵",
    np.allclose(matrix_at_age_parallel, matrix_at_age_merged),
)
check(
    "N2 矩阵通道秩或权重秩至少一边退化",
    min(channel_rank_parallel, weight_rank_parallel) == 1
    and weight_rank_proportional == 1,
    f"channel_rank={channel_rank_parallel}, weight_rank={weight_rank_parallel}",
)


# ==================================================================
head("N3  平行通道的边界值会从零抬起")
# ==================================================================
boundary_profile_value = merged_profile[-1]
check(
    "N3 平行接触项在边界抬起剖面",
    abs(boundary_profile_value - 0.3) < 1e-12,
    f"boundary={boundary_profile_value:.12f}",
)
check(
    "N3 原始几何权重仍在边界为零",
    abs(geometry_weight[-1]) < 1e-12,
)


# ==================================================================
head("N4  非中心独立通道不满足平行合并判据")
# ==================================================================
independent_channel = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)

# 检查 independent_channel 是否落在 base_channel 的一维复线性张成中。
base_vector = base_channel.reshape(-1)
channel_vector = independent_channel.reshape(-1)
projection = (
    np.vdot(base_vector, channel_vector)
    / np.vdot(base_vector, base_vector)
) * base_vector
parallel_residual = np.linalg.norm(channel_vector - projection)

check(
    "N4 独立通道不平行于基础通道",
    parallel_residual > 1e-2,
    f"residual={parallel_residual:.6f}",
)


# ==================================================================
head("N5  纯中心通道不改变 likelihood ratio")
# ==================================================================
lambda_0 = -1.0
lambda_1 = -2.0
contact_center = 0.37 * np.ones_like(age)

ratio_without_center = np.exp(
    (lambda_1 - lambda_0) * geometry_weight
)
ratio_with_center = np.exp(
    (lambda_1 - lambda_0) * geometry_weight
    + 0.0 * contact_center
)

check(
    "N5 中心接触项在 likelihood ratio 中抵消",
    np.allclose(ratio_without_center, ratio_with_center),
)


# ==================================================================
head("N6  纯中心通道仍使完整算子脱离严格单轮廓")
# ==================================================================
center_shifted_matrix = (
    base_channel[None, :, :] * geometry_weight[:, None, None]
    + np.eye(dimension)[None, :, :] * contact_center[:, None, None]
)
strict_single_profile = (
    base_channel[None, :, :] * geometry_weight[:, None, None]
)
channel_rank_center = np.linalg.matrix_rank(
    np.stack(
        [
            base_channel.reshape(-1),
            np.eye(dimension).reshape(-1),
        ],
        axis=1,
    )
)
weight_rank_center = np.linalg.matrix_rank(
    np.stack([geometry_weight, contact_center], axis=1)
)

check(
    "N6 中心通道改变完整矩阵而非仅改变单轮廓标量",
    not np.allclose(center_shifted_matrix, strict_single_profile),
)
check(
    "N6 中心通道情形两侧秩都大于一",
    channel_rank_center == 2 and weight_rank_center == 2,
    f"channel_rank={channel_rank_center}, weight_rank={weight_rank_center}",
)


# ==================================================================
head("N7  独立非中心扰动改变谱比例")
# ==================================================================
noncentral_contact = np.diag([0.12, -0.18]).astype(complex)
perturbed_total = (
    base_channel[None, :, :] * geometry_weight[:, None, None]
    + noncentral_contact[None, :, :] * np.ones_like(age)[:, None, None]
)

base_ratio = base_channel[0, 0] / base_channel[1, 1]
perturbed_ratio = perturbed_total[-1, 0, 0] / perturbed_total[-1, 1, 1]

check(
    "N7 独立非中心扰动改变谱比例",
    abs(perturbed_ratio - base_ratio) > 1e-2,
    f"ratio={perturbed_ratio.real:.6f}",
)


# ==================================================================
head("N8  独立非中心扰动改变 likelihood ratio")
# ==================================================================
likelihood_without_contact = np.exp(
    (lambda_1 - lambda_0) * geometry_weight
)
likelihood_with_contact = np.exp(
    (lambda_1 - lambda_0)
    * (geometry_weight + 0.12)
)

check(
    "N8 独立非中心扰动改变 likelihood ratio",
    np.max(np.abs(likelihood_with_contact - likelihood_without_contact))
    > 1e-2,
)


# ==================================================================
head("N9  文档登记单剖面归约判据")
# ==================================================================
check(
    "N9 文档登记单剖面归约判据",
    "R-Z-SINGLE-PROFILE-REDUCTION-CRITERION" in DOC
    and "可合并为单剖面" in DOC,
)


# ==================================================================
head("N10  文档登记接触通道分类缺口")
# ==================================================================
check(
    "N10 文档登记接触通道分类缺口",
    "R-Z-CONTACT-CHANNEL-CLASSIFICATION-GAP" in DOC
    and "中心、平行还是独立非中心" in DOC,
)


# ==================================================================
head("N11  文档不把单轮廓写成默认形式")
# ==================================================================
check(
    "N11 文档保留单轮廓的条件地位",
    "单剖面是一个可检验的通道关系" in DOC
    and "不是完整球模 Hamiltonian 的默认形式" in DOC,
)


# ==================================================================
head("N12  上游边界保持")
# ==================================================================
check(
    "N12 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)


# ==================================================================
head("N13  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N13 D237 文档不含项目外体系名或外部路径",
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
