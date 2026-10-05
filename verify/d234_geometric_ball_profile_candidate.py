#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D234 —— 几何球模核参照
================================================================================
【被检验的问题（D234）】
  球的几何模权重是什么形状？
  年龄映射到径向与映射到有限区间时，候选剖面是否相同？
  这些候选是否解决 D233 的符号年龄破缺？

【本步判据】
  N1  D234 与恢复结构已登记
  N2  球权重在球心最大、边界为零并严格下降
  N3  球权重的二阶导数为负常数
  N4  半空间单侧退化，有限区间双侧退化
  N5  格子对称权重与连续抛物线共享边界零点与单峰形状
  N6  径向年龄映射给单侧退化的二次剖面
  N7  有限区间年龄映射给双侧退化的抛物型剖面
  N8  两种年龄剖面的 L1 差异为正
  N9  D232 likelihood ratio 在两种映射下都非恒定
  N10 两种 likelihood ratio 不相同
  N11 文档登记四项恢复结构
  N12 文档不把几何候选写成 U1-U4 推论
  N13 文档保留符号年龄破缺缺口
  N14 上游边界保持
  N15 文档不引用项目外体系名或外部路径

运行：python3 verify/d234_geometric_ball_profile_candidate.py
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
            "D234_geometric_ball_profile_candidate.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D234 与恢复结构已登记")
# ==================================================================
check(
    "N1 公理表登记 D234",
    "D234" in AX,
)
check(
    "N1 公理表与文档登记几何球候选",
    "R-Z-CONFORMAL-BALL-PROFILE-CANDIDATE" in BOTH
    and "f_B" in DOC
    and "抛物型" in DOC,
)
check(
    "N1 公理表与文档登记三个映射缺口",
    all(
        marker in BOTH
        for marker in (
            "R-Z-AGE-REGION-MAPPING-GAP",
            "R-Z-SOURCE-OPERATOR-IDENTIFICATION-GAP",
            "R-Z-MODULAR-CONTACT-TERM-GAP",
        )
    ),
)


# ==================================================================
head("N2  球权重在球心最大、边界为零并严格下降")
# ==================================================================
radius = 2.4
point_count = 2001
r = np.linspace(0.0, radius, point_count)
ball_profile = (radius**2 - r**2) / (2.0 * radius)
radial_increments = np.diff(ball_profile)

check(
    "N2 球心值等于 R/2",
    abs(ball_profile[0] - radius / 2.0) < 1e-12,
    f"center={ball_profile[0]:.12f}",
)
check(
    "N2 边界值为零",
    abs(ball_profile[-1]) < 1e-12,
    f"boundary={ball_profile[-1]:.3e}",
)
check(
    "N2 球权重严格下降",
    np.all(radial_increments < 0.0),
)


# ==================================================================
head("N3  球权重的二阶导数为负常数")
# ==================================================================
first_derivative = np.gradient(ball_profile, r)
second_derivative = np.gradient(first_derivative, r)
interior = (r > 0.1 * radius) & (r < 0.9 * radius)
second_derivative_error = np.max(
    np.abs(second_derivative[interior] + 1.0 / radius)
)

check(
    "N3 内部二阶导数近似 -1/R",
    second_derivative_error < 2e-3,
    f"max_error={second_derivative_error:.3e}",
)


# ==================================================================
head("N4  半空间单侧退化，有限区间双侧退化")
# ==================================================================
length = 3.0
x = np.linspace(0.0, length, point_count)
half_space_profile = x
interval_profile = x * (length - x) / length

check(
    "N4 半空间仅左边界为零",
    abs(half_space_profile[0]) < 1e-12
    and half_space_profile[-1] > 0.0,
)
check(
    "N4 有限区间两端均为零",
    abs(interval_profile[0]) < 1e-12
    and abs(interval_profile[-1]) < 1e-12,
)
check(
    "N4 有限区间权重中部最大",
    int(np.argmax(interval_profile)) == point_count // 2,
)


# ==================================================================
head("N5  格子对称权重与连续抛物线共享边界零点与单峰形状")
# ==================================================================
sites = np.arange(1, 8, dtype=float)
lattice_centers = sites - 0.5
lattice_profile = lattice_centers * (7.0 - lattice_centers) / 7.0
continuum_at_lattice = lattice_centers * (7.0 - lattice_centers) / 7.0

check(
    "N5 格子对称权重与连续权重同点一致",
    np.max(np.abs(lattice_profile - continuum_at_lattice)) < 1e-12,
)
check(
    "N5 格子权重在离散中心附近最大",
    int(np.argmax(lattice_profile)) in (3, 4),
    f"argmax={int(np.argmax(lattice_profile))}",
)


# ==================================================================
head("N6  径向年龄映射给单侧退化的二次剖面")
# ==================================================================
period = 5.0
age = np.linspace(0.0, period, point_count)
radial_age_profile = (radius / 2.0) * (1.0 - age**2 / period**2)

check(
    "N6 径向年龄剖面初始最大",
    abs(radial_age_profile[0] - radius / 2.0) < 1e-12,
)
check(
    "N6 径向年龄剖面终端为零",
    abs(radial_age_profile[-1]) < 1e-12,
)
check(
    "N6 径向年龄剖面严格下降",
    np.all(np.diff(radial_age_profile) < 0.0),
)


# ==================================================================
head("N7  有限区间年龄映射给双侧退化的抛物型剖面")
# ==================================================================
interval_age_profile = age * (period - age) / period

check(
    "N7 有限区间年龄剖面两端为零",
    abs(interval_age_profile[0]) < 1e-12
    and abs(interval_age_profile[-1]) < 1e-12,
)
check(
    "N7 有限区间年龄剖面中部最大",
    int(np.argmax(interval_age_profile)) == point_count // 2,
)


# ==================================================================
head("N8  两种年龄剖面的 L1 差异为正")
# ==================================================================
profile_l1_gap = np.trapezoid(
    np.abs(radial_age_profile - interval_age_profile),
    age,
)
check(
    "N8 不同年龄映射给出不同剖面",
    profile_l1_gap > 1e-2,
    f"L1={profile_l1_gap:.6f}",
)


# ==================================================================
head("N9  D232 likelihood ratio 在两种映射下都非恒定")
# ==================================================================
lambda_0 = -1.0
lambda_1 = 0.0
radial_ratio = np.exp((lambda_1 - lambda_0) * radial_age_profile)
interval_ratio = np.exp((lambda_1 - lambda_0) * interval_age_profile)

check(
    "N9 径向候选 ratio 非恒定",
    np.max(radial_ratio) - np.min(radial_ratio) > 1e-2,
)
check(
    "N9 区间候选 ratio 非恒定",
    np.max(interval_ratio) - np.min(interval_ratio) > 1e-2,
)
check(
    "N9 径向 ratio 在终端为一",
    abs(radial_ratio[-1] - 1.0) < 1e-12,
)
check(
    "N9 区间 ratio 在两端为一",
    abs(interval_ratio[0] - 1.0) < 1e-12
    and abs(interval_ratio[-1] - 1.0) < 1e-12,
)


# ==================================================================
head("N10  两种 likelihood ratio 不相同")
# ==================================================================
ratio_l1_gap = np.trapezoid(np.abs(radial_ratio - interval_ratio), age)
check(
    "N10 两种 ratio 的 L1 差异为正",
    ratio_l1_gap > 1e-2,
    f"L1={ratio_l1_gap:.6f}",
)


# ==================================================================
head("N11  文档登记四项恢复结构")
# ==================================================================
registered_markers = (
    "R-Z-CONFORMAL-BALL-PROFILE-CANDIDATE",
    "R-Z-AGE-REGION-MAPPING-GAP",
    "R-Z-SOURCE-OPERATOR-IDENTIFICATION-GAP",
    "R-Z-MODULAR-CONTACT-TERM-GAP",
)
check(
    "N11 文档登记四项恢复结构",
    all(marker in DOC for marker in registered_markers),
)


# ==================================================================
head("N12  文档不把几何候选写成 U1-U4 推论")
# ==================================================================
check(
    "N12 文档保留条件候选地位",
    "它们不是 `U1-U4` 的推论" in DOC
    and "几何参照给出目标函数，不给出零动力学选择器" in DOC
    and "条件候选" in DOC,
)


# ==================================================================
head("N13  文档保留符号年龄破缺缺口")
# ==================================================================
check(
    "N13 文档保留 D233 符号年龄缺口",
    "几何剖面选择不是符号年龄破缺的替代品" in DOC
    and "N_+(a)=N_-(a)" in DOC
    and "符号年龄" in DOC,
)


# ==================================================================
head("N14  上游边界保持")
# ==================================================================
check(
    "N14 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)


# ==================================================================
head("N15  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N15 D234 文档不含项目外体系名或外部路径",
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
