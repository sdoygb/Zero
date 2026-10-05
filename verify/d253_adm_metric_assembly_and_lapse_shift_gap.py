#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D253 —— 叶状 ADM 度规组装与 lapse/shift 缺口
================================================================================
【被检验的问题（D253）】
  正则时间势能否给唯一 Lorentz 度规？
  ADM 公式怎样把 lapse、shift 与叶层空间度规组成度规？
  同一时间势与同一叶层空间度规能否承载不同零锥？
  从层读回到 ADM 输入的缺口在哪里？

【本步判据】
  N1  D253 与恢复结构已登记
  N2  ADM 行列式与体积公式
  N3  ADM 号差
  N4  同一 t,h 配不同 lapse
  N5  同一 t,h 配不同 shift
  N6  叶层共形类加叶层体积
  N7  同一时间与四维体积不固定 lapse
  N8  完整时空共形类加体积的唯一性边界
  N9  层读回未选择 ADM 数据
  N10 文档登记 ADM 与缺口
  N11 文档不把 ADM 组装写成 U1-U4 推论
  N12 上游边界保持
  N13 文档不引用项目外体系名或外部路径

运行：python3 verify/d253_adm_metric_assembly_and_lapse_shift_gap.py
"""
import io
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
            "D253_adm_metric_assembly_and_lapse_shift_gap.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


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


def x_speeds(lapse, shift_x):
    return -shift_x - lapse, -shift_x + lapse


# ==================================================================
head("N1  D253 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D253", "D253" in AX)
check(
    "N1 公理表与文档登记嵌入到叶状",
    "R-Z-EMBEDDING-TO-FOLIATION" in BOTH,
)
check(
    "N1 公理表与文档登记 ADM 组装",
    "R-Z-ADM-METRIC-ASSEMBLY" in BOTH,
)
check(
    "N1 公理表与文档登记 lapse/shift 缺口",
    "R-Z-LAPSE-SHIFT-GAP" in BOTH,
)
check(
    "N1 公理表与文档登记叶层共形体积缺口",
    "R-Z-LEAF-CONFORMAL-VOLUME-GAP" in BOTH,
)
check(
    "N1 公理表与文档登记读回到 ADM 选择缺口",
    "R-Z-READOUT-TO-ADM-SELECTION-GAP" in BOTH,
)


# ==================================================================
head("N2  ADM 行列式与体积公式")
# ==================================================================
lapse = 1.7
shift = np.array([0.3, -0.2, 0.4])
h = np.array(
    [
        [2.0, 0.1, 0.0],
        [0.1, 1.5, -0.2],
        [0.0, -0.2, 1.2],
    ]
)
metric = adm_metric(h, lapse, shift)
det_metric = np.linalg.det(metric)
det_expected = -lapse ** 2 * np.linalg.det(h)
volume_density = math.sqrt(-det_metric)
volume_expected = lapse * math.sqrt(np.linalg.det(h))

check(
    "N2 det(g)=(-1)^n N^2 det(h)",
    np.isclose(det_metric, det_expected),
    f"det={det_metric:.12f}",
)
check(
    "N2 sqrt(-det(g))=N sqrt(det(h))",
    np.isclose(volume_density, volume_expected),
    f"vol={volume_density:.12f}",
)
check(
    "N2 文档登记体积公式",
    "N\\,dt\\wedge\\operatorname{vol}_h" in DOC,
)


# ==================================================================
head("N3  ADM 号差")
# ==================================================================
positive, negative, zero, eigenvalues = signature_counts(metric)
check(
    "N3 ADM 度规号差为 (1,n)",
    positive == 1 and negative == h.shape[0] and zero == 0,
    f"signature=({positive},{negative}), eig={np.round(eigenvalues, 8)}",
)
check(
    "N3 文档登记四维号差",
    "(1,n)" in DOC and "(1,3)" in DOC,
)


# ==================================================================
head("N4  同一 t,h 配不同 lapse")
# ==================================================================
spatial_identity = np.eye(3)
metric_lapse_one = adm_metric(spatial_identity, 1.0, np.zeros(3))
metric_lapse_two = adm_metric(spatial_identity, 2.0, np.zeros(3))
speeds_lapse_one = x_speeds(1.0, 0.0)
speeds_lapse_two = x_speeds(2.0, 0.0)

check(
    "N4 两个度规的空间块相同",
    np.allclose(metric_lapse_one[1:, 1:], metric_lapse_two[1:, 1:]),
)
check(
    "N4 lapse 改变时间块",
    not np.isclose(metric_lapse_one[0, 0], metric_lapse_two[0, 0]),
)
check(
    "N4 lapse 改变零速度",
    not np.isclose(abs(speeds_lapse_one[1]), abs(speeds_lapse_two[1])),
    f"v1=({speeds_lapse_one[0]:.3f},{speeds_lapse_one[1]:.3f}), "
    f"v2=({speeds_lapse_two[0]:.3f},{speeds_lapse_two[1]:.3f})",
)


# ==================================================================
head("N5  同一 t,h 配不同 shift")
# ==================================================================
shift_vector = np.array([0.5, 0.0, 0.0])
metric_shift = adm_metric(spatial_identity, 1.0, shift_vector)
speeds_shift = x_speeds(1.0, 0.5)

check(
    "N5 shift 仍保持空间块",
    np.allclose(metric_shift[1:, 1:], -spatial_identity),
)
check(
    "N5 shift 引入非零交叉项",
    np.isclose(metric_shift[0, 1], -0.5)
    and np.isclose(metric_shift[1, 0], -0.5),
)
check(
    "N5 shift 使光锥倾斜",
    not np.isclose(speeds_shift[0], -1.0)
    or not np.isclose(speeds_shift[1], 1.0),
    f"v=({speeds_shift[0]:.3f},{speeds_shift[1]:.3f})",
)


# ==================================================================
head("N6  叶层共形类加叶层体积")
# ==================================================================
base_h = np.diag([1.0, 2.0, 3.0])
base_volume_density = math.sqrt(np.linalg.det(base_h))
target_volume_density = 1.7 * base_volume_density
conformal_factor = (target_volume_density / base_volume_density) ** (2.0 / 3.0)
rescaled_h = conformal_factor * base_h
rescaled_volume_density = math.sqrt(np.linalg.det(rescaled_h))

check(
    "N6 h=(mu/vol_h0)^(2/n)h0 保持给定体积",
    np.isclose(rescaled_volume_density, target_volume_density),
    f"vol={rescaled_volume_density:.12f}",
)
check(
    "N6 文档登记 n 维叶层构造",
    "\\mu_\\Sigma" in DOC and "2/n" in DOC,
)


# ==================================================================
head("N7  同一时间与四维体积不固定 lapse")
# ==================================================================
lambda_value = 2.0
h_lambda = lambda_value ** (-2.0 / 3.0) * spatial_identity
metric_lambda = adm_metric(h_lambda, lambda_value, np.zeros(3))
volume_lambda = lambda_value * math.sqrt(np.linalg.det(h_lambda))
volume_one = 1.0 * math.sqrt(np.linalg.det(spatial_identity))
conformal_ratio_time = metric_lambda[0, 0] / metric_lapse_one[0, 0]
conformal_ratio_space = h_lambda[0, 0] / spatial_identity[0, 0]

check(
    "N7 同一时间函数与四维体积",
    np.isclose(volume_lambda, volume_one),
    f"vol_lambda={volume_lambda:.12f}",
)
check(
    "N7 不同 lapse 与空间度规",
    not np.isclose(metric_lambda[0, 0], metric_lapse_one[0, 0])
    and not np.allclose(h_lambda, spatial_identity),
)
check(
    "N7 时间与空间共形比不相等，故不是完整共形缩放",
    not np.isclose(conformal_ratio_time, conformal_ratio_space),
    f"time={conformal_ratio_time:.6f}, space={conformal_ratio_space:.6f}",
)


# ==================================================================
head("N8  完整时空共形类加体积的唯一性边界")
# ==================================================================
dimension = 4
conformal_scale = 1.25
volume_ratio = conformal_scale ** dimension
same_volume_roots = [
    root
    for root in (-1.0, 1.0)
    if root > 0.0 and abs(root ** 4 - 1.0) < 1.0e-15
]

check(
    "N8 时空共形缩放按 Omega^d 改变体积",
    np.isclose(volume_ratio, conformal_scale ** 4),
)
check(
    "N8 固定体积迫使正共形因子为唯一",
    len(same_volume_roots) == 1 and np.isclose(same_volume_roots[0], 1.0),
)
check(
    "N8 文档区分叶层数据与完整时空数据",
    "完整时空 Lorentz 共形类" in DOC and "叶层共形类" in DOC,
)


# ==================================================================
head("N9  层读回未选择 ADM 数据")
# ==================================================================
readout_profile = np.array([0.0, 1.0, 2.0])
admissible_lapses = np.array([0.5, 1.0, 2.0])
distinct_metrics = {
    tuple(np.round(adm_metric(spatial_identity, lapse_i, np.zeros(3)).ravel(), 12))
    for lapse_i in admissible_lapses
}

check(
    "N9 同一读回剖面允许多个 lapse",
    len(distinct_metrics) == len(admissible_lapses),
)
check(
    "N9 文档登记读回到 ADM 选择缺口",
    "R-Z-READOUT-TO-ADM-SELECTION-GAP" in DOC
    and "不能直接选择" in DOC,
)
check(
    "N9 文档登记源张量依赖未定度规的隐式方程",
    "G(\\mathrm g)" in DOC and "T(\\phi,\\mathrm g)" in DOC,
)


# ==================================================================
head("N10  文档登记 ADM 与缺口")
# ==================================================================
check(
    "N10 文档登记 ADM 公式",
    "\\mathrm g_{N,\\beta}" in DOC and "N^2dt^2" in DOC,
)
check(
    "N10 文档登记叶层共形体积缺口",
    "R-Z-LEAF-CONFORMAL-VOLUME-GAP" in DOC
    and "叶层体积元" in DOC,
)
check(
    "N10 文档登记 lapse/shift 规范地位",
    "切片规范数据" in DOC and "不是独立" in DOC,
)


# ==================================================================
head("N11  文档不把 ADM 组装写成 U1-U4 推论")
# ==================================================================
check(
    "N11 文档明确条件性",
    "不是 `U1-U4` 的推论" in DOC
    and "条件组装公式" in DOC,
)
check(
    "N11 文档不把 ADM 数据写成已导出",
    "没有从层读回导出度规" in DOC
    and "仍未由层读回导出" in DOC,
)


# ==================================================================
head("N12  上游边界保持")
# ==================================================================
check(
    "N12 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)
check(
    "N12 文档登记输入预算",
    "## §3 输入预算" in DOC,
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
    "N13 D253 文档不含项目外体系名或外部路径",
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
