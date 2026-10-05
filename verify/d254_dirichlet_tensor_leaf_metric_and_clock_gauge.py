#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D254 —— Dirichlet 张量到叶层度规与读回时钟规范
================================================================================
【被检验的问题（D254）】
  标量读回与单位 lapse、零 shift 规范能否写成 g=dt^2-h？
  连续 Dirichlet 张量何时唯一给叶层空间度规？
  二维与三维为什么不同？
  离散交换权重能否直接恢复连续张量？

【本步判据】
  N1  D254 与恢复结构已登记
  N2  时钟规范给出 g=dt^2-h
  N3  同一 t,h 的其他 lapse/shift 不被自动排除
  N4  三维 Dirichlet 张量唯一反解 h
  N5  二维 Dirichlet 张量只给共形类
  N6  测度与度规体积一致性
  N7  独立测度的共形重标度歧义
  N8  离散交换权重不能直接给连续张量
  N9  文档登记 ADM 缺口更新
  N10 文档不把时钟规范写成上游推论
  N11 上游边界保持
  N12 文档不引用项目外体系名或外部路径

运行：python3 verify/d254_dirichlet_tensor_leaf_metric_and_clock_gauge.py
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
            "D254_dirichlet_tensor_leaf_metric_and_clock_gauge.md",
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


# ==================================================================
head("N1  D254 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D254", "D254" in AX)
check(
    "N1 公理表与文档登记读回时钟规范",
    "R-Z-READOUT-CLOCK-GAUGE" in BOTH,
)
check(
    "N1 公理表与文档登记 Dirichlet 张量到度规",
    "R-Z-DIRICHLET-TENSOR-METRIC" in BOTH,
)
check(
    "N1 公理表与文档登记度规体积一致性",
    "R-Z-METRIC-VOLUME-CONSISTENCY" in BOTH,
)
check(
    "N1 公理表与文档登记时钟规范缺口",
    "R-Z-CLOCK-GAUGE-SELECTION-GAP" in BOTH,
)
check(
    "N1 公理表与文档登记离散张量缺口",
    "R-Z-DISCRETE-DIRICHLET-TENSOR-GAP" in BOTH,
)


# ==================================================================
head("N2  时钟规范给出 g=dt^2-h")
# ==================================================================
h_sample = np.array(
    [
        [2.0, 0.1, 0.0],
        [0.1, 1.5, -0.2],
        [0.0, -0.2, 1.2],
    ]
)
metric_clock = adm_metric(h_sample, 1.0, np.zeros(3))
positive, negative, zero, eigenvalues = signature_counts(metric_clock)

check(
    "N2 N=1, beta=0 给时间块为一",
    np.isclose(metric_clock[0, 0], 1.0),
)
check(
    "N2 N=1, beta=0 给空间块为 -h",
    np.allclose(metric_clock[1:, 1:], -h_sample),
)
check(
    "N2 时钟规范度规号差为 (1,n)",
    positive == 1 and negative == 3 and zero == 0,
    f"eig={np.round(eigenvalues, 8)}",
)
check(
    "N2 文档登记 g=dt^2-h",
    "dt^2-h_{ij}dx^i dx^j" in DOC
    or "\\mathrm g=dt^2-h" in DOC,
)


# ==================================================================
head("N3  同一 t,h 的其他 lapse/shift 不被自动排除")
# ==================================================================
metric_lapse = adm_metric(h_sample, 2.0, np.zeros(3))
metric_shift = adm_metric(h_sample, 1.0, np.array([0.5, 0.0, 0.0]))

check(
    "N3 lapse 变体与时钟规范不同",
    not np.allclose(metric_lapse, metric_clock),
)
check(
    "N3 shift 变体与时钟规范不同",
    not np.allclose(metric_shift, metric_clock),
)
check(
    "N3 shift 变体仍保持空间块",
    np.allclose(metric_shift[1:, 1:], -h_sample),
)


# ==================================================================
head("N4  三维 Dirichlet 张量唯一反解 h")
# ==================================================================
det_h = float(np.linalg.det(h_sample))
Q_three = math.sqrt(det_h) * np.linalg.inv(h_sample)
h_recovered = float(np.linalg.det(Q_three)) * np.linalg.inv(Q_three)
det_relation_three = math.sqrt(det_h)

check(
    "N4 det(Q)=sqrt(det h) 在 n=3 成立",
    np.isclose(np.linalg.det(Q_three), det_relation_three),
    f"detQ={np.linalg.det(Q_three):.12f}",
)
check(
    "N4 h=(det Q)Q^{-1} 精确恢复",
    np.allclose(h_recovered, h_sample),
    f"max_error={np.max(np.abs(h_recovered - h_sample)):.3e}",
)
check(
    "N4 文档登记三维唯一反解",
    "(\\det Q)\\,Q^{-1}" in DOC,
)


# ==================================================================
head("N5  二维 Dirichlet 张量只给共形类")
# ==================================================================
h_two = np.array([[2.0, 0.3], [0.3, 1.4]])
det_h_two = float(np.linalg.det(h_two))
Q_two = math.sqrt(det_h_two) * np.linalg.inv(h_two)
scaled_h_two = 1.7 * h_two
Q_two_scaled = math.sqrt(float(np.linalg.det(scaled_h_two))) * np.linalg.inv(scaled_h_two)

check(
    "N5 二维 det(Q)=1",
    np.isclose(np.linalg.det(Q_two), 1.0),
    f"detQ={np.linalg.det(Q_two):.12f}",
)
check(
    "N5 二维共形缩放保持 Q",
    np.allclose(Q_two_scaled, Q_two),
)
check(
    "N5 二维共形缩放改变 h",
    not np.allclose(scaled_h_two, h_two),
)


# ==================================================================
head("N6  测度与度规体积一致性")
# ==================================================================
omega_values = np.array([0.5, 0.75, 1.0, 1.25, 1.5, 2.0])
consistent_roots_three = [
    value
    for value in omega_values
    if np.isclose(value ** 3, value ** 2)
]

check(
    "N6 三维一致性条件 Omega^(n-2)=1 只有唯一正解",
    len(consistent_roots_three) == 1
    and np.isclose(consistent_roots_three[0], 1.0),
)
check(
    "N6 二维一致性条件对全部正 Omega 成立",
    all(np.isclose(value ** 2, value ** 2) for value in omega_values),
)
check(
    "N6 文档登记度规体积一致性",
    "\\mu=\\operatorname{vol}_h" in DOC
    or "\\mu=\\operatorname{vol}_h" in DOC,
)


# ==================================================================
head("N7  独立测度的共形重标度歧义")
# ==================================================================
omega = 1.6
h_scaled = omega ** 2 * h_sample
mu_base = math.sqrt(det_h)
mu_scaled = omega ** 2 * mu_base
Q_base = mu_base * np.linalg.inv(h_sample)
Q_scaled = mu_scaled * np.linalg.inv(h_scaled)

check(
    "N7 独立测度与 h 同步缩放保持 Q",
    np.allclose(Q_base, Q_scaled),
)
check(
    "N7 若测度不绑定度规体积则存在共形歧义",
    not np.isclose(omega, 1.0)
    and np.allclose(Q_base, Q_scaled),
)


# ==================================================================
head("N8  离散交换权重不能直接给连续张量")
# ==================================================================
edge_weight = 2.0
spacing_a = 0.1
spacing_b = 0.2
coefficient_a = edge_weight * spacing_a
coefficient_b = edge_weight * spacing_b

check(
    "N8 同一图权重配不同边长给不同连续系数",
    not np.isclose(coefficient_a, coefficient_b),
    f"A_a={coefficient_a:.6f}, A_b={coefficient_b:.6f}",
)
check(
    "N8 文档登记离散到连续张量缺口",
    "R-Z-DISCRETE-DIRICHLET-TENSOR-GAP" in DOC
    and "边长" in DOC,
)


# ==================================================================
head("N9  文档登记 ADM 缺口更新")
# ==================================================================
check(
    "N9 文档登记条件度规组装",
    "\\mathrm g=dt^2-h" in DOC,
)
check(
    "N9 文档登记时钟规范选择缺口",
    "R-Z-CLOCK-GAUGE-SELECTION-GAP" in DOC
    and "不是读回标量的自动推论" in DOC,
)
check(
    "N9 文档登记三维空间侧反解",
    "R-Z-DIRICHLET-TENSOR-METRIC" in DOC
    and "唯一反解" in DOC,
)


# ==================================================================
head("N10  文档不把时钟规范写成上游推论")
# ==================================================================
check(
    "N10 文档明确条件性",
    "不是 `U1-U4` 的推论" in DOC
    and "条件结构" in DOC,
)
check(
    "N10 文档不把读回时钟写成已选择",
    "物理时间解释仍是输入" in DOC
    and "仍未导出" in DOC,
)


# ==================================================================
head("N11  上游边界保持")
# ==================================================================
check(
    "N11 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)
check(
    "N11 文档登记输入预算",
    "## §3 输入预算" in DOC,
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
    "N12 D254 文档不含项目外体系名或外部路径",
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
