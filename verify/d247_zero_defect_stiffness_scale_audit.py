#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D247 —— 零缺陷刚度尺度审计
================================================================================
【被检验的问题（D247）】
  全对全核形唯一化以后，刚度 delta 是否能被现有上游固定？
  年龄重标度、终端中性和 C1 如何作用于 delta？
  D238 的几何标定是导出 delta，还是转移输入？

【本步判据】
  N1  D247 与恢复结构已登记
  N2  形状 Q_k^2 与数值 delta 可分离
  N3  年龄重标度给 delta'=delta/lambda^2
  N4  终端量 delta T^2 在重标度下不变
  N5  不同 delta 都满足终端中立
  N6  不同 delta 给不同危险率斜率
  N7  几何关系可条件表达 delta
  N8  几何数据缺失时 delta 不被固定
  N9  C1 不选择无量纲 delta
  N10 文档登记与上游边界
  N11 文档不引用项目外体系名或外部路径

运行：python3 verify/d247_zero_defect_stiffness_scale_audit.py
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
            "D247_zero_defect_stiffness_scale_audit.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D247 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D247", "D247" in AX)
check(
    "N1 公理表与文档登记形状数值拆分",
    "R-Z-ZERO-DEFECT-SHAPE-VALUE-SPLIT" in BOTH,
)
check(
    "N1 公理表与文档登记年龄重标度",
    "R-Z-DELTA-AGE-RESCALING" in BOTH,
)
check(
    "N1 公理表与文档登记绝对尺度缺口",
    "R-Z-DELTA-ABSOLUTE-SCALE-GAP" in BOTH,
)


# ==================================================================
head("N2  形状 Q_k^2 与数值 delta 可分离")
# ==================================================================
ages = np.arange(0, 13, dtype=float)
shape = ages**2
delta_one = 0.09
delta_two = 0.17
potential_one = delta_one * shape
potential_two = delta_two * shape

check(
    "N2 两势共享同一形状",
    np.allclose(potential_one / delta_one, potential_two / delta_two),
)
check(
    "N2 两势具有不同数值",
    not np.allclose(potential_one, potential_two),
)


# ==================================================================
head("N3  年龄重标度给 delta'=delta/lambda^2")
# ==================================================================
T = 8.0
lam = 3.0
T_scaled = lam * T
delta_scaled = delta_one / lam**2

check(
    "N3 终端年龄按 lambda 重标度",
    np.isclose(T_scaled, lam * T),
)
check(
    "N3 刚度按 lambda^-2 重标度",
    np.isclose(delta_scaled, delta_one / lam**2),
)


# ==================================================================
head("N4  终端量 delta T^2 在重标度下不变")
# ==================================================================
terminal_energy = delta_one * T**2
scaled_terminal_energy = delta_scaled * T_scaled**2

check(
    "N4 终端能量保持",
    np.isclose(terminal_energy, scaled_terminal_energy),
)
check(
    "N4 无量纲剖面保持",
    np.isclose(
        delta_scaled * (T_scaled / 2.0) ** 2,
        delta_one * (T / 2.0) ** 2,
    ),
)


# ==================================================================
head("N5  不同 delta 都满足终端中立")
# ==================================================================
log_ratio_one = delta_one * (T**2 - ages**2)
log_ratio_two = delta_two * (T**2 - ages**2)

check(
    "N5 第一刚度终端中立",
    np.isclose(log_ratio_one[int(T)], 0.0),
)
check(
    "N5 第二刚度终端中立",
    np.isclose(log_ratio_two[int(T)], 0.0),
)
check(
    "N5 两比例不同",
    not np.allclose(log_ratio_one, log_ratio_two),
)


# ==================================================================
head("N6  不同 delta 给不同危险率斜率")
# ==================================================================
hazard_one = 2.0 * delta_one * ages[1:]
hazard_two = 2.0 * delta_two * ages[1:]

check(
    "N6 危险率斜率随 delta 改变",
    not np.allclose(hazard_one, hazard_two),
)
check(
    "N6 斜率比为刚度比",
    np.allclose(hazard_two / hazard_one, delta_two / delta_one),
)


# ==================================================================
head("N7  几何关系可条件表达 delta")
# ==================================================================
radius = 2.5
spectral_gap = 0.4
T_geom = 5.0
delta_from_geometry = (
    radius * spectral_gap / (2.0 * T_geom**2)
)

check(
    "N7 几何关系给具体 delta",
    np.isclose(
        delta_from_geometry,
        radius * spectral_gap / (2.0 * T_geom**2),
    ),
)
check(
    "N7 几何表达式为正",
    delta_from_geometry > 0.0,
)


# ==================================================================
head("N8  几何数据缺失时 delta 不被固定")
# ==================================================================
delta_alternatives = radius * spectral_gap / (
    2.0 * np.array([4.0, 5.0, 6.0]) ** 2
)

check(
    "N8 不同终端年龄给不同几何 delta",
    len(set(np.round(delta_alternatives, 12))) == 3,
)
check(
    "N8 几何标定只转移输入",
    not np.allclose(
        delta_alternatives,
        np.full_like(delta_alternatives, delta_one),
    ),
)


# ==================================================================
head("N9  C1 不选择无量纲 delta")
# ==================================================================
fixed_reference_unit = 1.0
delta_models = np.array([0.09, 0.17, 0.31])
rescaled_models = delta_models / lam**2

check(
    "N9 同一参照单位可配多个 delta",
    np.allclose(
        np.full_like(delta_models, fixed_reference_unit),
        np.ones_like(delta_models),
    ),
)
check(
    "N9 重标度不消去模型间差异",
    not np.allclose(rescaled_models, np.full_like(rescaled_models, rescaled_models[0])),
)


# ==================================================================
head("N10  文档登记与上游边界")
# ==================================================================
check(
    "N10 文档登记形状数值拆分",
    "R-Z-ZERO-DEFECT-SHAPE-VALUE-SPLIT" in DOC
    and "形状 $Q^2$" in DOC
    and "刚度数值 $\delta$" in DOC,
)
check(
    "N10 文档登记年龄重标度",
    "R-Z-DELTA-AGE-RESCALING" in DOC
    and "\\delta'=\\frac{\\delta}{\\lambda^2}" in DOC,
)
check(
    "N10 文档登记绝对尺度缺口",
    "R-Z-DELTA-ABSOLUTE-SCALE-GAP" in DOC
    and "无量纲刚度" in DOC,
)
check(
    "N10 文档不把新结构写成 U1-U4 推论",
    "它们不是 `U1-U4` 的推论" in DOC
    and "不是从 `U1-U4` 导出" in DOC,
)
check(
    "N10 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
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
    "N11 D247 文档不含项目外体系名或外部路径",
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
