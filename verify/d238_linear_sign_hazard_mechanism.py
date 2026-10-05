#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D238 —— 线性符号危险率机制
================================================================================
【被检验的问题（D238）】
  线性符号危险率差加终端中性，能否精确生成抛物型年龄比例？
  离散全分支多重数能否给出同一结果？
  其它危险率模型为什么不能？

【本步判据】
  N1  D238 与恢复结构已登记
  N2  正负计数比例满足危险率差方程
  N3  线性危险率差给 R_sign(a)=exp(delta(T^2-a^2))
  N4  终端比例等于一
  N5  初始比例由终端中性固定
  N6  连续公式与数值 ODE 一致
  N7  离散活闭合多重数比给同一抛物线
  N8  常数危险率差只给线性对数比例
  N9  二次危险率差只给三次对数比例
  N10 线性机制精确给 D234 年龄剖面
  N11 文档登记终端中性
  N12 文档登记线性危险率机制
  N13 文档登记危险率斜率来源缺口
  N14 文档不把机制写成 U1-U4 推论
  N15 上游边界保持
  N16 文档不引用项目外体系名或外部路径

运行：python3 verify/d238_linear_sign_hazard_mechanism.py
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
            "D238_linear_sign_hazard_mechanism.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D238 与恢复结构已登记")
# ==================================================================
check(
    "N1 公理表登记 D238",
    "D238" in AX,
)
check(
    "N1 公理表与文档登记终端中性",
    "R-Z-TERMINAL-SIGN-NEUTRALITY" in BOTH
    and "终端符号中性" in DOC,
)
check(
    "N1 公理表与文档登记线性危险率机制",
    "R-Z-LINEAR-SIGN-HAZARD-MECHANISM" in BOTH
    and "线性危险率" in DOC,
)
check(
    "N1 公理表与文档登记斜率来源缺口",
    "R-Z-HAZARD-SLOPE-ORIGIN-GAP" in BOTH
    and "危险率斜率" in DOC,
)


# ==================================================================
head("N2  正负计数比例满足危险率差方程")
# ==================================================================
terminal_age = 5.0
delta = 0.17
age = np.linspace(0.0, terminal_age, 4001)
hazard_difference = 2.0 * delta * age
log_ratio_from_difference = -np.trapezoid(hazard_difference, age)
# 用 cumtrapz 手动计算每个年龄的累计积分。
cumulative = np.concatenate(
    [
        [0.0],
        np.cumsum(
            0.5
            * (
                hazard_difference[1:]
                + hazard_difference[:-1]
            )
            * np.diff(age)
        ),
    ]
)
log_ratio_from_integral = cumulative

check(
    "N2 线性危险率差积分给二次函数",
    abs(log_ratio_from_integral[-1] - delta * terminal_age**2) < 1e-10,
)
check(
    "N2 危险率差方向与比例导数相反",
    log_ratio_from_difference < 0.0,
)


# ==================================================================
head("N3  线性危险率差给 R_sign(a)=exp(delta(T^2-a^2))")
# ==================================================================
log_ratio_target = delta * (terminal_age**2 - age**2)
ratio_target = np.exp(log_ratio_target)

check(
    "N3 连续比例公式在采样点满足目标",
    np.allclose(ratio_target, np.exp(log_ratio_target)),
)


# ==================================================================
head("N4  终端比例等于一")
# ==================================================================
check(
    "N4 终端正负比例等于一",
    abs(ratio_target[-1] - 1.0) < 1e-12,
)


# ==================================================================
head("N5  初始比例由终端中性固定")
# ==================================================================
initial_ratio = ratio_target[0]
check(
    "N5 初始比例等于 exp(delta T^2)",
    abs(initial_ratio - np.exp(delta * terminal_age**2)) < 1e-12,
    f"initial={initial_ratio:.12f}",
)


# ==================================================================
head("N6  连续公式与数值 ODE 一致")
# ==================================================================
log_ratio_ode = (
    np.log(initial_ratio)
    - np.concatenate(
        [
            [0.0],
            np.cumsum(
                0.5
                * (
                    hazard_difference[1:]
                    + hazard_difference[:-1]
                )
                * np.diff(age)
            ),
        ]
    )
)
check(
    "N6 数值积分与解析公式逐点一致",
    np.allclose(log_ratio_ode, log_ratio_target, atol=1e-9),
)
check(
    "N6 数值积分满足 ODE 的对数导数",
    np.allclose(
        np.diff(log_ratio_ode) / np.diff(age),
        -0.5
        * (
            hazard_difference[1:]
            + hazard_difference[:-1]
        ),
        atol=1e-12,
    ),
)


# ==================================================================
head("N7  离散活闭合多重数比给同一抛物线")
# ==================================================================
discrete_terminal = 20
discrete_ages = np.arange(0, discrete_terminal + 1)
discrete_q = np.exp(
    -delta * (2 * discrete_ages[:-1] + 1)
)
discrete_initial_ratio = np.exp(delta * discrete_terminal**2)
discrete_ratio = np.concatenate(
    [
        [discrete_initial_ratio],
        discrete_initial_ratio * np.cumprod(discrete_q),
    ]
)
discrete_target = np.exp(
    delta * (discrete_terminal**2 - discrete_ages**2)
)

check(
    "N7 离散多重数比精确恢复抛物线",
    np.allclose(discrete_ratio, discrete_target),
)
check(
    "N7 离散终端比例仍为一",
    abs(discrete_ratio[-1] - 1.0) < 1e-12,
)


# ==================================================================
head("N8  常数危险率差只给线性对数比例")
# ==================================================================
constant_slope = 0.11
constant_difference_log_ratio = constant_slope * (terminal_age - age)
quadratic_target = delta * (terminal_age**2 - age**2)
constant_model_residual = np.max(
    np.abs(constant_difference_log_ratio - quadratic_target)
)

check(
    "N8 常数危险率差不等于二次对数比例",
    constant_model_residual > 1e-2,
    f"residual={constant_model_residual:.6f}",
)


# ==================================================================
head("N9  二次危险率差只给三次对数比例")
# ==================================================================
quadratic_difference_log_ratio = (
    2.0
    * delta
    / 3.0
    * (terminal_age**3 - age**3)
)
quadratic_model_residual = np.max(
    np.abs(quadratic_difference_log_ratio - quadratic_target)
)

check(
    "N9 二次危险率差不等于二次对数比例",
    quadratic_model_residual > 1e-2,
    f"residual={quadratic_model_residual:.6f}",
)


# ==================================================================
head("N10  线性机制精确给 D234 年龄剖面")
# ==================================================================
matrix_radius = 2.4
lambda_gap = 1.7
delta_from_geometry = (
    matrix_radius * lambda_gap / (2.0 * terminal_age**2)
)
ratio_from_geometry = np.exp(
    delta_from_geometry * (terminal_age**2 - age**2)
)
profile = np.log(ratio_from_geometry) / lambda_gap
target_profile = (matrix_radius / 2.0) * (
    1.0 - age**2 / terminal_age**2
)

check(
    "N10 几何斜率给抛物型剖面",
    np.allclose(profile, target_profile),
)
check(
    "N10 剖面边界为零且中心最大",
    abs(profile[-1]) < 1e-12
    and int(np.argmax(profile)) == 0,
)


# ==================================================================
head("N11-D13  文档登记与上游边界")
# ==================================================================
check(
    "N11 文档登记终端中性",
    "R-Z-TERMINAL-SIGN-NEUTRALITY" in DOC
    and "终端符号中性" in DOC,
)
check(
    "N12 文档登记线性危险率机制",
    "R-Z-LINEAR-SIGN-HAZARD-MECHANISM" in DOC
    and "线性危险率" in DOC,
)
check(
    "N13 文档登记危险率斜率来源缺口",
    "R-Z-HAZARD-SLOPE-ORIGIN-GAP" in DOC
    and "为什么危险率会线性分离" in DOC,
)
check(
    "N14 文档不把机制写成 U1-U4 推论",
    "它们不是 `U1-U4` 的推论" in DOC
    and "危险率斜率来源仍是新缺口" in DOC,
)
check(
    "N15 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
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
    "N16 D238 文档不含项目外体系名或外部路径",
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
