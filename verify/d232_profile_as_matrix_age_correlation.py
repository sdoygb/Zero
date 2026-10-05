#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D232 —— 模密度剖面与矩阵年龄相关性
================================================================================
【被检验的问题（D232）】
  模密度剖面 f 是否等价于矩阵相位的年龄 likelihood ratio？
  乘积态是否只能给恒定剖面？
  非恒定剖面的相关性是否可测？

【本步判据】
  N1  D232 与两个恢复结构已登记
  N2  常量剖面给出年龄无关矩阵比例
  N3  线性剖面给出年龄相关矩阵比例
  N4  二次剖面给出年龄相关矩阵比例
  N5  likelihood ratio 可恢复剖面
  N6  常量剖面协方差为零
  N7  非恒定剖面协方差非零
  N8  D224 乘积态只给恒定剖面
  N9  文档登记剖面与相关性等价
  N10 文档登记矩阵年龄 likelihood 缺口
  N11 文档不把相关性来源写成已导出
  N12 上游边界保持
  N13 文档不引用项目外体系或外部路径

运行：python3 verify/d232_profile_as_matrix_age_correlation.py
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
            "D232_profile_as_matrix_age_correlation.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D232 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记剖面相关性改写",
    "D232" in AX
    and "R-Z-PROFILE-AS-MATRIX-AGE-CORRELATION" in BOTH
    and "矩阵相位" in DOC,
)
check(
    "N1 公理表与文档登记矩阵年龄 likelihood 缺口",
    "R-Z-MATRIX-AGE-LIKELIHOOD-GAP" in BOTH
    and "likelihood ratio" in DOC,
)


# ==================================================================
head("N2  常量剖面给出年龄无关矩阵比例")
# ==================================================================
age = np.linspace(0.0, 1.0, 301)
r_matrix = 0.75
lambda_0 = np.log(r_matrix)
lambda_1 = np.log(1.0 - r_matrix)

f_constant = np.ones_like(age)
f_linear = 1.0 - age
f_quadratic = 1.0 - age**2

profiles = {
    "constant": f_constant,
    "linear": f_linear,
    "quadratic": f_quadratic,
}


def conditional_probabilities(profile):
    weight_0 = np.exp(-lambda_0 * profile)
    weight_1 = np.exp(-lambda_1 * profile)
    normalization = weight_0 + weight_1
    return weight_0 / normalization, weight_1 / normalization


p0_constant, p1_constant = conditional_probabilities(f_constant)
check(
    "N2 常量剖面矩阵比例为常数",
    np.max(np.abs(p0_constant - p0_constant[0])) < 1e-12
    and np.max(np.abs(p1_constant - p1_constant[0])) < 1e-12,
)


# ==================================================================
head("N3  线性剖面给出年龄相关矩阵比例")
# ==================================================================
p0_linear, p1_linear = conditional_probabilities(f_linear)
check(
    "N3 线性剖面矩阵比例随年龄变化",
    np.linalg.norm(p0_linear - p0_linear[0]) > 1e-3,
)


# ==================================================================
head("N4  二次剖面给出年龄相关矩阵比例")
# ==================================================================
p0_quadratic, p1_quadratic = conditional_probabilities(f_quadratic)
check(
    "N4 二次剖面矩阵比例随年龄变化",
    np.linalg.norm(p0_quadratic - p0_quadratic[0]) > 1e-3,
)


# ==================================================================
head("N5  likelihood ratio 可恢复剖面")
# ==================================================================
profile_error = {}
for name, probabilities in {
    "constant": (p0_constant, p1_constant),
    "linear": (p0_linear, p1_linear),
    "quadratic": (p0_quadratic, p1_quadratic),
}.items():
    p0, p1 = probabilities
    recovered = np.log(p0 / p1) / (lambda_1 - lambda_0)
    profile_error[name] = float(np.max(np.abs(recovered - profiles[name])))

check(
    "N5 likelihood ratio 精确恢复三种剖面",
    all(error < 1e-10 for error in profile_error.values()),
    f"max_errors={ {name: float(f'{error:.3e}') for name, error in profile_error.items()} }",
)


# ==================================================================
head("N6  常量剖面协方差为零")
# ==================================================================
age_centered = age - np.mean(age)
matrix_imbalance = p1_constant - p0_constant
covariance_constant = float(np.mean(age_centered * matrix_imbalance))
check(
    "N6 常量剖面协方差为零",
    abs(covariance_constant) < 1e-12,
    f"cov={covariance_constant:.3e}",
)


# ==================================================================
head("N7  非恒定剖面协方差非零")
# ==================================================================
covariance_linear = float(
    np.mean(age_centered * (p1_linear - p0_linear))
)
covariance_quadratic = float(
    np.mean(age_centered * (p1_quadratic - p0_quadratic))
)
check(
    "N7 非恒定剖面协方差非零",
    abs(covariance_linear) > 1e-3 and abs(covariance_quadratic) > 1e-3,
    f"cov=({covariance_linear:.6f},{covariance_quadratic:.6f})",
)


# ==================================================================
head("N8  D224 乘积态只给恒定剖面")
# ==================================================================
product_p0 = 0.75 * np.ones_like(age)
product_p1 = 0.25 * np.ones_like(age)
product_profile = np.log(product_p0 / product_p1) / (lambda_1 - lambda_0)
check(
    "N8 乘积态恢复剖面为常数",
    np.max(np.abs(product_profile - product_profile[0])) < 1e-12,
    f"profile={product_profile[0]:.6f}",
)


# ==================================================================
head("N9  文档登记剖面与相关性等价")
# ==================================================================
check(
    "N9 文档登记等价关系",
    "R-Z-PROFILE-AS-MATRIX-AGE-CORRELATION" in DOC
    and "模密度剖面等价于矩阵相位的年龄 likelihood ratio" in DOC,
)


# ==================================================================
head("N10  文档登记矩阵年龄 likelihood 缺口")
# ==================================================================
check(
    "N10 文档登记 likelihood 缺口",
    "R-Z-MATRIX-AGE-LIKELIHOOD-GAP" in DOC
    and "年龄条件矩阵比例" in DOC,
)


# ==================================================================
head("N11  文档不把相关性来源写成已导出")
# ==================================================================
check(
    "N11 文档保留未解来源",
    "本文没有导出其中任何一项" in DOC
    and "仍未导出" in DOC,
)


# ==================================================================
head("N12  上游边界保持")
# ==================================================================
check(
    "N12 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)


# ==================================================================
head("N13  文档不引用项目外体系或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N13 D232 文档不含项目外体系名或外部路径",
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
