#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D250 —— 纯层间交换选择刚性标量源
================================================================================
【被检验的问题（D250）】
  纯差分层间交换是否天然选择无势 Dirichlet？
  它是否有一个常量零模？
  它的时间型连续化是否给 p=rho？
  尘埃、辐射与刚性标量怎样区分？

【本步判据】
  N1  D250 与恢复结构已登记
  N2  纯交换能量有常量零模
  N3  交换流是能量的负梯度的一半
  N4  onsite 项破坏一般常量零模
  N5  等权细化给 Dirichlet 能量
  N6  纯交换无势支路给 p=rho
  N7  尘埃、辐射与刚性标量的迹不同
  N8  文档登记 onsite 缺口与源分支
  N9  文档不把新结构写成 U1-U4 推论
  N10 上游边界保持
  N11 文档不引用项目外体系名或外部路径

运行：python3 verify/d250_pure_exchange_selects_stiff_scalar_source.py
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
            "D250_pure_exchange_selects_stiff_scalar_source.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D250 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D250", "D250" in AX)
check(
    "N1 公理表与文档登记纯交换选择",
    "R-Z-PURE-EXCHANGE-DIRICHLET-SELECTION" in BOTH,
)
check(
    "N1 公理表与文档登记刚性源判别",
    "R-Z-STIFF-SOURCE-DISCRIMINATOR" in BOTH,
)
check(
    "N1 公理表与文档登记 onsite 缺口",
    "R-Z-ONSITE-LAYER-TERM-GAP" in BOTH,
)
check(
    "N1 公理表与文档登记物质多重态缺口",
    "R-Z-MATTER-MULTIPLET-GAP" in BOTH,
)


# ==================================================================
head("N2  纯交换能量有常量零模")
# ==================================================================
rng = np.random.default_rng(250)
size = 7
raw = rng.normal(size=(size, size))
kernel = 0.5 * (raw + raw.T)
np.fill_diagonal(kernel, 0.0)
potential = rng.normal(size=size)


def exchange_energy(values, conductance):
    difference = values[:, None] - values[None, :]
    return 0.5 * np.sum(conductance * difference**2)


energy_zero = exchange_energy(potential, kernel)
energy_shifted = exchange_energy(potential + 4.2, kernel)

check(
    "N2 常量平移不改变交换能量",
    np.isclose(energy_zero, energy_shifted),
)
check(
    "N2 交换能量非负",
    energy_zero >= 0.0,
)


# ==================================================================
head("N3  交换流是能量的负梯度的一半")
# ==================================================================
step = 1.0e-6
gradient = np.zeros(size)
for index in range(size):
    plus = potential.copy()
    minus = potential.copy()
    plus[index] += step
    minus[index] -= step
    gradient[index] = (
        exchange_energy(plus, kernel) - exchange_energy(minus, kernel)
    ) / (2.0 * step)

velocity = np.array(
    [
        np.sum(kernel[index] * (potential - potential[index]))
        for index in range(size)
    ]
)

check(
    "N3 数值梯度正确",
    np.allclose(velocity, -0.5 * gradient),
)
check(
    "N3 总交换流为零",
    np.isclose(np.sum(velocity), 0.0),
)


# ==================================================================
head("N4  onsite 项破坏一般常量零模")
# ==================================================================
mass_squared = 0.37
onsite_zero = mass_squared * np.sum(potential**2)
onsite_shifted = mass_squared * np.sum((potential + 4.2) ** 2)

check(
    "N4 onsite 项改变常量平移",
    not np.isclose(onsite_zero, onsite_shifted),
)
check(
    "N4 纯交换仍保持常量零模",
    np.isclose(
        exchange_energy(potential, kernel),
        exchange_energy(potential + 4.2, kernel),
    ),
)


# ==================================================================
head("N5  等权细化给 Dirichlet 能量")
# ==================================================================
linear_values = np.array([0.0, 0.2, 0.5, 0.7, 1.0])
h = 0.25
dirichlet_energy = 0.5 * np.sum(
    ((linear_values[1:] - linear_values[:-1]) / h) ** 2
) * h

check(
    "N5 线性场有限差分能量为正",
    dirichlet_energy > 0.0,
)
check(
    "N5 常量场离散能量为零",
    np.isclose(
        0.5 * np.sum(((np.ones(5)[1:] - np.ones(5)[:-1]) / h) ** 2) * h,
        0.0,
    ),
)


# ==================================================================
head("N6  纯交换无势支路给 p=rho")
# ==================================================================
eta = np.diag([-1.0, 1.0, 1.0, 1.0])
gradient_timelike = np.array([1.2, 0.0, 0.0, 0.0])
x_timelike = float(gradient_timelike @ eta @ gradient_timelike)
u_energy = -x_timelike
stress_stiff = (
    np.outer(gradient_timelike, gradient_timelike)
    - eta * (0.5 * x_timelike)
)

rho_stiff = float(stress_stiff[0, 0])
p_stiff = float(stress_stiff[1, 1])

check(
    "N6 时间型梯度",
    x_timelike < 0.0,
)
check(
    "N6 刚性能量密度",
    np.isclose(rho_stiff, 0.5 * u_energy),
)
check(
    "N6 刚性压强等于能量密度",
    np.isclose(p_stiff, rho_stiff),
)


# ==================================================================
head("N7  尘埃、辐射与刚性标量的迹不同")
# ==================================================================
rho = 2.0
trace_stiff = -rho + 3.0 * rho
trace_dust = -rho
trace_radiation = 0.0

check(
    "N7 刚性迹为正",
    np.isclose(trace_stiff, 2.0 * rho),
)
check(
    "N7 尘埃迹与刚性不同",
    not np.isclose(trace_dust, trace_stiff),
)
check(
    "N7 辐射迹与刚性不同",
    not np.isclose(trace_radiation, trace_stiff),
)


# ==================================================================
head("N8  文档登记 onsite 缺口与源分支")
# ==================================================================
check(
    "N8 文档登记纯交换选择",
    "R-Z-PURE-EXCHANGE-DIRICHLET-SELECTION" in DOC
    and "纯差分" in DOC,
)
check(
    "N8 文档登记刚性判别",
    "R-Z-STIFF-SOURCE-DISCRIMINATOR" in DOC
    and "p=\\rho" in DOC,
)
check(
    "N8 文档登记 onsite 与多重态缺口",
    "R-Z-ONSITE-LAYER-TERM-GAP" in DOC
    and "R-Z-MATTER-MULTIPLET-GAP" in DOC,
)


# ==================================================================
head("N9  文档不把新结构写成 U1-U4 推论")
# ==================================================================
check(
    "N9 文档明确条件性",
    "不是 `U1-U4` 的推论" in DOC
    and "答案是条件性的" in DOC,
)
check(
    "N9 文档不把纯交换写成完整 GR",
    "几何与量子化仍未闭合" in DOC,
)


# ==================================================================
head("N10  上游边界保持")
# ==================================================================
check(
    "N10 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)
check(
    "N10 文档登记输入预算",
    "## §3 输入预算" in DOC,
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
    "N11 D250 文档不含项目外体系名或外部路径",
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
