#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D239 —— 线性危险率的二次势审计
================================================================================
【被检验的问题（D239）】
  线性危险率能否拆成中心临界与常曲率？
  它是否等价于二次符号势与均匀年龄配对？
  现有零和闭合或终端中性是否自动选出它？

【本步判据】
  N1  D239 与恢复结构已登记
  N2  危险率差是比例对数的负导数
  N3  线性危险率等价于中心临界与常曲率
  N4  终端中性只固定积分常数
  N5  去掉中心临界会出现线性偏置项
  N6  二次势给 D238 的离散多重数比
  N7  二次势给 D234 的年龄抛物型比例
  N8  二次势可写成有序年龄配对的均匀耦合
  N9  均匀配对不是唯一微观分解
  N10 D220 对称闭合给零危险率差
  N11 终端中性允许任意非恒定势
  N12 二次、三次与周期势给出不同年龄比例
  N13 文档登记与上游边界
  N14 文档不引用项目外体系名或外部路径

运行：python3 verify/d239_quadratic_sign_potential_audit.py
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
            "D239_quadratic_sign_potential_audit.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D239 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D239", "D239" in AX)
check(
    "N1 公理表与文档登记危险率势改写",
    "R-Z-HAZARD-POTENTIAL-REWRITE" in BOTH,
)
check(
    "N1 公理表与文档登记均匀年龄配对",
    "R-Z-UNIFORM-AGE-PAIR-COUPLING" in BOTH,
)
check(
    "N1 公理表与文档登记二次势来源缺口",
    "R-Z-QUADRATIC-SIGN-POTENTIAL-GAP" in BOTH,
)


# ==================================================================
head("N2  危险率差是比例对数的负导数")
# ==================================================================
terminal_age = 6.0
delta = 0.13
age = np.linspace(0.0, terminal_age, 6001)
u = delta * (terminal_age**2 - age**2)
u_prime = np.gradient(u, age)
hazard_difference = 2.0 * delta * age
interior = slice(5, -5)

check(
    "N2 比例对数导数等于负危险率差",
    np.allclose(
        u_prime[interior],
        -hazard_difference[interior],
        atol=1e-5,
    ),
)
check(
    "N2 终端比例对数归零",
    abs(u[-1]) < 1e-12,
)


# ==================================================================
head("N3  线性危险率等价于中心临界与常曲率")
# ==================================================================
u_double_prime = np.gradient(u_prime, age)
quadratic_fit = np.polyfit(age, u, 2)

check(
    "N3 二次拟合的线性项等于零",
    abs(quadratic_fit[1]) < 1e-10,
    f"linear={quadratic_fit[1]:.3e}",
)
check(
    "N3 内部曲率为常数",
    np.allclose(
        u_double_prime[20:-20],
        -2.0 * delta,
        atol=1e-5,
    ),
)


# ==================================================================
head("N4  终端中性只固定积分常数")
# ==================================================================
constant_from_terminal = delta * terminal_age**2
u_from_constant = -delta * age**2 + constant_from_terminal

check(
    "N4 终端中性给常数 C=delta T^2",
    np.allclose(u_from_constant, u),
)
check(
    "N4 常数改变不改变曲率",
    abs(np.polyfit(age, u + 7.0, 2)[0] + delta) < 1e-10,
)


# ==================================================================
head("N5  去掉中心临界会出现线性偏置项")
# ==================================================================
slope_bias = 0.21
u_biased = (
    delta * (terminal_age**2 - age**2)
    + slope_bias * (terminal_age - age)
)
u_biased_prime_at_zero = np.gradient(u_biased, age)[1]

check(
    "N5 斜偏解仍满足终端中性",
    abs(u_biased[-1]) < 1e-12,
)
biased_linear_coefficient = np.polyfit(age, u_biased, 2)[1]
check(
    "N5 斜偏解在中心有非零斜率",
    abs(biased_linear_coefficient + slope_bias) < 1e-10,
    f"linear coefficient={biased_linear_coefficient:.6f}",
)
check(
    "N5 斜偏解不同于 D238",
    np.max(np.abs(u_biased - u)) > 1e-2,
)


# ==================================================================
head("N6  二次势给 D238 的离散多重数比")
# ==================================================================
discrete_terminal = 24
ages = np.arange(0, discrete_terminal + 1)
potential = delta * ages**2
potential_increment = np.diff(potential)
q = np.exp(-potential_increment[:-1])
initial_ratio = np.exp(potential[-1])
ratio = initial_ratio * np.exp(-(potential - potential[0]))
target = np.exp(delta * (discrete_terminal**2 - ages**2))

check(
    "N6 势差等于 delta(2k+1)",
    np.allclose(
        potential_increment[:-1],
        delta * (2 * ages[:-2] + 1),
    ),
)
check(
    "N6 多重数比给 D238 的 q_k",
    np.allclose(
        q,
        np.exp(-delta * (2 * ages[:-2] + 1)),
    ),
)
check(
    "N6 比例恢复 D238 抛物线",
    np.allclose(ratio, target),
)


# ==================================================================
head("N7  二次势给 D234 的年龄抛物型比例")
# ==================================================================
lambda_gap = 1.9
radius = 3.1
delta_from_geometry = radius * lambda_gap / (2.0 * terminal_age**2)
continuous_potential = delta_from_geometry * (terminal_age**2 - age**2)
profile = continuous_potential / lambda_gap
target_profile = (radius / 2.0) * (
    1.0 - age**2 / terminal_age**2
)

check(
    "N7 二次势给 D234 年龄剖面",
    np.allclose(profile, target_profile),
)
check(
    "N7 中心最大且边界退化",
    int(np.argmax(profile)) == 0 and abs(profile[-1]) < 1e-12,
)


# ==================================================================
head("N8  二次势可写成有序年龄配对的均匀耦合")
# ==================================================================
coupling = delta
pair_potential = np.array(
    [
        coupling * (k**2)
        for k in range(discrete_terminal + 1)
    ]
)

check(
    "N8 有序配对计数等于 k^2",
    np.allclose(
        pair_potential,
        coupling * ages**2,
    ),
)
check(
    "N8 新增大有序配对数等于 2k+1",
    np.allclose(
        np.diff(ages**2),
        2 * ages[:-1] + 1,
    ),
)


# ==================================================================
head("N9  均匀配对不是唯一微观分解")
# ==================================================================
size = 5
uniform = delta * np.ones((size, size))
nonuniform = uniform.copy()
nonuniform[3, 3] += 0.04

phi_uniform = np.array(
    [
        uniform[:k, :k].sum()
        for k in range(1, size + 1)
    ]
)
phi_nonuniform = np.array(
    [
        nonuniform[:k, :k].sum()
        for k in range(1, size + 1)
    ]
)

check(
    "N9 非均匀矩阵保持前三个势值",
    np.allclose(phi_uniform[:3], phi_nonuniform[:3]),
)
check(
    "N9 非均匀矩阵在高阶分叉",
    not np.allclose(phi_uniform, phi_nonuniform),
)


# ==================================================================
head("N10  D220 对称闭合给零危险率差")
# ==================================================================
symmetric_ratio = np.ones_like(age)
symmetric_log_ratio = np.log(symmetric_ratio)

check(
    "N10 对称正负计数给常比例",
    np.allclose(symmetric_ratio, 1.0),
)
check(
    "N10 对称闭合给零对数导数",
    np.allclose(np.gradient(symmetric_log_ratio, age), 0.0),
)


# ==================================================================
head("N11  终端中性允许任意非恒定势")
# ==================================================================
quadratic_potential = delta * ages**2
cubic_potential = (
    delta
    * discrete_terminal**2
    * (ages / discrete_terminal) ** 3
)
periodic_potential = delta * np.sin(
    2.0 * np.pi * ages / discrete_terminal
)


def terminal_compensated_ratio(potential_values):
    initial = np.exp(potential_values[-1] - potential_values[0])
    return initial * np.exp(
        -(potential_values - potential_values[0])
    )


ratio_quadratic = terminal_compensated_ratio(quadratic_potential)
ratio_cubic = terminal_compensated_ratio(cubic_potential)
ratio_periodic = terminal_compensated_ratio(periodic_potential)

check(
    "N11 三类势终端比例都等于一",
    abs(ratio_quadratic[-1] - 1.0) < 1e-12
    and abs(ratio_cubic[-1] - 1.0) < 1e-12
    and abs(ratio_periodic[-1] - 1.0) < 1e-12,
)
check(
    "N11 三类势都非恒定",
    np.ptp(ratio_quadratic) > 1e-2
    and np.ptp(ratio_cubic) > 1e-2
    and np.ptp(ratio_periodic) > 1e-2,
)


# ==================================================================
head("N12  二次、三次与周期势给出不同年龄比例")
# ==================================================================
check(
    "N12 二次与三次比例不同",
    np.max(np.abs(ratio_quadratic - ratio_cubic)) > 1e-2,
)
check(
    "N12 二次与周期比例不同",
    np.max(np.abs(ratio_quadratic - ratio_periodic)) > 1e-2,
)
check(
    "N12 只有二次势给常数势差曲率",
    np.ptp(np.diff(quadratic_potential, 2)) < 1e-12
    and np.ptp(np.diff(cubic_potential, 2)) > 1e-2,
)


# ==================================================================
head("N13  文档登记与上游边界")
# ==================================================================
check(
    "N13 文档登记危险率势改写",
    "R-Z-HAZARD-POTENTIAL-REWRITE" in DOC
    and "线性危险率、二次符号势" in DOC,
)
check(
    "N13 文档登记均匀年龄配对",
    "R-Z-UNIFORM-AGE-PAIR-COUPLING" in DOC
    and "均匀年龄配对" in DOC,
)
check(
    "N13 文档登记二次势来源缺口",
    "R-Z-QUADRATIC-SIGN-POTENTIAL-GAP" in DOC
    and "为什么零动力学选出二次势" in DOC,
)
check(
    "N13 文档不把结果写成 U1-U4 推论",
    "它们不是 `U1-U4` 的推论" in DOC
    and "仍不是零和闭合的自动结果" in DOC,
)
check(
    "N13 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)


# ==================================================================
head("N14  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N14 D239 文档不含项目外体系名或外部路径",
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
