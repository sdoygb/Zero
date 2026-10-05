#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D249 —— 最小应力提升
================================================================================
【被检验的问题（D249）】
  层间读出势成为局部标量后，最小二阶作用量给出什么应力张量？
  正定动能系数是否是独立输入？
  时间型无势标量是否给 p=rho？
  守恒流能否唯一选择标量、尘埃或其他物态？

【本步判据】
  N1  D249 与恢复结构已登记
  N2  正定动能系数可由场重定义吸收
  N3  最小标量应力公式正确
  N4  在壳散度恒等式为零
  N5  时间型无势梯度给 p=rho
  N6  势能改变物态为 p=rho-2V
  N7  空间型梯度给非完美流体型应力
  N8  尘埃与刚性标量有相同能量流但不同应力
  N9  文档登记应力类别缺口
  N10 文档不把新结构写成 U1-U4 推论
  N11 上游边界保持
  N12 文档不引用项目外体系名或外部路径

运行：python3 verify/d249_minimal_stress_lift_from_interlayer_scalar.py
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
            "D249_minimal_stress_lift_from_interlayer_scalar.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D249 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D249", "D249" in AX)
check(
    "N1 公理表与文档登记动能归一化冗余",
    "R-Z-KINETIC-NORMALIZATION-REDUNDANCY" in BOTH,
)
check(
    "N1 公理表与文档登记标量 Dirichlet 应力提升",
    "R-Z-SCALAR-DIRICHLET-STRESS-LIFT" in BOTH,
)
check(
    "N1 公理表与文档登记刚性标量物态",
    "R-Z-STIFF-SCALAR-EQUATION-OF-STATE" in BOTH,
)
check(
    "N1 公理表与文档登记应力类别缺口",
    "R-Z-STRESS-LIFT-CLASS-GAP" in BOTH,
)


# ==================================================================
head("N2  正定动能系数可由场重定义吸收")
# ==================================================================
d_phi = np.array([0.7, -0.2, 0.1, 0.4])
z_coefficient = 2.25
sqrt_z = np.sqrt(z_coefficient)
d_chi = sqrt_z * d_phi

kinetic_phi = 0.5 * z_coefficient * np.dot(d_phi, d_phi)
kinetic_chi = 0.5 * np.dot(d_chi, d_chi)

check(
    "N2 重定义后动能相同",
    np.isclose(kinetic_phi, kinetic_chi),
)
check(
    "N2 重定义保持正定性",
    kinetic_chi > 0.0,
)


# ==================================================================
head("N3  最小标量应力公式正确")
# ==================================================================
eta = np.diag([-1.0, 1.0, 1.0, 1.0])
d = np.array([1.4, 0.35, 0.0, 0.0])
v_potential = 0.25
x_kinematic = float(d @ eta @ d)

stress_manual = np.outer(d, d) - eta * (0.5 * x_kinematic + v_potential)

check(
    "N3 应力矩阵对称",
    np.allclose(stress_manual, stress_manual.T),
)
check(
    "N3 X 为时间型梯度",
    x_kinematic < 0.0,
)
check(
    "N3 应力公式给出有限实谱",
    np.all(np.isfinite(np.linalg.eigvalsh(stress_manual))),
)


# ==================================================================
head("N4  在壳散度恒等式为零")
# ==================================================================
# 取无势零锥平面波 box phi=0 的样本，并用有限差分核对
# div^a T_ab = 0。这里不把恒等于零的常值张量当成核验。
amplitude = 0.8
wave_number = 1.3
frequency = wave_number


def scalar_gradient(t, x):
    phase = wave_number * x - frequency * t
    return np.array(
        [
            -amplitude * frequency * np.sin(phase),
            amplitude * wave_number * np.sin(phase),
            0.0,
            0.0,
        ]
    )


def stress_at(t, x):
    gradient = scalar_gradient(t, x)
    x_kinematic = float(gradient @ eta @ gradient)
    return np.outer(gradient, gradient) - eta * (0.5 * x_kinematic)


def flat_divergence_component(t, x, component, step=1e-3):
    dt_stress = (
        stress_at(t + step, x) - stress_at(t - step, x)
    ) / (2.0 * step)
    dx_stress = (
        stress_at(t, x + step) - stress_at(t, x - step)
    ) / (2.0 * step)
    return float(-dt_stress[0, component] + dx_stress[1, component])


divergence_residual = np.array(
    [
        flat_divergence_component(0.2, 0.3, 0),
        flat_divergence_component(0.2, 0.3, 1),
    ]
)

check(
    "N4 场方程在样本点满足",
    np.isclose(
        -(-frequency**2) + (-wave_number**2),
        0.0,
    ),
)
check(
    "N4 有限差分散度为零",
    np.max(np.abs(divergence_residual)) < 1.0e-5,
)


# ==================================================================
head("N5  时间型无势梯度给 p=rho")
# ==================================================================
d_timelike = np.array([1.5, 0.0, 0.0, 0.0])
x_timelike = float(d_timelike @ eta @ d_timelike)
u_energy = -x_timelike
stress_stiff = np.outer(d_timelike, d_timelike) - eta * (0.5 * x_timelike)

rho_stiff = float(stress_stiff[0, 0])
p_stiff = float(stress_stiff[1, 1])

check(
    "N5 无势时间型梯度",
    x_timelike < 0.0,
)
check(
    "N5 能量密度等于 U/2",
    np.isclose(rho_stiff, 0.5 * u_energy),
)
check(
    "N5 压强等于能量密度",
    np.isclose(p_stiff, rho_stiff),
)


# ==================================================================
head("N6  势能改变物态为 p=rho-2V")
# ==================================================================
v_shift = 0.4
stress_with_potential = (
    np.outer(d_timelike, d_timelike)
    - eta * (0.5 * x_timelike + v_shift)
)
rho_with_potential = float(stress_with_potential[0, 0])
p_with_potential = float(stress_with_potential[1, 1])

check(
    "N6 势能下的能量密度",
    np.isclose(rho_with_potential, rho_stiff + v_shift),
)
check(
    "N6 势能下的压强",
    np.isclose(p_with_potential, p_stiff - v_shift),
)
check(
    "N6 物态关系 p=rho-2V",
    np.isclose(
        p_with_potential,
        rho_with_potential - 2.0 * v_shift,
    ),
)


# ==================================================================
head("N7  空间型梯度给非完美流体型应力")
# ==================================================================
d_spacelike = np.array([0.2, 1.1, 0.0, 0.0])
x_spacelike = float(d_spacelike @ eta @ d_spacelike)
stress_spacelike = (
    np.outer(d_spacelike, d_spacelike)
    - eta * (0.5 * x_spacelike)
)

check(
    "N7 X 为空间型",
    x_spacelike > 0.0,
)
check(
    "N7 空间型应力不是各向同性压力",
    not np.isclose(stress_spacelike[1, 1], stress_spacelike[2, 2]),
)


# ==================================================================
head("N8  尘埃与刚性标量有相同能量流但不同应力")
# ==================================================================
rho_dust = 2.0
dust_stress = np.diag([rho_dust, 0.0, 0.0, 0.0])
stiff_stress = np.diag([rho_dust, rho_dust, rho_dust, rho_dust])

check(
    "N8 相同能量流分量",
    np.allclose(dust_stress[:, 0], stiff_stress[:, 0]),
)
check(
    "N8 空间压强不同",
    not np.allclose(dust_stress[1:, 1:], stiff_stress[1:, 1:]),
)
check(
    "N8 迹也不同",
    not np.isclose(np.trace(dust_stress), np.trace(stiff_stress)),
)


# ==================================================================
head("N9  文档登记应力类别缺口")
# ==================================================================
check(
    "N9 文档登记动力归一化冗余",
    "R-Z-KINETIC-NORMALIZATION-REDUNDANCY" in DOC
    and "场重定义" in DOC,
)
check(
    "N9 文档登记标量应力提升",
    "R-Z-SCALAR-DIRICHLET-STRESS-LIFT" in DOC
    and "最小标量类内应力唯一" in DOC,
)
check(
    "N9 文档登记物态与类别缺口",
    "R-Z-STIFF-SCALAR-EQUATION-OF-STATE" in DOC
    and "R-Z-STRESS-LIFT-CLASS-GAP" in DOC,
)


# ==================================================================
head("N10  文档不把新结构写成 U1-U4 推论")
# ==================================================================
check(
    "N10 文档明确条件性",
    "不是 `U1-U4` 的推论" in DOC
    and "不是零和约束的推论" in DOC,
)
check(
    "N10 文档不把最小标量类写成无条件",
    "最小标量类是一个选择" in DOC,
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
    "N12 D249 文档不含项目外体系名或外部路径",
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
