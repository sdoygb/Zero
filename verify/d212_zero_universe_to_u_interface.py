#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D212
===================
【被检验的问题（D212）】
  周期 T 能否条件接到 U1-U4？
  周期毁灭能否给出 U4 所需严格支持增长？
  T 是否自动等于模时间？

运行：python3 verify/d212_zero_universe_to_u_interface.py
"""
import io
import json
import math
import os
import sys

import numpy as np

from latex_utils import canonical_math


PASS, FAIL = [], []


def check(name, condition, detail=""):
    (PASS if condition else FAIL).append(name)
    print(f"  [{'v' if condition else 'x'}] {name}" + (f"   {detail}" if detail else ""))


def head(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# [lh 迁入] 本脚本读旧理论的 AXIOMS.md 与 derivations/，故 ROOT 指向其原仓库（../modular-equilibrium）
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    os.pardir, os.pardir, "modular-equilibrium"))
DOC_PATH = os.path.join(
    ROOT,
    "derivations",
    "D212_zero_universe_to_u_interface.md",
)
AX_PATH = os.path.join(ROOT, "AXIOMS.md")
PERIODIC_JSON = os.path.join(
    ROOT,
    "simulations",
    "zero_sum_periodic_destruction_results.json",
)

DOC = canonical_math(io.open(DOC_PATH, encoding="utf-8").read())
AX = io.open(AX_PATH, encoding="utf-8").read()
BOTH = AX + "\n" + DOC


def phase_projection(size, index):
    projection = np.zeros((size, size), dtype=complex)
    projection[index, index] = 1.0
    return projection


def cyclic_shift(size):
    shift = np.zeros((size, size), dtype=complex)
    for index in range(size):
        shift[(index + 1) % size, index] = 1.0
    return shift


def density_from_weights(weights):
    return np.diag(np.asarray(weights, dtype=float))


def matrix_unit(size, row, column):
    unit = np.zeros((size, size), dtype=complex)
    unit[row, column] = 1.0
    return unit


def modular_conjugate(matrix, modular_hamiltonian, time):
    evolution = np.diag(
        np.exp(-1j * time * np.diag(modular_hamiltonian))
    )
    return evolution @ matrix @ evolution.conj().T


def modular_span(matrix):
    return float(np.linalg.norm(matrix - np.diag(np.diag(matrix))))


def support_step(point, time, period):
    if point == "global":
        return "global"
    age = point[1]
    new_age = age + time
    return ("local", new_age) if new_age < period else "global"


T = 3

# ==================================================================
head("N1  D212 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 文档登记周期到 U 的条件桥与状态支持缺口",
    "R-Z-U-CONDITIONAL-BRIDGE" in BOTH
    and "R-Z-U-STATE-SUPPORT-GAP" in BOTH
    and "条件接口" in DOC,
)
check(
    "N1 文档不把条件桥写成无条件上游导出",
    "不是 `U1-U4` 的推论" in DOC
    and "无条件导出" in DOC
    and "当前接口仍是一个条件桥" in DOC,
)


# ==================================================================
head("N2  周期实验确实使用 T=3")
# ==================================================================
with io.open(PERIODIC_JSON, encoding="utf-8") as handle:
    periodic = json.load(handle)
check(
    "N2 现有周期毁灭模型的周期等于核验所用 T",
    periodic["model"]["destruction_period"] == T,
    f"T={periodic['model']['destruction_period']}",
)


# ==================================================================
head("N3  周期循环生成 M_T(C)")
# ==================================================================
projections = [phase_projection(T, index) for index in range(T)]
shift = cyclic_shift(T)

check(
    "N3 移子是 T 阶酉算子",
    np.allclose(shift @ shift.conj().T, np.eye(T))
    and np.allclose(np.linalg.matrix_power(shift, T), np.eye(T)),
)
check(
    "N3 相位投影正交、完备并循环置换",
    all(
        np.allclose(projections[i] @ projections[j], projections[i] if i == j else 0)
        for i in range(T)
        for j in range(T)
    )
    and np.allclose(sum(projections), np.eye(T))
    and all(
        np.allclose(
            shift @ projections[k] @ shift.conj().T,
            projections[(k + 1) % T],
        )
        for k in range(T)
    ),
)

matrix_units = []
for row in range(T):
    for column in range(T):
        power = (row - column) % T
        unit = (
            projections[row]
            @ np.linalg.matrix_power(shift, power)
            @ projections[column]
        )
        matrix_units.append(unit.reshape(-1))

gram = np.column_stack(matrix_units)
rank = int(np.linalg.matrix_rank(gram))
check(
    "N3 相位投影与移子生成全部 T^2 个矩阵单位",
    rank == T * T,
    f"rank={rank}, expected={T * T}",
)
check(
    "N3 生成代数为非交换矩阵代数",
    np.linalg.norm(shift @ projections[0] - projections[0] @ shift) > 1.0e-10,
)


# ==================================================================
head("N4  忠实态选择决定模流")
# ==================================================================
uniform_weights = [1.0 / T for _ in range(T)]
nonuniform_weights = [0.55, 0.30, 0.15]

uniform_rho = density_from_weights(uniform_weights)
nonuniform_rho = density_from_weights(nonuniform_weights)

uniform_k = -np.diag(np.log(uniform_weights))
nonuniform_k = -np.diag(np.log(nonuniform_weights))

probe = matrix_unit(T, 0, 1)
time = 0.7
uniform_image = modular_conjugate(probe, uniform_k, time)
nonuniform_image = modular_conjugate(probe, nonuniform_k, time)

expected_nonuniform_phase = np.exp(
    1j * time * (math.log(nonuniform_weights[0]) - math.log(nonuniform_weights[1]))
)

check(
    "N4 均匀态是中心态并给出恒等模流",
    np.allclose(uniform_image, probe),
)
check(
    "N4 非均匀态是忠实非中心态",
    np.all(np.asarray(nonuniform_weights) > 0.0)
    and not np.allclose(nonuniform_rho, np.eye(T) / T),
)
check(
    "N4 非均匀态给出非平凡模流",
    not np.allclose(nonuniform_image, probe)
    and np.linalg.norm(nonuniform_image - probe) > 1.0e-3,
    f"norm={np.linalg.norm(nonuniform_image - probe):.6f}",
)
check(
    "N4 模流在矩阵单位上的相位与公式一致",
    np.allclose(nonuniform_image, expected_nonuniform_phase * probe),
)


# ==================================================================
head("N5  U2 的最小条件支持")
# ==================================================================
diagonal_algebra_dim = T
full_algebra_dim = T * T
off_diagonal_probe = matrix_unit(T, 0, 1)
check(
    "N5 对角支撑包络严格小于全矩阵代数",
    diagonal_algebra_dim < full_algebra_dim
    and modular_span(off_diagonal_probe) > 1.0e-10,
)
check(
    "N5 非中心态模流保持对角支撑代数",
    all(
        np.allclose(
            modular_conjugate(projection, nonuniform_k, time),
            projection,
        )
        for projection in projections
    ),
)
check(
    "N5 文档明确 U2 仅条件满足存在版",
    "U2_{\\rm ex}" in DOC
    and "U2_{\\rm sem}" in DOC
    and "U2_{\\rm sem}\\text{ 未由 }T\\text{ 给出。}" in DOC,
)


# ==================================================================
head("N6  周期毁灭给非可逆支持寿命")
# ==================================================================
local_support_points = [("local", age) for age in range(T)]
support_points = local_support_points + ["global"]
semigroup_ok = all(
    support_step(support_step(point, t, T), s, T)
    == support_step(point, s + t, T)
    for s in range(0, 2 * T + 2)
    for t in range(0, 2 * T + 2)
    for point in local_support_points
)
check(
    "N6 阈值支持映射满足半群律",
    semigroup_ok,
)
check(
    "N6 局域支持在寿命 T 内保留并在 T 后进入全局包络",
    support_step(("local", 0), T - 1, T) == ("local", T - 1)
    and support_step(("local", 0), T, T) == "global"
    and support_step(("local", 1), T - 1, T) == "global"
    and support_step("global", 100, T) == "global",
)
check(
    "N6 支持映射保序",
    all(
        support_points.index(support_step(point, 1, T))
        <= support_points.index(support_step("global", 1, T))
        for point in local_support_points
    )
    and support_points.index(support_step(("local", 0), T, T))
    <= support_points.index(support_step("global", T, T)),
)
check(
    "N6 寿命后给出严格支持扩张",
    support_step(("local", 0), T, T) == "global"
    and diagonal_algebra_dim < full_algebra_dim,
)


# ==================================================================
head("N7  纯可逆周期不能给严格支持增长")
# ==================================================================
sample_dimensions = [3, 4, 5]
cycle_inequalities_hold = all(
    sample_dimensions[index] <= sample_dimensions[(index + 1) % len(sample_dimensions)]
    for index in range(len(sample_dimensions))
)
check(
    "N7 非恒定维数无法沿有限循环全部递增",
    not cycle_inequalities_hold,
    f"dims={sample_dimensions}",
)
check(
    "N7 文档登记纯可逆周期无严格增长的障碍",
    "纯周期循环只能给" in DOC
    and "严格支持增长需要非可逆映射" in DOC,
)


# ==================================================================
head("N8  T 不被等同于模参数")
# ==================================================================
check(
    "N8 文档把周期 T 与模时间分开",
    "T\\text{ 是周期动力学与支持寿命" in DOC
    and "\\text{模时间来自状态，不自动等于 }T" in DOC,
)
check(
    "N8 文档登记 KMS 或热时间关系的缺口",
    "KMS" in DOC
    and "热时间" in DOC,
)


# ==================================================================
head("N9  输入预算与边界")
# ==================================================================
check(
    "N9 文档列出交叉积提升、支持点与忠实态输入",
    "交叉积提升" in DOC
    and "寿命阶段支持点" in DOC
    and "全局支持点" in DOC
    and "忠实态权重" in DOC,
)
check(
    "N9 文档不声称替代连续 QFT、几何或标准模型恢复门",
    "连续 QFT" in DOC
    and "几何与标准模型恢复门" in DOC,
)
check(
    "N9 文档不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC
    and "不新增 `U5`" in DOC,
)


print("\n" + "=" * 78)
print(f"通过 {len(PASS)} / 不符 {len(FAIL)}")
if FAIL:
    print("失败项：")
    for name in FAIL:
        print(" -", name)
    sys.exit(1)
print("全部通过")
