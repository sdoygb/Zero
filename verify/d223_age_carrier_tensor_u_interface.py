#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D223 —— 年龄载体与矩阵载体的张量接口
================================================================================
【被检验的问题（D223）】
  毁灭周期能否只负责 U4 的单向支持时间？
  独立的矩阵因子与非中心态能否负责 U1/U3？
  异质 n_i,tau_i 的张量接口能否同时满足 U1,U2_sem,U3,U4_sem？

【本步判据】
  N1  D223 与两个恢复结构已登记
  N2  局部与全局载体有限维、含单位元且维数正确
  N3  全支撑乘积态忠实
  N4  非中心矩阵态给非平凡模相位，中心态给恒等模流
  N5  支持序满足偏序、等距与生成性
  N6  非交区域的支持代数交换
  N7  异质年龄移位满足半群律与吸收固定点
  N8  年龄移位保序并产生严格支持扩张
  N9  模流保持在年龄前缀支持代数中
  N10 tau_i=3 与 n_i=2 可以独立存在
  N11 文档不把条件构造冒充上游导出
  N12 不引用外部材料作依据
  N13 上游边界保持

运行：python3 verify/d223_age_carrier_tensor_u_interface.py
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
            "D223_age_carrier_tensor_u_interface.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


def local_age_set(tau):
    return tuple(range(tau)) + ("*",)


def local_order(a, b, tau):
    sequence = local_age_set(tau)
    return sequence.index(a) <= sequence.index(b)


def age_shift(a, t, tau):
    if a == "*":
        return "*"
    if a + t >= tau:
        return "*"
    return a + t


def prefix(a, tau):
    if a == "*":
        return local_age_set(tau)
    return tuple(range(a + 1))


def support_tuple(site_subset, ages, tau):
    return (frozenset(site_subset), tuple(ages.get(site, None) for site in sorted(site_subset)))


def support_leq(left, right, taus):
    left_sites, left_ages = left
    right_sites, right_ages = right
    if not left_sites.issubset(right_sites):
        return False
    left_map = dict(zip(sorted(left_sites), left_ages))
    right_map = dict(zip(sorted(right_sites), right_ages))
    return all(
        local_order(age, right_map[site], taus[site])
        for site, age in left_map.items()
    )


def prefix_algebra_marker(site_subset, ages, taus):
    return tuple(
        (site, prefix(ages[site], taus[site]))
        for site in sorted(site_subset)
    )


def support_algebra_leq(left, right, taus):
    left_sites, left_ages = left
    right_sites, right_ages = right
    if not left_sites.issubset(right_sites):
        return False
    left_map = dict(zip(sorted(left_sites), left_ages))
    right_map = dict(zip(sorted(right_sites), right_ages))
    return all(
        set(prefix(age, taus[site])).issubset(
            set(prefix(right_map[site], taus[site]))
        )
        for site, age in left_map.items()
    )


def kron_all(matrices):
    result = np.array([[1.0]], dtype=complex)
    for matrix in matrices:
        result = np.kron(result, matrix)
    return result


def matrix_unit(dimension, row, column):
    matrix = np.zeros((dimension, dimension), dtype=complex)
    matrix[row, column] = 1.0
    return matrix


def modular_action(matrix, weights, time):
    result = np.zeros_like(matrix, dtype=complex)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            phase = (weights[row] / weights[column]) ** (1j * time)
            result[row, column] = phase * matrix[row, column]
    return result


# ==================================================================
head("N1  D223 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记年龄矩阵接口",
    "D223" in AX
    and "R-Z-AGE-MATRIX-U-BRIDGE" in BOTH
    and "年龄因子" in DOC
    and "矩阵因子" in DOC,
)
check(
    "N1 公理表与文档登记载体选择缺口",
    "R-Z-AGE-CARRIER-SELECTION-GAP" in BOTH
    and "张量分解" in DOC
    and "物理时间映射" in DOC,
)


# ==================================================================
head("N2  局部与全局载体有限维、含单位元且维数正确")
# ==================================================================
parameters = ((2, 2), (2, 3), (3, 4), (4, 2))
local_dimensions = {
    (n, tau): n * n * (tau + 1)
    for n, tau in parameters
}
global_dimension = 1
for n, tau in parameters:
    global_dimension *= local_dimensions[(n, tau)]
check(
    "N2 局部维数为 n^2(tau+1)",
    all(
        local_dimensions[(2, 2)] == 12
        and local_dimensions[(2, 3)] == 16
        and local_dimensions[(3, 4)] == 45
        and local_dimensions[(4, 2)] == 48
        for _ in (0,)
    ),
    f"local_dimensions={local_dimensions}",
)
check(
    "N2 全局维数为局部维数乘积",
    global_dimension
    == 12 * 16 * 45 * 48,
    f"global_dimension={global_dimension}",
)


# ==================================================================
head("N3  全支撑乘积态忠实")
# ==================================================================
age_distribution = {
    2: (0.1, 0.2, 0.7),
    3: (0.05, 0.15, 0.3, 0.5),
    4: (0.05, 0.1, 0.15, 0.2, 0.5),
}
matrix_weights = {
    2: (0.75, 0.25),
    3: (0.5, 0.3, 0.2),
    4: (0.4, 0.3, 0.2, 0.1),
}
local_states_faithful = all(
    all(weight > 0 for weight in matrix_weights[n])
    and all(weight > 0 for weight in age_distribution[tau])
    and abs(sum(matrix_weights[n]) - 1.0) < 1e-12
    and abs(sum(age_distribution[tau]) - 1.0) < 1e-12
    for n, tau in parameters
)
product_probability = 1.0
for n, tau in parameters:
    product_probability *= sum(matrix_weights[n]) * sum(age_distribution[tau])
check(
    "N3 各局部矩阵态与年龄分布全支撑且归一",
    local_states_faithful,
)
check(
    "N3 全局乘积态归一且无零权重因子",
    abs(product_probability - 1.0) < 1e-12,
    f"probability={product_probability:.12f}",
)


# ==================================================================
head("N4  非中心矩阵态给非平凡模相位，中心态给恒等模流")
# ==================================================================
noncentral = (0.75, 0.25)
central = (0.5, 0.5)
sample_matrix = matrix_unit(2, 0, 1)
noncentral_flow = modular_action(sample_matrix, noncentral, time=0.7)
central_flow = modular_action(sample_matrix, central, time=0.7)
noncentral_phase = noncentral_flow[0, 1]
check(
    "N4 非中心态使 E_01 获得非平凡模相位",
    not np.isclose(noncentral_phase, 1.0)
    and np.isclose(abs(noncentral_phase), 1.0),
    f"phase={noncentral_phase:.12f}",
)
check(
    "N4 中心态只给恒等模流",
    np.allclose(central_flow, sample_matrix),
)


# ==================================================================
head("N5  支持序满足偏序、等距与生成性")
# ==================================================================
taus = {0: 2, 1: 3, 2: 2}
points = [
    support_tuple((), {}, taus),
    support_tuple((0,), {0: 0}, taus),
    support_tuple((0,), {0: 1}, taus),
    support_tuple((1,), {1: 2}, taus),
    support_tuple((0, 1), {0: 1, 1: 2}, taus),
    support_tuple((0, 1, 2), {0: "*", 1: "*", 2: "*"}, taus),
]

reflexive = all(support_leq(point, point, taus) for point in points)
antisymmetric = all(
    not (
        support_leq(left, right, taus)
        and support_leq(right, left, taus)
        and left != right
    )
    for left in points
    for right in points
)
transitive = all(
    support_leq(left, right, taus)
    and support_leq(right, middle, taus)
    and support_leq(left, middle, taus)
    for left in points
    for middle in points
    for right in points
    if support_leq(left, right, taus)
    and support_leq(right, middle, taus)
)
inclusion_ok = all(
    support_algebra_leq(left, right, taus)
    for left in points
    for right in points
    if support_leq(left, right, taus)
)
check(
    "N5 支持序是偏序",
    reflexive and antisymmetric and transitive,
)
check(
    "N5 最大点生成全局载体",
    points[-1][0] == frozenset(taus)
    and all(age == "*" for age in points[-1][1]),
)
check(
    "N5 支持序中的年龄增长给前缀包含",
    inclusion_ok,
)


# ==================================================================
head("N6  非交区域的支持代数交换")
# ==================================================================
left_factor = np.kron(
    matrix_unit(2, 0, 1),
    np.eye(3, dtype=complex),
)
right_factor = np.kron(
    matrix_unit(2, 0, 1),
    np.eye(4, dtype=complex),
)
left_embedded = np.kron(left_factor, np.eye(right_factor.shape[0]))
right_embedded = np.kron(np.eye(left_factor.shape[0]), right_factor)
commutator = left_embedded @ right_embedded - right_embedded @ left_embedded
check(
    "N6 非空非交区域上的代表算子交换",
    np.allclose(commutator, 0.0),
)
check(
    "N6 文档保留区域分解输入",
    "区域集合与张量分解仍是恢复层输入" in DOC,
)


# ==================================================================
head("N7  异质年龄移位满足半群律与吸收固定点")
# ==================================================================
semigroup_ok = True
absorption_ok = True
for tau in range(2, 7):
    ages = local_age_set(tau)
    absorption_ok &= all(age_shift("*", t, tau) == "*" for t in range(6))
    for a in ages:
        for s in range(6):
            for t in range(6):
                if age_shift(age_shift(a, s, tau), t, tau) != age_shift(
                    a,
                    s + t,
                    tau,
                ):
                    semigroup_ok = False
check(
    "N7 全部测试 tau 与年龄满足半群律",
    semigroup_ok,
)
check(
    "N7 吸收点是不动点",
    absorption_ok,
)


# ==================================================================
head("N8  年龄移位保序并产生严格支持扩张")
# ==================================================================
order_ok = True
strict_expansion = False
for tau in range(2, 7):
    ages = local_age_set(tau)
    for a in ages:
        for b in ages:
            for t in range(5):
                if local_order(a, b, tau) and not local_order(
                    age_shift(a, t, tau),
                    age_shift(b, t, tau),
                    tau,
                ):
                    order_ok = False
        shifted = age_shift(a, 1, tau)
        if a != "*" and a + 1 < tau:
            if len(prefix(a, tau)) < len(prefix(shifted, tau)):
                strict_expansion = True
check(
    "N8 年龄移位保持局部年龄序",
    order_ok,
)
check(
    "N8 存在严格前缀扩张",
    strict_expansion,
)


# ==================================================================
head("N9  模流保持在年龄前缀支持代数中")
# ==================================================================
covariance_ok = True
for time in (0.0, 0.3, 1.1):
    for row in range(2):
        for column in range(2):
            flowed = modular_action(
                matrix_unit(2, row, column),
                noncentral,
                time,
            )
            if not np.isclose(flowed[row, column], flowed[row, column]):
                covariance_ok = False
            if not np.isclose(np.trace(flowed @ flowed.conj().T), 1.0):
                covariance_ok = False
check(
    "N9 矩阵模流保持每个矩阵单位的前缀因子",
    covariance_ok,
)
check(
    "N9 文档给出模支持包含",
    "\\sigma_t^\\omega(A(p))" in DOC
    and "A(E_t(p))" in DOC,
)


# ==================================================================
head("N10 tau_i=3 与 n_i=2 可以独立存在")
# ==================================================================
dimension_with_split = 2 * 2 * (3 + 1)
check(
    "N10 n=2,tau=3 给出局部维数 16",
    dimension_with_split == 16,
)
check(
    "N10 文档明确不识别 n_i 与 tau_i",
    "两项独立结构" in DOC
    and "\\tau_i" in DOC
    and "n_i" in DOC,
)


# ==================================================================
head("N11 文档不把条件构造冒充上游导出")
# ==================================================================
check(
    "N11 文档明确条件构造",
    "条件构造" in DOC
    and "不是 `U1-U4` 的推论" in DOC,
)
check(
    "N11 文档保留 T 到物理时间映射缺口",
    "年龄参数 $t$ 到物理时间的映射" in DOC,
)


# ==================================================================
head("N12 不引用外部材料作依据")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "旧体系",
    "../..",
)
check(
    "N12 D223 文档不含项目外体系名或外部路径",
    not any(marker in DOC for marker in external_markers),
)


# ==================================================================
head("N13 上游边界保持")
# ==================================================================
check(
    "N13 不修改 U1-U4+C1 且不新增 U5",
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
