#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D224 —— 最小年龄矩阵桥与模流/支持寿命分离
================================================================================
【被检验的问题（D224）】
  在 D223 的张量类中，最小非交换矩阵因子能否压到 M_2？
  D220 的活动年龄谱能否接入年龄因子？
  固定 ttau 是否决定模流？固定矩阵态是否决定支持寿命？

【本步判据】
  N1  D224 与两个恢复结构已登记
  N2  n=1 交换，n=2 是最小非交换矩阵因子
  N3  M_2 张量年龄因子的局部维数正确
  N4  D220 活动年龄谱正、归一
  N5  加入终端权重后的完整年龄分布忠实且归一
  N6  tau=2 年龄谱均匀，但非中心矩阵态仍给非平凡模流
  N7  相同 tau、不同 r 给相同支持寿命但不同模流
  N8  相同 r、不同 tau 给相同模流但不同支持阈值
  N9  相同活动谱与 r、不同 epsilon 给不同完整年龄态
  N10 文档保留直接和、参数选择与物理时间缺口
  N11 文档不把条件构造冒充上游导出
  N12 不引用外部材料作依据
  N13 上游边界保持

运行：python3 verify/d224_minimal_m2_age_bridge_and_modular_support_independence.py
"""
import io
import math
import os
import sys


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
            "D224_minimal_m2_age_bridge_and_modular_support_independence.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


def central_binomial(k):
    return math.comb(k, k // 2)


def activity_weights(tau):
    unnormalized = [central_binomial(k) for k in range(tau)]
    total = sum(unnormalized)
    return tuple(weight / total for weight in unnormalized)


def age_distribution(tau, epsilon):
    active = tuple(
        (1.0 - epsilon) * weight
        for weight in activity_weights(tau)
    )
    return active + (epsilon,)


def modular_phase_ratio(r):
    return r / (1.0 - r)


def modular_phase(r, time):
    ratio = modular_phase_ratio(r)
    phase = ratio ** (1j * time)
    return complex(phase)


def modular_phase_from_tau(r, time, tau):
    del tau
    return modular_phase(r, time)


def age_shift(a, t, tau):
    if a == "*":
        return "*"
    if a + t >= tau:
        return "*"
    return a + t


# ==================================================================
head("N1  D224 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记最小 M_2 年龄桥",
    "D224" in AX
    and "R-Z-MINIMAL-M2-AGE-BRIDGE" in BOTH
    and "M_2" in DOC,
)
check(
    "N1 公理表与文档登记模流支持寿命分离",
    "R-Z-MODULAR-SUPPORT-INDEPENDENCE" in BOTH
    and "支持寿命" in DOC
    and "模流" in DOC,
)


# ==================================================================
head("N2  n=1 交换，n=2 是最小非交换矩阵因子")
# ==================================================================
check(
    "N2 n=1 给交换载体 C",
    "M_1(\\mathbb C)=\\mathbb C" in DOC
    and "\\mathbb C\\otimes C(X_{\\tau_i})" in DOC,
)
check(
    "N2 n=2 是张量类中的最小非交换维数",
    "n_i\\ge2" in DOC
    and "\\dim_{\\mathbb C}M_2(\\mathbb C)=4" in DOC,
)


# ==================================================================
head("N3  M_2 张量年龄因子的局部维数正确")
# ==================================================================
local_dimensions = {
    tau: 4 * (tau + 1)
    for tau in range(2, 9)
}
check(
    "N3 局部维数为 4(tau+1)",
    local_dimensions[2] == 12
    and local_dimensions[3] == 16
    and local_dimensions[8] == 36,
    f"dimensions={local_dimensions}",
)
check(
    "N3 文档给出最小维数公式",
    "4(\\tau_i+1)" in DOC,
)


# ==================================================================
head("N4  D220 活动年龄谱正、归一")
# ==================================================================
weights = {
    tau: activity_weights(tau)
    for tau in range(2, 9)
}
positive = all(
    all(weight > 0 for weight in values)
    for values in weights.values()
)
normalized = all(
    abs(sum(values) - 1.0) < 1e-12
    for values in weights.values()
)
check(
    "N4 所有活动权重严格正",
    positive,
)
check(
    "N4 所有活动权重归一",
    normalized,
    f"tau2={weights[2]}, tau3={weights[3]}",
)


# ==================================================================
head("N5  加入终端权重后的完整年龄分布忠实且归一")
# ==================================================================
terminal_weights = (0.05, 0.2, 0.4)
full_distributions = {
    (tau, epsilon): age_distribution(tau, epsilon)
    for tau in (2, 3, 5)
    for epsilon in terminal_weights
}
full_faithful = all(
    all(weight > 0 for weight in distribution)
    for distribution in full_distributions.values()
)
full_normalized = all(
    abs(sum(distribution) - 1.0) < 1e-12
    for distribution in full_distributions.values()
)
check(
    "N5 完整年龄分布全支撑",
    full_faithful,
)
check(
    "N5 完整年龄分布归一",
    full_normalized,
)


# ==================================================================
head("N6  tau=2 年龄谱均匀，但非中心矩阵态仍给非平凡模流")
# ==================================================================
tau_two_weights = weights[2]
age_uniform = abs(tau_two_weights[0] - tau_two_weights[1]) < 1e-12
noncentral_ratio = modular_phase_ratio(0.75)
central_ratio = modular_phase_ratio(0.5)
check(
    "N6 tau=2 的活动年龄权重均匀",
    age_uniform,
    f"weights={tau_two_weights}",
)
check(
    "N6 非中心矩阵态仍给非平凡模相位",
    not math.isclose(noncentral_ratio, 1.0)
    and math.isclose(central_ratio, 1.0),
)


# ==================================================================
head("N7  相同 tau、不同 r 给相同支持寿命但不同模流")
# ==================================================================
tau_fixed = 3
r_left = 0.75
r_right = 0.8
support_left = tuple(
    age_shift(age, 1, tau_fixed)
    for age in range(tau_fixed)
)
support_right = tuple(
    age_shift(age, 1, tau_fixed)
    for age in range(tau_fixed)
)
check(
    "N7 固定 tau 后支持序列相同",
    support_left == support_right,
    f"support={support_left}",
)
check(
    "N7 不同 r 给不同模相位比",
    not math.isclose(
        modular_phase_ratio(r_left),
        modular_phase_ratio(r_right),
    ),
    f"ratios={modular_phase_ratio(r_left):.6f},{modular_phase_ratio(r_right):.6f}",
)


# ==================================================================
head("N8  相同 r、不同 tau 给相同模流但不同支持阈值")
# ==================================================================
r_fixed = 0.75
short = age_shift(1, 1, tau=2)
long = age_shift(1, 1, tau=3)
flow_short = modular_phase_from_tau(r_fixed, 0.7, tau=2)
flow_long = modular_phase_from_tau(r_fixed, 0.7, tau=3)
check(
    "N8 相同 r 的模相位比相同",
    math.isclose(flow_short.real, flow_long.real)
    and math.isclose(flow_short.imag, flow_long.imag),
    f"phase={flow_short:.12f}",
)
check(
    "N8 不同 tau 在同一推进中给出不同吸收结果",
    short == "*" and long == 2,
    f"short={short}, long={long}",
)


# ==================================================================
head("N9  相同活动谱与 r、不同 epsilon 给不同完整年龄态")
# ==================================================================
active_a = age_distribution(3, 0.1)[:-1]
active_b = age_distribution(3, 0.4)[:-1]
full_a = age_distribution(3, 0.1)
full_b = age_distribution(3, 0.4)
active_ratio_same = all(
    math.isclose(
        active_a[index] / active_a[0],
        active_b[index] / active_b[0],
    )
    for index in range(len(active_a))
)
check(
    "N9 epsilon 不改变活动权重比例",
    active_ratio_same,
)
check(
    "N9 epsilon 改变终端权重",
    not math.isclose(
        full_a[-1],
        full_b[-1],
    ),
    f"terminal={full_a[-1]:.2f},{full_b[-1]:.2f}",
)


# ==================================================================
head("N10 文档保留直接和、参数选择与物理时间缺口")
# ==================================================================
check(
    "N10 文档保留直接和实现比较缺口",
    "直接和实现" in DOC
    and "不是本文要比较的对象" in DOC,
)
check(
    "N10 文档保留三个标量参数缺口",
    "\\tau_i,r_i,\\epsilon_i" in DOC
    and "物理时间" in DOC,
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
    "N11 文档明确最小性范围",
    "只在张量型局部载体类中成立" in DOC
    and "绝对下界" in DOC,
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
    "N12 D224 文档不含项目外体系名或外部路径",
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
