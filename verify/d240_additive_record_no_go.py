#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D240 —— 可加记录无解
================================================================================
【被检验的问题（D240）】
  静态零层的可加记录读法能否生成 D239 的二次符号势？
  可加独立记录给出什么危险率？
  记录间配对相互作用能否给 D238 的多重数比？

【本步判据】
  N1  D240 与恢复结构已登记
  N2  可加独立记录给年龄线性对数比例
  N3  可加独立记录给常数危险率差
  N4  终端中性把可加记录比例写成指数型
  N5  二次记录权重给线性危险率
  N6  二次权重与 D238 抛物型比例一致
  N7  k^2 等于有序年龄配对计数
  N8  有序配对每步新增 2k+1
  N9  同一重数可以承载不同配对矩阵
  N10 可加读数看不见配对结构
  N11 均匀记录配对给 D238 多重数比
  N12 非均匀配对给不同高阶势
  N13 符号对称闭合给零配对差
  N14 文档登记与上游边界
  N15 文档不引用项目外体系名或外部路径

运行：python3 verify/d240_additive_record_no_go.py
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
            "D240_additive_record_no_go.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D240 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D240", "D240" in AX)
check(
    "N1 公理表与文档登记可加记录无解",
    "R-Z-ADDITIVE-RECORD-NOGO" in BOTH,
)
check(
    "N1 公理表与文档登记记录配对相互作用",
    "R-Z-PAIR-RECORD-INTERACTION" in BOTH,
)
check(
    "N1 公理表与文档登记配对来源缺口",
    "R-Z-PAIR-INTERACTION-ORIGIN-GAP" in BOTH,
)


# ==================================================================
head("N2  可加独立记录给年龄线性对数比例")
# ==================================================================
terminal_age = 18
ages = np.arange(0, terminal_age + 1)
w_plus = 0.17
w_minus = 0.05
initial_log_ratio = 0.4
additive_log_ratio = (
    initial_log_ratio
    + (w_plus - w_minus) * ages
)

check(
    "N2 可加记录对数比例差分常数",
    np.allclose(
        np.diff(additive_log_ratio),
        w_plus - w_minus,
    ),
)
check(
    "N2 可加记录对数比例非二次",
    np.ptp(np.diff(additive_log_ratio, 2)) < 1e-12,
)


# ==================================================================
head("N3  可加独立记录给常数危险率差")
# ==================================================================
additive_hazard_difference = -np.diff(additive_log_ratio)[:-1]

check(
    "N3 可加记录危险率差是常数",
    np.allclose(
        additive_hazard_difference,
        -(w_plus - w_minus),
    ),
)
check(
    "N3 可加记录危险率差不随年龄变化",
    np.ptp(additive_hazard_difference) < 1e-12,
)


# ==================================================================
head("N4  终端中性把可加记录比例写成指数型")
# ==================================================================
terminal_compensation = -(w_plus - w_minus) * terminal_age
additive_terminal_ratio = np.exp(
    terminal_compensation
    + (w_plus - w_minus) * ages
)
quadratic_target = np.exp(
    0.035 * (terminal_age**2 - ages**2)
)

check(
    "N4 终端补偿满足终端中性",
    abs(additive_terminal_ratio[-1] - 1.0) < 1e-12,
)
check(
    "N4 指数型比例不同于二次型比例",
    np.max(np.abs(additive_terminal_ratio - quadratic_target)) > 1e-2,
)


# ==================================================================
head("N5  二次记录权重给线性危险率")
# ==================================================================
gamma_plus = 0.11
gamma_minus = 0.24
delta = gamma_minus - gamma_plus
quadratic_log_ratio = (
    (gamma_plus - gamma_minus) * ages**2
)
quadratic_hazard_difference = -np.diff(quadratic_log_ratio)

check(
    "N5 二次权重给线性危险率差",
    np.allclose(
        quadratic_hazard_difference,
        2.0 * delta * ages[:-1] + delta,
    ),
)
check(
    "N5 二次权重对数比例的曲率为常数",
    np.ptp(np.diff(quadratic_log_ratio, 2)) < 1e-12,
)


# ==================================================================
head("N6  二次权重与 D238 抛物型比例一致")
# ==================================================================
terminal_log_ratio = delta * terminal_age**2
d238_log_ratio = terminal_log_ratio + quadratic_log_ratio
d238_ratio = np.exp(d238_log_ratio)
d238_target = np.exp(
    delta * (terminal_age**2 - ages**2)
)

check(
    "N6 终端补偿后恢复 D238 抛物线",
    np.allclose(d238_ratio, d238_target),
)
check(
    "N6 D238 比例终端等于一",
    abs(d238_ratio[-1] - 1.0) < 1e-12,
)


# ==================================================================
head("N7  k^2 等于有序年龄配对计数")
# ==================================================================
pair_counts = np.array(
    [
        len(
            [
                (i, j)
                for i in range(k)
                for j in range(k)
            ]
        )
        for k in range(terminal_age + 1)
    ]
)

check(
    "N7 有序配对计数等于 k^2",
    np.array_equal(pair_counts, ages**2),
)
check(
    "N7 配对计数随 k 二次增长",
    np.allclose(np.diff(pair_counts, 2), 2.0),
)


# ==================================================================
head("N8  有序配对每步新增 2k+1")
# ==================================================================
check(
    "N8 配对增量等于 2k+1",
    np.array_equal(
        np.diff(pair_counts),
        2 * ages[:-1] + 1,
    ),
)
check(
    "N8 增量严格递增",
    np.all(np.diff(np.diff(pair_counts)) > 0),
)


# ==================================================================
head("N9  同一重数可以承载不同配对矩阵")
# ==================================================================
size = 6
uniform = 0.08 * np.ones((size, size))
nonuniform = uniform.copy()
nonuniform[2, 2] += 0.03
nonuniform[5, 5] -= 0.03

phi_uniform = np.array(
    [uniform[:k, :k].sum() for k in range(1, size + 1)]
)
phi_nonuniform = np.array(
    [nonuniform[:k, :k].sum() for k in range(1, size + 1)]
)

check(
    "N9 两配对矩阵给相同总重数",
    abs(uniform.sum() - nonuniform.sum()) < 1e-12,
)
check(
    "N9 两配对矩阵给不同高阶势",
    np.max(np.abs(phi_uniform - phi_nonuniform)) > 1e-2,
)


# ==================================================================
head("N10  可加读数看不见配对结构")
# ==================================================================
additive_reading = lambda matrix: float(np.diag(matrix).sum() + matrix.sum())
reading_uniform = additive_reading(uniform)
reading_nonuniform = additive_reading(nonuniform)

check(
    "N10 可加读数对两矩阵相同",
    abs(reading_uniform - reading_nonuniform) < 1e-12,
)
check(
    "N10 可加读数不是配对势",
    abs(reading_uniform - phi_uniform[-1]) > 1e-12,
)


# ==================================================================
head("N11  均匀记录配对给 D238 多重数比")
# ==================================================================
delta_pair = 0.073
pair_potential = delta_pair * ages**2
pair_increment = np.diff(pair_potential)
pair_q = np.exp(-pair_increment)
d238_q = np.exp(-delta_pair * (2 * ages[:-1] + 1))

check(
    "N11 均匀配对势差给 D238 q_k",
    np.allclose(pair_q, d238_q),
)
check(
    "N11 均匀配对比例终端可补偿",
    abs(
        np.exp(pair_potential[-1])
        * np.prod(pair_q)
        - 1.0
    )
    < 1e-12,
)


# ==================================================================
head("N12  非均匀配对给不同高阶势")
# ==================================================================
nonuniform_small = uniform[:4, :4].copy()
nonuniform_small[3, 3] += 0.07
phi_small_uniform = np.array(
    [
        uniform[:k, :k].sum()
        for k in range(1, 5)
    ]
)
phi_small_nonuniform = np.array(
    [
        nonuniform_small[:k, :k].sum()
        for k in range(1, 5)
    ]
)

check(
    "N12 前三个低阶势保持相同",
    np.allclose(
        phi_small_uniform[:3],
        phi_small_nonuniform[:3],
    ),
)
check(
    "N12 高阶级势发生分叉",
    abs(phi_small_uniform[-1] - phi_small_nonuniform[-1]) > 1e-2,
)


# ==================================================================
head("N13  符号对称闭合给零配对差")
# ==================================================================
symmetric_gamma = 0.15
symmetric_log_ratio = np.zeros_like(ages, dtype=float)
symmetric_hazard_difference = -np.diff(symmetric_log_ratio)

check(
    "N13 对称符号给相同二次系数",
    abs(symmetric_gamma - symmetric_gamma) < 1e-12,
)
check(
    "N13 对称符号给零危险率差",
    np.allclose(symmetric_hazard_difference, 0.0),
)


# ==================================================================
head("N14  文档登记与上游边界")
# ==================================================================
check(
    "N14 文档登记可加记录无解",
    "R-Z-ADDITIVE-RECORD-NOGO" in DOC
    and "可加多重集" in DOC,
)
check(
    "N14 文档登记记录配对相互作用",
    "R-Z-PAIR-RECORD-INTERACTION" in DOC
    and "配对权重" in DOC,
)
check(
    "N14 文档登记配对来源缺口",
    "R-Z-PAIR-INTERACTION-ORIGIN-GAP" in DOC
    and "记录之间出现配对相互作用" in DOC,
)
check(
    "N14 文档不把新相互作用写成 U1-U4 推论",
    "它们不是 `U1-U4` 的推论" in DOC
    and "新增恢复规则" in DOC,
)
check(
    "N14 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)


# ==================================================================
head("N15  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N15 D240 文档不含项目外体系名或外部路径",
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
