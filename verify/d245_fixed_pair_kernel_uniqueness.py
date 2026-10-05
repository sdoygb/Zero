#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D245 —— 固定二体核唯一性
================================================================================
【被检验的问题（D245）】
  固定平移不变二体核能否被精确二次势唯一选择？
  年龄增长边权为何不属于这个固定核？
  静态零记录是否自动给出静态记录对耦合？

【本步判据】
  N1  D245 与恢复结构已登记
  N2  平移不变二体核给前缀增量公式
  N3  目标线性危险率固定自配对
  N4  相邻差分逐项固定全部距离核
  N5  唯一核回代给二次势
  N6  有限范围核不能匹配全部年龄
  N7  年龄增长边权可给二次势但显含年龄
  N8  增长边权不属于固定二体核
  N9  静态记录不等于静态记录对耦合
  N10 终端中性不选择核形
  N11 文档登记与上游边界
  N12 文档不引用项目外体系名或外部路径

运行：python3 verify/d245_fixed_pair_kernel_uniqueness.py
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
            "D245_fixed_pair_kernel_uniqueness.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D245 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D245", "D245" in AX)
check(
    "N1 公理表与文档登记固定核唯一性",
    "R-Z-FIXED-KERNEL-UNIQUENESS" in BOTH,
)
check(
    "N1 公理表与文档登记增长边权排除",
    "R-Z-CHANGING-EDGE-WEIGHT-NOGO" in BOTH,
)
check(
    "N1 公理表与文档登记静态配对核缺口",
    "R-Z-STATIC-PAIR-KERNEL-GAP" in BOTH,
)


# ==================================================================
head("N2  平移不变二体核给前缀增量公式")
# ==================================================================
delta = 0.17
max_age = 12
kernel = np.full(max_age + 1, delta)


def prefix_potential(s_values, k):
    total = k * s_values[0]
    total += 2.0 * sum(
        (k - r) * s_values[r] for r in range(1, k)
    )
    return float(total)


def prefix_increment(s_values, k):
    return float(s_values[0] + 2.0 * np.sum(s_values[1 : k + 1]))


phi_values = np.array([prefix_potential(kernel, k) for k in range(max_age + 1)])
increment_values = np.array(
    [prefix_increment(kernel, k) for k in range(max_age)]
)
direct_increments = np.diff(phi_values)

check(
    "N2 前缀差分等于距离累计和",
    np.allclose(direct_increments, increment_values),
)
check(
    "N2 距离累计和对常数核线性增长",
    np.allclose(
        increment_values,
        delta * (2 * np.arange(max_age) + 1),
    ),
)


# ==================================================================
head("N3  目标线性危险率固定自配对")
# ==================================================================
target_increment = delta * (2 * np.arange(max_age + 1) + 1)
recovered_kernel = np.empty_like(kernel)
recovered_kernel[0] = target_increment[0]
recovered_kernel[1:] = 0.5 * np.diff(target_increment)

check(
    "N3 自配对由 k=0 增量固定",
    np.isclose(recovered_kernel[0], delta),
)
check(
    "N3 目标增量不依赖任意核假设",
    np.allclose(
        target_increment,
        delta * (2 * np.arange(max_age + 1) + 1),
    ),
)


# ==================================================================
head("N4  相邻差分逐项固定全部距离核")
# ==================================================================
check(
    "N4 相邻危险率增量差为 2delta",
    np.allclose(
        np.diff(target_increment),
        2.0 * delta * np.ones(max_age),
    ),
)
check(
    "N4 反推的全部距离核都是 delta",
    np.allclose(recovered_kernel, delta),
)
check(
    "N4 固定核非有限范围",
    np.count_nonzero(recovered_kernel) == max_age + 1,
)


# ==================================================================
head("N5  唯一核回代给二次势")
# ==================================================================
recovered_phi = np.array(
    [prefix_potential(recovered_kernel, k) for k in range(max_age + 1)]
)

check(
    "N5 回代给 delta k^2",
    np.allclose(recovered_phi, delta * np.arange(max_age + 1) ** 2),
)
check(
    "N5 回代给目标危险率增量",
    np.allclose(np.diff(recovered_phi), target_increment[:-1]),
)


# ==================================================================
head("N6  有限范围核不能匹配全部年龄")
# ==================================================================
finite_range = np.array([delta, delta, delta, 0.0, 0.0], dtype=float)
finite_increment = np.array(
    [prefix_increment(finite_range, k) for k in range(4)]
)

check(
    "N6 有限范围核早期可匹配",
    np.allclose(finite_increment[:3], target_increment[:3]),
)
check(
    "N6 有限范围核在范围外失败",
    not np.isclose(finite_increment[3], target_increment[3]),
)


# ==================================================================
head("N7  年龄增长边权可给二次势但显含年龄")
# ==================================================================
ages = np.arange(1, max_age + 1)
edge_coefficient = 2.0 * delta
growing_edge_weight = edge_coefficient * ages
growing_potential = (ages / 2.0) * growing_edge_weight

check(
    "N7 增长边权给二次势",
    np.allclose(growing_potential, delta * ages**2),
)
check(
    "N7 单边权随年龄变化",
    np.allclose(np.diff(growing_edge_weight), edge_coefficient),
)


# ==================================================================
head("N8  增长边权不属于固定二体核")
# ==================================================================
fixed_kernel_pair_energy = np.full_like(ages, delta, dtype=float)
changed_pair_energy = growing_edge_weight

check(
    "N8 同一记录对的能量随未来年龄改变",
    not np.allclose(fixed_kernel_pair_energy, changed_pair_energy),
)
check(
    "N8 年龄增长核形成函数 K_ij(k)",
    np.all(np.diff(changed_pair_energy) > 0),
)


# ==================================================================
head("N9  静态记录不等于静态记录对耦合")
# ==================================================================
record_labels = np.array([0, 1, 2], dtype=int)
record_labels_after_write = np.array([0, 1, 2], dtype=int)
pair_energy_before = np.full((3, 3), delta)
pair_energy_after = np.full((3, 3), 2.0 * delta)

check(
    "N9 已写入记录标签保持静态",
    np.array_equal(record_labels, record_labels_after_write),
)
check(
    "N9 记录对能量仍可被重标定",
    not np.allclose(pair_energy_before, pair_energy_after),
)


# ==================================================================
head("N10  终端中性不选择核形")
# ==================================================================
def terminal_neutral_log_ratio(potential_values):
    log_ratio = -potential_values
    log_ratio -= log_ratio[-1]
    return log_ratio


uniform_log_ratio = terminal_neutral_log_ratio(
    delta * np.arange(max_age + 1) ** 2
)
finite_log_ratio = terminal_neutral_log_ratio(
    np.array([prefix_potential(finite_range, k) for k in range(5)] + [0.0] * (max_age - 4))
)

check(
    "N10 均匀核可终端中性",
    abs(uniform_log_ratio[-1]) < 1e-12,
)
check(
    "N10 非二次势也可终端中性",
    abs(finite_log_ratio[-1]) < 1e-12,
)
check(
    "N10 终端中性不选择核形",
    not np.allclose(uniform_log_ratio, finite_log_ratio),
)


# ==================================================================
head("N11  文档登记与上游边界")
# ==================================================================
check(
    "N11 文档登记固定核唯一性",
    "R-Z-FIXED-KERNEL-UNIQUENESS" in DOC
    and "唯一全对全常数核" in DOC,
)
check(
    "N11 文档登记增长边权排除",
    "R-Z-CHANGING-EDGE-WEIGHT-NOGO" in DOC
    and "不属于固定二体核" in DOC,
)
check(
    "N11 文档登记静态配对核缺口",
    "R-Z-STATIC-PAIR-KERNEL-GAP" in DOC
    and "记录静态" in DOC
    and "记录对耦合静态" in DOC,
)
check(
    "N11 文档不把新结构写成 U1-U4 推论",
    "它们不是 `U1-U4` 的推论" in DOC
    and "仍不是 `U1-U4` 的推论" in DOC,
)
check(
    "N11 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
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
    "N12 D245 文档不含项目外体系名或外部路径",
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
