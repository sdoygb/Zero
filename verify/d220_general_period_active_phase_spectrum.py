#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D220 —— 一般周期的活动相位谱
================================================================================
【被检验的问题（D220）】
  D219 的 (1,1,2) 是否只属于 T=3？
  闭环重播种后的活动相位计数是否由中心二项系数给出？
  一般 T>=3 能否给显式非中心忠实 U3 态族？

【本步判据】
  N1  D220 与一般谱结构已登记
  N2  中心二项系数满足正游走总计数
  N3  单历史在多种周期下都等于 2 c_k
  N4  多种闭合历史形状给出相同相位谱
  N5  历史集合大小只给整体倍率
  N6  不同区域给相同归一化权重
  N7  T=3 回收 D219 的 (1,1,2)
  N8  T=2 给中心态，T>=3 给非中心忠实态
  N9  一般 T 的末端相位计数为 A c_T
  N10 D215 接口的一般 T 权重与模相位正确
  N11 文献不把活动谱冒充物理态或连续 QFT
  N12 不引用旧体系或外部路径
  N13 不修改上游、不新增 U5，失败时非零退出

运行：python3 verify/d220_general_period_active_phase_spectrum.py
"""
import io
import itertools
import math
import os
import sys

import numpy as np

# [lh 迁入] 本脚本读旧理论的 AXIOMS.md 与 derivations/，故 ROOT 指向其原仓库（../modular-equilibrium）
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    os.pardir, os.pardir, "modular-equilibrium"))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from latex_utils import canonical_math
from simulations.zero_sum_periodic_destruction import (
    evolve_one_step,
    reseed_from_history,
)


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
            "D220_general_period_active_phase_spectrum.md",
        ),
        encoding="utf-8",
    ).read()
)
D215 = canonical_math(
    io.open(
        os.path.join(
            ROOT,
            "derivations",
            "D215_cycle_local_semantic_u_interface.md",
        ),
        encoding="utf-8",
    ).read()
)
D219 = canonical_math(
    io.open(
        os.path.join(
            ROOT,
            "derivations",
            "D219_active_phase_seed_independence.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


def is_closed(word):
    return bool(word) and sum(word) == 0


def closed_words_upto(max_length):
    words = []
    for length in range(2, max_length + 1, 2):
        for plus_positions in itertools.combinations(
            range(length),
            length // 2,
        ):
            word = [-1] * length
            for position in plus_positions:
                word[position] = 1
            words.append(tuple(word))
    return words


def c(k):
    return math.comb(k, k // 2)


def run_activity_spectrum(history_sets, period, cycles=1):
    histories = [set(book) for book in history_sets]
    zero_layer = set()
    cycle_spectra = []

    for _ in range(cycles):
        active = reseed_from_history(histories)
        layers = [active]
        for _step in range(period):
            active, _ = evolve_one_step(active, histories, zero_layer)
            layers.append(active)
        cycle_spectra.append(
            [
                [len(paths) for paths in layer]
                for layer in layers
            ]
        )

    return cycle_spectra


def normalized_spectrum(period):
    values = np.asarray([c(k) for k in range(period)], dtype=float)
    return values / float(np.sum(values))


closed_words = closed_words_upto(8)
probe_seeds = (
    (1, -1),
    (1, 1, -1, -1),
    (1, -1, 1, -1),
    (1, 1, 1, -1, -1, -1),
)


# ==================================================================
head("N1  D220 与一般谱结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记一般活动相位谱",
    "D220" in AX
    and "R-Z-U3-ACTIVE-PHASE-SPECTRUM" in BOTH
    and "中心二项系数" in DOC,
)
check(
    "N1 文档保留活动占用与物理局域性缺口",
    "R-Z-U3-ACTIVE-OCCUPANCY-GAP" in BOTH
    and "R-Z-U-PHYSICAL-LOCALITY-GAP" in BOTH
    and "物理量子态" in DOC,
)


# ==================================================================
head("N2  中心二项系数满足正游走总计数")
# ==================================================================
distribution = {1: 1}
recurrence_totals = []
for _step in range(15):
    recurrence_totals.append(sum(distribution.values()))
    next_distribution = {}
    for height, count in distribution.items():
        if height == 1:
            next_distribution[2] = next_distribution.get(2, 0) + count
        else:
            next_distribution[height - 1] = (
                next_distribution.get(height - 1, 0) + count
            )
            next_distribution[height + 1] = (
                next_distribution.get(height + 1, 0) + count
            )
    distribution = next_distribution
expected_totals = [c(k) for k in range(15)]
check(
    "N2 正游走递推总数等于中心二项系数",
    recurrence_totals == expected_totals,
    f"sequence={recurrence_totals[:10]}",
)
check(
    "N2 文档给出偶奇中心二项形式",
    "c_{2m}" in DOC
    and "\\binom{2m}{m}" in DOC
    and "c_{2m+1}" in DOC
    and "\\binom{2m+1}{m}" in DOC,
)


# ==================================================================
head("N3  单历史在多种周期下都等于 2 c_k")
# ==================================================================
periods = tuple(range(1, 11))
single_ok = True
single_examples = {}
for period in periods:
    for seed in probe_seeds:
        spectrum = run_activity_spectrum([{seed}], period, cycles=1)[0]
        counts = [layer[0] for layer in spectrum]
        expected = [2 * c(k) for k in range(period + 1)]
        if counts != expected:
            single_ok = False
            break
    single_examples[period] = expected
check(
    "N3 每种单历史都给出 A c_k，A=2",
    single_ok,
    f"periods={periods}",
)
check(
    "N3 文档列出一般相位谱前项",
    "1,1,2,3,6,10,20,35,\\ldots" in DOC,
)


# ==================================================================
head("N4  多种闭合历史形状给出相同相位谱")
# ==================================================================
shape_runs = [
    run_activity_spectrum([{seed}], 8, cycles=1)[0]
    for seed in closed_words[:20]
]
shape_spectra = [
    tuple(layer[0] for layer in spectrum)
    for spectrum in shape_runs
]
check(
    "N4 前二十个闭合词的相位谱完全一致",
    len(set(shape_spectra)) == 1
    and shape_spectra[0] == tuple(2 * c(k) for k in range(9)),
    f"distinct_spectra={len(set(shape_spectra))}",
)
check(
    "N4 文档说明词形作商只留绝对值",
    "只留下 }n=|\\sum x_j|" in DOC
    or "只留下 $n=|\\sum x_j|$" in DOC,
)


# ==================================================================
head("N5  历史集合大小只给整体倍率")
# ==================================================================
history_counts = (1, 2, 3, 5, 10, 20)
scaling_runs = [
    run_activity_spectrum([set(closed_words[:count])], 6, cycles=1)[0]
    for count in history_counts
]
scaling_ok = True
for run, history_count in zip(scaling_runs, history_counts):
    counts = [layer[0] for layer in run]
    expected = [2 * history_count * c(k) for k in range(7)]
    if counts != expected:
        scaling_ok = False
        break
check(
    "N5 相位谱按 2|H| 精确缩放",
    scaling_ok,
    f"history_counts={history_counts}",
)


# ==================================================================
head("N6  不同区域给相同归一化权重")
# ==================================================================
region_history_sets = [
    set(closed_words[:1]),
    set(closed_words[1:4]),
    set(closed_words[4:9]),
]
region_run = run_activity_spectrum(region_history_sets, 7, cycles=1)[0]
region_weights = []
for site in range(len(region_history_sets)):
    counts = np.asarray(
        [layer[site] for layer in region_run[:7]],
        dtype=float,
    )
    region_weights.append(counts / float(np.sum(counts)))
check(
    "N6 不同区域的历史大小不同但归一化谱相同",
    all(np.allclose(weight, region_weights[0]) for weight in region_weights),
    f"weights={[weight.tolist() for weight in region_weights]}",
)


# ==================================================================
head("N7  T=3 回收 D219 的 (1,1,2)")
# ==================================================================
t3_spectrum = run_activity_spectrum([{probe_seeds[0]}], 3, cycles=1)[0]
t3_counts = [layer[0] for layer in t3_spectrum]
check(
    "N7 T=3 的相位计数为 (2,2,4,6)",
    t3_counts == [2, 2, 4, 6],
    f"t3_counts={t3_counts}",
)
check(
    "N7 D219 与 D220 在文档中相互衔接",
    "D219" in DOC
    and "D220" in D219
    and "T=3" in DOC,
)


# ==================================================================
head("N8  T=2 给中心态，T>=3 给非中心忠实态")
# ==================================================================
t2_weights = normalized_spectrum(2)
t3_weights = normalized_spectrum(3)
t4_weights = normalized_spectrum(4)
check(
    "N8 T=2 权重均匀且态中心",
    np.allclose(t2_weights, (0.5, 0.5))
    and np.allclose(t2_weights, np.full(2, 0.5)),
)
check(
    "N8 T=3 与 T=4 权重严格正且非均匀",
    np.all(t3_weights > 0.0)
    and np.all(t4_weights > 0.0)
    and not np.allclose(t3_weights, np.full(3, 1.0 / 3.0))
    and not np.allclose(t4_weights, np.full(4, 0.25)),
    f"T3={t3_weights.tolist()}, T4={t4_weights.tolist()}",
)
check(
    "N8 文档给出 T=3,4,5 权重样本",
    "(1/4,1/4,1/2)" in DOC.replace(" ", "")
    and "(1/7,1/7,2/7,3/7)" in DOC.replace(" ", "")
    and "(1/13,1/13,2/13,3/13,6/13)" in DOC.replace(" ", ""),
)


# ==================================================================
head("N9  一般 T 的末端相位计数为 A c_T")
# ==================================================================
terminal_ok = True
for period in periods:
    spectrum = run_activity_spectrum([{probe_seeds[1]}], period, cycles=1)[0]
    terminal_count = spectrum[-1][0]
    if terminal_count != 2 * c(period):
        terminal_ok = False
        break
check(
    "N9 末端计数精确等于 A c_T",
    terminal_ok,
    f"periods={periods}",
)
check(
    "N9 文档区分末端相位与局部相位",
    "a_{i,T}=A_i c_T" in DOC
    and "不是局部相位投影" in DOC,
)


# ==================================================================
head("N10 D215 接口的一般 T 权重与模相位正确")
# ==================================================================
modular_ok = True
for period in (3, 4, 5):
    weights = normalized_spectrum(period)
    modular_phase = np.exp(
        1j * 0.7 * (math.log(weights[0]) - math.log(weights[period - 1]))
    )
    if abs(modular_phase - 1.0) <= 1.0e-6:
        modular_ok = False
check(
    "N10 T>=3 的首末模相位非平凡",
    modular_ok,
)
check(
    "N10 D215 接口已更新为一般 T 态族",
    "D219" in D215
    and "D220" in D215
    and "q_{v,k}" in D215
    and "\\rho_{\\rm cyc}^{(T)}" in D215,
)


# ==================================================================
head("N11 文献不把活动谱冒充物理态或连续 QFT")
# ==================================================================
check(
    "N11 文档保留物理态映射缺口",
    "活动占用到物理量子态" in DOC
    and "连续 QFT" in DOC
    and "物理时间" in DOC,
)
check(
    "N11 文档只声称条件结构定理",
    "条件结构定理" in DOC
    and "不是物理量子引力构造" in DOC,
)


# ==================================================================
head("N12 不引用旧体系或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "旧体系",
    "外部路径",
    "../..",
)
check(
    "N12 D220 文档不含旧体系名或外部路径",
    not any(marker in DOC for marker in external_markers),
)


# ==================================================================
head("N13 上游边界")
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
