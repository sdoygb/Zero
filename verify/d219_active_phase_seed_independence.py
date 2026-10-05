#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D219 —— 活动相位结构的种子无关性
================================================================================
【被检验的问题（D219）】
  D218 的 (1,1,2) 是否依赖特定初始词或历史集合？
  任意非空闭合历史集合在 T=3 周期末重播种后是否都给 A(1,1,2)？
  该结构比例是否仍不消除活动占用到量子态的识别输入？

【本步判据】
  N1  D219 与结构比例已登记
  N2  D219 只研究局部相位比例，不混同总增长率
  N3  闭合历史都是偶数长度，重播种后词和为 ±1
  N4  穷举的最小闭合词都给出 (2,2,4,6)
  N5  任意多个闭合历史只给整体等比例放大
  N6  不同区域可独立放大，但局部比例仍为 (1,1,2)
  N7  相位 0 到 1 恰好一半闭合、一半存活
  N8  相位 1 到 2 两个延拓都不闭合
  N9  末端相位 3 计数为 3A 且不并回局部相位
  N10 D218 当前模型的第二周期后满足结构引理
  N11 归一化态为 (1/4,1/4,1/2)，忠实且非中心
  N12 文献不把结构引理冒充物理态或连续 QFT
  N13 不引用旧体系或外部路径
  N14 不修改上游、不新增 U5，失败时非零退出

运行：python3 verify/d219_active_phase_seed_independence.py
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
    CYCLES,
    DESTRUCTION_PERIOD,
    INITIAL_WORDS,
    evolve_one_step,
    reseed_from_history,
    word_from_string,
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
            "D219_active_phase_seed_independence.md",
        ),
        encoding="utf-8",
    ).read()
)
D218 = canonical_math(
    io.open(
        os.path.join(
            ROOT,
            "derivations",
            "D218_active_phase_occupancy_and_d_boundary_state.md",
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


def run_cycles(history_sets, cycles=1):
    histories = [set(book) for book in history_sets]
    zero_layer = set()
    cycle_layers = []

    for _ in range(cycles):
        active = reseed_from_history(histories)
        layers = [active]
        for _step in range(DESTRUCTION_PERIOD):
            active, _ = evolve_one_step(active, histories, zero_layer)
            layers.append(active)
        cycle_layers.append(layers)
        active = [set() for _ in histories]

    return cycle_layers


def phase_counts(layers):
    return [
        [len(paths) for paths in layer]
        for layer in layers
    ]


def run_d218_windows(cycles):
    active = [{word_from_string(text)} for text in INITIAL_WORDS]
    histories = [set() for _ in INITIAL_WORDS]
    zero_layer = set()
    windows = []

    for _ in range(cycles):
        phases = [[len(paths) for paths in active]]
        for _step in range(DESTRUCTION_PERIOD):
            active, _ = evolve_one_step(active, histories, zero_layer)
            phases.append([len(paths) for paths in active])
        windows.append(phases)
        active = [set() for _ in active]
        active = reseed_from_history(histories)

    return windows


closed_words = closed_words_upto(8)


# ==================================================================
head("N1  D219 与结构比例已登记")
# ==================================================================
check(
    "N1 文档与公理表登记活动相位结构比例",
    "D219" in AX
    and "R-Z-U3-ACTIVE-PHASE-RATIO" in BOTH
    and "R-Z-U3-ACTIVE-PHASE-INTERFACE" in BOTH
    and "种子无关性" in DOC,
)
check(
    "N1 文档保留活动占用到量子态的缺口",
    "R-Z-U3-ACTIVE-OCCUPANCY-GAP" in BOTH
    and "物理量子态" in DOC
    and "识别输入" in DOC,
)


# ==================================================================
head("N2  D219 只研究局部相位比例，不混同总增长率")
# ==================================================================
check(
    "N2 文档明确分开相位比例与周期总增长",
    "局部相位比例" in DOC
    and "周期总增长速率" in DOC
    and "总规模由历史层中新闭合词的数量决定" in DOC,
)
check(
    "N2 D218 已把结构定理与初始样本分开",
    "模型读数而不是种子无关性定理" in D218
    and "D219" in D218
    and "(a_0,a_1,a_2)=A(1,1,2)" in D218.replace(" ", ""),
)


# ==================================================================
head("N3  闭合历史都是偶数长度，重播种后词和为 ±1")
# ==================================================================
check(
    "N3 枚举闭合词长度全为偶数",
    len(closed_words) == 98
    and all(is_closed(word) and len(word) % 2 == 0 for word in closed_words),
    f"closed_words={len(closed_words)}",
)
reseeded_signs = []
for word in closed_words:
    reseeded_signs.extend((sum(word + (1,)), sum(word + (-1,))))
check(
    "N3 闭合历史的两个初始延拓词和都为 ±1",
    all(abs(value) == 1 for value in reseeded_signs),
)
check(
    "N3 文档证明偶数长度与 ±1 重播种",
    "|w|=r+s=2r" in DOC
    and "\\sum(w_+)=+1" in DOC
    and "\\sum(w_-)=-1" in DOC,
)


# ==================================================================
head("N4  穷举的最小闭合词都给出 (2,2,4,6)")
# ==================================================================
single_seed_runs = [
    run_cycles([{word}], cycles=2)
    for word in closed_words
]
single_seed_ok = True
single_seed_first_cycle_ok = True
for runs in single_seed_runs:
    first_counts = [site[0] for site in phase_counts(runs[0])[:4]]
    if first_counts != [2, 2, 4, 6]:
        single_seed_first_cycle_ok = False
    for cycle in runs:
        counts = phase_counts(cycle)[:4]
        site_counts = [site[0] for site in counts[:3]]
        if not (
            site_counts[1] == site_counts[0]
            and site_counts[2] == 2 * site_counts[0]
        ):
            single_seed_ok = False
            break
check(
    "N4 所有平衡词的首个重播种周期都给出 (2,2,4,6)",
    single_seed_first_cycle_ok,
    f"tested={len(single_seed_runs)}",
)
check(
    "N4 所有平衡词后续周期都保持局部 (1,1,2)",
    single_seed_ok,
    f"tested={len(single_seed_runs)}",
)
check(
    "N4 文档给出单历史基准 (2,2,4,6)",
    "(2,2,4,6)" in DOC.replace(" ", ""),
)


# ==================================================================
head("N5  任意多个闭合历史只给整体等比例放大")
# ==================================================================
history_counts = (1, 2, 3, 5, 10, 20)
history_sets = [set(closed_words[:count]) for count in history_counts]
multi_runs = [
    run_cycles([history_set], cycles=1)[0]
    for history_set in history_sets
]
multi_ok = True
for index, count in enumerate(history_counts):
    counts = [
        layer[0]
        for layer in phase_counts(multi_runs[index])[:4]
    ]
    if counts != [2 * count, 2 * count, 4 * count, 6 * count]:
        multi_ok = False
        break
check(
    "N5 历史集合大小只决定共同倍率",
    multi_ok,
    f"history_counts={history_counts}",
)


# ==================================================================
head("N6  不同区域可独立放大，但局部比例仍为 (1,1,2)")
# ==================================================================
region_history_sets = [
    set(closed_words[:1]),
    set(closed_words[1:3]),
    set(closed_words[3:7]),
]
region_run = run_cycles(region_history_sets, cycles=1)[0]
region_phase_counts = phase_counts(region_run)
region_ratios = []
for site, count in enumerate((1, 2, 4)):
    site_counts = [layer[site] for layer in region_phase_counts[:3]]
    expected = [2 * count, 2 * count, 4 * count]
    region_ratios.append(
        tuple(float(value) / site_counts[0] for value in site_counts)
    )
    if site_counts != expected:
        region_ratios[-1] = None
check(
    "N6 每区的相位 0,1,2 比例都为 (1,1,2)",
    all(
        ratio is not None
        and np.allclose(ratio, (1.0, 1.0, 2.0))
        for ratio in region_ratios
    ),
    f"region_ratios={region_ratios}",
)


# ==================================================================
head("N7  相位 0 到 1 恰好一半闭合、一半存活")
# ==================================================================
probe_run = run_cycles([{closed_words[3]}], cycles=1)[0]
phase0_signs = [sum(word) for word in probe_run[0][0]]
phase1_signs = [sum(word) for word in probe_run[1][0]]
phase2_signs = [sum(word) for word in probe_run[2][0]]
check(
    "N7 相位 0 和相位 1 的词和层结构正确",
    all(abs(value) == 1 for value in phase0_signs)
    and all(abs(value) == 2 for value in phase1_signs),
    f"phase0={phase0_signs}, phase1={phase1_signs}",
)
check(
    "N7 相位 2 的词和只取 ±1 或 ±3",
    all(abs(value) in (1, 3) for value in phase2_signs),
    f"phase2={phase2_signs}",
)
check(
    "N7 文档给出 0→1 存活保持的证明",
    "相位 }0\\text{ 到相位 }1\\text{ 的存活计数保持不变" in DOC
    and "恰好产生" in DOC,
)


# ==================================================================
head("N8  相位 1 到 2 两个延拓都不闭合")
# ==================================================================
check(
    "N8 相位 2 计数精确为相位 1 的两倍",
    len(probe_run[2][0]) == 2 * len(probe_run[1][0]),
    f"phase1={len(probe_run[1][0])}, phase2={len(probe_run[2][0])}",
)
check(
    "N8 文档用 ±2 词和排除闭合",
    "\\sum_jy_j=\\pm2" in DOC
    and "都不为零，因此都不闭合" in DOC,
)


# ==================================================================
head("N9  末端相位 3 计数为 3A 且不并回局部相位")
# ==================================================================
probe_counts = [len(layer[0]) for layer in probe_run]
check(
    "N9 单历史的末端计数为 3A",
    probe_counts[:4] == [2, 2, 4, 6],
    f"probe_counts={probe_counts}",
)
check(
    "N9 文档拒绝把相位 3 并回局部相位",
    "a_{i,3}=3a_{i,1}=3A_i" in DOC
    and "不并入" in DOC
    and "第四个局部相位" in DOC,
)


# ==================================================================
head("N10 D218 当前模型的第二周期后满足结构引理")
# ==================================================================
d218_windows = run_d218_windows(6)
d218_late_windows = d218_windows[1:]
d218_ok = True
for window in d218_late_windows:
    for site_counts in range(len(INITIAL_WORDS)):
        counts = [snapshot[site_counts] for snapshot in window[:3]]
        if counts[1] != counts[0] or counts[2] != 2 * counts[0]:
            d218_ok = False
            break
check(
    "N10 D218 的第二周期起所有区域都保持 (1,1,2)",
    d218_ok,
    f"cycles={len(d218_windows)}",
)
check(
    "N10 文档说明 D218 数值表与 D219 定理的区别",
    "D219" in D218
    and "种子无关性定理" in D218,
)


# ==================================================================
head("N11 归一化态为 (1/4,1/4,1/2)，忠实且非中心")
# ==================================================================
weights = np.asarray([1.0, 1.0, 2.0])
weights = weights / float(np.sum(weights))
phase = np.exp(
    1j
    * 0.7
    * (math.log(weights[0]) - math.log(weights[2]))
)
check(
    "N11 权重严格正且非均匀",
    np.allclose(weights, (0.25, 0.25, 0.5))
    and np.all(weights > 0.0)
    and not np.allclose(weights, np.full(3, 1.0 / 3.0)),
)
check(
    "N11 模相位非平凡",
    abs(phase - 1.0) > 1.0e-6,
    f"phase={phase}",
)
check(
    "N11 文档给出同一归一化态",
    "\\rho_{\\rm act}" in DOC
    and "\\frac14P_0+\\frac14P_1+\\frac12P_2" in DOC,
)


# ==================================================================
head("N12 文献不把结构引理冒充物理态或连续 QFT")
# ==================================================================
check(
    "N12 文档保留物理态映射缺口",
    "活动占用到量子态" in DOC
    and "物理量子态" in DOC
    and "连续 QFT" in DOC,
)
check(
    "N12 文档只声称条件结构引理",
    "条件结构引理" in DOC
    and "不等于完成" in DOC
    and "没有被消除" in DOC,
)


# ==================================================================
head("N13 不引用旧体系或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "旧体系",
    "外部路径",
    "../..",
)
check(
    "N13 D219 文档不含旧体系名或外部路径",
    not any(marker in DOC for marker in external_markers),
)
check(
    "N13 D218 的 D219 交叉引用仍自包含",
    "D219" in D218
    and not any(marker in D218 for marker in external_markers),
)


# ==================================================================
head("N14 上游边界")
# ==================================================================
check(
    "N14 不修改 U1-U4+C1 且不新增 U5",
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
