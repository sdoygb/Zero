#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D218 —— 活动相位占用与 D 层边界
================================================================================
【被检验的问题（D218）】
  周期末端的 D 边界能否让活动相位占用稳定非均匀？
  最小 T=3 模型是否给出持久非中心忠实态？
  它是否仍然依赖活动占用到量子态的识别输入？

【本步判据】
  N1  D218 与两个恢复结构已登记
  N2  周期相位与末端 D 边界分离
  N3  活动相位计数稳定为 (1,1,2)
  N4  跨周期比例持续保持
  N5  活动占用态忠实且非中心
  N6  模流出现持续非平凡相位
  N7  历史累计计数会均匀化，本文不混用两者
  N8  跨区域稳定比例一致
  N9  活动占用到量子态的识别仍是输入
  N10 不把活动占用冒充物理态或连续 QFT
  N11 不引用旧体系或外部路径
  N12 不修改上游、不新增 U5，失败时非零退出

运行：python3 verify/d218_active_phase_occupancy_and_d_boundary_state.py
"""
import io
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
            "D218_active_phase_occupancy_and_d_boundary_state.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


def run_active_phase_windows(cycles):
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


def normalized_occupancy(counts):
    vector = np.asarray(counts, dtype=float)
    return vector / float(np.sum(vector))


# ==================================================================
head("N1  D218 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 文档登记 D 边界相位态与活动占用缺口",
    "D218" in AX
    and "R-Z-U3-D-BOUNDARY-PHASE-STATE" in BOTH
    and "R-Z-U3-ACTIVE-OCCUPANCY-GAP" in BOTH
    and "活动相位占用" in DOC,
)
check(
    "N1 文档不把活动占用直接写成物理量子态",
    "不能声称" in DOC
    and "物理量子态" in DOC
    and "识别输入" in DOC,
)


# ==================================================================
head("N2  周期相位与末端 D 边界分离")
# ==================================================================
check(
    "N2 局域相位数等于 T=3 且末端相位单独吸收",
    DESTRUCTION_PERIOD == 3
    and "P_0,P_1,P_2" in DOC
    and "E_i\\longrightarrow D_i" in DOC,
)
check(
    "N2 文档明确末端不回写局部相位",
    "不并回局域相位" in DOC
    and "D_i" in DOC,
)


# ==================================================================
head("N3  活动相位计数稳定为 (1,1,2)")
# ==================================================================
windows = run_active_phase_windows(6)
site_zero_windows = [[snapshot[0] for snapshot in window] for window in windows]
check(
    "N3 首个局域相位的活动计数比例稳定",
    all(
        counts[1] == counts[0]
        and counts[2] == 2 * counts[0]
        for counts in site_zero_windows
    ),
    f"windows={site_zero_windows}",
)
check(
    "N3 文档列出同样的稳定比例",
    "(1,1,2)" in DOC.replace(" ", "")
    and "\\frac{a_{0,2}}{a_{0,1}}" in DOC,
)


# ==================================================================
head("N4  跨周期比例持续保持")
# ==================================================================
ratios = [
    tuple(
        float(value) / float(counts[0])
        for value in counts[:3]
    )
    for counts in site_zero_windows
]
check(
    "N4 每周期比例都保持 (1,1,2)",
    all(np.allclose(ratio, (1.0, 1.0, 2.0)) for ratio in ratios),
    f"ratios={ratios}",
)
check(
    "N4 文档给出多个周期样本",
    "(4,4,8)" in DOC.replace(" ", "")
    and "(340,340,680)" in DOC.replace(" ", ""),
)


# ==================================================================
head("N5  活动占用态忠实且非中心")
# ==================================================================
final_counts = site_zero_windows[-1][:3]
weights = normalized_occupancy(final_counts)
check(
    "N5 稳定权重严格正",
    np.all(weights > 0.0),
    f"weights={weights.tolist()}",
)
check(
    "N5 稳定权重非均匀",
    np.allclose(weights, np.asarray([0.25, 0.25, 0.5]))
    and not np.allclose(weights, np.full(3, 1.0 / 3.0)),
)


# ==================================================================
head("N6  模流出现持续非平凡相位")
# ==================================================================
phase = np.exp(
    1j
    * 0.7
    * (math.log(weights[0]) - math.log(weights[2]))
)
check(
    "N6 q0/q2 给非平凡模相位",
    abs(phase - 1.0) > 1.0e-6,
    f"phase={phase}",
)
check(
    "N6 文档给出 E_02 的非平凡相位",
    "2^{-it}E_{02}" in DOC
    and "持续非平凡模流候选" in DOC,
)


# ==================================================================
head("N7  历史累计计数与活动占用分开")
# ==================================================================
check(
    "N7 文档区分活动占用和闭合历史计数",
    "历史累计计数会均匀化" in DOC
    and "即时活动层" in DOC
    and "D217" in DOC,
)
check(
    "N7 活动占用计数不是闭合历史长度计数",
    all(
        len(counts) == DESTRUCTION_PERIOD + 1
        and counts[-1] == 3 * counts[0]
        for counts in site_zero_windows
    )
    and "闭合长度只能为偶数" in DOC,
)


# ==================================================================
head("N8  跨区域稳定比例一致")
# ==================================================================
last_cycle = windows[-1]
site_ratios = {}
for site in range(len(INITIAL_WORDS)):
    counts = [snapshot[site] for snapshot in last_cycle][:3]
    site_ratios[site] = tuple(
        float(value) / float(counts[0])
        for value in counts
    )
check(
    "N8 末周期所有区域的相位比例回到 (1,1,2)",
    all(np.allclose(ratio, (1.0, 1.0, 2.0)) for ratio in site_ratios.values()),
    f"site_ratios={site_ratios}",
)


# ==================================================================
head("N9  活动占用到量子态的识别仍是输入")
# ==================================================================
check(
    "N9 文档登记活动占用态输入",
    "活动路径数" in DOC
    and "相位投影" in DOC
    and "作为态权重" in DOC,
)
check(
    "N9 文档登记即时窗口与区域合并缺口",
    "即时周期活动占用" in DOC
    and "区域合并" in DOC
    and "R-Z-U3-ACTIVE-OCCUPANCY-GAP" in BOTH,
)


# ==================================================================
head("N10 不把活动占用冒充物理态或连续 QFT")
# ==================================================================
check(
    "N10 文档拒绝唯一状态映射",
    "唯一状态映射" in DOC
    and "连续 QFT" in DOC,
)
check(
    "N10 文档只登记条件构造",
    "条件构造" in DOC
    and "未解选择器" in DOC,
)


# ==================================================================
head("N11 不引用旧体系或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "旧体系",
    "外部路径",
    "../..",
)
check(
    "N11 文档不含旧体系名或外部路径",
    not any(marker in DOC for marker in external_markers),
)


# ==================================================================
head("N12 上游边界")
# ==================================================================
check(
    "N12 不修改 U1-U4+C1 且不新增 U5",
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
