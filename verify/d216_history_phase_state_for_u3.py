#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D216 —— 精确历史相位计数与 U3 非中心态
================================================================================
【被检验的问题（D216）】
  全局零层为什么会自然给出均匀中心态？
  精确局部历史能否给出非中心忠实态并启动非平凡模流？
  历史相位态是否仍然依赖窗口、相位解释与区域选择？

【本步判据】
  N1  D216 与两个恢复结构已登记
  N2  循环不变的全局零层只能给均匀中心态
  N3  局部精确历史保留相位锚点
  N4  最小周期模型的历史窗口可计算相位计数
  N5  四个周期后历史相位计数全正
  N6  历史相位计数非均匀
  N7  历史相位态忠实且非中心
  N8  模流出现非平凡相位
  N9  历史窗口依赖仍被登记
  N10 相位解释与均匀加权仍是输入
  N11 T 未被等同于模参数
  N12 不引用旧体系或外部路径
  N13 不修改上游、不新增 U5，失败时非零退出

运行：python3 verify/d216_history_phase_state_for_u3.py
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
            "D216_history_phase_state_for_u3.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


def history_phase_counts(histories, site, period):
    counts = [0] * period
    for word in histories[site]:
        counts[len(word) % period] += 1
    return counts


def run_periodic_histories():
    active = [{word_from_string(text)} for text in INITIAL_WORDS]
    histories = [set() for _ in INITIAL_WORDS]
    zero_layer = set()
    snapshots = []

    for time in range(1, DESTRUCTION_PERIOD * CYCLES + 1):
        active, _ = evolve_one_step(active, histories, zero_layer)
        if time % DESTRUCTION_PERIOD == 0:
            snapshots.append(
                [set(book) for book in histories]
            )
            active = [set() for _ in active]
            active = reseed_from_history(histories)
    return snapshots


def modular_phase(weights, source, target, time=0.7):
    return np.exp(1j * time * (math.log(weights[source]) - math.log(weights[target])))


# ==================================================================
head("N1  D216 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 文档登记历史相位态与窗口缺口",
    "D216" in AX
    and "R-Z-U3-HISTORY-PHASE-STATE" in BOTH
    and "R-Z-U3-HISTORY-WINDOW-GAP" in BOTH
    and "历史相位态" in DOC,
)
check(
    "N1 文档不把历史计数写成唯一态",
    "不能声称" in DOC
    and "唯一非中心态" in DOC
    and "条件选择器" in DOC,
)


# ==================================================================
head("N2  循环不变的全局零层只能给均匀中心态")
# ==================================================================
period = DESTRUCTION_PERIOD
uniform_weights = np.full(period, 1.0 / period)
shifted_weights = np.roll(uniform_weights, 1)
check(
    "N2 均匀权重在循环平移下不变",
    np.allclose(uniform_weights, shifted_weights),
)
check(
    "N2 均匀态是中心态并给恒等模流",
    np.allclose(np.diag(uniform_weights), np.eye(period) / period)
    and "\\sigma_t^{\\rho_{\\rm unif}}" in DOC
    and "\\operatorname{id}" in DOC,
)


# ==================================================================
head("N3  局部精确历史保留相位锚点")
# ==================================================================
sample_words = (
    (1, -1, 1, -1),
    (-1, 1, -1, 1),
)
canonical_forms = {
    min(word[index:] + word[:index] for index in range(len(word)))
    for word in sample_words
}
check(
    "N3 同一循环类的两个代表可在精确历史层中不同",
    len(sample_words) == 2
    and len(canonical_forms) == 1
    and sample_words[0] != sample_words[1],
)
check(
    "N3 文档区分全局零层与精确局部历史",
    "P_i" in DOC
    and "循环等价类" in DOC
    and "精确词" in DOC,
)


# ==================================================================
head("N4  最小周期模型的历史窗口可计算相位计数")
# ==================================================================
snapshots = run_periodic_histories()
check(
    "N4 周期模型产生每周期历史快照",
    len(snapshots) == CYCLES,
    f"snapshots={len(snapshots)}",
)
final_histories = snapshots[-1]
counts_site_0 = history_phase_counts(final_histories, 0, period)
check(
    "N4 末窗口可按长度模 T 计数",
    sum(counts_site_0) == len(final_histories[0])
    and sum(counts_site_0) > 0,
    f"counts={counts_site_0}, total={sum(counts_site_0)}",
)


# ==================================================================
head("N5  四个周期后历史相位计数全正")
# ==================================================================
check(
    "N5 末窗口全部相位都有正计数",
    all(count > 0 for count in counts_site_0),
    f"counts={counts_site_0}",
)
check(
    "N5 文档给出忠实条件",
    "N_{i,k}(m)>0" in DOC
    and "忠实" in DOC,
)


# ==================================================================
head("N6  历史相位计数非均匀")
# ==================================================================
expected_counts = [60, 55, 55]
check(
    "N6 最小周期见证的历史计数与文档一致",
    counts_site_0 == expected_counts,
    f"computed={counts_site_0}, expected={expected_counts}",
)
check(
    "N6 历史计数非均匀",
    len(set(counts_site_0)) > 1
    and "60,55,55" in DOC.replace(" ", ""),
)


# ==================================================================
head("N7  历史相位态忠实且非中心")
# ==================================================================
weights = np.asarray(counts_site_0, dtype=float) / float(sum(counts_site_0))
check(
    "N7 权重全正因而态忠实",
    np.all(weights > 0.0),
    f"weights={weights.tolist()}",
)
check(
    "N7 权重不全相等因而态非中心",
    not np.allclose(weights, np.full(period, 1.0 / period))
    and np.linalg.norm(weights - np.full(period, 1.0 / period)) > 1.0e-3,
    f"L1-distance={np.linalg.norm(weights - np.full(period, 1.0 / period), ord=1):.6f}",
)


# ==================================================================
head("N8  模流出现非平凡相位")
# ==================================================================
phase = modular_phase(weights, 0, 1)
check(
    "N8 q0/q1 给非平凡模相位",
    not np.isclose(phase, 1.0)
    and abs(phase - 1.0) > 1.0e-6,
    f"phase={phase}",
)
check(
    "N8 文档给出矩阵单位相位公式",
    "\\sigma_t^{\\rho_i}(E_{jk})" in DOC
    and "60}{55" in DOC,
)


# ==================================================================
head("N9  历史窗口依赖仍被登记")
# ==================================================================
forbidden = all(count > 0 for count in counts_site_0)
earlier_counts = [
    history_phase_counts(snapshot, 0, period)
    for snapshot in snapshots
]
check(
    "N9 历史窗口不同计数不同",
    len({tuple(counts) for counts in earlier_counts}) > 1,
    f"windows={earlier_counts}",
)
check(
    "N9 文档登记窗口截断缺口",
    "历史窗口" in DOC
    and "R-Z-U3-HISTORY-WINDOW-GAP" in BOTH
    and forbidden,
)


# ==================================================================
head("N10 相位解释与均匀加权仍是输入")
# ==================================================================
check(
    "N10 文档登记长度相位解释",
    "|w|\\equiv k\\pmod T" in DOC
    and "长度模" in DOC,
)
check(
    "N10 文档登记精确历史均匀加权",
    "不偏好任一精确历史" in DOC
    and "均匀计数" in DOC,
)


# ==================================================================
head("N11 T 未被等同于模参数")
# ==================================================================
check(
    "N11 文档分离 T 与模时间",
    "t_{\\rm mod}=T" in DOC
    and "模时间仍来自态的模谱" in DOC,
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
    "N12 文档不含旧体系名或外部路径",
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
