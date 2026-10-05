#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D217 —— 历史相位态的渐近障碍
================================================================================
【被检验的问题（D217）】
  闭合历史计数态的非中心性能否长期保留？
  偶数周期是否因奇偶屏障而不忠实？
  奇数周期的最小 T=3 模型是否渐近回到均匀中心态？

【本步判据】
  N1  D217 与渐近障碍结构已登记
  N2  闭合词长度必为偶数
  N3  偶数 T 只占据半数相位
  N4  奇数 T 的偶数位移生成全部相位
  N5  T=3 递推矩阵的显式形式
  N6  T=3 非均匀模式衰减率严格小于 1
  N7  数值窗口偏离均匀量下降
  N8  渐近均匀态是中心态且模流恒等
  N9  U3 字面接口与持久非平凡模时间被分开
  N10 持久非中心性需要额外相位破缺源
  N11 不引用旧体系或外部路径
  N12 不修改上游、不新增 U5，失败时非零退出

运行：python3 verify/d217_history_phase_asymptotic_no_go.py
"""
import io
import itertools
import math
import os
import sys

import numpy as np

from latex_utils import canonical_math


PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))


def head(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# [lh 迁入] 本脚本读旧理论的 AXIOMS.md 与 derivations/，故 ROOT 指向其原仓库（../modular-equilibrium）
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    os.pardir, os.pardir, "modular-equilibrium"))
AX = io.open(os.path.join(ROOT, "AXIOMS.md"), encoding="utf-8").read()
DOC = canonical_math(
    io.open(
        os.path.join(
            ROOT,
            "derivations",
            "D217_history_phase_asymptotic_no_go.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


def is_closed_word(word):
    return bool(word) and sum(word) == 0


def phase_deviation(counts):
    vector = np.asarray(counts, dtype=float)
    probabilities = vector / float(np.sum(vector))
    uniform = np.full(len(vector), 1.0 / len(vector))
    return float(np.linalg.norm(probabilities - uniform, ord=1))


def recurrence_eigenvalues(period):
    shift = np.roll(np.eye(period), 1, axis=0)
    matrix = np.eye(period) + 4.0 * np.linalg.matrix_power(shift, 2)
    return np.linalg.eigvals(matrix)


# ==================================================================
head("N1  D217 与渐近障碍结构已登记")
# ==================================================================
check(
    "N1 文档登记历史相位渐近障碍",
    "D217" in AX
    and "R-Z-U3-PHASE-ASYMPTOTIC-NOGO" in BOTH
    and "渐近障碍" in DOC,
)
check(
    "N1 文档不否定 U3 字面接口但否定稳定非平凡模流",
    "U3\\text{ 字面要求" in DOC
    and "持久非平凡模时间" in DOC
    and "没有" in DOC,
)


# ==================================================================
head("N2  闭合词长度必为偶数")
# ==================================================================
sample_lengths = list(range(1, 12))
even_closed_words = []
for length in sample_lengths:
    for word in itertools.product((1, -1), repeat=length):
        if is_closed_word(word):
            even_closed_words.append(word)
            if length % 2 != 0:
                break
check(
    "N2 所有枚举到的闭合词长度都为偶数",
    all(len(word) % 2 == 0 for word in even_closed_words)
    and "|w|=r+s=2r" in DOC,
    f"closed_words={len(even_closed_words)}",
)


# ==================================================================
head("N3  偶数 T 只占据半数相位")
# ==================================================================
even_periods = (2, 4, 6, 8)
odd_phase_residues_empty = all(
    all((2 * r) % period != k for r in range(period) for k in range(period) if k % 2 == 1)
    for period in even_periods
)
check(
    "N3 偶数周期的奇数相位不可达",
    odd_phase_residues_empty,
)
check(
    "N3 文档登记偶数 T 的不忠实结论",
    "T\\text{ 偶}" in DOC
    and "不满足 }U3\\text{ 的忠实性" in DOC,
)


# ==================================================================
head("N4  奇数 T 的偶数位移生成全部相位")
# ==================================================================
odd_periods = (3, 5, 7, 9)
reach_ok = True
for period in odd_periods:
    reachable = {(2 * r) % period for r in range(period)}
    if reachable != set(range(period)):
        reach_ok = False
check(
    "N4 奇数周期由偶数位移生成整个相位环",
    reach_ok,
)
check(
    "N4 文档区分奇偶屏障与相位混合",
    "gcd(2,T)=1" in DOC
    and "相位混合" in DOC,
)


# ==================================================================
head("N5  T=3 递推矩阵的显式形式")
# ==================================================================
period = 3
shift = np.roll(np.eye(period), 1, axis=0)
matrix = np.eye(period) + 4.0 * np.linalg.matrix_power(shift, 2)
expected_matrix = np.eye(period)
for phase in range(period):
    expected_matrix[(phase + 2) % period, phase] += 4.0
check(
    "N5 T=3 递推矩阵为 I+4S^2",
    np.allclose(matrix, expected_matrix)
    and "(I+4S^2)N_n" in DOC,
)


# ==================================================================
head("N6  T=3 非均匀模式衰减率严格小于 1")
# ==================================================================
eigenvalues = recurrence_eigenvalues(period)
dominant = max(abs(value) for value in eigenvalues)
nonuniform = sorted(abs(value) for value in eigenvalues if abs(abs(value) - 5.0) > 1.0e-9)
decay_ratio = nonuniform[0] / dominant
check(
    "N6 主导特征值为 5",
    np.isclose(dominant, 5.0),
    f"dominant={dominant:.12f}",
)
check(
    "N6 非均匀模式比率为 sqrt(13)/5",
    np.allclose(nonuniform, [math.sqrt(13.0), math.sqrt(13.0)])
    and np.isclose(decay_ratio, math.sqrt(13.0) / 5.0),
    f"ratio={decay_ratio:.12f}",
)


# ==================================================================
head("N7  数值窗口偏离均匀量下降")
# ==================================================================
count_windows = (
    [0, 1, 1],
    [4, 3, 3],
    [12, 15, 15],
    [60, 55, 55],
    [220, 231, 231],
    [924, 903, 903],
    [3612, 3655, 3655],
    [14620, 14535, 14535],
)
deviations = [phase_deviation(counts) for counts in count_windows]
check(
    "N7 偏离均匀量整体下降且末期很小",
    deviations[-1] < deviations[1]
    and deviations[-1] < 0.01,
    f"deviations={['%.6f' % value for value in deviations]}",
)
check(
    "N7 文档列出同一组窗口计数",
    "(220,231,231)" in DOC.replace(" ", "")
    and "(14620,14535,14535)" in DOC.replace(" ", ""),
)


# ==================================================================
head("N8  渐近均匀态是中心态且模流恒等")
# ==================================================================
uniform = np.full(period, 1.0 / period)
check(
    "N8 均匀态为 I/T 且模流恒等",
    np.allclose(uniform, np.ones(period) / period)
    and "均匀中心态" in DOC
    and "模流恒等" in DOC,
)
check(
    "N8 文档给出渐近均匀结论",
    "q_{i,k}(m)" in DOC
    and "\\frac1T" in DOC,
)


# ==================================================================
head("N9  U3 字面接口与持久非平凡模时间被分开")
# ==================================================================
check(
    "N9 文档分开字面 U3 与物理非平凡时间",
    "U3\\text{ 字面要求}" in DOC
    and "非平凡物理时间" in DOC
    and "上游字面接口仍在" in DOC
    and "物理时间需要额外破缺源" in DOC,
)
check(
    "N9 文档不推翻 D215 的 U1-U4 字面接口",
    "字面接口没有被推翻" in DOC
    and "D215" in DOC,
)


# ==================================================================
head("N10 持久非中心性需要额外相位破缺源")
# ==================================================================
check(
    "N10 文档列出非均匀相位源与初始分支",
    "非均匀相位源" in DOC
    and "初始相位分支" in DOC
    and "跨区域关联" in DOC,
)
check(
    "N10 文档拒绝把历史计数当成稳定时间源",
    "暂态非中心模流" in DOC
    and "不能得到" in DOC,
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
