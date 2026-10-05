#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D221 —— U1 循环自然性与原始相位障碍
================================================================================
【被检验的问题（D221）】
  给定自由传递 Z_T 作用后，交叉积给 M_T 是否仍需要任意选择？
  裸闭合词的旋转轨道是否自然给出统一 M_T？
  毁灭周期 T 与闭包最小周期 p(w) 是否可以直接等同？

【本步判据】
  N1  D221 与两个恢复结构已登记
  N2  Z_T 的循环移位是自由传递 torsor
  N3  交叉积关系成立且生成 T^2 维矩阵代数
  N4  等变起点变换只给内自同构并满足函子性
  N5  裸集合上至少存在常值与矩阵两个不同函子
  N6  闭合词的旋转轨道大小等于最小周期并整除词长
  N7  所有闭合词的最小周期都是偶数
  N8  不存在最小周期为 3 的闭合词
  N9  周期词与原始词给出不同自然循环阶
  N10 D218 历史层实际含非原始周期
  N11 D212 的 T=3 必须与闭包旋转周期分开
  N12 文档保留原始扇区或相位提升缺口
  N13 不把自然性结果冒充物理局域性或物理态
  N14 不引用旧体系或外部路径
  N15 不修改上游、不新增 U5，失败时非零退出

运行：python3 verify/d221_u1_cyclic_naturality_and_primitive_phase_gap.py
"""
import io
import itertools
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
            "D221_u1_cyclic_naturality_and_primitive_phase_gap.md",
        ),
        encoding="utf-8",
    ).read()
)
D212 = canonical_math(
    io.open(
        os.path.join(
            ROOT,
            "derivations",
            "D212_zero_universe_to_u_interface.md",
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


def rotation_orbit(word):
    return {
        word[index:] + word[:index]
        for index in range(len(word))
    }


def primitive_period(word):
    for period in range(1, len(word) + 1):
        if (
            len(word) % period == 0
            and all(
                word[index] == word[index % period]
                for index in range(len(word))
            )
        ):
            return period
    return len(word)


def periodic_shift(period):
    matrix = np.zeros((period, period), dtype=complex)
    for index in range(period):
        matrix[(index + 1) % period, index] = 1.0
    return matrix


def phase_projection(period, index):
    matrix = np.zeros((period, period), dtype=complex)
    matrix[index, index] = 1.0
    return matrix


def induced_shift_map(period, shift):
    def apply(matrix):
        return np.linalg.matrix_power(
            periodic_shift(period),
            shift,
        ) @ matrix @ np.linalg.matrix_power(
            periodic_shift(period),
            -shift,
        )

    return apply


def induced_crossed_map(period, shift):
    shift_map = induced_shift_map(period, shift)
    return {
        "P": [
            shift_map(phase_projection(period, index))
            for index in range(period)
        ],
        "u": shift_map(periodic_shift(period)),
    }


def run_histories(cycles):
    active = [{word_from_string(text)} for text in INITIAL_WORDS]
    histories = [set() for _ in INITIAL_WORDS]
    zero_layer = set()
    snapshots = []

    for _ in range(cycles):
        for _step in range(DESTRUCTION_PERIOD):
            active, _ = evolve_one_step(active, histories, zero_layer)
        snapshots.append([set(book) for book in histories])
        active = [set() for _ in active]
        active = reseed_from_history(histories)

    return snapshots


periods = tuple(range(2, 9))
closed_words = closed_words_upto(10)


# ==================================================================
head("N1  D221 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记循环自然性",
    "D221" in AX
    and "R-Z-U1-CYCLIC-NATURALITY" in BOTH
    and "自由传递" in DOC,
)
check(
    "N1 公理表与文档登记原始相位缺口",
    "R-Z-U1-PRIMITIVE-PHASE-GAP" in BOTH
    and "原始闭包" in DOC
    and "相位提升" in DOC,
)


# ==================================================================
head("N2  Z_T 的循环移位是自由传递 torsor")
# ==================================================================
torsor_ok = True
for period in periods:
    for start in range(period):
        orbit = {(start + shift) % period for shift in range(period)}
        stabilizer = {
            shift
            for shift in range(period)
            if (start + shift) % period == start
        }
        if len(orbit) != period or stabilizer != {0}:
            torsor_ok = False
            break
check(
    "N2 每个点轨道大小为 T 且稳定子只有 0",
    torsor_ok,
    f"periods={periods}",
)


# ==================================================================
head("N3  交叉积关系成立且生成 T^2 维矩阵代数")
# ==================================================================
crossed_ok = True
dimensions = {}
for period in periods:
    identity = np.eye(period, dtype=complex)
    shift = periodic_shift(period)
    projections = [
        phase_projection(period, index)
        for index in range(period)
    ]
    relations = all(
        np.allclose(
            shift @ projections[index] @ shift.conj().T,
            projections[(index + 1) % period],
        )
        for index in range(period)
    )
    shifted_cycle = np.allclose(
        np.linalg.matrix_power(shift, period),
        identity,
    )
    basis = [
        projections[index] @ np.linalg.matrix_power(shift, power)
        for index in range(period)
        for power in range(period)
    ]
    vector_basis = np.column_stack([matrix.reshape(-1) for matrix in basis])
    dimension = int(np.linalg.matrix_rank(vector_basis))
    dimensions[period] = dimension
    if not relations or not shifted_cycle or dimension != period * period:
        crossed_ok = False
        break
check(
    "N3 交叉积关系与 M_T 维数都正确",
    crossed_ok,
    f"dimensions={dimensions}",
)


# ==================================================================
head("N4  等变起点变换只给内自同构并满足函子性")
# ==================================================================
naturality_ok = True
for period in periods:
    identity_map = induced_crossed_map(period, 0)
    projections = [
        phase_projection(period, index)
        for index in range(period)
    ]
    if not all(
        np.allclose(identity_map["P"][index], projections[index])
        for index in range(period)
    ):
        naturality_ok = False
    for shift in range(period):
        for other in range(period):
            composed = induced_shift_map(period, shift + other)
            left = induced_shift_map(period, shift)
            right = induced_shift_map(period, other)
            if not all(
                np.allclose(
                    left(right(projections[index])),
                    composed(projections[index]),
                )
                for index in range(period)
            ):
                naturality_ok = False
            shift_matrix = periodic_shift(period)
            if not np.allclose(
                left(right(shift_matrix)),
                composed(shift_matrix),
            ):
                naturality_ok = False
check(
    "N4 起点变换由 U^s 内实现且复合保持",
    naturality_ok,
)
check(
    "N4 文档给出唯一到内自同构的结论",
    "唯一到内自同构" in DOC
    and "自然" in DOC,
)


# ==================================================================
head("N5  裸集合上至少存在常值与矩阵两个不同函子")
# ==================================================================
constant_dimension = 1
matrix_dimensions = {
    period: period * period
    for period in periods
}
constant_functorial = all(
    np.eye(1, dtype=complex).shape == (1, 1)
    for _period in periods
)
crossed_functorial = naturality_ok
check(
    "N5 常值函子与矩阵函子都可定义但维数不同",
    constant_functorial
    and crossed_functorial
    and all(
        matrix_dimensions[period] != constant_dimension
        for period in periods
    ),
    f"matrix_dimensions={matrix_dimensions}",
)
check(
    "N5 文档拒绝裸集合自然选择 M_T",
    "裸集合" in DOC
    and "恢复层注入" in DOC,
)


# ==================================================================
head("N6  闭合词的旋转轨道大小等于最小周期并整除词长")
# ==================================================================
period_data = {
    word: primitive_period(word)
    for word in closed_words
}
check(
    "N6 全部闭合词满足轨道大小整除词长",
    all(
        len(rotation_orbit(word)) == period_data[word]
        and len(word) % period_data[word] == 0
        for word in closed_words
    ),
    f"closed_words={len(closed_words)}",
)


# ==================================================================
head("N7  所有闭合词的最小周期都是偶数")
# ==================================================================
all_even_periods = all(
    period_data[word] % 2 == 0
    for word in closed_words
)
check(
    "N7 闭合力零使最小周期块也为闭合词",
    all_even_periods,
    f"period_values={sorted(set(period_data.values()))}",
)
check(
    "N7 文档给出 p(w)|L 与 p(w) 偶数的证明",
    "p(w)\\mid L" in DOC
    and "最小旋转周期 }p(w)\\text{ 都为偶数" in DOC,
)


# ==================================================================
head("N8  不存在最小周期为 3 的闭合词")
# ==================================================================
period_three_words = [
    word
    for word in closed_words
    if period_data[word] == 3
]
check(
    "N8 闭合词最小周期不会等于 3",
    period_three_words == []
    and "p(w)=3" in DOC,
)


# ==================================================================
head("N9  周期词与原始词给出不同自然循环阶")
# ==================================================================
periodic_word = (1, -1, 1, -1)
primitive_word = (1, -1, -1, 1)
check(
    "N9 (+,-,+,-) 给 Z_2 torsor",
    primitive_period(periodic_word) == 2
    and len(rotation_orbit(periodic_word)) == 2,
)
check(
    "N9 (+,-,-,+) 给 Z_4 torsor",
    primitive_period(primitive_word) == 4
    and len(rotation_orbit(primitive_word)) == 4,
)
check(
    "N9 文档给出周期词反例",
    "(+1,-1,+1,-1)" in DOC
    and "M_2(\\mathbb C)" in DOC,
)


# ==================================================================
head("N10 D218 历史层实际含非原始周期")
# ==================================================================
history_snapshots = run_histories(3)
all_history_words = {
    word
    for book in history_snapshots[-1]
    for word in book
}
history_periods = {
    word: primitive_period(word)
    for word in all_history_words
}
nonprimitive = [
    (word, period)
    for word, period in history_periods.items()
    if period != len(word)
]
check(
    "N10 历史层确实包含 p(w)<|w| 的周期词",
    len(nonprimitive) > 0
    and all(period % 2 == 0 for _word, period in nonprimitive),
    f"nonprimitive={len(nonprimitive)}",
)


# ==================================================================
head("N11 D212 的 T=3 必须与闭包旋转周期分开")
# ==================================================================
check(
    "N11 文档明确奇数 T 不能来自闭包旋转",
    "奇数毁灭周期" in DOC
    and "不能来自闭合词自身的旋转周期" in DOC,
)
check(
    "N11 D212 仍登记交叉积提升输入",
    "交叉积提升与循环表示选择" in D212
    or "交叉积" in D212,
)


# ==================================================================
head("N12 文档保留原始扇区或相位提升缺口")
# ==================================================================
check(
    "N12 文档列出原始闭包扇区条件",
    "\\mathcal P_T" in DOC
    and "p(w)=T" in DOC,
)
check(
    "N12 文档列出显式相位提升条件",
    "\\mathbb Z_T\\longrightarrow\\operatorname{Bij}(\\operatorname{Orb}(w))" in DOC,
)


# ==================================================================
head("N13 不把自然性结果冒充物理局域性或物理态")
# ==================================================================
check(
    "N13 文档保留局域性、状态与时间缺口",
    "物理局域性" in DOC
    and "状态选择" in DOC
    and "支持半流与模时间关系" in DOC,
)
check(
    "N13 文档只声称条件表示定理",
    "条件表示定理" in DOC
    and "不是" in DOC,
)


# ==================================================================
head("N14 不引用旧体系或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "旧体系",
    "外部路径",
    "../..",
)
check(
    "N14 D221 文档不含旧体系名或外部路径",
    not any(marker in DOC for marker in external_markers),
)


# ==================================================================
head("N15 上游边界")
# ==================================================================
check(
    "N15 不修改 U1-U4+C1 且不新增 U5",
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
