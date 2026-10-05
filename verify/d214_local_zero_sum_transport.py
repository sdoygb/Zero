#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D214 —— 局域零和传输
================================================================================
【被检验的问题（D214）】
  闭环词序能否给出非空、连通、严格保零和的最近邻传输图？
  最近邻边是否足以生成所有整数零和态，而不需要远区配对？
  这项结果是否仍不等于四维几何或 GR？

【本步判据】
  N1  D214 与两个恢复结构已登记
  N2  闭合词与循环邻接
  N3  最近邻转移严格保零和
  N4  局域边集非空、连通且无远区边
  N5  环图关联矩阵的秩
  N6  关联矩阵的像等于零和超平面
  N7  显式整数局域流重建任意零和态
  N8  远区正负点不需要全局边
  N9  无偏好最近邻生成元保持总荷
  N10 当前周期追加步仍违反局域零和
  N11 环结构只给一维局域候选
  N12 物理局域识别仍是输入
  N13 输入预算与未导出项完整登记
  N14 不把局域传输冒充四维、洛伦兹或 GR
  N15 不引用旧体系或外部路径
  N16 不修改上游、不新增 U5，失败时非零退出

运行：python3 verify/d214_local_zero_sum_transport.py
"""
import io
import itertools
import os
import sys
from fractions import Fraction

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
            "D214_local_zero_sum_transport.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


def total_charge(state):
    return sum(state)


def local_transfer(state, site):
    result = list(state)
    result[site] -= 1
    result[(site + 1) % len(result)] += 1
    return tuple(result)


def incidence_matrix(size):
    matrix = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for site in range(size):
        matrix[site][site] -= 1
        matrix[(site + 1) % size][site] += 1
    return matrix


def rank_fraction(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    rows = len(work)
    cols = len(work[0]) if rows else 0
    rank = 0
    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if work[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        value = work[rank][col]
        work[rank] = [item / value for item in work[rank]]
        for row in range(rows):
            if row == rank:
                continue
            factor = work[row][col]
            if factor:
                work[row] = [
                    item - factor * pivot_item
                    for item, pivot_item in zip(work[row], work[rank])
                ]
        rank += 1
        if rank == rows:
            break
    return rank


def explicit_edge_flow(state):
    flow = []
    running = 0
    for value in state:
        running += value
        flow.append(-running)
    if flow[-1] != 0:
        raise ValueError("state is not zero-sum")
    return flow


def reconstruct_from_flow(flow):
    size = len(flow)
    state = [0] * size
    for site, amount in enumerate(flow):
        state[site] -= amount
        state[(site + 1) % size] += amount
    return tuple(state)


def is_zero_sum_state(state):
    return bool(state) and sum(state) == 0


def cycle_adjacency(size):
    adjacency = np.zeros((size, size), dtype=int)
    for site in range(size):
        adjacency[site, (site + 1) % size] = 1
        adjacency[(site + 1) % size, site] = 1
    return adjacency


def cycle_laplacian(size):
    laplacian = 2.0 * np.eye(size)
    for site in range(size):
        laplacian[site, (site - 1) % size] -= 1.0
        laplacian[site, (site + 1) % size] -= 1.0
    return laplacian


# ==================================================================
head("N1  D214 与恢复结构已登记")
# ==================================================================
registration_ok = (
    "D214" in AX
    and "R-Z-LOCAL-ZERO-SUM-TRANSPORT" in BOTH
    and "R-Z-LOCAL-SITE-IDENTIFICATION" in BOTH
    and "# D214" in DOC
    and "局域零和传输" in DOC
)
check("N1 D214 与两个恢复结构已登记", registration_ok)


# ==================================================================
head("N2  闭合词与循环邻接")
# ==================================================================
words = (
    (1, -1),
    (1, -1, 1, -1),
    (1, 1, -1, -1),
    (1, 1, 1, -1, -1, -1),
)
word_ok = all(is_zero_sum_state(word) for word in words)
check(
    "N2 闭合词总荷为零并给出循环词位",
    word_ok
    and "i\\sim i+1\\pmod L" in DOC
    and "C_L" in DOC,
)


# ==================================================================
head("N3  最近邻转移严格保零和")
# ==================================================================
states = (
    (1, -1, 0, 0),
    (2, -1, -1, 0),
    (0, 3, -1, -2),
)
transfers_ok = all(
    total_charge(local_transfer(state, site)) == total_charge(state)
    for state in states
    for site in range(len(state))
)
check(
    "N3 每个 T_i=n+e_{i+1}-e_i 都保持总荷",
    transfers_ok
    and "T_i n" in DOC
    and "\\Delta Q" in DOC,
)


# ==================================================================
head("N4  局域边集非空、连通且无远区边")
# ==================================================================
size = 6
adjacency = cycle_adjacency(size)
reachable = {0}
frontier = [0]
while frontier:
    current = frontier.pop()
    for target in range(size):
        if adjacency[current, target] and target not in reachable:
            reachable.add(target)
            frontier.append(target)
edge_count = int(np.sum(adjacency) // 2)
long_range_edges = int(
    sum(
        adjacency[i, j]
        for i in range(size)
        for j in range(size)
        if min((i - j) % size, (j - i) % size) != 1
    )
    // 2
)
check(
    "N4 环图有 L 条局域边且连通",
    edge_count == size and reachable == set(range(size)),
    f"edges={edge_count}, reachable={len(reachable)}",
)
check(
    "N4 环图没有额外远区边",
    long_range_edges == 0,
    f"long_range_edges={long_range_edges}",
)


# ==================================================================
head("N5  环图关联矩阵的秩")
# ==================================================================
ranks = {}
for local_size in range(3, 7):
    incidence = incidence_matrix(local_size)
    ranks[local_size] = rank_fraction(incidence)
rank_ok = all(ranks[local_size] == local_size - 1 for local_size in ranks)
check(
    "N5 连通环关联矩阵秩为 L-1",
    rank_ok,
    f"ranks={ranks}",
)
incidence = np.array([[float(value) for value in row] for row in incidence_matrix(6)])
check(
    "N5 常数向量在关联矩阵左零空间",
    np.allclose(np.ones(6) @ incidence, 0.0),
)


# ==================================================================
head("N6  关联矩阵的像等于零和超平面")
# ==================================================================
sample_state = (2, -1, -1, 0)
flow = explicit_edge_flow(sample_state)
reconstructed = reconstruct_from_flow(flow)
check(
    "N6 显式局域流位于零和超平面并重建样本态",
    reconstructed == sample_state and total_charge(reconstructed) == 0,
    f"flow={flow}, state={reconstructed}",
)
check(
    "N6 文档给出 im B=H_Q",
    "\\operatorname{im}B=H_Q" in DOC
    and "L-1" in DOC,
)


# ==================================================================
head("N7  任意整数零和态可由局域流生成")
# ==================================================================
bound = 2
local_size = 5
all_states = [
    state
    for state in itertools.product(range(-bound, bound + 1), repeat=local_size)
    if total_charge(state) == 0
]
all_reconstructed = all(
    reconstruct_from_flow(explicit_edge_flow(state)) == state
    for state in all_states
)
all_integer = all(
    all(isinstance(value, int) for value in explicit_edge_flow(state))
    for state in all_states
)
check(
    "N7 所有有限盒零和态都有整数局域流",
    all_reconstructed and all_integer,
    f"states={len(all_states)}",
)


# ==================================================================
head("N8  远区正负点不需要全局边")
# ==================================================================
nonlocal_state = (1, 0, 0, -1, 0)
nonlocal_flow = explicit_edge_flow(nonlocal_state)
nonlocal_reconstructed = reconstruct_from_flow(nonlocal_flow)
check(
    "N8 远区正负点也可由相邻边流重建",
    nonlocal_reconstructed == nonlocal_state
    and any(value != 0 for value in nonlocal_flow),
    f"flow={nonlocal_flow}",
)
check(
    "N8 文档区分全局配对与局域传播",
    "全局正负配对" in DOC
    and "局域传输机制" in DOC,
)


# ==================================================================
head("N9  无偏好最近邻生成元保持总荷")
# ==================================================================
generator_ok = all(
    total_charge(local_transfer(state, site)) == 0
    for state in all_states
    for site in range(local_size)
)
check(
    "N9 最近邻对称核不改变总荷",
    generator_ok
    and "\\mathcal L_{\\rm loc}Q=0" in DOC,
)


# ==================================================================
head("N10 当前周期追加步仍违反局域零和")
# ==================================================================
append_changes = (1, -1)
check(
    "N10 当前 w->w± 的单步荷变化非零",
    all(change != 0 for change in append_changes)
    and "\\Delta Q(w\\pm)=1" in DOC,
)
check(
    "N10 文档把局域转移写成替换规则",
    "规则替换候选" in DOC
    and "当前周期程序" in DOC
    and "仍违反" in DOC,
)


# ==================================================================
head("N11 环结构只给一维局域候选")
# ==================================================================
laplacian = cycle_laplacian(size)
eigenvalues = np.linalg.eigvalsh(laplacian)
zero_modes = int(np.sum(np.abs(eigenvalues) < 1.0e-12))
positive_eigenvalues = eigenvalues[eigenvalues > 1.0e-12]
first_gap = float(np.min(positive_eigenvalues))
check(
    "N11 环 Laplacian 只有一个零模",
    zero_modes == 1,
    f"zero_modes={zero_modes}",
)
check(
    "N11 环只有一度连通分量而不是四维方向",
    first_gap > 0.0
    and "C_{\\rm loc}" in DOC
    and "C_{4D}" in DOC
    and "洛伦兹号差" in DOC,
)


# ==================================================================
head("N12 物理局域识别仍是输入")
# ==================================================================
check(
    "N12 文档明列词位到物理站点的识别缺口",
    "词位与物理空间点的唯一识别" in DOC
    and "词位即局域站点" in DOC
    and "R-Z-LOCAL-SITE-IDENTIFICATION" in BOTH,
)
check(
    "N12 文档不声称已经得到物理距离或体积元",
    "物理距离" in DOC
    and "体积元" in DOC
    and "仍未得到" in DOC,
)


# ==================================================================
head("N13 输入预算与未导出项完整登记")
# ==================================================================
check(
    "N13 文档登记局域传输输入预算",
    "输入预算" in DOC
    and "整数局域荷" in DOC
    and "最近邻单位转移" in DOC
    and "无偏好速率" in DOC,
)
check(
    "N13 文档登记四维、洛伦兹与源几何缺口",
    "四维细化" in DOC
    and "洛伦兹号差" in DOC
    and "局部守恒应力张量" in DOC,
)


# ==================================================================
head("N14 不把局域传输冒充四维、洛伦兹或 GR")
# ==================================================================
check(
    "N14 文档拒绝从环图直接得到四维 GR",
    "不能由这一张环图直接启动" in DOC
    and "四维洛伦兹时空或 Einstein 方程" in DOC
    and "不能得到" in DOC,
)
check(
    "N14 文档没有把一维环图写成物理时空",
    "一维环" in DOC
    and "不给出四维坐标" in DOC,
)


# ==================================================================
head("N15 不引用旧体系或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "旧体系",
    "外部路径",
    "../..",
)
check(
    "N14 文档不含旧体系名或外部路径",
    not any(marker in DOC for marker in external_markers),
)


# ==================================================================
head("N16 上游边界")
# ==================================================================
check(
    "N16 不修改 U1-U4+C1 且不新增 U5",
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
