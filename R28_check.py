#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R28_check.py -- 核验相位身份簇归约、纯规范 no-go 与张满条件。

对应 R28_phase_identity_cluster_reduction.md。
"""

from __future__ import annotations

import io
import os
import sys
from fractions import Fraction
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def head(title: str) -> None:
    print("")
    print("=" * 72)
    print(title)
    print("=" * 72)


def check(name: str, cond: bool, detail: str = "") -> None:
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def read(name: str) -> str:
    with io.open(os.path.join(HERE, name), encoding="utf-8") as handle:
        return handle.read()


def rank_matrix(rows: list[list[Fraction]]) -> int:
    if not rows:
        return 0
    a = [row[:] for row in rows]
    n_rows = len(a)
    n_cols = len(a[0])
    rank = 0
    col = 0
    while rank < n_rows and col < n_cols:
        pivot = None
        for r in range(rank, n_rows):
            if a[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            col += 1
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [x / scale for x in a[rank]]
        for r in range(n_rows):
            if r != rank and a[r][col] != 0:
                factor = a[r][col]
                a[r] = [a[r][c] - factor * a[rank][c] for c in range(n_cols)]
        rank += 1
        col += 1
    return rank


def transpose(a: list[list[Fraction]]) -> list[list[Fraction]]:
    return [list(row) for row in zip(*a)]


def incidence_edges(m: int) -> tuple[list[tuple[int, int]], list[list[Fraction]]]:
    edges = [(i, j) for i in range(m) for j in range(i + 1, m)]
    d_cols: list[list[Fraction]] = []
    for vertex in range(m):
        column = []
        for i, j in edges:
            value = Fraction(0)
            if i == vertex:
                value -= 1
            if j == vertex:
                value += 1
            column.append(value)
        d_cols.append(column)
    return edges, d_cols


R28 = read("R28_phase_identity_cluster_reduction.md")
R27 = read("R27_phase_cochain_quotient_and_dimension_dictionary.md")
R29 = read("R29_full_support_ledger_factorization_no_go.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")

head("F1  文档范围与结论")

check("标题与三重条件在位",
      "相位身份簇归约" in R28
      and "PHASE-IDENTITY-DER" in R28
      and "EDGE-CONNECTION" in R28
      and "HOLONOMY-FULL-SPAN" in R28
      and "INHERITANCE-IDENTITY" in R28)
check("核心 no-go 没有被写成正面定理",
      "PHASE-1-COCHAIN+PAIR-ID-QUOTIENT" in R28
      and "不能自动给出" in R28
      and "没有证明" in R28)
check("没有推出四维 GR",
      "没有由四维标签推出 Lorentz" in R28
      and "Lovelock" in R28
      and "GR" in R28)

head("F2  商维数与记录身份子空间")

for D in range(2, 30):
    m = D + 1
    edges, d_cols = incidence_edges(m)
    rank_d = rank_matrix(d_cols)
    quotient_dim = len(edges) - rank_d
    check("D=%d 时 rank(d)=m-1 且 quot=C(D,2)" % D,
          rank_d == m - 1 and quotient_dim == comb(D, 2),
          "rank_d=%d quotient=%d" % (rank_d, quotient_dim))

edges4, d_cols4 = incidence_edges(5)
D4 = 4
identity_rows = [[Fraction(1) if i == j else Fraction(0) for j in range(len(edges4))]
                 for i in range(len(edges4))]
rank_full = rank_matrix(d_cols4 + identity_rows)
check("D=4 时全部边记录张满 C^1",
      rank_full == len(edges4),
      "rank=%d" % rank_full)
check("D=4 时全部边记录在商中张满 C(D,2)",
      rank_full - rank_matrix(d_cols4) == comb(D4, 2))

head("F3  纯规范 no-go")

theta_rows = [
    [Fraction(1, 1), Fraction(2, 1), Fraction(3, 1), Fraction(4, 1), Fraction(5, 1)],
    [Fraction(0, 1), Fraction(1, 1), Fraction(-1, 1), Fraction(2, 1), Fraction(-2, 1)],
]
gauge_records = []
for theta in theta_rows:
    out = []
    for i, j in edges4:
        out.append(theta[j] - theta[i])
    gauge_records.append(out)

gauge_plus_d = gauge_records + d_cols4
rank_gauge_plus_d = rank_matrix(gauge_plus_d)
check("纯规范记录加到 dC^0 后不增加秩",
      rank_gauge_plus_d == rank_matrix(d_cols4),
      "rank=%d rank_d=%d" % (rank_gauge_plus_d, rank_matrix(d_cols4)))
check("纯规范记录的商身份维数为零",
      rank_gauge_plus_d - rank_matrix(d_cols4) == 0)

single_nongauge = [Fraction(1) if i == 0 else Fraction(0) for i in range(len(edges4))]
single_plus_d = [single_nongauge] + d_cols4
check("单条非纯规范边记录只给低维身份",
      rank_matrix(single_plus_d) - rank_matrix(d_cols4) == 1)

head("F4  张满判据与旧输入不足")

proper_subspace = [
    [Fraction(1) if i == 0 else Fraction(0) for i in range(len(edges4))],
]
proper_plus_d = proper_subspace + d_cols4
check("真子空间不满足 I_rec+dC^0=C^1",
      rank_matrix(proper_plus_d) < len(edges4))
check("文档写明维数相等不等于身份已经出现",
      "\\dim Q=\\binom D2" in R28
      and "\\dim H_{\\rm rec}=\\binom D2" in R28
      and "not\\Longrightarrow" in R28)
check("文档明确 R27 两条输入不足",
      "命题 R28.4" in R28
      and "PHASE-1-COCHAIN+PAIR-ID-QUOTIENT" in R28)

head("F5  与代价和代际问题正交")

check("身份秩与代价分开",
      "身份秩与身份代价是两个正交问题" in R28
      and "FULL-SUPPORT-LEDGER" in R28
      and "HOLONOMY-FULL-SPAN" in R28)
check("R28 已把代价簇交给 R29 四项分解",
      "R29_full_support_ledger_factorization_no_go.md" in R28
      and "LEDGER-FACTORIZATION" in R28
      and "DIR-SUPPORT-D" in R28
      and "RECORD-FAMILY-D" in R28
      and "PRODUCT-LEDGER" in R28
      and "SAME-Q" in R28)
check("代际账本仍独立",
      "WIPE-RESET-LEDGER" in R28
      and "L=4" in R28
      and "共同代际账本" in R28)

head("F6  全项目口径同步")

check("R27 已指向 R28",
      "R28_phase_identity_cluster_reduction.md" in R27
      and "PHASE-IDENTITY-DER" in R27)
check("R28 已指向 R29",
      "FULL-SUPPORT-LEDGER" in R29
      and "q^D" in R29
      and "R28-30" in R28)
check("STATUS 已登记 R28",
      "### 2.28" in STATUS
      and "PHASE-IDENTITY-DER" in STATUS
      and "HOLONOMY-FULL-SPAN" in STATUS
      and "R28_check.py" in STATUS)
check("STATUS 已登记 R29",
      "### 2.29" in STATUS
      and "LEDGER-FACTORIZATION" in STATUS
      and "R29_check.py" in STATUS)
check("INDEX 已收录 R28",
      "R28_phase_identity_cluster_reduction.md" in INDEX)
check("INDEX 已收录 R29",
      "R29_full_support_ledger_factorization_no_go.md" in INDEX)

head("F7  诚实边界")

check("R28 不关闭主要上游",
      "没有关闭 `DIM-SECTOR`、`EVO-NORM` 或 `SURV4-GLOBAL`" in R28)
check("R28 不把归约写成原生构造",
      "仍然开放" in R28
      or "仍是开放输入" in R28
      or "新增开放输入" in R28)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
