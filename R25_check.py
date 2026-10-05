#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R25_check.py -- 核验原生账本单方向 no-go 与成对连接四维峰。

对应 R25_native_pair_cost_and_four_dim_peak.md。
"""

from __future__ import annotations

import io
import json
import math
import os
import sys
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb, exp

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


def canonical(word: tuple[int, ...]) -> tuple[int, ...]:
    return min(word[i:] + word[:i] for i in range(len(word)))


def ledger_q(L: int) -> Fraction:
    words = [w for w in product((1, -1), repeat=L) if sum(w) == 0]
    classes = Counter(canonical(w) for w in words)
    return sum(Fraction(size, len(words)) ** 2 for size in classes.values())


def argmax(values: dict[int, Fraction]) -> list[int]:
    maximum = max(values.values())
    return [D for D, value in values.items() if value == maximum]


def single(q: Fraction, dmax: int = 40) -> dict[int, Fraction]:
    return {D: D * q**D for D in range(1, dmax + 1)}


def pair(q: Fraction, dmax: int = 40) -> dict[int, Fraction]:
    return {D: comb(D, 2) * q**D for D in range(2, dmax + 1)}


R25 = read("R25_native_pair_cost_and_four_dim_peak.md")
R23 = read("R23_dimension_descendant_selection.md")
R24 = read("R24_global_four_survival_gate.md")
R26 = read("R26_pair_carrier_reduction_no_go.md")
R3 = read("R3_dimension_selection.md")
R0 = read("R0_publication_theorem.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")

with io.open(os.path.join(HERE, "R25_native_pair_cost_results.json"), encoding="utf-8") as handle:
    RESULTS = json.load(handle)

head("F1  文档边界与正面候选")

check("标题与对象在位",
      "成对维度增益" in R25 and "四维全局峰" in R25 and "PAIR-CARRIER" in R25)
check("明确没有使用 GR 或 D>=4",
      "不用 GR" in R25 and "不用 `D≥4`" in R25)
check("明确仍未推出四维 GR",
      "没有" in R25 and "四维 Lorentz" in R25 and "Einstein 方程" in R25)
check("PAIR-CARRIER 登记为具名结构输入",
      "\\texttt{PAIR-CARRIER}" in R25 and "具名结构输入" in R25)
check("LEDGER-ROT 登记为具名读出",
      "\\texttt{LEDGER-ROT}" in R25 and "具名读出" in R25)

head("F2  原生旋转类账本纯度")

for L in (4, 6, 8, 10):
    q = ledger_q(L)
    bound = Fraction(L, comb(L, L // 2))
    check("L=%d 的 q 可由旋转类独立复算且 q<=L/N" % L,
          q <= bound,
          "q=%s, L/N=%s" % (q, bound))
    check("L=%d 的 L/N < 3/4" % L, bound < Fraction(3, 4), str(bound))

check("q4=5/9 精确",
      ledger_q(4) == Fraction(5, 9) and "q_4" in R25 and "\\frac59" in R25)
check("单方向 no-go 引理在位",
      "引理 R25.2" in R25 and "q_L\\le\\frac{L}{N_L}\\le\\frac23<\\frac34" in R25)

head("F3  单方向路线 no-go")

q4 = ledger_q(4)
check("q4 单方向唯一峰为 D=2",
      argmax(single(q4)) == [2],
      str(argmax(single(q4))))
check("q4 单方向 F4<F3",
      single(q4)[4] < single(q4)[3],
      "F3=%s, F4=%s" % (single(q4)[3], single(q4)[4]))
check("单方向 no-go 与 R23 窗口区分",
      "旧窗口" not in R25
      and "单方向模型" in R25
      and "\\frac34<q<\\frac45" in R25
      and "\\frac12<q<\\frac35" in R25)

head("F4  成对连接路线")

check("成对窗口左右边界正确",
      argmax(pair(Fraction(1, 2))) == [3, 4]
      and argmax(pair(Fraction(3, 5))) == [4, 5],
      "%s / %s" % (argmax(pair(Fraction(1, 2))), argmax(pair(Fraction(3, 5)))))
check("q=5/9 在窗口内且唯一峰 D=4",
      Fraction(1, 2) < q4 < Fraction(3, 5)
      and argmax(pair(q4)) == [4],
      str(argmax(pair(q4))))
check("q=5/9 的 pair 前六项与文档一致",
      pair(q4)[2] == Fraction(25, 81)
      and pair(q4)[3] == Fraction(125, 243)
      and pair(q4)[4] == Fraction(1250, 2187)
      and pair(q4)[5] == Fraction(31250, 59049)
      and pair(q4)[6] == Fraction(78125, 177147)
      and "\\frac{25}{81}" in R25
      and "\\frac{125}{243}" in R25
      and "\\frac{1250}{2187}" in R25
      and "\\frac{31250}{59049}" in R25
      and "\\frac{78125}{177147}" in R25)

head("F5  L=4 特殊性与连续代价负对照")

check("L=6 的原生 q 给出 pair 峰 D=2",
      argmax(pair(ledger_q(6))) == [2],
      str(argmax(pair(ledger_q(6)))))
check("L>=6 的原生 q 小于 1/2",
      all(ledger_q(L) < Fraction(1, 2) for L in range(6, 21, 2)))
check("L=2 不提供有限四维峰",
      ledger_q(2) == 1 and "L=2" in R25)

cont = exp(-0.25)
cont_values = {D: comb(D, 2) * cont**D for D in range(2, 41)}
cont_max = max(cont_values.values())
cont_winners = [D for D, value in cont_values.items() if math.isclose(value, cont_max, rel_tol=1e-14)]
check("q=exp(-1/4) 的 pair 负对照不选四维",
      4 not in cont_winners,
      "argmax=%s" % cont_winners)
check("文档明确 q=exp(-1/L) 与离散记录过程不一致",
      "记录率 }1/L" in R25 and "插值错误" in R25 and "e^{-1/4}" in R25)

head("F6  独立探针 JSON")

rows = {row["L"]: row for row in RESULTS["rows"]}
check("探针 JSON 覆盖 L=2..20",
      sorted(rows) == list(range(2, 21, 2)))
check("探针 q4 与核验脚本一致",
      Fraction(rows[4]["q"]) == q4)
check("探针单方向 q4 峰为 2",
      rows[4]["single_argmax"] == [2])
check("探针成对 q4 峰为 4",
      rows[4]["pair_argmax"] == [4])
check("探针边界与文档一致",
      RESULTS["pair_boundaries"]["q=1/2"] == [3, 4]
      and RESULTS["pair_boundaries"]["q=3/5"] == [4, 5])

head("F7  全项目口径同步")

check("R24 指向 R25 且登记 PAIR-CARRIER",
      "R25_native_pair_cost_and_four_dim_peak.md" in R24
      and "PAIR-CARRIER" in R24)
check("R23 指向 R25 并区分单方向与成对窗口",
      "R25_native_pair_cost_and_four_dim_peak.md" in R23
      and "PAIR-CARRIER" in R23
      and "1/2<q<3/5" in R23)
check("R3 与 R0 已登记 R25",
      "R25_native_pair_cost_and_four_dim_peak.md" in R3
      and "R25_native_pair_cost_and_four_dim_peak.md" in R0)
check("STATUS 登记 R25 且不把条件候选写成已证",
      "R25" in STATUS
      and "PAIR-CARRIER" in STATUS
      and "条件证成" in STATUS
      and "PAIR-CARRIER-DER" in STATUS)
check("R26 已把 PAIR-CARRIER-DER 拆为六项且保留条件边界",
      "R26_pair_carrier_reduction_no_go.md" in R25
      and "PAIR-GRAPH-KD" in R25
      and "PAIR-COST-FACTORIZATION" in R25
      and "PAIR-NO-EXTRA-MULT" in R25
      and "D=3" in R25
      and "六个可分别审计的接口" in R26)
check("INDEX 自动收录 R25",
      "R25_native_pair_cost_and_four_dim_peak.md" in INDEX)

head("F8  反过度主张")

check("四维峰没有升级成四维 GR",
      "仍只是候选维数标签" in R25
      and "四维 Lorentz" in R25
      and "没有" in R25)
check("PAIR-CARRIER-DER 仍未证",
      "PAIR-CARRIER-DER" in R25
      and "**开放**" in R25)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
