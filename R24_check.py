#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R24_check.py -- 核验全局四维生存峰门槛与 GR 临时脚手架边界。

对应 R24_global_four_survival_gate.md。
"""

from __future__ import annotations

import io
import math
import os
import sys

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


def argmax(values: dict[int, float], rel_tol: float = 1.0e-12) -> list[int]:
    maximum = max(values.values())
    return [
        key
        for key, value in values.items()
        if math.isclose(value, maximum, rel_tol=rel_tol, abs_tol=0.0)
    ]


R24 = read("R24_global_four_survival_gate.md")
R23 = read("R23_dimension_descendant_selection.md")
R26 = read("R26_pair_carrier_reduction_no_go.md")
R3 = read("R3_dimension_selection.md")
R0 = read("R0_publication_theorem.md")
STATUS = read("STATUS.md")
DIMS = list(range(1, 25))

head("F1  GR 临时脚手架边界")

check("GR-LB 明确写成 TEMP-GR 而非最终证明",
      "TEMP-GR" in R24
      and "临时脚手架" in R24
      and "最终必须撤掉" in R24)
check("最终目标 SURV4-GLOBAL 使用全域候选",
      "SURV4-GLOBAL" in R24
      and "\\mathcal D=\\{1,2,3,4,\\ldots\\}" in R24
      and "不得使用 (R24-1)" in R24)
check("GR 内排序单独登记为脚手架结果",
      "SURV4-GR-SCAFFOLD" in R24
      and "不是最终目标" in R24)

head("F2  乘积存活率全局 no-go")

for a in (0.10, 0.50, 0.90, 0.99):
    winners = argmax({D: a**D for D in DIMS})
    check("a=%.2f 时全局峰唯一在 D=1" % a, winners == [1])

check("乘积 no-go 与不含 GR 的目标区分清楚",
      "\\arg\\max_{D\\ge1}a_L^D=\\{1\\}" in R24
      and "不能证明 (R24-2)" in R24
      and "只能与 `GR-LB` 合起来证明 (R24-3)" in R24)

head("F3  不使用 GR 的全局四维窗口")

for q in (0.76, 0.7788, 0.79):
    winners = argmax({D: D * q**D for D in DIMS})
    check("q=%.4f 时全局唯一峰在 D=4" % q, winners == [4])

check("q=3/4 与 4/5 是并列边界",
      argmax({D: D * 0.75**D for D in DIMS}) == [3, 4]
      and argmax({D: D * 0.80**D for D in DIMS}) == [4, 5])
check("窗口公式与 R23 定理一致",
      "\\frac34<q<\\frac45" in R24
      and "命题 R24.2" in R24
      and "不使用 GR" in R24)

head("F4  缺失桥与循环性边界")

check("DIM-COST-Q 被登记为开放且不得按四维调参",
      "DIM-COST-Q" in R24
      and "**开放**" in R24
      and "定义不使用 }D=4" in R24
      and "不能只对 $D=4$ 调参" in R24)
check("DIM-INTERACT 不是把答案写进选择器",
      "DIM-INTERACT" in R24
      and "那只是把答案写进选择器，判为循环" in R24)
check("R25 把 DIM-INTERACT 收窄为 PAIR-CARRIER 且保留开放边界",
      "R25_native_pair_cost_and_four_dim_peak.md" in R24
      and "PAIR-CARRIER" in R24
      and "1/2<q<3/5" in R24
      and "PAIR-CARRIER-DER" in R24
      and "**开放**" in R24)
check("R26 把 PAIR-CARRIER-DER 继续拆为六项",
      "R26_pair_carrier_reduction_no_go.md" in R24
      and "PAIR-GRAPH-KD" in R24
      and "PAIR-COST-FACTORIZATION" in R24
      and "PAIR-NO-EXTRA-MULT" in R24
      and "D=3" in R24
      and "R26_pair_carrier_reduction_no_go.md" in R26)
check("几何探针的反向信号在位",
      "no stable four-dimensional plateau detected" in R24
      and "不能把 `DIM-DESC` 的 $q$ 窗口冒充 Zero 原生选维" in R24)

head("F5  全项目口径同步")

check("R23 指向 R24 并声明 GR 只是临时脚手架",
      "R24_global_four_survival_gate.md" in R23
      and "GR-LB" in R23
      and "临时" in R23)
check("R3 与 R0 都已登记 R24",
      "R24_global_four_survival_gate.md" in R3
      and "R24" in R0)
check("STATUS 登记 R24 且不把目标说成已证",
      "R24" in STATUS
      and "全局四维生存峰" in STATUS
      and "未证" in STATUS)

head("F6  反过度主张")

check("没有把条件窗口说成 Zero 原生四维选择",
      "不能声称“四维生存率较高”已从 Zero 证成" in R24)
check("没有把四维标签等同于四维 GR",
      "四维度规、Lorentz、GR" in R24 and "仍未由此推出" in R24)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
