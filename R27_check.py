#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R27_check.py -- 核验相位边 1-上链商、字典汇流与周期归一化风险。

对应 R27_phase_cochain_quotient_and_dimension_dictionary.md。
"""

from __future__ import annotations

import io
import math
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


def argmax_fraction(values: dict[int, Fraction]) -> list[int]:
    maximum = max(values.values())
    return [key for key, value in values.items() if value == maximum]


def pair_mode_q(q: Fraction, offset: int = 0, dmax: int = 60) -> dict[int, Fraction]:
    return {D: comb(D, 2) * q ** (D + offset) for D in range(2, dmax + 1)}


R27 = read("R27_phase_cochain_quotient_and_dimension_dictionary.md")
R26 = read("R26_pair_carrier_reduction_no_go.md")
R25 = read("R25_native_pair_cost_and_four_dim_peak.md")
R23 = read("R23_dimension_descendant_selection.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")

head("F1  文档范围与结论")

check("标题和核心恒等式在位",
      "相位上链商" in R27
      and "C(m,2)" in R27
      and "\\binom m2-(m-1)=\\binom D2" in R27)
check("条件结果没有被写成无条件定理",
      "这四项仍不是 Z0 的无条件推论" in R27
      and "条件每代四维峰" in R27)
check("采用共同毁灭代际账本且未保留旧并行周期条款",
      "WIPE-RESET-LEDGER" in R27
      and "共同毁灭代际" in R27
      and "PARALLEL-PERIOD" not in R27)
check("R27 已把身份簇交给 R28 补满秩",
      "R28_phase_identity_cluster_reduction.md" in R27
      and "PHASE-IDENTITY-DER" in R27
      and "HOLONOMY-FULL-SPAN" in R27)
check("R27 已把代价簇交给 R29 分解",
      "R29_full_support_ledger_factorization_no_go.md" in R27
      and "LEDGER-FACTORIZATION" in R27
      and "PRODUCT-LEDGER" in R27
      and "SAME-Q" in R27)
check("没有推出四维 GR",
      "没有由四维标签推出 Lorentz" in R27
      and "Lovelock" in R27
      and "GR" in R27)

head("F2  零和秩与相位商维数")

for D in range(2, 40):
    m = D + 1
    raw_edges = comb(m, 2)
    gauge_rank = m - 1
    quotient = raw_edges - gauge_rank
    check("D=%d 时 C(m,2)-D=C(D,2)" % D,
          quotient == comb(D, 2),
          "m=%d raw=%d gauge=%d quotient=%d" % (m, raw_edges, gauge_rank, quotient))

check("梯度核维数为 1 的证明锚点在位",
      "\\ker d=\\mathbb R\\mathbf 1" in R27
      and "\\operatorname{rank}d" in R27)
check("圈秩改写无算术错误",
      "|E|-|V|+1" in R27
      and "\\binom m2-m+1" in R27
      and "\\frac{(m-1)(m-2)}2" in R27)

head("F3  原始边读法与商类读法分叉")

for D in range(2, 20):
    check("D=%d 时原始边数为 C(D+1,2)=D+C(D,2)" % D,
          comb(D + 1, 2) == D + comb(D, 2),
          "%d=%d+%d" % (comb(D + 1, 2), D, comb(D, 2)))

check("原始边身份仍保留 R26 的字典 no-go",
      "R26 的字典 no-go 没有被取消" in R27
      and "原始边身份" in R27
      and "商类身份" in R27)
check("R27 明确采用 alpha 而不是 beta",
      "D:=m-1" in R27
      and "DIR-DICT-ALPHA" in R27
      and "不再要求 `β`" in R27)

head("F4  q=5/9 的每代峰")

q = Fraction(5, 9)

for offset in range(0, 12):
    check("固定偏移 %d 下唯一峰仍为 D=4" % offset,
          argmax_fraction(pair_mode_q(q, offset)) == [4],
          str(argmax_fraction(pair_mode_q(q, offset))))

check("成对窗口为 1/2<q<3/5",
      argmax_fraction(pair_mode_q(Fraction(1, 2))) == [3, 4]
      and argmax_fraction(pair_mode_q(Fraction(3, 5))) == [4, 5])
check("q=5/9 在窗口内",
      Fraction(1, 2) < q < Fraction(3, 5))

head("F5  代价重数反例仍受控")

support = {D: comb(D, 2) * q**2 for D in range(2, 60)}
edge = {D: comb(D, 2) * q ** comb(D, 2) for D in range(2, 60)}
check("只付两维支撑代价时无有限峰",
      all(support[D + 1] > support[D] for D in range(2, 59)))
check("每条边独立付代价时唯一峰为 D=2",
      argmax_fraction(edge) == [2],
      str(argmax_fraction(edge)))
check("全支撑账本未被写成自动推论",
      "FULL-SUPPORT-LEDGER" in R27 and "开放具名输入" in R27)

head("F6  共同毁灭代际下的长期比较")

generation_factor = pair_mode_q(q)
check("每代乘法 F_D 的唯一峰仍为 D=4",
      argmax_fraction(generation_factor) == [4],
      str(argmax_fraction(generation_factor)))

lambda_generation = {
    D: math.log(float(value))
    for D, value in generation_factor.items()
}
max_generation = max(lambda_generation.values())
generation_winners = [
    D for D, value in lambda_generation.items()
    if math.isclose(value, max_generation, rel_tol=1e-14)
]
check("log F_D 的代际增长峰仍为 D=4",
      generation_winners == [4],
      "argmax=%s" % generation_winners)

population = {D: Fraction(1) for D in range(2, 40)}
for _ in range(20):
    population = {D: population[D] * generation_factor[D] for D in population}
check("二十个共同毁灭代际后的谱系峰仍为 D=4",
      argmax_fraction(population) == [4],
      str(argmax_fraction(population)))

lambda_serial = {
    D: math.log(float(comb(D, 2) * q**D)) / D
    for D in range(2, 40)
}
max_serial = max(lambda_serial.values())
serial_winners = [
    D for D, value in lambda_serial.items()
    if math.isclose(value, max_serial, rel_tol=1e-14)
]
check("自由运行串行模型只作边界且不再伪装成演化层定律",
      4 not in serial_winners and "不再适用于 D211" in R27,
      "argmax=%s" % serial_winners)
check("文档登记共同代际输入、局部异步边界与代际重复",
      "WIPE-RESET-LEDGER" in R27
      and "D211" in R27
      and "D222" in R27
      and "活动层全部毁灭" in R27
      and "历史层只保留最高两层亚层" in R27
      and "全部保留" in R27
      and "LOCAL-GENERATION-LEDGER" in R27
      and "GENERATION-MULTIPLICITY" in R27
      and "\\tau D" in R27)
check("代际相对峰没有被写成绝对演化层占比",
      "绝对演化层占比" in R27
      and "EVO-NORM" in R27
      and "不能替代 `EVO-NORM`" in R27)

head("F7  全项目口径同步")

check("R27 明确 R26 六项接口的重新组织",
      "与 R26 六项接口的对照" in R27
      and "PHASE-1-COCHAIN" in R27
      and "PAIR-ID-QUOTIENT" in R27
      and "FULL-SUPPORT-LEDGER" in R27)
check("R27 只把身份簇更新为 R28 的条件接口",
      "PHASE-IDENTITY-DER" in R27
      and "EDGE-CONNECTION" in R27
      and "INHERITANCE-IDENTITY" in R27)
check("R27 只把代价簇更新为 R29 的四项接口",
      "LEDGER-FACTORIZATION" in R27
      and "DIR-SUPPORT-D" in R27
      and "RECORD-FAMILY-D" in R27
      and "PRODUCT-LEDGER" in R27
      and "SAME-Q" in R27)
check("R26 仍保持原始边身份 no-go",
      "连通性不推出 `K_D`" in R26
      and "C(D+1,2)" in R26
      and "唯一峰是三维" in R26)
check("R25 原始窗口仍在文档中",
      "1/2<q<3/5" in R25.replace(" ", "")
      or "\\frac12<q<\\frac35" in R25)
check("STATUS 已登记 R27",
      "### 2.27" in STATUS and "R27" in STATUS)
check("STATUS 已登记 R29",
      "### 2.29" in STATUS and "LEDGER-FACTORIZATION" in STATUS)
check("INDEX 已收录 R27",
      "R27_phase_cochain_quotient_and_dimension_dictionary.md" in INDEX)

head("F8  反过度主张")

check("R27 不关闭 SURV4-GLOBAL",
      "没有关闭 `DIM-SECTOR`、`EVO-NORM` 或 `SURV4-GLOBAL`" in R27)
check("R27 不把条件四维峰升级成 GR",
      "条件每代四维峰" in R27 and "四维 Lorentz" in R27)
check("R27 保留四项原生物理桥（周期项已改为共同代际账本）",
      "PHASE-IDENTITY-DER" in R27
      and "FULL-SUPPORT-LEDGER" in R27
      and "WIPE-RESET-LEDGER" in R27)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
