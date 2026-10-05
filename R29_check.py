#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R29_check.py -- 核验全支撑账本归约、乘积记录 no-go 与最小四项输入。

对应 R29_full_support_ledger_factorization_no_go.md。
"""

from __future__ import annotations

import io
import os
import sys
from fractions import Fraction
from math import comb, sqrt

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


def argmax_d(values: list[Fraction], dmin: int = 2) -> list[int]:
    maximum = max(values)
    return [dmin + index for index, value in enumerate(values) if value == maximum]


R29 = read("R29_full_support_ledger_factorization_no_go.md")
R28 = read("R28_phase_identity_cluster_reduction.md")
R27 = read("R27_phase_cochain_quotient_and_dimension_dictionary.md")
R26 = read("R26_pair_carrier_reduction_no_go.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")

head("F1  文档范围与核心结论")

check("标题与四项归约在位",
      "全支撑账本归约" in R29
      and "DIR-SUPPORT-D" in R29
      and "RECORD-FAMILY-D" in R29
      and "PRODUCT-LEDGER" in R29
      and "SAME-Q" in R29)
check("全支撑没有写成乘积记录",
      "FULL-SUPPORT-LEDGER" in R29
      and "not\\Longrightarrow" in R29
      and "PRODUCT-LEDGER" in R29
      and "不能直接写" in R29)
check("没有把条件结果写成无条件定理",
      "条件因子化" in R29
      and "仍没有从 Zero 原生导出" in R29
      and "不关闭 O3" in R29)
check("没有推出四维 GR",
      "没有由四维标签推出 Lorentz" in R29
      and "Lovelock" in R29
      and "GR" in R29)

head("F2  全支撑—非乘积反例")

for q_num, q_den in [(5, 9), (1, 3), (1, 2), (3, 5), (4, 5)]:
    q = q_num / q_den
    discriminant = 1.0 - q * q
    a_sq = (1.0 + sqrt(discriminant)) / 2.0
    b_sq = (1.0 - sqrt(discriminant)) / 2.0
    overlap = 2.0 * sqrt(a_sq * b_sq)
    check("q=%d/%d 时有正 a,b 且重叠保持 q" % (q_num, q_den),
          a_sq > 0 and b_sq > 0 and abs(overlap - q) < 1e-12,
          "a^2=%.12f b^2=%.12f overlap=%.12f" % (a_sq, b_sq, overlap))

check("两态构造只给张量秩至多 2",
      "联合张量秩至多为 2" in R29
      and "不是 `D` 个独立记录因子的张量积" in R29)
check("文档明确全支撑与乘积分离",
      "\\texttt{DIR-SUPPORT-D}" in R29
      and "\\texttt{PRODUCT-LEDGER}" in R29
      and "A_D\\ne q^D" in R29)

head("F3  共同记录反例与有限峰")

q = Fraction(5, 9)
values = [comb(D, 2) * q for D in range(2, 80)]
check("共同记录模型随 D 单调增长",
      all(values[i + 1] > values[i] for i in range(len(values) - 1)))
check("共同记录模型没有有限峰",
      argmax_d(values) == [79])
check("二项系数比值正确",
      all(Fraction(comb(D + 1, 2), comb(D, 2)) == Fraction(D + 1, D - 1)
          for D in range(2, 40)))
check("文档给出共同记录模型",
      "共同记录合并模型" in R29
      and "B\\binom D2q" in R29
      and "没有有限维峰" in R29)

head("F4  条件因子化与四维峰")

pair_values = [comb(D, 2) * q ** D for D in range(2, 80)]
check("q=5/9 的成对乘积模型唯一峰为 D=4",
      argmax_d(pair_values) == [4],
      "argmax=%s" % argmax_d(pair_values))
check("四维峰只在四项条件同时采用后给出",
      "PHASE-IDENTITY-DER" in R29
      and "LEDGER-FACTORIZATION" in R29
      and "WIPE-RESET-LEDGER" in R29
      and "L=4" in R29
      and "同时成立时" in R29)
check("身份与代价分开记账",
      "M_D\\text{ 管身份" in R29
      and "A_D\\text{ 管代价" in R29
      and "两个独立问题" in R29)

head("F5  旧理论源审计与边界")

check("本仓库源审计覆盖 G71/G72/R25/R26/R27/R28",
      all(name in R29 for name in ["G71", "G72", "R25", "R26", "R27", "R28"]))
check("旧理论源审计覆盖关键反例和边界",
      all(name in R29 for name in ["D224", "D225", "D125", "D127", "D122", "D183", "D194", "D186"]))
check("cosmos-construct 没有被写成能补乘积记录",
      "Q673" in R29
      and "Q549" in R29
      and "不能产生乘积记录结构" in R29)
check("源审计限定于已审计材料",
      "在被审计的材料中" in R29
      and "没有发现" in R29)

head("F6  五项清单与全项目口径")

check("五项清单被重新记账",
      "PHASE-1-COCHAIN" in R29
      and "PAIR-ID-QUOTIENT" in R29
      and "PARALLEL-PERIOD" in R29
      and "L=4" in R29)
check("PARALLEL-PERIOD 不再作为当前必要项",
      "不适用；当前采用 D211／D222 的共同毁灭代际口径" in R29
      and "否" in R29)
check("R28 已指向 R29",
      "R29_full_support_ledger_factorization_no_go.md" in R28
      and "LEDGER-FACTORIZATION" in R28)
check("R27 已指向 R29",
      "R29_full_support_ledger_factorization_no_go.md" in R27
      and "PRODUCT-LEDGER" in R27)
check("R26 成对代价 no-go 仍被引用",
      "重数不蕴含 `q^D` 代价" in R26
      and "PAIR-COST-FACTORIZATION" in R26)
check("STATUS 已登记 R29",
      "### 2.29" in STATUS
      and "R29_full_support_ledger_factorization_no_go.md" in STATUS
      and "DIR-SUPPORT-D" in STATUS
      and "PRODUCT-LEDGER" in STATUS
      and "R29_check.py" in STATUS)
check("INDEX 已收录 R29",
      "R29_full_support_ledger_factorization_no_go.md" in INDEX)

head("F7  诚实边界")

check("R29 不关闭主要上游",
      "没有关闭 `DIM-SECTOR`、`EVO-NORM` 或 `SURV4-GLOBAL`" in R29)
check("R29 不排除未来共同构造",
      "没有排除未来可构造一个共同的新结构一次推出这四项" in R29)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
