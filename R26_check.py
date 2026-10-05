#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R26_check.py -- 核验成对载体约化的图论 no-go、维数字典 no-go 与代价因子化。

对应 R26_pair_carrier_reduction_no_go.md。
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


def argmax(values: dict[int, Fraction]) -> list[int]:
    maximum = max(values.values())
    return [key for key, value in values.items() if value == maximum]


def alpha_pair(q: Fraction, dmax: int = 40) -> dict[int, Fraction]:
    return {D: comb(D + 1, 2) * q**D for D in range(1, dmax + 1)}


def beta_pair(q: Fraction, dmax: int = 40) -> dict[int, Fraction]:
    return {D: comb(D, 2) * q**D for D in range(2, dmax + 1)}


def support_cost(q: Fraction, dmax: int = 40) -> dict[int, Fraction]:
    return {D: comb(D, 2) * q**2 for D in range(2, dmax + 1)}


def edge_cost(q: Fraction, dmax: int = 40) -> dict[int, Fraction]:
    return {D: comb(D, 2) * q ** comb(D, 2) for D in range(2, dmax + 1)}


def replica(q: Fraction, c: Fraction, dmax: int = 80) -> dict[int, Fraction]:
    return {D: comb(D, 2) * c ** comb(D, 2) * q**D for D in range(2, dmax + 1)}


R26 = read("R26_pair_carrier_reduction_no_go.md")
R25 = read("R25_native_pair_cost_and_four_dim_peak.md")
R24 = read("R24_global_four_survival_gate.md")
R23 = read("R23_dimension_descendant_selection.md")
R3 = read("R3_dimension_selection.md")
R0 = read("R0_publication_theorem.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")

head("F1  文档范围与结论")

check("标题与核心判定在位",
      "成对载体约化" in R26
      and "C(D,2)" in R26
      and "不能由连通性直接得到" in R26)
check("没有把条件计数写成 Zero 原生推导",
      "条件计数定理，不是从 Z0／Z 条款的推导" in R26)
check("没有把四维峰等同于四维 GR",
      "没有由四维峰推出 Lorentz" in R26 and "Lovelock" in R26 and "GR" in R26)
check("R26 明确 PAIR-CARRIER-DER 已被分解",
      "不是一个原子命题" in R26
      and "六个可分别审计的接口" in R26
      and "PAIR-COST-FACTORIZATION" in R26)

head("F2  连通图与边轨道计数")

for D in range(3, 31):
    lower = D - 1
    upper = comb(D, 2)
    check("D=%d 时连通图边数界成立" % D,
          lower <= D <= upper,
          "cycle=%d, tree=%d, K=%d" % (D, lower, upper))

check("D=4 的环图只给 4 条边而不是 6 条",
      4 < comb(4, 2) and "|E(C_D)|=D" in R26)
check("身份数等于边数的定理在位",
      "定理 R26.1（身份数等于边数）" in R26
      and "N_{\\rm id}(\\Gamma_D)=|E_D|" in R26)
check("连通性 no-go 与 G89 合法环图模型相连",
      "G89" in R26 and "命题 1" in R26 and "环图 `C_m`" in R26)
check("全局读出与局部继承图分开",
      "Z-E1" in R26
      and "D243" in R26
      and "全局可达性不是局部方向对" in R26)

head("F3  带条件的 C(D,2) 计数")

for D in range(2, 25):
    check("D=%d 的 K_D 边数为 C(D,2)" % D,
          comb(D, 2) == D * (D - 1) // 2,
          str(comb(D, 2)))

check("四项条件计数定理在位",
      "定理 R26.4" in R26
      and "DIR-DICT-BETA" in R26
      and "PAIR-GRAPH-KD" in R26
      and "PAIR-ID-EDGE" in R26
      and "PAIR-COUNT-1" in R26)

head("F4  路线 alpha 的字典 no-go")

q = Fraction(5, 9)
alpha = alpha_pair(q)
check("q=5/9 时 C(D+1,2)q^D 的唯一峰为 D=3",
      argmax(alpha) == [3],
      str(argmax(alpha)))
check("alpha 模型在 D=4 对 D=3 的比值为 25/27",
      alpha[4] / alpha[3] == Fraction(25, 27),
      str(alpha[4] / alpha[3]))
check("alpha 四维窗口为 3/5<q<2/3",
      argmax(alpha_pair(Fraction(3, 5))) == [3, 4]
      and argmax(alpha_pair(Fraction(2, 3))) == [4, 5]
      and "\\frac35<q<\\frac23" in R26)
check("D194 的标签骨架与字典冲突在位",
      "D194" in R26
      and "C(m,2)" in R26
      and "D=m-1" in R26
      and "唯一峰是三维" in R26)
check("beta 模型仍给 q=5/9 的四维峰",
      argmax(beta_pair(q)) == [4],
      str(argmax(beta_pair(q))))

head("F5  成对重数不等于代价")

support = support_cost(q)
check("只付两维支撑代价时序列严格增长",
      all(support[D + 1] > support[D] for D in range(2, 40)))
edge_ratio = {
    D: edge_cost(q)[D + 1] / edge_cost(q)[D]
    for D in range(2, 39)
}
check("每条边独立付代价时比值从 D=2 起小于 1",
      edge_ratio[2] == Fraction(25, 27),
      str(edge_ratio[2]))
check("每条边独立付代价时比值继续下降",
      all(edge_ratio[D + 1] < edge_ratio[D] for D in range(2, 38)))
check("每条边独立付代价时唯一峰为 D=2",
      argmax(edge_cost(q)) == [2],
      str(argmax(edge_cost(q))))
check("q^D 代价必须单独登记为开放输入",
      "PAIR-COST-FACTORIZATION" in R26
      and "开放具名输入" in R26
      and "重数不蕴含 `q^D` 代价" in R26)

rep_bad = replica(q, Fraction(11, 10))
check("按边额外复制 c>1 时尾部重新增长",
      rep_bad[80] > rep_bad[79] > rep_bad[78],
      "D=78<79<80 的比值重新大于 1")
check("额外复制必须单独登记为开放输入",
      "PAIR-NO-EXTRA-MULT" in R26 and "不存在有限全局峰" in R26)

head("F6  旧理论与当前材料审计")

for marker in ["D194", "D244", "D25", "D26", "G80", "G85", "G86", "D243"]:
    check("已审计材料 %s 在位" % marker, marker in R26)
check("D194 只作组合骨架",
      "组合骨架" in R26 and "不是“每个标签自动产生一个可继承后代”的动力学法则" in R26)
check("当前结论没有冒充一般不可能性",
      "没有证明任何形式的成对载体都不存在" in R26)

head("F7  全项目口径同步")

check("R25 已指向 R26 并登记代价因子化",
      "R26_pair_carrier_reduction_no_go.md" in R25
      and "PAIR-COST-FACTORIZATION" in R25
      and "PAIR-NO-EXTRA-MULT" in R25)
check("R25 不再把 PAIR-CARRIER-DER 写成唯一未拆分桥",
      "当前唯一关键物理桥" not in R25
      and "不是一个原子命题" in R25)
check("R24 已指向 R26",
      "R26_pair_carrier_reduction_no_go.md" in R24
      and "PAIR-COST-FACTORIZATION" in R24)
check("R23 已指向 R26",
      "R26_pair_carrier_reduction_no_go.md" in R23)
check("R3 已指向 R26",
      "R26_pair_carrier_reduction_no_go.md" in R3)
check("R0 已登记 R26 且保持 O3 未关闭",
      "R26" in R0
      and "PAIR-COST-FACTORIZATION" in R0
      and "O3 仍未关闭" in R0)
check("STATUS 新增 R26 当前状态",
      "### 2.26 成对载体约化与两条 no-go（`R26`）当前状态" in STATUS
      and "PAIR-GRAPH-KD" in STATUS
      and "PAIR-COST-FACTORIZATION" in STATUS
      and "q=5/9" in STATUS
      and "R26_check.py" in STATUS)
check("INDEX 已收录 R26",
      "R26_pair_carrier_reduction_no_go.md" in INDEX)

head("F8  反过度主张")

check("R26 不关闭 SURV4-GLOBAL",
      "没有关闭 `DIM-SECTOR`、`EVO-NORM` 或 `SURV4-GLOBAL`" in R26)
check("R26 不把条件四维峰升级成 GR",
      "R25 的正结果是条件模型内的精确结果" in R26
      and "四维 GR" in R26)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
