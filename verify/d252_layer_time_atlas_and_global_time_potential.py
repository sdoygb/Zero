#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D252 —— 分层时间图册与全局时间势
================================================================================
【被检验的问题（D252）】
  分层局部时间能否拼成全局时间势？
  正环是否给出精确 no-go？
  不可通约周期是否阻止严格全局周期？
  全局时间势能否自动给 Lorentz 光锥？

【本步判据】
  N1  D252 与恢复结构已登记
  N2  无正环约束系统可行
  N3  正环约束系统不可行
  N4  不可通约周期不存在共同精确周期
  N5  同一时间势可配不同空间度规
  N6  不同空间度规给不同光锥速度
  N7  事件序定向不能自动给号差
  N8  文档登记时间势与光锥缺口
  N9  文档不把时间势写成 U1-U4 推论
  N10 上游边界保持
  N11 文档不引用项目外体系名或外部路径

运行：python3 verify/d252_layer_time_atlas_and_global_time_potential.py
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
if os.path.join(ROOT, "verify") not in sys.path:
    sys.path.insert(0, os.path.join(ROOT, "verify"))

from latex_utils import canonical_math


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
            "D252_layer_time_atlas_and_global_time_potential.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D252 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D252", "D252" in AX)
check(
    "N1 公理表与文档登记层时间图册",
    "R-Z-LAYER-TIME-ATLAS" in BOTH,
)
check(
    "N1 公理表与文档登记时间势可行性",
    "R-Z-TIME-POTENTIAL-FEASIBILITY" in BOTH,
)
check(
    "N1 公理表与文档登记光锥缺口",
    "R-Z-TIME-TO-LORENTZ-CONE-GAP" in BOTH,
)
check(
    "N1 沿用层时间子级缺口",
    "R-Z-LAYER-TIME-SUBLEVEL-GAP" in BOTH,
)


# ==================================================================
head("N2  无正环约束系统可行")
# ==================================================================
def constraints_satisfied(time_values, constraints):
    for source, target, lower_bound in constraints:
        if time_values[target] - time_values[source] < lower_bound - 1.0e-12:
            return False
    return True


feasible_times = {"x": 0.0, "y": 2.0, "z": 5.0}
feasible_constraints = [
    ("x", "y", 1.0),
    ("y", "z", 2.0),
    ("x", "z", 3.0),
]

check(
    "N2 候选时间势满足全部约束",
    constraints_satisfied(feasible_times, feasible_constraints),
)
check(
    "N2 候选时间势有足够间隔",
    feasible_times["z"] - feasible_times["x"] >= 3.0,
)


# ==================================================================
head("N3  正环约束系统不可行")
# ==================================================================
positive_cycle = [
    ("x", "y", 1.0),
    ("y", "z", 1.0),
    ("z", "x", 1.0),
]
cycle_sum = sum(bound for _, _, bound in positive_cycle)


def brute_force_feasible(constraints, bound=3.0, step=0.25):
    grid = np.arange(-bound, bound + step, step)
    for tx in grid:
        for ty in grid:
            for tz in grid:
                candidate = {"x": tx, "y": ty, "z": tz}
                if constraints_satisfied(candidate, constraints):
                    return True
    return False


check(
    "N3 正环总下界为正",
    cycle_sum > 0.0,
)
check(
    "N3 正环在有限网格搜索中不可行",
    not brute_force_feasible(positive_cycle),
)


# ==================================================================
head("N4  不可通约周期不存在共同精确周期")
# ==================================================================
period_x = 1.0
period_y = math.sqrt(2.0)
common_period_found = False

for integer_x in range(1, 2001):
    for integer_y in range(1, 2001):
        if integer_x * integer_x == 2 * integer_y * integer_y:
            common_period_found = True
            break
    if common_period_found:
        break

check(
    "N4 周期比为 sqrt(2) 的倒数",
    np.isclose(period_x / period_y, 1.0 / math.sqrt(2.0)),
)
check(
    "N4 整数方程 m^2=2n^2 无解",
    not common_period_found,
)


# ==================================================================
head("N5  同一时间势可配不同空间度规")
# ==================================================================
time_function = lambda t, x: t
same_time_at_event = (
    time_function(2.0, 0.2)
    == time_function(2.0, 3.7)
)
metric_one = np.diag([-1.0, 1.0, 1.0])
metric_two = np.diag([-1.0, 4.0, 4.0])

check(
    "N5 两模型共享时间势",
    same_time_at_event,
)
check(
    "N5 空间度规不同",
    not np.allclose(metric_one[1:, 1:], metric_two[1:, 1:]),
)


# ==================================================================
head("N6  不同空间度规给不同光锥速度")
# ==================================================================
speed_one = 1.0 / math.sqrt(metric_one[1, 1])
speed_two = 1.0 / math.sqrt(metric_two[1, 1])

check(
    "N6 第一光锥速度为一",
    np.isclose(speed_one, 1.0),
)
check(
    "N6 第二光锥速度不同",
    not np.isclose(speed_two, speed_one),
)
check(
    "N6 同一时间势仍兼容不同光锥",
    same_time_at_event and not np.isclose(speed_two, speed_one),
)


# ==================================================================
head("N7  事件序定向不能自动给号差")
# ==================================================================
event_order = [(0, 1), (1, 2), (0, 2)]
reversed_order = [(2, 1), (1, 0), (2, 0)]
undirected_original = {frozenset(edge) for edge in event_order}
undirected_reversed = {frozenset(edge) for edge in reversed_order}

check(
    "N7 时间反演保持无向可比图",
    undirected_original == undirected_reversed,
)
check(
    "N7 时间反演改变定向",
    set(event_order) != set(reversed_order),
)


# ==================================================================
head("N8  文档登记时间势与光锥缺口")
# ==================================================================
check(
    "N8 文档登记时间图册",
    "R-Z-LAYER-TIME-ATLAS" in DOC
    and "局部时间图册" in DOC,
)
check(
    "N8 文档登记正环判据",
    "正环" in DOC
    and "全局时间势存在" in DOC,
)
check(
    "N8 文档登记光锥缺口",
    "R-Z-TIME-TO-LORENTZ-CONE-GAP" in DOC
    and "不同光锥" in DOC,
)


# ==================================================================
head("N9  文档不把时间势写成 U1-U4 推论")
# ==================================================================
check(
    "N9 文档明确条件性",
    "不是 `U1-U4` 的推论" in DOC
    and "条件给出全局时间势" in DOC,
)
check(
    "N9 文档不把时间势写成完整 Lorentz",
    "不是完整 Lorentz 时空" in DOC,
)


# ==================================================================
head("N10  上游边界保持")
# ==================================================================
check(
    "N10 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)
check(
    "N10 文档登记输入预算",
    "## §3 输入预算" in DOC,
)


# ==================================================================
head("N11  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N11 D252 文档不含项目外体系名或外部路径",
    not any(marker in DOC for marker in external_markers),
)


print("\n" + "=" * 78)
print(f"通过 {len(PASS)} / 不符 {len(FAIL)}")
print("=" * 78)
if FAIL:
    for name in FAIL:
        print(f"失败：{name}")
    sys.exit(1)
print("全部通过")
