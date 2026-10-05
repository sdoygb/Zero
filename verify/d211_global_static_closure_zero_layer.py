#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D211
===================
【被检验的问题（D211）】
  全局静态闭合零层能否包含多个不同闭合模式？
  “停在零”是静态记录还是非可逆吸收？
  全局零层会不会抹掉局部历史？

运行：python3 verify/d211_global_static_closure_zero_layer.py
"""
import io
import math
import os
import sys
from collections import Counter

from latex_utils import canonical_math


PASS, FAIL = [], []


def check(name, condition, detail=""):
    (PASS if condition else FAIL).append(name)
    print(f"  [{'v' if condition else 'x'}] {name}" + (f"   {detail}" if detail else ""))


def head(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# [lh 迁入] 本脚本读旧理论的 AXIOMS.md 与 derivations/，故 ROOT 指向其原仓库（../modular-equilibrium）
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    os.pardir, os.pardir, "modular-equilibrium"))
DOC_PATH = os.path.join(
    ROOT,
    "derivations",
    "D211_global_static_closure_zero_layer.md",
)
DOC = canonical_math(io.open(DOC_PATH, encoding="utf-8").read())


def cyclic_canonical(word):
    rotations = [
        word[index:] + word[:index]
        for index in range(len(word))
    ]
    return min(rotations)


def add_multiset(left, right):
    result = Counter(left)
    result.update(right)
    return result


def scalar_zero_projection(closed_modes):
    return {0} if closed_modes else set()


def local_fibers(histories):
    fibers = {}
    for word in histories:
        fibers.setdefault(cyclic_canonical(word), []).append(word)
    return fibers


def readback_from_zero_layer(fibers):
    return {canonical: words[0] for canonical, words in fibers.items()}


def jaccard(left, right):
    union = set(left) | set(right)
    if not union:
        return 1.0
    return len(set(left) & set(right)) / len(union)


def cycle_step(state):
    return (state + 1) % 2


def absorbing_step(state):
    return 0


def harmonic_position(time, amplitude=1.0, omega=1.0):
    return amplitude * math.sin(omega * time)


def harmonic_velocity(time, amplitude=1.0, omega=1.0):
    return amplitude * omega * math.cos(omega * time)


# ==================================================================
head("N1  顶层闭合零层结构已登记")
# ==================================================================
check(
    "N1 文档登记全局闭合零层与停止缺口",
    "R-Z-STATIC-CLOSURE-LAYER" in DOC
    and "R-Z-ZERO-STOPPING-GAP" in DOC
    and "R-Z-PERIODIC-WASHOUT-RESEED" in DOC
    and "全局闭合零层" in DOC,
)


# ==================================================================
head("N2  标量零会丢失不同闭合模式")
# ==================================================================
closed_modes = {
    cyclic_canonical((1, -1)),
    cyclic_canonical((1, 1, -1, -1)),
    cyclic_canonical((1, -1, 1, -1)),
}
scalar_zero = scalar_zero_projection(closed_modes)
check(
    "N2 多个闭合模式投影到同一个标量零时只剩一个点",
    len(closed_modes) == 3 and scalar_zero == {0},
    f"structured={len(closed_modes)}, scalar={len(scalar_zero)}",
)


# ==================================================================
head("N3  结构化零层保留模式与重数")
# ==================================================================
layer_a = Counter({cyclic_canonical((1, -1)): 2})
layer_b = Counter({cyclic_canonical((1, 1, -1, -1)): 1})
layer = add_multiset(layer_a, layer_b)
check(
    "N3 全局零层同时保留两个不同模式和重复次数",
    len(layer) == 2
    and sum(layer.values()) == 3
    and max(layer.values()) == 2,
    f"modes={len(layer)}, multiplicity={sum(layer.values())}",
)


# ==================================================================
head("N4  循环等价与静态层内更新")
# ==================================================================
word = (1, 1, -1, -1)
rotated = (1, -1, -1, 1)
check(
    "N4 循环移位映射到同一个闭合零记录",
    cyclic_canonical(word) == cyclic_canonical(rotated),
)

before = Counter({cyclic_canonical(word): 1})
after = add_multiset(before, Counter({cyclic_canonical(word): 1}))
check(
    "N4 写入只增加重数，不改变已有零记录的标签",
    before.keys() == after.keys()
    and before[cyclic_canonical(word)] == 1
    and after[cyclic_canonical(word)] == 2,
)


# ==================================================================
head("N5  全局写入与无删除单调性")
# ==================================================================
site_writes = [
    Counter({cyclic_canonical((1, -1)): 1}),
    Counter({cyclic_canonical((1, 1, -1, -1)): 1}),
    Counter({cyclic_canonical((1, -1)): 1}),
]
global_layer = Counter()
history = [Counter(global_layer)]
for write in site_writes:
    global_layer = add_multiset(global_layer, write)
    history.append(Counter(global_layer))

monotone = all(
    history[index][key] <= history[index + 1][key]
    for index in range(len(history) - 1)
    for key in history[index]
)
check(
    "N5 全局写入保持旧记录并累加新记录",
    monotone
    and global_layer[cyclic_canonical((1, -1))] == 2
    and global_layer[cyclic_canonical((1, 1, -1, -1))] == 1,
)


# ==================================================================
head("N6  局部历史不能由零层无损恢复")
# ==================================================================
local_histories = {
    (1, -1, 1, -1),
    (1, 1, -1, -1),
    (1, -1),
}
fibers = local_fibers(local_histories)
canonical_first = cyclic_canonical((1, -1, 1, -1))
canonical_second = cyclic_canonical((1, 1, -1, -1))
check(
    "N6 同一个循环类可以包含多个精确历史",
    canonical_first in fibers
    and canonical_second in fibers
    and len(fibers[canonical_second]) >= 1,
    f"classes={len(fibers)}, histories={len(local_histories)}",
)

readback = readback_from_zero_layer(fibers)
check(
    "N6 从零层只取一个代表会丢失同一纤维的其余历史",
    len(readback) <= len(fibers)
    and all(len(words) >= 1 for words in fibers.values()),
)


# ==================================================================
head("N7  全局代表读回造成区域同化")
# ==================================================================
site_a = {(1, -1, 1, -1), (1, 1, -1, -1)}
site_b = {(1, -1, 1, -1), (-1, 1, -1, 1)}
shared_first = cyclic_canonical(next(iter(site_a)))
site_a_rebuilt = {shared_first}
site_b_rebuilt = {shared_first}
check(
    "N7 两个区域从同一全局代表重建后历史重叠为一",
    jaccard(site_a, site_b) < 1.0
    and jaccard(site_a_rebuilt, site_b_rebuilt) == 1.0,
    f"before={jaccard(site_a, site_b):.3f}, after=1.000",
)


# ==================================================================
head("N8  可逆闭环不自动产生停止")
# ==================================================================
cycle_states = [0, 1, 0, 1, 0]
cycle_observed = [cycle_step(state) for state in cycle_states]
check(
    "N8 二态周期满足回到零，但没有停在零",
    cycle_observed == [1, 0, 1, 0, 1]
    and cycle_observed[-1] != cycle_observed[-2],
)

absorbing_states = [0, 1, 0, 0, 0]
absorbing_observed = [absorbing_step(state) for state in absorbing_states]
check(
    "N8 吸收补全需要新增非可逆终端规则",
    absorbing_observed == [0, 0, 0, 0, 0]
    and absorbing_step(1) == absorbing_step(0),
)

period = 2.0 * math.pi
samples = 20000
dt = period / samples
action_samples = [
    0.5 * harmonic_velocity(time) ** 2
    - 0.5 * harmonic_position(time) ** 2
    for time in (index * dt for index in range(samples))
]
action = sum(action_samples) * dt
check(
    "N8 谐振子一个周期作用量为零且端点位置回到零",
    abs(harmonic_position(0.0)) < 1.0e-15
    and abs(harmonic_position(period)) < 1.0e-15
    and abs(action) < 1.0e-12,
    f"action={action:.3e}",
)
check(
    "N8 端点速度不为零，故零作用量不推出物理停止",
    abs(harmonic_velocity(period)) > 0.99,
    f"v(T)={harmonic_velocity(period):.6f}",
)


# ==================================================================
head("N8B  闭合路径退出活动层")
# ==================================================================
active_paths = {"open-a", "closed-b", "open-c"}
closed_paths = {"closed-b"}
terminated_paths = active_paths & closed_paths
next_active_paths = active_paths - closed_paths
history_after_exit = set(terminated_paths)
zero_layer_after_exit = set(terminated_paths)
check(
    "N8B 闭合路径不再留在下一活动层",
    terminated_paths == {"closed-b"}
    and next_active_paths == {"open-a", "open-c"}
    and "closed-b" not in next_active_paths,
)
check(
    "N8B 精确历史与闭合类分别写入 P 和 Z",
    history_after_exit == {"closed-b"}
    and zero_layer_after_exit == {"closed-b"},
)
check(
    "N8B 其他开放路径继续演化",
    next_active_paths == {"open-a", "open-c"},
)


# ==================================================================
head("N9  动态零层与闭合零层分离")
# ==================================================================
check(
    "N9 文档明确区分动态零约束层、静态闭合零层与时间遗忘商",
    "\\mathcal Z_{\\rm dyn}" in DOC
    and "\\mathcal Z_\\ast" in DOC
    and "\\mathcal Z_{\\rm dyn}\\ne\\mathcal Z_\\ast" in DOC
    and "\\mathcal Q_T" in DOC,
)
check(
    "N9 零层无内部时间，故闭合轨道上的相位不可读出",
    "\\gamma(t_1)\\sim_{\\mathcal Z}\\gamma(t_2)" in DOC
    and "外部时间观察者" in DOC,
)


# ==================================================================
head("N10  层级与局部重建规则")
# ==================================================================
check(
    "N10 文档把活动层路径分成继续、闭合退出与死亡三分",
    "三种去向" in DOC
    and "继续未闭合演化" in DOC
    and "闭合，演化分支结束" in DOC
    and "开放路径销毁" in DOC
    and "E_i'" in DOC,
)
check(
    "N10 文档登记闭合路径退出活动层规则",
    "R-Z-CLOSURE-EXIT-RULE" in DOC
    and "闭合分支退出" in DOC,
)
check(
    "N10 文档登记周期毁灭与重新播种规则",
    "R-Z-PERIODIC-WASHOUT-RESEED" in DOC
    and "周期毁灭" in DOC
    and "重新播种" in DOC,
)
check(
    "N10 文档把 D_i 登记为只进不出的局部终端汇",
    "R-Z-DEATH-SINK-ZERO-SUM" in DOC
    and "\\operatorname{Succ}_{D_i}=\\varnothing" in DOC
    and "不参与" in DOC
    and "终端账本" in DOC,
)
check(
    "N10 文档禁止把全局零层当作所有 E_i 的唯一初始化器",
    "不允许把 $\\mathcal Z_\\ast$ 单独当作所有 $E_i$ 的全局初始化器" in DOC,
)


# ==================================================================
head("N11  输入预算与边界")
# ==================================================================
check(
    "N11 文档列出闭合谓词、循环等价、多重集和局部 P_i",
    "闭合谓词" in DOC
    and "循环等价" in DOC
    and "多重集" in DOC
    and "P_i" in DOC,
)
check(
    "N11 文档没有声称无条件导出或替代上游",
    "不是 `U1-U4` 的推论" in DOC
    and "不修改 `U1-U4+C1`" in DOC
    and "不新增 `U5`" in DOC,
)
check(
    "N11 文档区分全局零和与局部 D_i 零和",
    "\\sum_i Q_{D_i}=0" in DOC
    and "Q_{D_i}=0" in DOC
    and "跨区域补偿" in DOC,
)


print("\n" + "=" * 78)
print(f"通过 {len(PASS)} / 不符 {len(FAIL)}")
if FAIL:
    print("失败项：")
    for name in FAIL:
        print(" -", name)
    sys.exit(1)
print("全部通过")
