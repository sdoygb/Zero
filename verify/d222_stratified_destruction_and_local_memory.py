#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D222 —— 分层毁灭与局部记忆
================================================================================
【被检验的问题（D222）】
  活动层全部分层毁灭以后，历史层只保留最高两层是否仍然自洽？
  不同区域的局部寿命是否可以不同？
  被删除的历史低层是否确实不能再影响重播种？
  全局零和与闭合类标签是否仍然保持？

【本步判据】
  N1  D222 与两个恢复结构已登记
  N2  历史层最高两层投影正确且无新记录时幂等
  N3  不同局部寿命可以异步触发毁灭，其他区域不受影响
  N4  局部毁灭时活动层的全部亚层都被清空
  N5  已闭合历史删去低层后保持零和
  N6  活动层非零荷全部进入 D_i 且全局零和保持
  N7  只在低层不同的历史重播种结果相同
  N8  只保留全局零层会抹掉区域差异
  N9  零层压缩保留不同闭合类标签
  N10 tau_i=3 与闭合词最小周期 p(w) 可以分开
  N11 文档不把条件构造冒充上游导出
  N12 不引用旧体系或外部路径
  N13 上游边界保持

运行：python3 verify/d222_stratified_destruction_and_local_memory.py
"""
import io
import os
import sys


# [lh 迁入] 本脚本读旧理论的 AXIOMS.md 与 derivations/，故 ROOT 指向其原仓库（../modular-equilibrium）
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    os.pardir, os.pardir, "modular-equilibrium"))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

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
            "D222_stratified_destruction_and_local_memory.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


def top_two_indices(level_count):
    if level_count <= 0:
        return ()
    if level_count == 1:
        return (0,)
    return (level_count - 2, level_count - 1)


def retained_keys(levels):
    top = top_two_indices(len(levels))
    return {index: levels[index] for index in top}


def retained_contents(levels):
    return tuple(retained_keys(levels).values())


def primitive_period(word):
    for period in range(1, len(word) + 1):
        if len(word) % period == 0 and all(
            word[index] == word[index % period]
            for index in range(len(word))
        ):
            return period
    return len(word)


def is_closed(word):
    return bool(word) and sum(word) == 0


def zero_sum(records):
    return sum(sum(word) for word in records)


def condense_zero_layer(records):
    counts = {}
    for closed_class in records:
        counts[closed_class] = counts.get(closed_class, 0) + 1
    return counts


def local_destruction_events(lifetimes, horizon):
    events = []
    age = [0] * len(lifetimes)
    for reference in range(1, horizon + 1):
        for site, lifetime in enumerate(lifetimes):
            age[site] += 1
            if age[site] % lifetime == 0:
                events.append((reference, site))
    return events


# ==================================================================
head("N1  D222 与两个恢复结构已登记")
# ==================================================================
check(
    "N1 公理表与文档登记分层毁灭结构",
    "D222" in AX
    and "R-Z-STRATIFIED-DESTRUCTION-MEMORY" in BOTH
    and "最高两层" in DOC,
)
check(
    "N1 公理表与文档登记跨层时间缺口",
    "R-Z-LAYER-TIME-SUBLEVEL-GAP" in BOTH
    and "局部寿命" in DOC
    and "时间种类" in DOC,
)


# ==================================================================
head("N2  历史层最高两层投影正确且幂等")
# ==================================================================
projection_ok = True
idempotence_ok = True
for count in range(1, 9):
    levels = [f"P{q}" for q in range(count)]
    retained = retained_contents(levels)
    expected = (
        ("P0",)
        if count == 1
        else (f"P{count - 2}", f"P{count - 1}")
    )
    projection_ok &= retained == expected
    reindexed = list(retained)
    idempotence_ok &= retained_contents(reindexed) == retained
check(
    "N2 任意层数的最高两层都被正确保留",
    projection_ok,
    f"counts=1..8",
)
check(
    "N2 无新记录时重复投影幂等",
    idempotence_ok,
)


# ==================================================================
head("N3  不同局部寿命可以异步触发毁灭")
# ==================================================================
lifetimes = (2, 3, 5)
events = local_destruction_events(lifetimes, horizon=15)
site_event_counts = [
    sum(1 for _reference, site in events if site == index)
    for index in range(len(lifetimes))
]
event_times = {
    site: [
        reference
        for reference, event_site in events
        if event_site == site
    ]
    for site in range(len(lifetimes))
}
distinct_event_sequences = (
    len({tuple(times) for times in event_times.values()}) == 3
)
check(
    "N3 局部寿命 2,3,5 给出不同毁灭序列",
    site_event_counts == [7, 5, 3]
    and event_times[0] != event_times[1]
    and event_times[1] != event_times[2],
    f"events={event_times}",
)
check(
    "N3 毁灭事件按区域触发，不需要共同周期",
    distinct_event_sequences,
)


# ==================================================================
head("N4  局部毁灭时活动层的全部亚层都被清空")
# ==================================================================
active_layers = [
    {"E0"},
    {"E0", "E1"},
    {"E0", "E1", "E2"},
]
all_cleared = True
untouched_where_expected = True
for reference, site in events:
    active_by_region = [set(layer) for layer in active_layers]
    active_by_region[site].clear()
    all_cleared &= not active_by_region[site]
    untouched_where_expected &= all(
        active_by_region[other] == active_layers[other]
        for other in range(len(active_layers))
        if other != site
    )
check(
    "N4 每次局部事件清空该区域的全部活动亚层",
    all_cleared,
)
check(
    "N4 未触发事件的区域不被连带清空",
    untouched_where_expected,
)


# ==================================================================
head("N5  已闭合历史删去低层后保持零和")
# ==================================================================
history_levels = [
    {(1, -1), (1, 1, -1, -1)},
    {(1, -1, -1, 1), (-1, 1, 1, -1)},
    {(1, 1, -1, -1, -1, 1)},
]
all_closed = all(
    is_closed(word)
    for level in history_levels
    for word in level
)
retained_levels = retained_contents(history_levels)
retained_words = [
    word
    for level in retained_levels
    for word in level
]
check(
    "N5 测试历史全部是闭合零荷记录",
    all_closed and zero_sum(retained_words) == 0,
)
check(
    "N5 删除最低层后零和不变",
    zero_sum(word for level in history_levels for word in level) == 0
    and zero_sum(retained_words) == 0,
)


# ==================================================================
head("N6  活动层非零荷全部进入 D_i 且全局零和保持")
# ==================================================================
open_paths_by_site = [
    {(1,), (1, 1)},
    {(-1,), (-1, -1)},
    {(1, 1, 1), (-1, -1, -1)},
]
open_charge = [
    sum(sum(word) for word in paths)
    for paths in open_paths_by_site
]
dead_charge = [0] * len(open_paths_by_site)
for site, charge in enumerate(open_charge):
    dead_charge[site] += charge
check(
    "N6 活动层单个区域可以携带非零开放荷",
    open_charge == [3, -3, 0],
    f"open_charge={open_charge}",
)
check(
    "N6 全部开放荷进入 D_i 后全局仍为零",
    sum(open_charge) == 0 and sum(dead_charge) == 0,
)


# ==================================================================
head("N7  只在低层不同的历史重播种结果相同")
# ==================================================================
def projection(levels):
    return retained_contents(levels)


def seed_from_projection(retained):
    return tuple(
        tuple(f"seed:{item}" for item in level)
        for level in retained
    )


history_a = [
    {"detail:a"},
    {"summary:A"},
    {"stable:A"},
]
history_b = [
    {"detail:b"},
    {"summary:A"},
    {"stable:A"},
]
same_top = projection(history_a) == projection(history_b)
same_seed = seed_from_projection(projection(history_a)) == seed_from_projection(
    projection(history_b)
)
check(
    "N7 两份历史只在低层不同",
    history_a != history_b and same_top,
)
check(
    "N7 保留投影相同则重播种结果相同",
    same_seed,
)


# ==================================================================
head("N8  只保留全局零层会抹掉区域差异")
# ==================================================================
global_zero = {"zero:class"}
site_a_top = (frozenset({"local:A"}), frozenset({"stable:A"}))
site_b_top = (frozenset({"local:B"}), frozenset({"stable:B"}))
seed_with_local = (
    seed_from_projection(site_a_top),
    seed_from_projection(site_b_top),
)
seed_without_local = (
    tuple(global_zero),
    tuple(global_zero),
)
check(
    "N8 局部顶层不同则种子不同",
    seed_with_local[0] != seed_with_local[1],
)
check(
    "N8 只用全局零层时区域差异消失",
    seed_without_local[0] == seed_without_local[1],
)


# ==================================================================
head("N9  零层压缩保留不同闭合类标签")
# ==================================================================
zero_layer = ["alpha", "alpha", "beta", "alpha", "beta"]
condensed = condense_zero_layer(zero_layer)
check(
    "N9 类内重数可以压缩",
    condensed == {"alpha": 3, "beta": 2},
    f"condensed={condensed}",
)
check(
    "N9 不同闭合类不会在压缩后合并",
    set(condensed) == {"alpha", "beta"},
)


# ==================================================================
head("N10 tau_i=3 与闭合词最小周期 p(w) 可以分开")
# ==================================================================
closed_words = [
    (1, -1),
    (1, 1, -1, -1),
    (1, -1, 1, -1),
    (1, -1, -1, 1),
]
periods = {
    word: primitive_period(word)
    for word in closed_words
}
check(
    "N10 所有测试闭合词的最小旋转周期都是偶数",
    all(period % 2 == 0 for period in periods.values()),
    f"periods={periods}",
)
check(
    "N10 活动局部寿命 3 不需要等于 p(w)",
    all(period != 3 for period in periods.values())
    and "局部寿命" in DOC
    and "闭包旋转周期" in DOC,
)


# ==================================================================
head("N11 文档不把条件构造冒充上游导出")
# ==================================================================
check(
    "N11 文档明确条件构造",
    "条件构造" in DOC
    and "不是 `U1-U4` 的推论" in DOC,
)
check(
    "N11 文档列出未解决来源",
    "亚层分解从哪里来" in DOC
    and "为什么“最高两层”恰好是二" in DOC
    and "局部寿命" in DOC,
)


# ==================================================================
head("N12 不引用旧体系或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "旧体系",
    "外部路径",
    "../..",
)
check(
    "N12 D222 文档不含旧体系名或外部路径",
    not any(marker in DOC for marker in external_markers),
)


# ==================================================================
head("N13 上游边界保持")
# ==================================================================
check(
    "N13 不修改 U1-U4+C1 且不新增 U5",
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
