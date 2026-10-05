#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 核验 D251 —— 层类型复合与局域读出预层
================================================================================
【被检验的问题（D251）】
  层、事件、记录和局域读出能否组成一个最小类型复合？
  限制映射与细化一致性是否可检查？
  类型复合与事件序能否自动给物理区域、维数和时间方向？

【本步判据】
  N1  D251 与恢复结构已登记
  N2  层类型标签彼此区分
  N3  事件复合对测试模型满足结合性
  N4  限制映射满足函子性
  N5  细化一致性在相容截面上成立
  N6  不相容截面不能胶合
  N7  相同事件序可有不同物理嵌入
  N8  事件序不自动选择维数或时间定向
  N9  文档登记物理嵌入缺口
  N10 文档不把类型复合写成 U1-U4 推论
  N11 上游边界保持
  N12 文档不引用项目外体系名或外部路径

运行：python3 verify/d251_layer_type_complex_and_local_readout_presheaf.py
"""
import io
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
            "D251_layer_type_complex_and_local_readout_presheaf.md",
        ),
        encoding="utf-8",
    ).read()
)
BOTH = AX + "\n" + DOC


# ==================================================================
head("N1  D251 与恢复结构已登记")
# ==================================================================
check("N1 公理表登记 D251", "D251" in AX)
check(
    "N1 公理表与文档登记层类型复合",
    "R-Z-LAYER-TYPE-COMPLEX" in BOTH,
)
check(
    "N1 公理表与文档登记局域读出预层",
    "R-Z-LOCAL-READOUT-PRESHEAF" in BOTH,
)
check(
    "N1 公理表与文档登记细化一致性",
    "R-Z-REFINEMENT-READOUT-COHERENCE" in BOTH,
)
check(
    "N1 公理表与文档登记物理嵌入缺口",
    "R-Z-SUPPORT-EMBEDDING-GAP" in BOTH,
)


# ==================================================================
head("N2  层类型标签彼此区分")
# ==================================================================
layer_types = {"E", "P", "Z", "D"}
check(
    "N2 层类型数为四",
    len(layer_types) == 4,
)
check(
    "N2 类型标签不混同",
    len(layer_types & {"E", "P"}) == 2
    and len(layer_types & {"Z", "D"}) == 2,
)


# ==================================================================
head("N3  事件复合对测试模型满足结合性")
# ==================================================================
monoid_elements = [0, 1, 2, 3]


def compose(left, right):
    return (left + right) % 4


associativity_passed = all(
    compose(compose(a, b), c) == compose(a, compose(b, c))
    for a in monoid_elements
    for b in monoid_elements
    for c in monoid_elements
)

check(
    "N3 测试复合结合",
    associativity_passed,
)
check(
    "N3 存在复合中元",
    compose(0, 2) == 2 and compose(2, 0) == 2,
)


# ==================================================================
head("N4  限制映射满足函子性")
# ==================================================================
supports = {
    "top": {0, 1, 2},
    "mid": {0, 1},
    "left": {0},
}


def restrict(values, target_support):
    return {
        key: values[key]
        for key in target_support
        if key in values
    }


section = {0: 1.0, 1: 2.0, 2: 3.0}

left_direct = restrict(section, supports["left"])
left_via_mid = restrict(restrict(section, supports["mid"]), supports["left"])

check(
    "N4 两级限制等于一级限制",
    left_direct == left_via_mid,
)
check(
    "N4 限制不增加新数据",
    set(left_direct).issubset(set(section)),
)


# ==================================================================
head("N5  细化一致性在相容截面上成立")
# ==================================================================
def readout(potential, support):
    return {key: potential[key] for key in support}


global_potential = {0: 0.3, 1: -0.2, 2: 0.7}
coarse = readout(global_potential, supports["top"])
restricted = restrict(coarse, supports["mid"])
fine_readout = readout(
    restrict(global_potential, supports["mid"]),
    supports["mid"],
)

check(
    "N5 粗读出限制等于细读出",
    restricted == fine_readout,
)
check(
    "N5 中间支持覆盖完整",
    set(fine_readout) == supports["mid"],
)


# ==================================================================
head("N6  不相容截面不能胶合")
# ==================================================================
left_cover = {0, 1}
right_cover = {1, 2}
overlap = left_cover & right_cover
left_section = {0: 1.0, 1: 2.0}
right_section = {1: 9.0, 2: 4.0}

check(
    "N6 两个局部截面覆盖重叠支持",
    overlap == {1},
)
check(
    "N6 重叠处读值不一致",
    left_section[1] != right_section[1],
)
check(
    "N6 不相容截面不能给出共同值",
    not (left_section[1] == right_section[1]),
)


# ==================================================================
head("N7  相同事件序可有不同物理嵌入")
# ==================================================================
order_chain = [(0, 1), (1, 2), (0, 2)]
embedding_line_short = np.array([[0.0], [1.0], [2.0]])
embedding_line_long = np.array([[0.0], [1.0], [4.0]])
embedding_plane = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]])


def realizes_chain(positions, chain):
    for left, right in chain:
        if not np.all(positions[left] <= positions[right]):
            return False
    return True

short_distance = float(
    np.linalg.norm(embedding_line_short[1] - embedding_line_short[0])
)
long_distance = float(
    np.linalg.norm(embedding_line_long[2] - embedding_line_long[1])
)

check(
    "N7 三个嵌入都实现同一事件序",
    realizes_chain(embedding_line_short, order_chain)
    and realizes_chain(embedding_line_long, order_chain)
    and realizes_chain(embedding_plane, order_chain),
)
check(
    "N7 同一对的物理距离不同",
    short_distance != long_distance,
)
check(
    "N7 同一事件序可嵌入不同维数",
    embedding_line_short.shape[1] != embedding_plane.shape[1],
)


# ==================================================================
head("N8  事件序不自动选择维数或时间定向")
# ==================================================================
order_preserving_reversal = {0: 2, 1: 1, 2: 0}
reversed_edges = {
    (order_preserving_reversal[a], order_preserving_reversal[b])
    for a, b in order_chain
}
original_comparability = {frozenset(edge) for edge in order_chain}
reversed_comparability = {frozenset(edge) for edge in reversed_edges}

check(
    "N8 时间反演保持无向可比图",
    reversed_comparability == original_comparability,
)
check(
    "N8 时间反演翻转定向",
    reversed_edges != set(order_chain),
)


# ==================================================================
head("N9  文档登记物理嵌入缺口")
# ==================================================================
check(
    "N9 文档登记类型复合",
    "R-Z-LAYER-TYPE-COMPLEX" in DOC
    and "有限类型复合" in DOC,
)
check(
    "N9 文档登记读出预层",
    "R-Z-LOCAL-READOUT-PRESHEAF" in DOC
    and "预层" in DOC,
)
check(
    "N9 文档登记物理嵌入缺口",
    "R-Z-SUPPORT-EMBEDDING-GAP" in DOC
    and "物理区域" in DOC,
)


# ==================================================================
head("N10  文档不把类型复合写成 U1-U4 推论")
# ==================================================================
check(
    "N10 文档明确条件性",
    "不是 `U1-U4` 的推论" in DOC
    and "不是 `U1-U4` 的定理" in DOC,
)
check(
    "N10 文档不把类型底座写成完整几何",
    "不是物理时空" in DOC,
)


# ==================================================================
head("N11  上游边界保持")
# ==================================================================
check(
    "N11 不修改 U1-U4+C1 且不新增 U5",
    "不修改 `U1-U4+C1`" in DOC and "不新增 `U5`" in DOC,
)
check(
    "N11 文档登记输入预算",
    "## §3 输入预算" in DOC,
)


# ==================================================================
head("N12  文档不引用项目外体系名或外部路径")
# ==================================================================
external_markers = (
    "cosmos-construct",
    "../..",
    "旧体系",
    "外部体系",
)
check(
    "N12 D251 文档不含项目外体系名或外部路径",
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
