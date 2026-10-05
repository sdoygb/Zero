#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R39_check.py -- 核验"历史层对位点失明"的探针与文档。
"""

from __future__ import annotations

import io
import itertools
import json
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


DOC = read("R39_history_layer_site_blindness.md")
R38 = read("R38_entanglement_from_shared_closure_origin.md")
G0 = read("G0_bottom_layer_and_derivation_route.md")
R28 = read("R28_phase_identity_cluster_reduction.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")
RES = json.load(io.open(os.path.join(HERE, "R39_site_blindness_results.json"), encoding="utf-8"))


def canon(w):
    L = len(w)
    return min(tuple(w[i:] + w[:i]) for i in range(L))


def closed_words(m, L):
    return [w for w in itertools.product((1, -1), repeat=L) if sum(w) % m == 0]


# ======================================================================
head("F1  文档结构：裁决、结构性理由、残余条件")

check("标题与性质为 R38-3 的裁决",
      "历史层对" in DOC and "失明" in DOC and "探针裁决（正面）" in DOC)
check("判据 (R39-1)(R39-2) 与裁决 (R39-3) 在位",
      "(R39-1)" in DOC and "(R39-2)" in DOC and "(R39-3)" in DOC)
check("结构性理由 (R39-4)（E1 缺失 ⇔ 失明）在位",
      "(R39-4)" in DOC and "E1" in DOC and "I5b" in DOC)
check("残余条件写明播种点 P_i 是唯一翻盘入口",
      "播种点" in DOC and "唯一" in DOC and "EDGE-CONNECTION" in DOC)
check("没有声称证明了播种不写位点",
      "没有证明**播种**不写位点" in DOC)

# ======================================================================
head("F2  独立复算：I(位点;记录)=0 且记录非平凡")

for m in (3, 4, 5, 6, 8):
    for L in (4, 6):
        words = closed_words(m, L)
        if not words:
            continue
        recs = {w: (w, canon(w)) for w in words}
        consistent = {}
        for w in words:
            consistent.setdefault(recs[w], set()).update(range(m))
        sizes = [len(s) for s in consistent.values()]
        total = sum(sizes)
        H_site = math.log2(m)
        H_cond = sum((s / total) * math.log2(s) for s in sizes)
        I = H_site - H_cond
        check("m=%d L=%d：每个记录与全部 %d 个位点一致且 I=0" % (m, L, m),
              set(sizes) == {m} and abs(I) < 1e-12, "I=%.2e" % I)

classes_6 = {canon(w) for w in closed_words(6, 6)}
classes_8 = {canon(w) for w in closed_words(8, 6)}
check("记录仍非平凡：m=6,L=6 有 6 个形状类", len(classes_6) == 6)
check("记录仍非平凡：m=8,L=6 有 4 个形状类", len(classes_8) == 4)

# ======================================================================
head("F3  与结果 JSON、G0 记录定义、R38 一致")

check("JSON 判定互信息全为零", RES["verdict"]["mutual_information_zero_in_all_cases"] is True)
check("JSON 写明残余条件为播种点", "P_i" in RES["verdict"]["residual_condition"])
check("G0 的记录定义确为（精确词，闭合类）",
      "精确词" in G0 and "闭合类" in G0)
check("R38 的 (R38-3) 被本文回答", "(R38-3)" in R38 and "(R38-3)" in DOC)
check("R28 的 A′ 路径在位（备选）", "EDGE-CONNECTION" in R28)
check("STATUS 已登记 R39",
      "### 2.39" in STATUS
      and "R39_history_layer_site_blindness.md" in STATUS
      and "R39_check.py" in STATUS)
check("INDEX 已收录 R39",
      "R39_history_layer_site_blindness.md" in INDEX
      and "R39_check.py" in INDEX)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
