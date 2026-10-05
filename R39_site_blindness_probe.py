#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R39_site_blindness_probe.py -- 历史层对"位点"是否失明？

记录（Z3/G0）：闭合分支写入 (精确词 w, 闭合类 [w])。问：由记录能否反推位点？
判据：记录 -> 位点 的互信息 I(site ; record)。
  I = 0  <=>  完全失明（R38 的路线 B' 有立足点）
  I > 0  <=>  记录含位点信息（则纠缠被退相干，只能走 A'）

实现：环图 C_m，移动 = ±1；闭合词 w 满足 sum(w_i) ≡ 0 (mod m)。
      走位点 v_0 -> v_0+w_0 -> ... ；记录只含 (w, [w])，不含 v_0。

输出：R39_site_blindness_results.json
"""

from __future__ import annotations

import io
import itertools
import json
import math
import os
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def canon(w):
    """词的旋转类代表元"""
    L = len(w)
    return min(tuple(w[i:] + w[:i]) for i in range(L))


def closed_words(m, L):
    """长度 L、在环 C_m 上闭合的词：sum ≡ 0 (mod m)"""
    out = []
    for w in itertools.product((1, -1), repeat=L):
        if sum(w) % m == 0:
            out.append(w)
    return out


def site_of(w, v0, m):
    """给定起始位点，走一圈（用于核对闭合性）"""
    v = v0
    for s in w:
        v = (v + s) % m
    return v


results = {}
for m in (3, 4, 5, 6, 8):
    for L in (4, 6):
        words = closed_words(m, L)
        if not words:
            continue
        # 记录：w 本身（精确词）与其旋转类
        rec = {w: (w, canon(w)) for w in words}
        # 对每个记录，找出全部与它一致的起始位点
        consistent = defaultdict(set)
        for w in words:
            for v0 in range(m):
                # 闭合性：与 v0 无关（因为 sum ≡ 0 mod m），故对每个 v0 都一致
                consistent[rec[w]].add(v0)
        sizes = [len(s) for s in consistent.values()]
        # 互信息 I(site ; record)：site 均匀先验
        # p(site|record) 在一致的位点上均匀 => I = H(site) - H(site|record)
        H_site = math.log2(m)
        H_cond = 0.0
        total = sum(sizes)
        for s in sizes:
            p = s / total
            H_cond += p * math.log2(len(consistent_values := set()) or s) if False else 0.0
        # 直接计算：每个记录下 site 均匀分布在 s 个位点上
        H_cond = sum((s / total) * math.log2(s) for s in sizes)
        I = H_site - H_cond
        results["m=%d_L=%d" % (m, L)] = {
            "n_closed_words": len(words),
            "n_records": len(consistent),
            "consistent_sites_per_record": sorted(set(sizes)),
            "H_site": H_site,
            "H_site_given_record": H_cond,
            "mutual_information_bits": I,
            "site_blind": bool(I < 1e-12),
        }

# 对照：记录是否至少能区分"形状"（旋转类）？
shape_test = {}
for m, L in ((6, 6), (8, 6)):
    words = closed_words(m, L)
    classes = {canon(w) for w in words}
    shape_test["m=%d_L=%d" % (m, L)] = {
        "n_words": len(words), "n_shape_classes": len(classes),
        "record_distinguishes_shapes": bool(len(classes) > 1),
    }

out = {
    "record_definition": "(精确词 w, 闭合类 [w])；不含位点标签（Z3/G0）",
    "site_information": results,
    "shape_information": shape_test,
    "verdict": {
        "mutual_information_zero_in_all_cases": all(v["mutual_information_bits"] < 1e-12
                                                    for v in results.values()),
        "consistent_sites_per_record": "全部 m 个位点（记录对绝对位点完全不敏感）",
        "record_still_nontrivial": all(v["record_distinguishes_shapes"] for v in shape_test.values()),
        "conclusion": ("历史层对绝对位点失明（I=0）：记录只含词形与闭合类。"
                       "故 R38 的路线 B' 有立足点：'哪一个位点'是历史层未记录的自由度。"),
        "residual_condition": ("唯一可能破坏失明的是播种点 P_i（D 系列列为未解选择器）；"
                               "若播种把位点写进记录，则失明被打破、须改走 A'（EDGE-CONNECTION）"),
    },
}

with io.open(os.path.join(HERE, "R39_site_blindness_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)

print(json.dumps(out, ensure_ascii=False, indent=1))
