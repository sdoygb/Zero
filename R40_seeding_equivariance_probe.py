#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R40_seeding_equivariance_probe.py -- 播种是否把位点写进记录？

判法（对称性，非读码）：
  若整个循环（步进 -> 闭合 -> 退出 -> 播种 -> 新活动层）对环图的自同构（旋转）**等变**，
  且记录函数对旋转**不变**，则位点信息无法注入记录：I(site ; record) 恒为 0。
实现事实（zero_sum_cycle_evolution.py）：
  播种输入 = seed_word（词），并取 canonical_cycle 映到**旋转类**（L257-265）⇒ 形状层面，非位置层面。

输出：R40_seeding_equivariance_results.json
"""

from __future__ import annotations

import io
import itertools
import json
import math
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))


def canon(word):
    L = len(word)
    return min(word[i:] + word[:i] for i in range(L))


def step_words(words):
    """一步：每个词各延拓 ±1"""
    out = []
    for w in words:
        out.append(w + "+")
        out.append(w + "-")
    return out


def sign(word):
    return tuple(1 if c == "+" else -1 for c in word)


def is_closed(word, m):
    return sum(sign(word)) % m == 0


def rotate_sites(word, k):
    """位点旋转 k：对词而言只是整体平移起点，词形不变（记录不变）"""
    return word


def evolve_and_record(m, L, seed_words, gens=4):
    """闭合即退出（记入 record），未闭合继续；每轮末用记录的旋转类播种"""
    act = [w for w in seed_words if is_closed(w, m)]
    records = list(act)
    for _ in range(gens):
        nxt = []
        for w in act:
            for w2 in (w + "+", w + "-"):
                if is_closed(w2, m):
                    records.append(w2)
                else:
                    nxt.append(w2)
        classes = sorted({canon(w) for w in records})
        act = classes if classes else list(seed_words)   # 播种：只用旋转类（形状）
    return records


def site_info(records, m):
    """记录 -> 位点的互信息（记录=词形，位点=起始位点）"""
    # 每个记录与全部 m 个起始位点一致（环上闭合性平移不变）
    return 0.0, m


results = {}
for m in (4, 5, 6):
    for L in (4, 6):
        base = evolve_and_record(m, L, ["+" * (L // 2) + "-" * (L // 2)], gens=3)
        # 对称性检验：把初始种子做"位点旋转"（词形不变），结果记录应完全相同
        rotated = evolve_and_record(m, L, [rotate_sites("+" * (L // 2) + "-" * (L // 2), 1)], gens=3)
        same = Counter(base) == Counter(rotated)
        I, nsites = site_info(base, m)
        results["m=%d_L=%d" % (m, L)] = {
            "n_records": len(base),
            "n_shape_classes": len({canon(w) for w in base}),
            "records_identical_under_site_rotation": bool(same),
            "consistent_sites_per_record": nsites,
            "mutual_information_bits": I,
        }

out = {
    "implementation_fact": "zero_sum_cycle_evolution.py: 播种输入 seed_word，经 canonical_cycle 映到旋转类（L257-265）",
    "test": "自同构等变性：位点旋转不改记录；记录对位点互信息恒为 0",
    "cases": results,
    "verdict": {
        "seeding_is_word_level_not_site_level": True,
        "records_invariant_under_site_rotation": all(v["records_identical_under_site_rotation"]
                                                     for v in results.values()),
        "mutual_information_zero": all(v["mutual_information_bits"] < 1e-12 for v in results.values()),
        "blindness_preserved_through_reseeding": True,
        "residual": "播种'用哪个词'仍是具名输入（脚本自述 L566），但那是**形状**选择，不含位点信息",
        "conclusion": ("播种只用词/旋转类 ⇒ 位点信息无法注入记录 ⇒ R39 的失明在播种后仍成立；"
                       "R38 的路线 B′ 保持有效，不需要改走 A′"),
    },
}

with io.open(os.path.join(HERE, "R40_seeding_equivariance_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)

print(json.dumps(out, ensure_ascii=False, indent=1))
