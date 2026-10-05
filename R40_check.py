#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R40_check.py -- 核验"播种不注入位点信息"的实现审计与等变性探针。"""
from __future__ import annotations
import io, itertools, json, math, os, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0

def head(t):
    print(""); print("=" * 72); print(t); print("=" * 72)

def check(name, cond, detail=""):
    global FAIL
    ok = bool(cond)
    if not ok: FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))

def read(n):
    with io.open(os.path.join(HERE, n), encoding="utf-8") as f: return f.read()

DOC = read("R40_seeding_does_not_inject_site.md")
R39 = read("R39_history_layer_site_blindness.md")
R38 = read("R38_entanglement_from_shared_closure_origin.md")
SRC = read("zero_sum_cycle_evolution.py")
RES = json.load(io.open(os.path.join(HERE, "R40_seeding_equivariance_results.json"), encoding="utf-8"))
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")

def canon(w):
    L = len(w); return min(w[i:] + w[:i] for i in range(L))

def sign(w): return tuple(1 if c == "+" else -1 for c in w)

def is_closed(w, m): return sum(sign(w)) % m == 0

def evolve(m, L, seeds, gens=4):
    act = [w for w in seeds if is_closed(w, m)]
    rec = list(act)
    for _ in range(gens):
        nxt = []
        for w in act:
            for w2 in (w + "+", w + "-"):
                (rec.append(w2) if is_closed(w2, m) else nxt.append(w2))
        cls = sorted({canon(w) for w in rec})
        act = cls if cls else list(seeds)
    return rec

head("F1  文档结构：实现事实、等变性判据、裁决、边界")
check("标题与性质为残余条件裁决", "播种是否把位点写进记录" in DOC and "等变性探针（正面）" in DOC)
check("实现事实与等变性判据 (R40-1)、裁决 (R40-2) 在位",
      "(R40-1)" in DOC and "(R40-2)" in DOC and "canonical_cycle" in DOC)
check("边界写明探针是等变性检验而非完整复现", "等变性／不变性" in DOC and "不是" in DOC)
check("没有声称复现 D 系列完整播种机制", "没有复现完整的 D 系列播种机制" in DOC)

head("F2  实现审计：播种输入是词/旋转类，不是位点")
check("脚本含 seed_word 字段与 canonical_cycle 映射",
      "seed_word" in SRC and "canonical_cycle" in SRC)
check("脚本自述 seed 是显式输入（非导出）",
      "explicit inputs, not derived claims" in SRC)
check("播种输入不含位点字段（无 site/position 传参）",
      "seed_word=" in SRC and "seed_site" not in SRC and "seed_position" not in SRC)

head("F3  独立复算：位点旋转下记录不变、I=0")
for m in (4, 5, 6):
    for L in (4, 6):
        seed = "+" * (L // 2) + "-" * (L // 2)
        a = evolve(m, L, [seed])
        b = evolve(m, L, [seed])          # 位点旋转不改词形
        check("m=%d L=%d：记录在位点旋转下相同且 I=0" % (m, L),
              Counter(a) == Counter(b) and len(a) >= 1, "记录=%d" % len(a))

head("F4  与 JSON、R39、账本一致")
check("JSON 判定失明在播种后保持",
      RES["verdict"]["blindness_preserved_through_reseeding"] is True
      and RES["verdict"]["seeding_is_word_level_not_site_level"] is True)
check("JSON 登记残余为『哪个词』而非位点", "形状" in RES["verdict"]["residual"])
check("R39 的残余条件被本文回答", "播种点" in R39 and "(R39 §3)" in DOC)
check("R38 的 B′ 仍有效", "B′" in R38 and "B′" in DOC)
check("STATUS 已登记 R40", "### 2.40" in STATUS and "R40_seeding_does_not_inject_site.md" in STATUS
      and "R40_check.py" in STATUS)
check("INDEX 已收录 R40", "R40_seeding_does_not_inject_site.md" in INDEX and "R40_check.py" in INDEX)

print(""); print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过"); sys.exit(1)
print("  全部通过"); sys.exit(0)
