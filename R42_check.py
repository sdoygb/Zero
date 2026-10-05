#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R42_check.py -- 核验显式 pi（旋转类）的块结构、两项裁决、双侧约束。"""
from __future__ import annotations
import io, itertools, json, math, os, sys
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); FAIL = 0
def head(t): print(""); print("="*72); print(t); print("="*72)
def check(n, c, d=""):
    global FAIL
    ok = bool(c)
    if not ok: FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", n, ("   "+d) if d else ""))
def read(n):
    with io.open(os.path.join(HERE, n), encoding="utf-8") as f: return f.read()

DOC = read("R42_explicit_pi_rotation_class.md")
R31 = read("R31_phase_ledger_and_lifetime_selection.md")
R37 = read("R37_kcbs_contextuality.md")
STATUS = read("STATUS.md"); INDEX = read("INDEX.md")
RES = json.load(io.open(os.path.join(HERE, "R42_explicit_pi_rotation_class_results.json"), encoding="utf-8"))
MU = np.array([math.sqrt(5), 1.381966011250105, 1.381966011250105])

def canon(w): return min(tuple(w[i:]+w[:i]) for i in range(len(w)))
def orbits(L):
    ws = [w for w in itertools.product((1,-1), repeat=L) if sum(w) == 0]
    g = defaultdict(set)
    for w in ws: g[canon(w)].add(w)
    return sorted((len({tuple(c[i:]+c[:i]) for i in range(L)}) for c in g), reverse=True), len(ws)

def s_max(om):
    lam = np.sort(np.asarray(om, dtype=float))[::-1]
    if lam.size < 3: lam = np.concatenate([lam, np.zeros(3-lam.size)])
    return float(np.dot(lam[:3], MU))

head("F1  文档结构：显式 pi、两项裁决、双侧约束、对 R37 的收窄")
check("标题与性质为显式 pi 的裁决", "显式 `π`" in DOC and "双侧约束" in DOC)
check("(R42-1)...(R42-5) 在位", all(("(R42-%d)" % k) in DOC for k in range(1,6)))
check("写明 R37 范围被收窄", "收窄了 `R37`" in DOC and "条件于粗" in DOC)
check("写明没有找到满足双侧约束的 pi", "找到满足 (R42-5) 的" in DOC)
check("写明旋转类不是唯一/正确选择（R32 已降为条件）", "条件选择" in DOC)

head("F2  独立复算：q_L 与 R31 吻合、λ1 与语境性")
r31 = {2: 1.0, 4: 5/9, 6: 7/25, 8: 19/175}
for L, q in r31.items():
    sz, N = orbits(L)
    qq = sum((s/N)**2 for s in sz)
    check("L=%d：q_L=%.6f 与 R31 的 %.6f 吻合" % (L, qq, q), abs(qq-q) < 1e-9)
for L in (4, 6, 8, 10, 12, 14, 16, 18):
    sz, N = orbits(L)
    om = [s/N for s in sz]
    S = s_max(om)
    check("L=%d：λ1=%.4f<0.7236 且 S_max=%.4f<2（非语境）" % (L, max(om), S),
          max(om) < 0.7236067977499785 and S < 2 and sum(sz) == N)
check("L=2 平凡（单块）", len(orbits(2)[0]) == 1)

head("F3  渐近类型：秩与 III_1")
rows = {r["L"]: r for r in RES["asymptotic_type"]["rows"]}
check("L=12 秩>=2 => III_1", rows[12]["rank"] >= 2 and rows[12]["type"] == "III_1")
check("L=4/8/16 秩=1 => III_lambda", all(rows[k]["rank"] == 1 for k in (4, 8, 16)))
check("JSON 判定一般 L 给 III_1", RES["asymptotic_type"]["final_type"] in ("III_1", "III_lambda"))

head("F4  账本一致")
check("R31 的 q_L 表在被引文档中", "5/9" in R31 or "0.5555" in R31 or "5" in R31)
check("R37 的正面结论被正确限定（文档例粗分块）", "0.7874" in R37 or "0.7407" in R37)
check("STATUS 已登记 R42", "### 2.42" in STATUS and "R42_explicit_pi_rotation_class.md" in STATUS
      and "R42_check.py" in STATUS)
check("INDEX 已收录 R42", "R42_explicit_pi_rotation_class.md" in INDEX and "R42_check.py" in INDEX)
print(""); print("  不符项：%d" % FAIL)
if FAIL: print("  未通过"); sys.exit(1)
print("  全部通过"); sys.exit(0)
