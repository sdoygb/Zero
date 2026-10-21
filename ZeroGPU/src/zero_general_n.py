#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""一般 N 的 L1 词基：定义检查 G1–G6（推导稿的承重定义）。

G1 |W_{L,N}| = C(L,N)
G2 半满极限 W_{L,L/2} = Z168 的平衡词集
G3 半满时 h/a/r 逐字回到 Z168（layer_of/word_area/word_runs）
G4 L=4,N=3 的表（h/a/r/z）
G5 粒子-空穴 w→-w 把 W_{L,N} 送到 W_{L,L-N}
G6 L=4,N=3 的激发胞腔（HF + 单激发；rank 0/1）

输出：results/zero_general_n.json
"""
from __future__ import annotations

import itertools
import json
import pathlib
import sys
from math import comb

HERE = pathlib.Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import zero_sim as zs  # noqa: E402

OUT = HERE.parents[1] / "results" / "zero_general_n.json"
OK, BAD, ROWS = [], [], []


def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def words_general(L, N):
    out = []
    for pos in itertools.combinations(range(L), N):
        w = [-1] * L
        for i in pos:
            w[i] = 1
        out.append(tuple(w))
    return out


def height(w):
    s = hi = lo = 0
    for x in w:
        s += x
        hi = max(hi, s)
        lo = min(lo, s)
    return max(hi, -lo)


def area(w):
    s = tot = 0
    for x in w:
        s += x
        tot += s
    return tot


def runs(w):
    return 1 + sum(1 for i in range(1, len(w)) if w[i] != w[i - 1])


def zc(w):
    s = 0
    n = 0
    for i, x in enumerate(w):
        s += x
        if 0 < i < len(w) - 1 and s == 0:
            n += 1
    return n


def main():
    print("=" * 72)
    print("zero_general_n.py — 一般 N 词基定义检查")
    print("=" * 72)
    # G1
    g1 = []
    for L, N in [(4, 2), (4, 3), (6, 3), (6, 2), (8, 4), (5, 2)]:
        n = len(words_general(L, N))
        g1.append(n == comb(L, N))
    ck("G1 |W_{L,N}| = C(L,N)", all(g1), str(g1))

    # G2
    g2 = []
    for L in [4, 6, 8]:
        eq = set(words_general(L, L // 2)) == set(zs.balanced_words(L))
        g2.append(eq)
    ck("G2 half-filled == Z168 balanced", all(g2), str(g2))

    # G3
    g3 = []
    for L in [4, 6, 8]:
        for w in zs.balanced_words(L):
            g3.append(height(w) == zs.layer_of(w) and area(w) == zs.word_area(w)
                      and runs(w) == zs.word_runs(w))
    ck("G3 features reduce to Z168 at N=L/2", all(g3), f"{sum(g3)}/{len(g3)} words")

    # G4
    tbl = {}
    for w in words_general(4, 3):
        tbl["".join("+" if x > 0 else "-" for x in w)] = [height(w), area(w), runs(w), zc(w)]
    exp = {"+++-": [3, 8, 2, 0], "++-+": [2, 6, 3, 0],
           "+-++": [2, 4, 3, 1], "-+++": [2, 2, 2, 1]}
    ck("G4 L=4,N=3 table", tbl == exp, str(tbl))

    # G5
    g5 = []
    for L, N in [(4, 3), (6, 2), (8, 3)]:
        ph = set(tuple(-x for x in w) for w in words_general(L, N))
        g5.append(ph == set(words_general(L, L - N)))
    ck("G5 particle-hole W_{L,N} -> W_{L,L-N}", all(g5), str(g5))

    # G6
    hf = (1, 1, 1, -1)
    ranks = {}
    for w in words_general(4, 3):
        k = sum(1 for a, b in zip(w, hf) if a != b) // 2
        ranks[k] = ranks.get(k, 0) + 1
    ck("G6 L=4,N=3 cells", ranks == {0: 1, 1: 3}, str(ranks))

    out = {"script": "zero_general_n.py", "checks": ROWS,
           "n_pass": len(OK), "n_fail": len(BAD),
           "table_L4_N3": tbl, "cells_L4_N3": ranks}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
