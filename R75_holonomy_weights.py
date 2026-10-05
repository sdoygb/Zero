#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R75 · III_1 的素数来源：D=4 单纯形上 6 个独立和乐类的权重

先更正 R73 的一处错误：
  R73 曾称"身份秩比通道对数短缺 D+1"。这是把**一般 K_n** 的圈秩公式
      dim H_1(K_n) = C(n,2) − n + 1 = C(n−1,2)
  误用于 D-单纯形。正确关系：
      D-单纯形（D+1 顶点）的圈秩 = C(D+1,2) − (D+1) + 1 = C(D,2)   ← 精确相等
  ⇒ 短缺的不是圈秩，而是"通道对数 = 树边 + 圈秩"这一分解：
      C(D+1,2) = (D+1) + C(D,2)
  ⇒ 6 个独立和乐类是对的（D=4）。

本探针：在 K_5（= 4-单纯形的 1-骨架）上
  A  找出 6 个独立和乐类（圈空间的一组基）
  B  用若干原生权重方案给它们赋权，测秩是否 = 5（III_1）
  C  检验素因子支撑
"""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}
D4 = 4
NVERT = D4 + 1                     # 5
VERTS = list(range(NVERT))
EDGES = list(combinations(VERTS, 2))


# ---------------------------------------------------------------- 基本工具
def fac(n):
    f, d = {}, 2
    n = int(n)
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def rank_of(vals):
    fs = [fac(v) for v in vals]
    primes = sorted({p for f in fs for p in f})
    rows = [[fs[i + 1].get(p, 0) - fs[i].get(p, 0) for p in primes]
            for i in range(len(fs) - 1)]
    if not rows:
        return 0, 0
    M = [[Fraction(x) for x in r] for r in rows]
    rk = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        pv = M[rk][c]
        M[rk] = [x / pv for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk, len(primes)


# ---------------------------------------------------------------- A 圈的基
def cycle_basis_K5():
    """返回 K_5 的 6 个独立圈（用顶点序列表示，规范到最小起点/方向）"""
    cycles = []
    # 方案：所有含 3 或 4 个顶点的圈，取自 {0,1,2,3} 的边指标化
    for (a, b) in combinations([0, 1, 2, 3], 2):
        # 每个 K_4 边 (a,b) 对应一个 4-圈：a→b→x→y→a，其中 {x,y} 是补集
        rest = [v for v in [0, 1, 2, 3] if v not in (a, b)]
        c = tuple([a, b, rest[0], rest[1]])
        cyc = canon_cycle(c)
        if cyc not in cycles:
            cycles.append(cyc)
    return cycles


def canon_cycle(c):
    n = len(c)
    rots = [c[i:] + c[:i] for i in range(n)]
    rev = tuple(reversed(c))
    rots += [rev[i:] + rev[:i] for i in range(n)]
    return min(rots)


def edges_of_cycle(c):
    n = len(c)
    return [tuple(sorted((c[i], c[(i + 1) % n]))) for i in range(n)]


# ---------------------------------------------------------------- B 权重方案
def closed_walks_upto(L, start):
    """长度 ≤ L 的闭合游走（回到起点）计数，按终点无关——返回所有闭游走"""
    walks = []
    def rec(path):
        if len(path) > 1 and len(path) <= L + 1:
            if path[-1] == start:
                walks.append(tuple(path))
        if len(path) >= L + 1:
            return
        for nb in VERTS:
            if nb != path[-1] and tuple(sorted((nb, path[-1]))) in EDGES:
                rec(path + [nb])
    rec([start])
    return walks


def weight_by_edge_traversal(cycles, L):
    """方案 1：类 c 的权重 = 所有长度 ≤ L 的闭游走穿越 c 的边次数之和"""
    w = defaultdict(int)
    for s in VERTS:
        for wk in closed_walks_upto(L, s):
            es = set(tuple(sorted((wk[i], wk[i + 1])))
                     for i in range(len(wk) - 1))
            for c in cycles:
                if set(edges_of_cycle(c)) <= es:
                    w[c] += 1
    return [w[c] for c in cycles]


def weight_by_cycle_length(cycles, L):
    """方案 2：权重 = 把该圈作为最短闭游走的实现数（长度 ≤ L）"""
    w = defaultdict(int)
    for s in VERTS:
        for wk in closed_walks_upto(L, s):
            c = canon_cycle(tuple(wk[:-1]))
            if c in cycles:
                w[c] += 1
    return [w[c] for c in cycles]


def weight_by_edges_only(cycles, L):
    """方案 3：权重 = 圈的边数（对照，应当全等 ⇒ 素贫乏）"""
    return [len(edges_of_cycle(c)) for c in cycles]


if __name__ == "__main__":
    print("=" * 74)
    print("更正：D-单纯形的圈秩 = C(D,2)（精确相等）")
    print("=" * 74)
    from math import comb
    for D in range(3, 8):
        e = comb(D + 1, 2)
        cr = e - (D + 1) + 1
        print("  D=%d: 边=%d 顶点=%d 圈秩=%d  C(D,2)=%d  %s"
              % (D, e, D + 1, cr, comb(D, 2),
                 "相等" if cr == comb(D, 2) else "不等"))
    print()
    print("  ⇒ R73 的\"短缺\"说法撤回：短缺的是分解 C(D+1,2) = (D+1) + C(D,2)")
    OUT["correction"] = "cycle rank of D-simplex = C(D,2) exactly"

    print()
    print("=" * 74)
    print("A · D=4 单纯形（K_5）的 6 个独立和乐类")
    print("=" * 74)
    cycles = cycle_basis_K5()
    print("  类数 = %d（应为 C(4,2)=6）" % len(cycles))
    for i, c in enumerate(cycles):
        print("    c%d = %s   边=%s" % (i, c, edges_of_cycle(c)))
    OUT["cycles"] = [list(c) for c in cycles]

    print()
    print("=" * 74)
    print("B · 各权重方案的秩（需 5 才能 III_1）")
    print("=" * 74)
    rows = []
    for L in [4, 5, 6, 7, 8]:
        w1 = weight_by_edge_traversal(cycles, L)
        w2 = weight_by_cycle_length(cycles, L)
        w3 = weight_by_edges_only(cycles, L)
        r1, p1 = rank_of(w1)
        r2, p2 = rank_of(w2)
        r3, p3 = rank_of(w3)
        print("  L=%d" % L)
        print("    方案1 边穿越: w=%s 秩=%d/%d 素数=%d" % (w1, r1, len(w1) - 1, p1))
        print("    方案2 圈实现: w=%s 秩=%d/%d 素数=%d" % (w2, r2, len(w2) - 1, p2))
        print("    方案3 仅边数: w=%s 秩=%d/%d 素数=%d" % (w3, r3, len(w3) - 1, p3))
        rows.append({"L": L, "w1": w1, "r1": r1, "w2": w2, "r2": r2,
                     "w3": w3, "r3": r3})
    OUT["weights"] = rows

    print()
    print("=" * 74)
    print("C · 素因子分析（方案 2 的各 L）")
    print("=" * 74)
    for r in rows:
        fs = [str(fac(x)) for x in r["w2"]]
        print("  L=%d  %s" % (r["L"], fs))

    with open(os.path.join(HERE, "R75_holonomy_weights_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R75_holonomy_weights_results.json")
