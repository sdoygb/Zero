#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R76 · III_1 的素数来源：D=4 单纯形的圈结构（修正版）

R75 的两处修正：
  (1) 圈基构造错误（只给了 2 个 4-圈，不是 6 维基）
  (2) 需区分"圈的**个数**（按长度分级）"与"6 个独立和乐类"

本探针：
  A  K_5 的圈空间维数与一组正确的基
  B  按长度分级的圈数 N(ℓ)（从闭游走枚举，与解析式对照）
  C  关键检验：整个圈结构的**素因子支撑**有多大？
  D  对称性障碍：K_5 是顶点传递的 ⇒ 任何对称权重给所有类相等值
"""
from __future__ import annotations

import json
import math
import os
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}
D4 = 4
NV = D4 + 1
VERTS = list(range(NV))
EDGES = set(combinations(VERTS, 2))


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
        return 0
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
    return rk


# ------------------------------------------------------------ A 圈基
def cycle_basis():
    """K_5 的圈空间基（维数 = 10 − 5 + 1 = 6）。
    用生成树 {0-i : i=1..4} 的基本圈：非树边 (i,j)，i,j>=1 给出三角形 (0,i,j)"""
    basis = []
    for i, j in combinations([1, 2, 3, 4], 2):
        basis.append((0, i, j))          # 三角形
    return basis


def cycle_rank_check():
    """用关联矩阵精确算圈秩"""
    idx = {v: i for i, v in enumerate(VERTS)}
    el = list(EDGES)
    B = [[0] * len(el) for _ in VERTS]
    for j, (a, b) in enumerate(el):
        B[idx[a]][j] = -1
        B[idx[b]][j] = 1
    # 精确秩
    M = [[Fraction(x) for x in row] for row in B]
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
    return len(el) - rk, rk, len(el)


# ------------------------------------------------------------ B 按长度的圈数
def simple_cycles_by_length(maxlen):
    """K_5 上所有简单圈，按长度分级（无向，规范到最小旋转）"""
    def canon(c):
        n = len(c)
        rots = [c[i:] + c[:i] for i in range(n)]
        rev = tuple(reversed(c))
        rots += [rev[i:] + rev[:i] for i in range(n)]
        return min(rots)

    out = defaultdict(set)
    for l in range(3, maxlen + 1):
        for perm in permutations(VERTS, l):
            ok = True
            for i in range(l):
                if tuple(sorted((perm[i], perm[(i + 1) % l]))) not in EDGES:
                    ok = False
                    break
            if not ok:
                continue
            out[l].add(canon(perm))
    return out


# ------------------------------------------------------------ C 素因子支撑
def prime_support(values):
    ps = set()
    for v in values:
        ps |= set(fac(v).keys())
    return ps


if __name__ == "__main__":
    print("=" * 74)
    print("A · K_5 的圈空间")
    print("=" * 74)
    cr, rk, ne = cycle_rank_check()
    print("  顶点=%d 边=%d rank(B)=%d  圈秩=%d" % (NV, ne, rk, cr))
    print("  C(4,2)=%d  ⇒ 相等: %s" % (math.comb(4, 2), cr == math.comb(4, 2)))
    basis = cycle_basis()
    print("  一组基（%d 个三角形）: %s" % (len(basis), basis))
    OUT["A"] = {"cycle_rank": cr, "basis": basis}

    print()
    print("=" * 74)
    print("B · 按长度分级的圈数 N(ℓ)")
    print("=" * 74)
    bylen = simple_cycles_by_length(5)
    print("  ℓ   圈数 N(ℓ)   分解")
    lens, counts = [], []
    for l in sorted(bylen):
        c = len(bylen[l])
        lens.append(l)
        counts.append(c)
        print("  %-3d %-11d %s" % (l, c, fac(c)))
    OUT["B"] = {"lengths": lens, "counts": counts}

    print()
    print("  解析对照：K_n 上长度 ℓ 的简单圈数 = C(n,ℓ)·(ℓ−1)!/2")
    for l in sorted(bylen):
        pred = math.comb(NV, l) * math.factorial(l - 1) // 2
        print("    ℓ=%d  实测=%d  解析=%d  %s"
              % (l, len(bylen[l]), pred,
                 "一致" if len(bylen[l]) == pred else "不一致"))

    print()
    print("=" * 74)
    print("C · 素因子支撑（III_1 的关键）")
    print("=" * 74)
    ps = prime_support(counts)
    print("  分级圈数的素因子支撑 = %s" % sorted(ps))
    print("  支撑大小 = %d ；需要的秩 = %d" % (len(ps), len(counts) - 1))
    print("  秩 = %d / 需 %d" % (rank_of(counts), len(counts) - 1))
    print()
    print("  逐项分解：")
    for l, c in zip(lens, counts):
        print("    ℓ=%d  N=%d = %s" % (l, c, fac(c)))
    OUT["C"] = {"primes": sorted(ps), "rank": rank_of(counts),
                "needed": len(counts) - 1}

    print()
    print("=" * 74)
    print("D · 对称性障碍")
    print("=" * 74)
    print("  K_5 的顶点传递性检验（自同构群在顶点/边上传递）：")
    vertex_transitive = True
    for a in VERTS:
        for b in VERTS:
            if a != b:
                found = False
                for perm in permutations(VERTS):
                    if perm[a] == b and all(
                            tuple(sorted((perm[x], perm[y]))) in EDGES
                            for x, y in EDGES):
                        found = True
                        break
                if not found:
                    vertex_transitive = False
    print("    顶点传递: %s" % vertex_transitive)
    print()
    print("  ⇒ 任何只依赖图结构的权重 w(c) 在所有 6 个基圈上相等 ⇒ 秩 = 0-1")
    print("  ⇒ 6 个独立和乐类**无法**通过图对称性区分为不同权重")
    print("  ⇒ III_1 需要**外来的等级结构**打破 S_4 对称")
    OUT["D"] = {"vertex_transitive": bool(vertex_transitive)}

    with open(os.path.join(HERE, "R76_cycle_structure_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R76_cycle_structure_results.json")
