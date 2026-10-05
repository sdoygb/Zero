#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R60 · 用局部配平（K11）打破 III_1 的素数支撑上限

问题（R59 之后）：
  初始段块权重 ω_h ∝ C(2h,h)/4^h，其素因子全部 ≤ 2h
  ⇒ 秩 = min(k-1, π(2h)) ≈ 2k/ln(2k) ≪ k-1
  ⇒ 不是 III_1 ⇒ L1 的 boost 没有容身之处。

本轮：按 K11（局部净电荷 + 全局配平）把路径分成 m 个局域站点，
      每站点有自己的壳层剖面与续接数，再聚合成块权重。
      站点标识（m 与站点的局部计数）会引入**与 2h 无关的素数**。

测：块权重的秩亏缺是否消失（亏缺/k → 0 即 III_1 方向）。
"""
from __future__ import annotations

import json
import math
import os
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}


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


def rank(vals):
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


# ---------------------------------------------------------------- 无站点（基线）
def baseline(L):
    """ω_h ∝ 高度为 h 的初始段权重和（无站点）"""
    b = Counter()
    for k in range(L + 1):
        for s in range(-k, k + 1, 2):
            rem = L - k
            if (rem - s) % 2:
                continue
            a = (rem - s) // 2
            if a < 0 or a > rem:
                continue
            b[abs(s)] += comb(rem, a)
    return {h: b[h] for h in sorted(b)}


# ---------------------------------------------------------------- 有站点（K11）
def with_sites(L, m, block_bits=None):
    """把 L 步分成 m 段（站点）。每段有自己的净电荷 c 与最大振幅 a。
    续接数由 DP 精确计算：w[k][c][a] = 长度 k、净电荷 c、最大|部分和| = a 的路径数。
    块 = 站点数 m 与 (c_a, a_a) 的组合；这里取块 = 全局高度 h，但权重含站点因子。"""
    if block_bits is None:
        block_bits = L // m
    # DP：w[k][c][a]
    w = [defaultdict(Counter) for _ in range(L + 1)]
    w[0][0][0] = 1
    for k in range(1, L + 1):
        for c, tbl in w[k - 1].items():
            for a, cnt in tbl.items():
                for st in (1, -1):
                    nc = c + st
                    na = max(a, abs(nc))
                    w[k][nc][na] += cnt
    # 每段端点位置
    cuts = [round(i * L / m) for i in range(m + 1)]
    # 分段：每段的 (净电荷, 最大振幅) —— 段内部分和以段起点为基准
    # 用 DP 分别算每段的剖面，再按站点求和
    site = Counter()
    for i in range(m):
        k0, k1 = cuts[i], cuts[i + 1]
        n = k1 - k0
        # 段内所有可能的 (净电荷 c, 振幅 a)，权重 = 该段路径数（朴素段模型）
        seg = Counter()
        for c, tbl in w[n].items():
            for a, cnt in tbl.items():
                seg[(c, a)] += cnt
        # 站点因子：站点的局部计数（把 (c,a) 组合聚合）
        for (c, a), v in seg.items():
            site[(i % 2, c, a)] += v      # 站点奇偶性作为标识
    return site, w


def rank_of(counter):
    vals = [v for _, v in sorted(counter.items())]
    return rank(vals)


if __name__ == "__main__":
    print("=" * 74)
    print("基线：无站点（块权重 ∝ C(2h,h)/4^h 的整数形式）")
    print("=" * 74)
    print(" h   k  秩  需 k-1  亏缺  亏缺/k")
    vals = []
    base_rows = []
    for h in range(0, 27):
        vals.append(comb(2 * h, h))
        r, npr = rank(vals)
        k = h + 1
        if h >= 3:
            d = (k - 1) - r
            base_rows.append({"h": h, "k": k, "rank": r, "need": k - 1,
                              "def": d, "ratio": round(d / k, 4)})
            print(" %-3d %-3d %-3d %-6d %-5d %.4f" % (h, k, r, k - 1, d, d / k))
    OUT["baseline"] = base_rows

    print()
    print("=" * 74)
    print("加入站点（K11）：块 = (站点奇偶, 净电荷, 振幅)")
    print("=" * 74)
    site_rows = []
    for L in [8, 10, 12, 14, 16]:
        for m in [2, 3, 4]:
            st, _ = with_sites(L, m)
            r, npr = rank_of(st)
            k = len(st)
            print("  L=%-3d m=%-2d 块数=%-4d 素数支撑=%-3d 秩=%-4d 需=%-4d 亏缺=%d (%.3f)"
                  % (L, m, k, npr, r, k - 1, (k - 1) - r,
                     ((k - 1) - r) / k if k else 0))
            site_rows.append({"L": L, "m": m, "blocks": k, "primes": npr,
                              "rank": r, "need": k - 1,
                              "def": (k - 1) - r})
    OUT["with_sites"] = site_rows

    print()
    print("=" * 74)
    print("判据：亏缺/k 是否随规模下降（→0 即 III_1 方向）")
    print("=" * 74)
    for L in [8, 10, 12, 14, 16]:
        rs = [x for x in site_rows if x["L"] == L]
        if rs:
            best = min(rs, key=lambda x: (x["need"] - x["rank"]) / max(x["blocks"], 1))
            print("  L=%-3d 最好: m=%d 亏缺/k=%.4f (块=%d 秩=%d 需=%d)"
                  % (L, best["m"], (best["need"] - best["rank"]) / best["blocks"],
                     best["blocks"], best["rank"], best["need"]))
    with open(os.path.join(HERE, "R60_site_prime_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("\n-> R60_site_prime_results.json")
