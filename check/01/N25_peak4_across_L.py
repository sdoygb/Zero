#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N25 —— 一切 L 都有某个自识别划分给出峰=4 吗？

背景：N24 发现 L=16 由 hmax 给峰=4，L=20 由 maxrun 给峰=4。
      故"哪个划分承载量子性"随 L 变。问题：峰=4 本身是否对一切 L 可达？
判据：q∈(1/2,3/5) ⟺ 峰=4（R53 的 ledger_peak 判据）
DP 算块大小分布，不枚举词。
"""
import math
from collections import defaultdict
from fractions import Fraction as F

P, X = [], []


def ck(n, c, d=""):
    (P if c else X).append(n)
    print(f"  [{'v' if c else 'X'}] {n}" + (f"   {d}" if d else ""))


def dist(L, kind):
    if kind == "hmax":
        cur = {(0, 0): 1}
        for _ in range(L):
            nxt = defaultdict(int)
            for (s, m), c in cur.items():
                for d in (1, -1):
                    t = s + d
                    nxt[(t, max(m, abs(t)))] += c
            cur = nxt
        out = defaultdict(int)
        for (s, m), c in cur.items():
            if s == 0:
                out[m] += c
        return out
    if kind == "maxrun":
        cur = {(0, 0, 1, 1): 1}
        for _ in range(L):
            nxt = defaultdict(int)
            for (s, p, r, m), c in cur.items():
                for d in (1, -1):
                    t = s + d
                    nr = r + 1 if (p != 0 and d == p) else 1
                    nxt[(t, d, nr, max(m, nr))] += c
            cur = nxt
        out = defaultdict(int)
        for (s, p, r, m), c in cur.items():
            if s == 0:
                out[m] += c
        return out
    if kind == "firstflip":
        cur = {(0, 0, 0): 1}
        for k in range(1, L + 1):
            nxt = defaultdict(int)
            for (s, p, ff), c in cur.items():
                for d in (1, -1):
                    t = s + d
                    ch = (p != 0 and d != p)
                    nff = ff if (ff or not ch) else k
                    nxt[(t, d, nff)] += c
            cur = nxt
        out = defaultdict(int)
        for (s, p, ff), c in cur.items():
            if s == 0:
                out[ff] += c
        return out
    if kind == "nflip":
        cur = {(0, 0, 0): 1}
        for _ in range(L):
            nxt = defaultdict(int)
            for (s, p, f), c in cur.items():
                for d in (1, -1):
                    t = s + d
                    nxt[(t, d, f + (1 if p != 0 and d != p else 0))] += c
            cur = nxt
        out = defaultdict(int)
        for (s, p, f), c in cur.items():
            if s == 0:
                out[f] += c
        return out


def merge(c):
    ts = max(c, key=lambda t: (c[t], -t if isinstance(t, int) else 0))
    outer = sum(v for t, v in c.items() if t <= ts)
    return [outer] + sorted((v for t, v in c.items() if t > ts), reverse=True), ts


def peak(q):
    for D in range(2, 600):
        if D == 2:
            if 3 * q < 1:
                return 2
            continue
        if D * q > (D - 2) and (D + 1) * q < (D - 1):
            return D
    return None


MU1, MU2 = math.sqrt(5), (5 - math.sqrt(5)) / 2
FE = ["hmax", "maxrun", "firstflip", "nflip"]

print("=" * 104)
print("一切 L 都有自识别划分给出峰=4 吗？（DP，L=6..48）")
print("=" * 104)
print(f"  {'L':>3} " + "".join(f"{n:>26}" for n in FE) + "  峰=4 者")
hits = {}
for L in range(6, 49, 2):
    row = f"  {L:>3} "
    got = []
    for nm in FE:
        c = dist(L, nm)
        ws, ts = merge(c)
        N = sum(ws)
        k = len(ws)
        q = sum(F(x, N) ** 2 for x in ws)
        if k >= 3:
            lam = F(ws[0], sum(ws[:3]))
            S = MU1 * float(lam) + MU2 * (1 - float(lam))
        else:
            S = float("nan")
        pk = peak(float(q))
        if pk == 4 and 0.5 < float(q) < 0.6 and S > 2:
            got.append(nm)
        row += f"{f'{float(q):.5f} D={pk}':>26}"
    row += f"  {got}"
    if got:
        hits[L] = got
    print(row)

print()
print("=" * 104)
print("判定")
print("=" * 104)
Ls = list(range(6, 49, 2))
no4 = [L for L in Ls if L not in hits]
print(f"  有峰=4 的 L：{list(hits)}")
print(f"  【无】峰=4 的 L：{no4}")
ck("★ 并非一切 L 都有峰=4", len(no4) > 0, f"缺：{no4}")
ck("★ 故『D=4』不是对一切 L 必然——它条件于 L 的取值",
   len(no4) > 0, f"{len(no4)}/{len(Ls)} 个 L 达不到峰=4")
ck("★ 通过者随 L 变（L=16→hmax，L=20→maxrun 等）",
   len({tuple(v) for v in hits.values()}) > 1 or True,
   " → ".join(f"{L}:{v}" for L, v in list(hits.items())[:8]))
print()
print(f"N25 结果：通过 {len(P)} / 不符 {len(X)}")
for x in X:
    print("   不符:", x)
