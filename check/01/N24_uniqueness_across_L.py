#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N24 —— 跨 L 的普查：能否承载量子性的自识别划分，是否唯一？

背景（U24 提出的可判问题）：
  L=16 上，过【量子三门 + 峰=4】的原生划分只有层高划分。
  问：换 L 呢？
    · 若对一切 L 都只有层高 ⟹ 那是【必然】
    · 若只在某些 L 有      ⟹ 那是【偶然】，且给出对 L 的约束

候选（自识别的原生划分，按 R53 合并规则）：
  hmax, firstflip, nflip, maxrun, area, orbit
判据：门1 q∈(0.5,0.6)；门2 S_max>2；门3 rank=k−1；峰=4
"""
import itertools
import math
from collections import defaultdict
from fractions import Fraction as F
from itertools import product

P, X = [], []


def ck(n, c, d=""):
    (P if c else X).append(n)
    print(f"  [{'v' if c else 'X'}] {n}" + (f"   {d}" if d else ""))


MU1, MU2 = math.sqrt(5), (5 - math.sqrt(5)) / 2


def pre(w):
    s, o = 0, []
    for x in w:
        s += x
        o.append(s)
    return o


def nflip(w):
    return sum(1 for i in range(1, len(w)) if w[i] != w[i - 1])


def firstflip(w):
    return next((i for i in range(1, len(w)) if w[i] != w[i - 1]), len(w))


def maxrun(w):
    return max(len(list(g)) for _, g in itertools.groupby(w))


def hmax(w):
    return max(abs(x) for x in pre(w))


def area(w):
    return sum(pre(w))


def orbit_size(w):
    L = len(w)
    return len({w[k:] + w[:k] for k in range(L)})


FE = {
    "hmax": hmax,
    "firstflip": firstflip,
    "nflip": nflip,
    "maxrun": maxrun,
    "area": area,
    "orbit": orbit_size,
}


def factorize(n):
    f, d, m = {}, 2, n
    while d * d <= m:
        while m % d == 0:
            f[d] = f.get(d, 0) + 1
            m //= d
        d += 1
    if m > 1:
        f[m] = f.get(m, 0) + 1
    return f


def rank_of(ws):
    exps = [factorize(v) for v in ws]
    primes = sorted({p for e in exps for p in e})
    if not primes:
        return 0
    M = [[F(exps[i].get(p, 0) - exps[i + 1].get(p, 0)) for p in primes]
         for i in range(len(exps) - 1)]
    r = 0
    for c in range(len(primes)):
        piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                g = M[i][c]
                M[i] = [x - g * y for x, y in zip(M[i], M[r])]
        r += 1
    return r


def ledger_peak(q):
    for D in range(2, 600):
        if D == 2:
            if 3 * q < 1:
                return 2
            continue
        if D * q > (D - 2) and (D + 1) * q < (D - 1):
            return D
    return None


def merge(c):
    """R53 合并规则：t* = 密度极大点；值 > t* 各自成块，值 <= t* 合并"""
    ts = max(c, key=lambda t: (c[t], -t if isinstance(t, int) else 0))
    outer = sum(v for t, v in c.items() if t <= ts)
    return [outer] + sorted((v for t, v in c.items() if t > ts), reverse=True), ts


def dist_dp(L, kind):
    """DP 算块大小分布（不枚举词）。kind ∈ {hmax, nflip, firstflip, maxrun}"""
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
    if kind == "maxrun":
        cur = {(0, 0, 1, 1): 1}       # (s, 上一步, 当前游程, 最大游程)
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
    return None


DPK = {"hmax": "hmax", "nflip": "nflip", "firstflip": "firstflip", "maxrun": "maxrun"}

print("=" * 100)
print("跨 L 普查：能承载量子性的自识别划分，是否唯一？")
print("=" * 100)
print(f"  {'L':>3} " + "".join(f"{n:>13}" for n in FE) + "   通过者")
summary = {}
for L in range(6, 25, 2):
    W = [w for w in product((1, -1), repeat=L) if sum(w) == 0]
    N = len(W)
    row = f"  {L:>3} "
    passed = []
    for nm, f in FE.items():
        if nm in DPK:
            c = dist_dp(L, DPK[nm])
        else:
            c = defaultdict(int)
            for w in W:
                c[f(w)] += 1
        ws, ts = merge(c)
        k = len(ws)
        q = sum(F(x, N) ** 2 for x in ws)
        if k >= 3:
            lam = F(ws[0], sum(ws[:3]))
            S = MU1 * float(lam) + MU2 * (1 - float(lam))
        else:
            S = float("nan")
        rk = rank_of(ws)
        pk = ledger_peak(float(q))
        ok = (0.5 < float(q) < 0.6) and (S > 2) and (rk == k - 1) and (pk == 4)
        if ok:
            passed.append(nm)
        mark = "★" if ok else " "
        row += f"{float(q):>12.6f}{mark}"
    row += f"   {passed}"
    print(row)
    summary[L] = passed

print()
print("=" * 100)
print("判定")
print("=" * 100)
hit = {L: p for L, p in summary.items() if p}
print(f"  有通过者的 L：{list(hit)}")
for L, p in hit.items():
    print(f"    L={L}: {p}")
ck("★ 存在通过量子门的 L", len(hit) > 0, f"{list(hit)}")
ck("★ 通过者【随 L 变化】：L=16 是 hmax，L=20 是 maxrun（意外，见 U26）",
   hit.get(16) == ["hmax"] and hit.get(20) == ["maxrun"], f"{hit}")
ck("★ 通过者唯一（每 L 恰一个）",
   all(len(p) == 1 for p in hit.values()),
   " → ".join(f"L={L}:{p}" for L, p in hit.items()))
ck("★ L=16 是其中一个（与 U24 一致）", 16 in hit, f"L=16 → {hit.get(16)}")
print()
print(f"N24 结果：通过 {len(P)} / 不符 {len(X)}")
for x in X:
    print("   不符:", x)
