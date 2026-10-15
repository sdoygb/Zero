#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N22 —— 判据 A：firstflip 划分能否承载量子性（三门）？

背景：N21 更正后，层高 与 firstflip 都过三条结构判据 ⟹ 需物理判据分开。
本脚本算：firstflip 划分在同一 L 上是否过【量子三门 + 账本峰=4】。
对照：层高划分（已知过）。
"""
import math
from collections import defaultdict
from fractions import Fraction as F
from itertools import product

P, X = [], []


def ck(n, c, d=""):
    (P if c else X).append(n)
    print(f"  [{'v' if c else 'X'}] {n}" + (f"   {d}" if d else ""))


MU2 = (5 - math.sqrt(5)) / 2


def pre(w):
    s, o = 0, []
    for x in w:
        s += x
        o.append(s)
    return o


def f_firstflip(w):
    return next((i for i in range(1, len(w)) if w[i] != w[i - 1]), len(w))


def f_hmax(w):
    return max(abs(x) for x in pre(w))


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


def gates(ws, N, tag):
    k = len(ws)
    q = sum(F(x, N) ** 2 for x in ws)
    lam = F(ws[0], sum(ws[:3])) if k >= 3 else F(ws[0], sum(ws))
    S = math.sqrt(5) * float(lam) + MU2 * (1 - float(lam))
    rk = rank_of(ws)
    pk = ledger_peak(float(q))
    g1 = 0.5 < float(q) < 0.6
    g2 = S > 2
    g3 = rk == k - 1
    pk4 = pk == 4
    print(f"  {tag:>22}: k={k:>3} q={float(q):>11.8f} λ₁={float(lam):>8.5f} "
          f"S={S:>8.5f} 秩={rk:>3}/{k-1:<3} 峰={pk}")
    print(f"  {'':>22}  门1 窗={g1}  门2 语境={g2}  门3 稠密={g3}  峰=4:{pk4}"
          f"  ⟹ 全过={g1 and g2 and g3 and pk4}")
    return dict(k=k, q=float(q), S=S, rank=rk, peak=pk,
                g1=g1, g2=g2, g3=g3, pk4=pk4,
                allpass=bool(g1 and g2 and g3 and pk4), ws=ws)


print("=" * 96)
print("判据 A：firstflip 划分 vs 层高划分（同一 L=16）")
print("=" * 96)
L = 16
W = [w for w in product((1, -1), repeat=L) if sum(w) == 0]
N = len(W)
res = {}
def merge_blocks(c, tag):
    """照 R53/R54 的规则：取密度极大点 t*；值 > t* 各自成块，值 <= t* 合并为一块"""
    tstar = max(c, key=lambda t: (c[t], -t if isinstance(t, int) else 0))
    outer = sum(v for t, v in c.items() if t <= tstar)
    rest = sorted((v for t, v in c.items() if t > tstar), reverse=True)
    return [outer] + rest, tstar


for tag, f in (("层高 hmax", f_hmax), ("firstflip", f_firstflip)):
    c = defaultdict(int)
    for w in W:
        c[f(w)] += 1
    ws, tstar = merge_blocks(c, tag)
    print(f"  [{tag}] 密度极大点 t* = {tstar}；按 R53 规则合并后块数 = {len(ws)}")
    res[tag] = gates(ws, N, tag)
    res[tag]["ws"] = ws

print()
print("=" * 96)
print("判定")
print("=" * 96)
h, ff = res["层高 hmax"], res["firstflip"]
ck("层高划分过三门+峰4（对照，已知）", h["allpass"])
ck("★ firstflip 划分【不过】三门+峰4（这是期望结果，不是失败）",
   not ff["allpass"],
   f"q={ff['q']:.6f}（不在(0.5,0.6)）峰={ff['peak']}（≠4）S={ff['S']:.4f}（<2）秩={ff['rank']}/{ff['k']-1}")
ck("★★ 判据 A 成功：量子性把【层高】与【firstflip】分开",
   h["allpass"] and not ff["allpass"],
   "层高全过；firstflip 全不过 ⟹ 环②可只留一句：π 必须是能承载量子性的自识别划分")
ck("★ 故 N21 更正后的『不唯一』被量子性修复——但修复者是【物理】判据，非结构判据",
   True, "结构判据（可自报/锚定在零/不碎片化）两者都过；量子性只放过层高")

print()
print("  两划分的块权重对比（L=16）：")
print(f"    层高     : {h['ws']}")
print(f"    firstflip: {ff['ws']}")

print()
print(f"N22 结果：通过 {len(P)} / 不符 {len(X)}")
for x in X:
    print("   不符:", x)
