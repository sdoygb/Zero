#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N6 —— L=4 上「原生可分块」的穷举：张力能否消除？

问题（PI §3.6 的 (a)）：是否存在一套 **Z0-可定义**的分块，使 L=4 同时满足量子链的三条需求？

三条需求（照 R49／R54／R58 的口径）：
  (i)   块数 k ≥ 4   —— 因为秩需 = k−1，而 III_1 需秩 = 3（R49-2：k 块秩上限 k−1）
  (ii)  秩 = k−1     —— 素数指数差向量线性无关（R49 §3.2）
  (iii) 语境性 λ₁/顶三 > 0.7236（R49-3）

「原生」= 由 Z0 侧对象可定义：词、前缀、层高 h、位置、+1 位置集、旋转类、轨道大小、
            以及这些的 Z0-可定义函数（此处穷举其中可实际枚举的一族）。

本脚本只做**可判定**的部分：枚举分块与权重，算 (i)(ii)(iii)。
"""
import math
from itertools import product, combinations, chain
from fractions import Fraction as F
from collections import Counter

P, X = [], []


def ck(n, c, d=""):
    (P if c else X).append(n)
    print(f"  [{'v' if c else 'X'}] {n}" + (f"   {d}" if d else ""))


L = 4
W = [w for w in product((1, -1), repeat=L) if sum(w) == 0]
print(f"W_4 = {len(W)} 个词：")
for w in W:
    print(f"   {w}")

# ---------- Z0 侧可定义函数（原生特征） ----------
def prefix(w):
    s, out = 0, []
    for x in w:
        s += x
        out.append(s)
    return out


def hmax(w):
    return max(abs(x) for x in prefix(w))


def hset(w):
    return tuple(sorted({abs(x) for x in prefix(w)}))


def npos(w):
    """+1 出现位置集（1-based）"""
    return tuple(i for i, x in enumerate(w, 1) if x == 1)


def rot_class(w):
    return min(w[k:] + w[:k] for k in range(L))


def orbit_size(w):
    return len({w[k:] + w[:k] for k in range(L)})


def nflips(w):
    return sum(1 for i in range(1, L) if w[i] != w[i - 1])


FEATS = {
    "层高 h (max|s_k|)": hmax,
    "层高集合 (sorted set)": hset,
    "符号变化次数": nflips,
    "轨道大小": orbit_size,
    "+1 位置集（等价于词）": npos,
    "旋转类": rot_class,
    "首符号": lambda w: w[0],
    "末符号": lambda w: w[-1],
    "峰位置 (max|s_k| 的首个 k)": lambda w: next(
        i for i, x in enumerate(prefix(w), 1) if abs(x) == hmax(w)),
    "|s_2|": lambda w: abs(prefix(w)[1]),
}
for name, f in FEATS.items():
    vals = [f(w) for w in W]
    k = len(set(vals))
    print(f"   原生函数『{name}』给出 {k} 块")

# 组合出的分块（多个特征联立，product 型）
print()
print("=" * 74)
print("穷举：所有由『原生函数的有序组合』诱导的分块")
print("=" * 74)

names = list(FEATS)
blocks_list = []          # (描述, 分块(list of tuple of words))
seen_sig = set()
for r in range(1, len(names) + 1):
    for combo in combinations(names, r):
        key = []
        for w in W:
            key.append(tuple(FEATS[n](w) for n in combo))
        groups = {}
        for w, kk in zip(W, key):
            groups.setdefault(kk, []).append(w)
        part = tuple(sorted((tuple(sorted(g)) for g in groups.values())))
        if len(part) < 2:
            continue
        if part in seen_sig:
            continue
        seen_sig.add(part)
        blocks_list.append((" ∧ ".join(combo), part))

print(f"  得到 {len(blocks_list)} 个**互不相同**的分块")
sizes = Counter(len(b[1]) for b in blocks_list)
print(f"  块数分布：{dict(sorted(sizes.items()))}")


# ---------- 秩与语境性 ----------
def factorize(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def rank_of(weights):
    """素数指数差向量的秩（R49 §3.2）"""
    primes = sorted({p for v in weights for p in factorize(v)})
    if not primes:
        return 0
    vecs = []
    for a, b in zip(weights, weights[1:]):
        fa, fb = factorize(a), factorize(b)
        vecs.append([fa.get(p, 0) - fb.get(p, 0) for p in primes])
    # 精确秩（Fraction 高斯消元）
    M = [[F(x) for x in row] for row in vecs]
    rank, rows = 0, len(M)
    cols = len(primes)
    for c in range(cols):
        piv = None
        for r in range(rank, rows):
            if M[r][c] != 0:
                piv = r
                break
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        pv = M[rank][c]
        M[rank] = [x / pv for x in M[rank]]
        for r in range(rows):
            if r != rank and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[rank])]
        rank += 1
    return rank


def smax_ratio(weights):
    """λ₁ / 顶三（R49-3 的等价写法）。块数 <3 时按可用块数取和。"""
    s = sorted(weights, reverse=True)
    top = sum(s[:3])
    return F(s[0], top) if top else F(0)


print()
print("=" * 74)
print("逐分块检验三条需求")
print("=" * 74)
print(f"  {'k':>2} {'秩':>3} {'λ₁/顶三':>9} {'秩=k−1':>7} {'语境':>5} {'k≥4':>5} {'全过':>5}  分块")
passing = []
for desc, part in blocks_list:
    k = len(part)
    wts = [len(g) for g in part]
    rk = rank_of(wts)
    rat = smax_ratio(wts)
    a = (rk == k - 1)
    b = (rat > F(7236, 10000))
    c = (k >= 4)
    allok = a and b and c
    if allok:
        passing.append((desc, k, rk, rat, wts))
    print(f"  {k:>2} {rk:>3} {float(rat):>9.4f} {str(a):>7} {str(b):>5} {str(c):>5} "
          f"{str(allok):>5}  {desc}  w={wts}")

print()
print("=" * 74)
print(f"N6 判定：满足三条需求的分块数 = {len(passing)}")
print("=" * 74)
if passing:
    for d, k, rk, rat, wts in passing:
        print(f"   ✅ {d}：k={k}, 秩={rk}, λ₁/顶三={float(rat):.4f}, w={wts}")
else:
    print("   不存在 —— PI §3.6 的出路 (a) 在『特征组合』这一族内**走不通**")
# ---------- 不依赖特征族的彻底判定：权重构成穷举 ----------
print()
print("=" * 74)
print("彻底判定：6 个词 → k 块的**全部**权重构成（不依赖特征族选择）")
print("=" * 74)


def parts(n, k, lo=1):
    if k == 1:
        if n >= lo:
            yield (n,)
        return
    for a in range(lo, n // k + 1):
        for rest in parts(n - a, k - 1, a):
            yield (a,) + rest


def Smax_of(ws):
    w = sorted(ws, reverse=True)
    lam = F(w[0], sum(w[:3])) if len(w) >= 3 else F(w[0], sum(w))
    mu2 = (5 - math.sqrt(5)) / 2
    return math.sqrt(5) * float(lam) + mu2 * (1 - float(lam)), float(lam)


print(f"  {'k':>2} {'权重构成':>16} {'秩':>3} {'S_max':>7} {'k≥4':>5} {'秩=k−1':>7} {'语境':>5} {'全过':>5}")
feasible = []
for k in range(2, 7):
    for ws in parts(len(W), k):
        rk = rank_of(list(ws))
        S, lam = Smax_of(list(ws))
        c1, c2, c3 = k >= 4, rk == k - 1, S > 2
        if c1 and c2 and c3:
            feasible.append((ws, rk, S))
        print(f"  {k:>2} {str(ws):>16} {rk:>3} {S:>7.4f} {str(c1):>5} {str(c2):>7} {str(c3):>5} "
              f"{str(c1 and c2 and c3):>5}")

print()
ck("L=4 上满足（k≥4 ∧ 秩=k−1 ∧ 语境）的权重构成数 = 0", len(feasible) == 0,
   f"实得 {len(feasible)}")
ck("k=4 的两种权重构成秩均 = 1（不同值个数 3 与 2，取小者−1）",
   rank_of([3, 1, 1, 1]) == 1 and rank_of([2, 2, 1, 1]) == 1,
   f"实得 {rank_of([3,1,1,1])}/{rank_of([2,2,1,1])}")
# ★ 精确机制：秩 = （权重序列中不同值的个数）− 1
def distinct_minus1(ws):
    return len(set(ws)) - 1
_bad = []
for k in range(2, 7):
    for ws in parts(len(W), k):
        if rank_of(list(ws)) != distinct_minus1(list(ws)):
            _bad.append((ws, rank_of(list(ws)), distinct_minus1(list(ws))))
ck("★ 秩 = （不同权重值的个数）− 1", not _bad, f"反例 {_bad}")
_ranks = {}
for k in range(2, 7):
    rs = [rank_of(list(ws)) for ws in parts(len(W), k)]
    _ranks[k] = max(rs) if rs else 0
print(f"  各块数下的最大可达秩：{_ranks}")
ck("★ 故 k ≥ 4 时秩 ≤ 1 < k−1（要求落空）",
   all(_ranks[k] <= 1 for k in (4, 5, 6)), f"实测 {_ranks}")
ck("★ 唯一的秩 2 出现在 k=3 的 (1,2,3)，但 k=3 不满足 k≥4",
   _ranks[3] == 2 and rank_of([1, 2, 3]) == 2)
ck("二者的 S_max 均 < 2（非语境）",
   Smax_of([3, 1, 1, 1])[0] < 2 and Smax_of([2, 2, 1, 1])[0] < 2)
ck("⇒ PI §3.6 出路 (a) 在 L=4 上**走不通**（两条独立理由）", True)

print()
print(f"N6 结果：通过 {len(P)} / 不符 {len(X)}")
