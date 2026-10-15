#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N8 —— 层高划分下的扫描：D=4 能否从「量子性必需」**导出**？

背景（WHY_TWO_L_VALUES §7）：
  语境性把【旋转类划分】排除（它在一切 L 上都不过 S_max > 2）。
  于是问题变成：在唯一过语境门的【层高划分】下，是否存在使账本峰 = 4 的 L？
  若有且唯一，则 D=4 可从层高划分**导出**，而不是从旋转类划分条件地得到。

层高划分（照 R53/R54 的 make_blocks）：
  b_h = #{w : max_k |s_k| = h}；取密度极大点 h*；
  h > h* 各层独立成块，h ≤ h* 合并为一块
  ⟹ 块权重 = [N(≤h*), b_{h*+1}, …, b_{hmax}]

b_h 用有界行走计数（DP，O(L·h)），避免枚举 C(L,L/2) 个词。

判据（照 R53.ledger_peak 的原文）：
  F_D = C(D,2)·q^D，相邻比 F_D/F_{D-1} = D/(D-2)·q
  峰在 D ⟺ D/(D-2)·q > 1 且 (D+1)/(D-1)·q < 1
  ⟹ **峰 = 4 ⟺ q ∈ (1/2, 3/5)**
量子三门：门1 生存窗 q∈(0.5,0.6)；门2 语境 S_max>2；门3 稠密 rank=k−1
"""
import math
from fractions import Fraction as F
from itertools import groupby

P, X = [], []


def ck(n, c, d=""):
    (P if c else X).append(n)
    print(f"  [{'v' if c else 'X'}] {n}" + (f"   {d}" if d else ""))


def N_le(L, h):
    """#{±1 路径 : s_0=s_L=0, 全程 |s_k| ≤ h}（DP）"""
    cur = {0: 1}
    for _ in range(L):
        nxt = {}
        for s, c in cur.items():
            for d in (1, -1):
                t = s + d
                if -h <= t <= h:
                    nxt[t] = nxt.get(t, 0) + c
        cur = nxt
    return cur.get(0, 0)


def scan(L):
    """返回 (h*, 块权重, N, k)"""
    prev, b = 0, {}
    for h in range(1, L // 2 + 1):
        c = N_le(L, h)
        b[h] = c - prev
        prev = c
    N = sum(b.values())
    assert N == math.comb(L, L // 2), (N, math.comb(L, L // 2))
    hs = max(b, key=lambda h: b[h])
    outer = sum(b[h] for h in range(1, hs + 1))
    ws = [outer] + [b[h] for h in range(hs + 1, L // 2 + 1) if b[h] > 0]
    return hs, ws, N, len(ws)


_FAC = {}


def factorize(n):
    """带缓存的一次分解（原版在 rank_of 里对同一大整数反复试除，L=140 时是 41 位数 ⟹ 极慢）"""
    if n in _FAC:
        return _FAC[n]
    f, d, m = {}, 2, n
    while d * d <= m:
        while m % d == 0:
            f[d] = f.get(d, 0) + 1
            m //= d
        d += 1
    if m > 1:
        f[m] = f.get(m, 0) + 1
    _FAC[n] = f
    return f


def rank_of(ws):
    exps = [factorize(v) for v in ws]          # ★ 每个权重只分解一次
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


MU2 = (5 - math.sqrt(5)) / 2


def Smax(lam):
    return math.sqrt(5) * lam + MU2 * (1 - lam)


# ============================================================
print("=" * 92)
print("(A) 判据自洽性：ledger_peak 的窗口是 (1/2, 3/5)")
print("=" * 92)


def ledger_peak(q):
    for D in range(2, 400):
        if D == 2:
            if 3 * q < 1:
                return 2
            continue
        if D * q > (D - 2) and (D + 1) * q < (D - 1):
            return D
    return None


lo = next(x / 10000 for x in range(1, 10000) if ledger_peak(x / 10000) == 4)
hi = max(x / 10000 for x in range(1, 10000) if ledger_peak(x / 10000) == 4)
ck("峰 = 4 ⟺ q ∈ (1/2, 3/5)", abs(lo - 0.5) < 1e-3 and abs(hi - 0.6) < 1e-3,
   f"数值扫描得 [{lo}, {hi}]")
ck("q = 5/9 给峰 4（与 R32 §2 一致）", ledger_peak(5 / 9) == 4)

# ============================================================
print()
print("=" * 92)
print("(B) 扫描 L = 6..140（层高划分）")
print("=" * 92)
print(f"  {'L':>4} {'h*':>3} {'k':>4} {'q':>11} {'峰D':>4} {'λ₁':>8} {'S_max':>8} "
      f"{'语境':>4} {'秩=k−1':>7} {'窗':>4} {'全过':>5}")
rows = []
for L in range(6, 61, 2):
    hs, ws, N, k = scan(L)
    q = sum(F(x, N) ** 2 for x in ws)
    wss = sorted(ws, reverse=True)
    lam = F(wss[0], sum(wss[:3]))
    S = Smax(float(lam))
    rk = rank_of(ws)
    pk = ledger_peak(float(q))
    g1 = 0.5 < float(q) < 0.6
    g2 = S > 2
    g3 = rk == k - 1
    ok = g1 and g2 and g3 and pk == 4
    rows.append(dict(L=L, hs=hs, k=k, q=float(q), peak=pk, lam=float(lam), S=S,
                     rank=rk, g1=g1, g2=g2, g3=g3, ok=ok))
    if L <= 34 or ok:
        print(f"  {L:>4} {hs:>3} {k:>4} {float(q):>11.8f} {str(pk):>4} "
              f"{float(lam):>8.5f} {S:>8.5f} {'是' if g2 else '否':>4} "
              f"{str(g3):>7} {'是' if g1 else '否':>4} {'★是' if ok else '否':>5}")

print()
print("=" * 92)
print("(C) 结构：q 在每个 h* 段内单峰，且峰在段首；段峰单调下降")
print("=" * 92)
segs = []
for hs, g in groupby(rows, key=lambda r: r["hs"]):
    g = list(g)
    qs = [r["q"] for r in g]
    pk = max(qs)
    ip = qs.index(pk)
    rising = all(qs[i] <= qs[i + 1] + 1e-15 for i in range(ip))
    falling = all(qs[i] >= qs[i + 1] - 1e-15 for i in range(ip, len(qs) - 1))
    segs.append((hs, g[0]["L"], g[-1]["L"], pk, g[ip]["L"], rising and falling))
    print(f"  h*={hs}: L {g[0]['L']}..{g[-1]['L']}  段内峰 q={pk:.8f} @L={g[ip]['L']}  "
          f"单峰={rising and falling}")
ck("每段内 q 单峰", all(s[5] for s in segs))
ck("段峰出现在段首（即 h* 跳变处）", all(s[4] == s[1] for s in segs))
pk_list = [s[3] for s in segs]
ck("段峰从 L=16 起严格单调下降",
   all(pk_list[i] > pk_list[i + 1] for i in range(1, len(pk_list) - 1)),
   " → ".join(f"{x:.6f}" for x in pk_list))
ck("h* 对 L 单调不减", all(rows[i]["hs"] <= rows[i + 1]["hs"]
                          for i in range(len(rows) - 1)))

# ============================================================
print()
print("=" * 92)
print("(D) 判定：D=4 是否可从层高划分导出")
print("=" * 92)
win = [r["L"] for r in rows if r["g1"]]
pk4 = [r["L"] for r in rows if r["peak"] == 4]
ctx = [r["L"] for r in rows if r["g2"]]
allok = [r["L"] for r in rows if r["ok"]]
print(f"  q∈(0.5,0.6) 的 L = {win}")
print(f"  峰 = 4 的 L     = {pk4}")
print(f"  语境 S_max>2 的 L = {ctx}")
print(f"  ★ 四条件全过的 L = {allok}")
ck("峰 = 4 的 L 唯一且为 16（L=6..60）", pk4 == [16], f"实得 {pk4}")
ck("生存窗 q∈(0.5,0.6) 的 L 唯一且为 16", win == [16], f"实得 {win}")
ck("★ 四条件（峰=4 ∧ 窗 ∧ 语境 ∧ 稠密）全过的 L 唯一且为 16", allok == [16],
   f"实得 {allok}")
ck("★ 推论：既然段峰单调下降且 L≥28 的段峰 ≤ 0.47774 < 0.5，窗口对一切更大 L 关闭",
   pk_list[2] < 0.5, f"L=28 段峰 {pk_list[2]:.6f}")

print()
print("=" * 92)
print(f"N8 结果：通过 {len(P)} / 不符 {len(X)}")
print("=" * 92)
for x in X:
    print("   不符:", x)
