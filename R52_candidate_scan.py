#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R52：候选 π 的穷尽筛选 —— 哪个"圆/高度"读法能堪进 Zero？

全部候选都在**同一批穷尽枚举的闭合 ±1 词**上测，保证可比。
闭合词：长度 L 偶，Σw=0，计数 = C(L, L/2)。

候选层指标 φ(w)（各自诱导一个 π：φ 的值 → 类）：
  C1 进度高度      h_max = k/L 峰值处（进度轴）
  C2 转向总量      T(w) = Σ|Δθ_i|   （方格上的总转角）
  C3 转向缺口      H(w) = T(w) − 2π  （"该转 2π，实际转了多少"）
  C4 浮出高度      max_k |S_k|      （部分和的最高点）
  C5 浮出层数      #不可分浮出      （词分解成几段独立浮出）
  C6 缠绕/净转向   净转向 W(w)      （有符号转向总量）
  C7 长度          L                （对照：已知失败）

三门：门1 q=Σω²∈(0.5,0.6)（R32 生存窗→D=4）
      门2 λ1(顶三归一)>0.723607（KCBS 语境性）
      门3 秩 = k−1（R49 稠密 ⇒ III₁）
"""
from __future__ import annotations

import json
import os
from collections import Counter
from fractions import Fraction
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
LAM = 0.723607
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
          53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]


# ---------------------------------------------------------------- 词生成
def closed_words(L):
    """长度 L 的平衡 ±1 词（Σ=0）。"""
    if L % 2:
        return []
    half = L // 2
    for pos in combinations(range(L), half):
        w = [-1] * L
        for p in pos:
            w[p] = 1
        yield tuple(w)


# ---------------------------------------------------------------- 候选层指标
def theta_steps(w):
    """方格转向：+1=右(0)，−1=上(π/2)。返回每步方向与相邻转角（单位 π/2）。"""
    dirs = [0 if x == 1 else 1 for x in w]          # 0=右, 1=上
    turns = []
    for i in range(len(dirs)):
        d = (dirs[(i + 1) % len(dirs)] - dirs[i]) % 4
        turns.append(0 if d == 0 else (1 if d == 1 else (2 if d == 2 else 3)))
    # 映射到 ±π/2：0→0, 1→+1, 2→2(掉头), 3→−1
    signed = {0: 0, 1: 1, 2: 2, 3: -1}
    return [signed[t] for t in turns]


def c1_progress(w):
    return len(w)                                   # 进度轴：每层被每词穿过一次

def c2_turning(w):
    return sum(abs(t) for t in theta_steps(w))

def c3_gap(w):
    """缺口 = |T − 4|（纯方格的 4 次直角 = 一圈 2π）"""
    return abs(c2_turning(w) - 4)

def c4_excursion_height(w):
    s = 0; mx = 0
    for x in w:
        s += x; mx = max(mx, abs(s))
    return mx

def c5_excursion_count(w):
    """分解成不可分浮出（只在末尾回到 0）的段数"""
    s = 0; cnt = 0
    for i, x in enumerate(w):
        s += x
        if s == 0:
            cnt += 1
    return cnt

def c6_net_turning(w):
    return sum(theta_steps(w))


# ---------------------------------------------------------------- 三门
def exp_vec(k):
    v, kk = [], k
    for p in PRIMES:
        e = 0
        while kk % p == 0:
            kk //= p; e += 1
        v.append(e)
    return None if kk != 1 else v


def exact_rank(vals):
    if len(vals) < 2:
        return 0
    vecs = []
    for i in range(len(vals) - 1):
        a, b = exp_vec(int(vals[i])), exp_vec(int(vals[i + 1]))
        if a is None or b is None:
            return None
        vecs.append([y - x for x, y in zip(a, b)])
    M = [[Fraction(x) for x in r] for r in vecs]
    rank = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(rank, len(M)) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        pv = M[rank][c]
        M[rank] = [x / pv for x in M[rank]]
        for i in range(len(M)):
            if i != rank and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[rank])]
        rank += 1
    return rank


def evaluate(cnt):
    """cnt: {层值: 计数} → 三门"""
    tot = sum(cnt.values())
    if tot == 0:
        return None
    ws = sorted((v / tot for v in cnt.values()), reverse=True)
    q = sum(x * x for x in ws)
    t3 = ws[:3]
    l1 = t3[0] / sum(t3) if len(t3) >= 3 else 1.0
    k = len(ws)
    vals = sorted(cnt.values(), reverse=True)
    rank = exact_rank(vals) if k <= 300 else None
    peak = None
    for D in range(2, 200):
        if D == 2:
            if 3 * q < 1: peak = 2
            continue
        if (D - 1) * q < (D + 1) and D * q > (D - 2):
            peak = D; break
    return {"blocks": k, "total": tot, "q": round(q, 6),
            "lam1": round(l1, 6), "rank": rank, "needed_rank": k - 1,
            "dense": (rank == k - 1) if rank is not None else None,
            "peak_D": peak,
            "gate1_survival": bool(0.5 < q < 0.6),
            "gate2_contextual": bool(l1 > LAM),
            "gate3_rank": (rank == k - 1) if rank is not None else None}


CANDS = [("C1 进度高度 h=k/L", c1_progress),
         ("C2 转向总量 T", c2_turning),
         ("C3 转向缺口 |T-4|", c3_gap),
         ("C4 浮出高度 max|S_k|", c4_excursion_height),
         ("C5 浮出层数", c5_excursion_count),
         ("C6 净转向 W", c6_net_turning)]

if __name__ == "__main__":
    OUT = {}
    for L in [10, 12, 14]:
        words = list(closed_words(L))
        print("=" * 78)
        print("L = %d   闭合词数 = %d" % (L, len(words)))
        print("=" * 78)
        for name, fn in CANDS:
            cnt = Counter(fn(w) for w in words)
            r = evaluate(cnt)
            OUT.setdefault(name, {})[str(L)] = r
            if r is None:
                continue
            print("  %-20s 层数=%-4d q=%-9.6f λ1=%-8.5f 峰D=%-3s "
                  "| 门1=%-5s 门2=%-5s 门3=%s"
                  % (name, r["blocks"], r["q"], r["lam1"], r["peak_D"],
                     r["gate1_survival"], r["gate2_contextual"], r["gate3_rank"]))
        print()
    # 汇总：哪些候选在任一 L 上同时过门1+门2
    print("=" * 78)
    print("汇总：同时过 门1(生存窗) 与 门2(语境性) 的候选")
    print("=" * 78)
    winners = []
    for name, byL in OUT.items():
        for L, r in byL.items():
            if r and r["gate1_survival"] and r["gate2_contextual"]:
                winners.append((name, L, r))
                print("  ✅ %s (L=%s): q=%.6f λ1=%.5f 秩=%s/%d"
                      % (name, L, r["q"], r["lam1"], r["rank"], r["needed_rank"]))
    if not winners:
        print("  （无）")
    # L=14 的层值分布（最有信息量）
    print()
    print("=" * 78)
    print("L=14 各候选的层值分布（前 12 项）")
    print("=" * 78)
    words = list(closed_words(14))
    for name, fn in CANDS:
        cnt = Counter(fn(w) for w in words)
        top = sorted(cnt.items(), key=lambda kv: -kv[1])[:12]
        print("  %-20s %s" % (name, top))
    OUT["_winners"] = [(n, L) for n, L, _ in winners]
    with open(os.path.join(HERE, "R52_candidate_scan_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R52_candidate_scan_results.json")
