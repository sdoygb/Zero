#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R51 探针：从原生 wipe_reseed 目标（L0/L2）抽出记录层（L1）的闭合词分布 π

直接复用 zero_sum_periodic_destruction.py 的演化与重播种规则：
    active: 每站点一组词 (tuple of ±1)
    一步：每个词派生两个子词 w±1；sum==0 者入 history 与 zero_layer，否则留 active
    重播种：每个 history 词 w 派生 w±1（非闭合者）作为新 active
即：**活动层全清 ⇒ 记录层（不可逆）留在 L1 ⇒ 从 L1 重播种**。

本探针测量两把尺子（R44 / R49）：
  G1  q = Σ ω_c²  是否落在 R32 生存窗口 (1/2, 3/5)
  G2  顶三归一 λ1 是否 > 0.723607（KCBS 语境性）
  G3  块权重尺度分布与素数差向量秩（R49 稠密判据）
"""
from __future__ import annotations

import json
import os
from collections import Counter
from fractions import Fraction
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
LAM_STAR = 0.723607
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
          53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]
OUT = {}


def word_from_string(text):
    return tuple(1 if c == "+" else -1 for c in text)


def is_closed(w):
    return bool(w) and sum(w) == 0


def canon(word):
    """闭合词的旋转等价类（R12 的 canonical cycle 口径）"""
    n = len(word)
    return min(word[i:] + word[:i] for i in range(n))


def evolve(active, histories, zero_layer):
    nxt = [set() for _ in active]
    closed = 0
    for site, paths in enumerate(active):
        for w in paths:
            for step in (1, -1):
                c = w + (step,)
                if is_closed(c):
                    histories[site].add(c)
                    zero_layer.add(canon(c))
                    closed += 1
                else:
                    nxt[site].add(c)
    return nxt, closed


def reseed(histories):
    active = [set() for _ in histories]
    for site, book in enumerate(histories):
        for w in book:
            for step in (1, -1):
                s = w + (step,)
                if not is_closed(s):
                    active[site].add(s)
    return active


def census(histories):
    """把 L1 的精确闭合词按 (站点无关的) 同构类聚合。

    三种读出候选（ω 的来源不同，必须并列报，见 R50 #5）：
      A 每个 distinct exact word 一次                    -> ω ∝ 1
      B 每个 distinct canonical cycle 一次                -> ω ∝ 1
      C 按精确词的出现重数（多重集，含跨站点重数）        -> ω ∝ count
    """
    exact_multi = Counter()
    exact_sets = Counter()
    for book in histories:
        for w in book:
            exact_multi[w] += 1
        for w in book:
            exact_sets[w] = 1
    cyc = Counter()
    for w in exact_multi:
        cyc[canon(w)] += 1
    return exact_multi, exact_sets, cyc


def q_of(counts):
    tot = sum(counts.values())
    if tot == 0:
        return None, 0
    return sum((c / tot) ** 2 for c in counts.values()), tot


def top3_l1(counts):
    tot = sum(counts.values())
    if tot == 0:
        return None
    s = sorted(counts.values(), reverse=True)[:3]
    return max(s) / sum(s)


def exp_vec_int(k):
    v, kk = [], k
    for p in PRIMES:
        e = 0
        while kk % p == 0:
            kk //= p
            e += 1
        v.append(e)
    return None if kk != 1 else v


def exact_rank(rows):
    M = [[Fraction(x) for x in r] for r in rows]
    if not M:
        return 0
    rank, ncol = 0, len(M[0])
    for c in range(ncol):
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
        if rank == len(M):
            break
    return rank


def rank_of_counts(counts):
    """R49 稠密判据：块权重 = 计数；秩 = k-1 ⇔ III_1。"""
    vals = sorted(counts.values(), reverse=True)
    if len(vals) < 2:
        return None, len(vals)
    nums = [int(v) for v in vals]
    vecs = []
    for i in range(len(nums) - 1):
        a, b = exp_vec_int(nums[i]), exp_vec_int(nums[i + 1])
        if a is None or b is None:
            return None, len(vals)
        vecs.append([y - x for x, y in zip(a, b)])
    return exact_rank(vecs), len(vals)


def run(generations, period=3, seed_words=("+", "-", "++", "--")):
    active = [{word_from_string(t)} for t in seed_words]
    histories = [set() for _ in seed_words]
    zero = set()
    trace = []
    for g in range(1, generations + 1):
        nxt, closed = evolve(active, histories, zero)
        active = nxt
        if g % period == 0:                       # 代际毁灭 + 重播种
            active = reseed(histories)
        em, es, cy = census(histories)
        trace.append({
            "gen": g,
            "exact_words": len(es),
            "canon_cycles": len(cy),
            "q_worduniform": round(q_of(es)[0], 6) if es else None,
            "q_cycleuniform": round(q_of(cy)[0], 6) if cy else None,
            "lam1_word": round(top3_l1(es), 5) if es else None,
            "lam1_cycle": round(top3_l1(cy), 5) if cy else None,
        })
    return histories, zero, trace


if __name__ == "__main__":
    N = int(os.environ.get("GENS", "9"))
    print("原生 wipe_reseed：生成 %d 代 ..." % N)
    hist, zero, trace = run(N)
    em, es, cy = census(hist)
    OUT["trace"] = trace
    OUT["final"] = {"distinct_exact_words": len(es),
                    "distinct_canon_cycles": len(cy),
                    "zero_layer": len(zero)}
    for nm, cnt in [("A_distinct_word_uniform", es),
                    ("B_distinct_cycle_uniform", cy),
                    ("C_exact_multiplicity", em)]:
        q, tot = q_of(cnt)
        l1 = top3_l1(cnt)
        r, k = rank_of_counts(cnt) if tot <= 4000 else (None, len(cnt))
        OUT[nm] = {"blocks": k, "total": tot,
                   "q": round(q, 6) if q else None,
                   "top3_l1": round(l1, 5) if l1 else None,
                   "rank": r, "dense_rank_needed": k - 1,
                   "in_survival_window": bool(q and 0.5 < q < 0.6) if q else False,
                   "contextual": bool(l1 and l1 > LAM_STAR)}
    print("\n--- 逐代 ---")
    print("  gen  exact  cycles   q(word)  q(cycle)  λ1(word) λ1(cycle)")
    for t in trace:
        print("  %3d  %5d  %6d   %-8s %-9s %-8s %s"
              % (t["gen"], t["exact_words"], t["canon_cycles"],
                 t["q_worduniform"], t["q_cycleuniform"],
                 t["lam1_word"], t["lam1_cycle"]))
    print("\n--- 末代读出（三种候选）---")
    for nm in ["A_distinct_word_uniform", "B_distinct_cycle_uniform",
               "C_exact_multiplicity"]:
        d = OUT[nm]
        print("  %-26s 块=%-6d q=%-9s λ1=%-8s 秩=%-6s 需要秩=%-5d 生存窗=%s 语境=%s"
              % (nm, d["blocks"], d["q"], d["top3_l1"], d["rank"],
                 d["dense_rank_needed"], d["in_survival_window"], d["contextual"]))
    with open(os.path.join(HERE, "R51_native_pi_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("\n-> R51_native_pi_results.json")
