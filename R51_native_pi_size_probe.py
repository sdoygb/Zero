#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R51 探针 II：把 L1 的闭合词按**长度类**聚合（= "主导块 + 幂律尾" 里的"类"）

上一探针按"每个 distinct 词/循环各算一块"读出，三种候选全部失败
（q→0，λ1→1/3）。本探针换读出：块 = 闭合词长度 2m 的等价类，
ω_m = 该类在 L1 中的质量。这才是"主导块 + 幂律尾"的自然含义。

测：q = Σ ω_m²；λ1（顶三归一）；长度直方图是否幂律；R49 秩（块权重 = 整数计数）。
"""
from __future__ import annotations

import json
import os
from collections import Counter
from fractions import Fraction
from math import log

HERE = os.path.dirname(os.path.abspath(__file__))
LAM_STAR = 0.723607
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
          53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109]
OUT = {}


def wfs(t):
    return tuple(1 if c == "+" else -1 for c in t)


def evolve(active, histories, zero):
    nxt = [set() for _ in active]
    for site, paths in enumerate(active):
        for w in paths:
            for step in (1, -1):
                c = w + (step,)
                if sum(c) == 0:
                    histories[site].add(c)
                    zero.add(min(c[i:] + c[:i] for i in range(len(c))))
                else:
                    nxt[site].add(c)
    return nxt


def reseed(histories):
    out = [set() for _ in histories]
    for site, book in enumerate(histories):
        for w in book:
            for step in (1, -1):
                s = w + (step,)
                if sum(s) != 0:
                    out[site].add(s)
    return out


def exp_vec(k):
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


def run(gens, period=3, seeds=("+", "-", "++", "--", "+++", "---")):
    active = [{wfs(t)} for t in seeds]
    hist = [set() for _ in seeds]
    zero = set()
    trace = []
    for g in range(1, gens + 1):
        active = evolve(active, hist, zero)
        if g % period == 0:
            active = reseed(hist)
        size_hist = Counter()
        for book in hist:
            for w in book:
                size_hist[len(w)] += 1
        tot = sum(size_hist.values())
        if tot:
            ws = sorted(size_hist.values(), reverse=True)
            q = sum((x / tot) ** 2 for x in ws)
            l1 = ws[0] / sum(ws[:3]) if len(ws) >= 3 else (
                ws[0] / sum(ws) if ws else None)
            trace.append({"gen": g, "words": tot, "size_classes": len(ws),
                          "q": round(q, 6), "lam1": round(l1, 5),
                          "top_share": round(ws[0] / tot, 5),
                          "hist": dict(sorted(size_hist.items()))})
    return hist, zero, trace


if __name__ == "__main__":
    G = int(os.environ.get("GENS", "10"))
    print("按长度类读出：%d 代 ..." % G)
    hist, zero, trace = run(G)
    for t in trace:
        print("  gen=%2d 词=%-6d 类=%-3d q=%.6f λ1=%.5f 首类占比=%.4f"
              % (t["gen"], t["words"], t["size_classes"], t["q"], t["lam1"],
                 t["top_share"]))
    last = trace[-1]
    print("\n末代长度直方图:", last["hist"])
    hs = last["hist"]
    ks = sorted(hs, key=lambda k: hs[k], reverse=True)
    print("最大类 m =", ks[0], " 计数 =", hs[ks[0]])
    # 幂律拟合（log-log）
    xs = [log(k) for k in ks if hs[k] > 0]
    ys = [log(hs[k]) for k in ks if hs[k] > 0]
    n = len(xs)
    if n >= 3:
        mx, my = sum(xs) / n, sum(ys) / n
        num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
        den = sum((a - mx) ** 2 for a in xs)
        slope = num / den if den else None
        print("log-log 斜率（幂律指数估计）= %.3f" % slope)
        OUT["powerlaw_slope"] = round(slope, 4)
    cnts = [hs[k] for k in ks]
    OUT["final"] = {"words": last["words"], "size_classes": last["size_classes"],
                    "q": last["q"], "lam1": last["lam1"],
                    "hist": last["hist"]}
    OUT["trace"] = [{k: v for k, v in t.items() if k != "hist"} for t in trace]
    OUT["gates"] = {
        "q_in_survival_window_0.5_0.6": bool(0.5 < last["q"] < 0.6),
        "contextual": bool(last["lam1"] > LAM_STAR),
        "rank": None, "needed_rank": last["size_classes"] - 1,
    }
    if len(cnts) <= 400:
        vecs = []
        ok = True
        for i in range(len(cnts) - 1):
            a, b = exp_vec(cnts[i]), exp_vec(cnts[i + 1])
            if a is None or b is None:
                ok = False
                break
            vecs.append([y - x for x, y in zip(a, b)])
        if ok:
            OUT["gates"]["rank"] = exact_rank(vecs)
    print("\n门 1 q∈(0.5,0.6) :", OUT["gates"]["q_in_survival_window_0.5_0.6"])
    print("门 2 λ1>0.7236   :", OUT["gates"]["contextual"])
    print("门 3 秩 = k-1    :", OUT["gates"]["rank"], "/",
          OUT["gates"]["needed_rank"])
    with open(os.path.join(HERE, "R51_native_pi_size_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("-> R51_native_pi_size_results.json")
