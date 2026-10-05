#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R51 探针 VI：记录层保留深度（D222 的 Top_2 概括分层）能否把 π 推向两道门？

D222：P_i 有"精确→概括"的亚层，毁灭时只保留**最概括的两层**，重播种只读保留层。
本探针把"保留深度"当作可调旋钮：毁灭时每个站点的历史只保留最概括的 keep 层，
其中"概括层级"按**长度桶**实现（越长的词 = 越精确，越短 = 越概括），并测：

  q = Σω²，λ1（顶三归一），有效类数，R49 秩 —— 随 keep 的变化

若某 keep 使 q∈(0.5,0.6) 且 λ1>0.7236，则"重播种影响粗粒化"确实能救回 π。
"""
from __future__ import annotations

import json
import os
from collections import Counter, defaultdict
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
LAM = 0.723607
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
          53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]
OUT = {}


def wfs(t): return tuple(1 if c == "+" else -1 for c in t)


def evolve(active, hist):
    nxt = [set() for _ in active]
    for site, paths in enumerate(active):
        for w in paths:
            for step in (1, -1):
                c = w + (step,)
                if sum(c) == 0:
                    hist[site].add(c)
                else:
                    nxt[site].add(c)
    return nxt


def keep_top_layers(hist, keep):
    """D222 口径：只保留最概括的 keep 层。
    以**长度桶**表示概括层级：长度越长越精确；保留 keep 个最长的长度桶。"""
    out = []
    for book in hist:
        if not book:
            out.append(set())
            continue
        lens = sorted({len(w) for w in book}, reverse=True)[:keep]
        out.append({w for w in book if len(w) in lens})
    return out


def reseed(hist):
    out = [set() for _ in hist]
    for site, book in enumerate(hist):
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


def exact_rank(vals):
    if len(vals) < 2:
        return None
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


def run(gens, keep, period=3, seeds=("+", "-", "++", "--", "+++", "---")):
    active = [{wfs(t)} for t in seeds]
    hist = [set() for _ in seeds]
    trace = []
    for g in range(1, gens + 1):
        active = evolve(active, hist)
        if g % period == 0:
            if keep is not None:
                hist = keep_top_layers(hist, keep)   # 先投影（D222）
            active = reseed(hist)                    # 重播种只读保留层
        cnt = Counter()
        for book in hist:
            for w in book:
                cnt[w] += 1
        if not cnt:
            break
        tot = sum(cnt.values())
        ws = sorted(cnt.values(), reverse=True)
        q = sum((x / tot) ** 2 for x in ws)
        l1 = ws[0] / sum(ws[:3]) if len(ws) >= 3 else 1.0
        trace.append({"gen": g, "classes": len(ws), "mass": tot,
                      "q": round(q, 6), "lam1": round(l1, 5)})
    return hist, trace


if __name__ == "__main__":
    print("保留深度 keep 扫描（D222 的 Top_keep 概括分层）")
    print(" keep   gen  类数    q          λ1        q∈(0.5,0.6)  λ1>0.7236")
    out = {}
    for keep in [1, 2, 3, 4, 6, 10, None]:
        hist, trace = run(24, keep)
        if not trace:
            continue
        f = trace[-1]
        out[str(keep)] = {"final": f, "trace": trace}
        print("  %-5s %3d  %-6d %.6f   %.5f   %-11s %s"
              % (keep, f["gen"], f["classes"], f["q"], f["lam1"],
                 0.5 < f["q"] < 0.6, f["lam1"] > LAM))
    # 找全局最优（最接近两道门）
    best = None
    for k, d in out.items():
        f = d["final"]
        score = (min(f["q"], 0.6) - 0.5) if f["q"] < 0.6 else -(f["q"] - 0.6)
        if best is None or score > best[0]:
            best = (score, k, f)
    print("\n最接近生存窗的 keep = %s : q=%.6f λ1=%.5f" % (best[1], best[2]["q"], best[2]["lam1"]))
    OUT["scan"] = {k: v["final"] for k, v in out.items()}
    with open(os.path.join(HERE, "R51_retention_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("-> R51_retention_results.json")
