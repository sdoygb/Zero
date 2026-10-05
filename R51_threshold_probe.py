#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R51 探针 IV：(R51-1) 精确 k 块语境性门限 + (R51-2) 代际修剪能否救回集中度

(R51-1) 设 k 块（降序）w，T=顶三之和，λ1=w1/T。则
        λ1 ≤ ( (k-1) - sqrt((k-1)(k-2)(1-q_k)) ) / k ,  q_k=Σw_i²
   等号在 前三个之外的所有块都等于 w_k 时达到。
   推论：语境性（λ1>λ*=0.723607）要求
        q_k ≥ q*_k = (1 + (k-1)λ*²)/k - ((k-1)(k-2)/k)·(λ*²/(1-λ*))
   k=3 -> 0.6（= R44）；k=4 -> 0.7753；k=9 -> 0.5962；k=10 -> 0.5916；k>=11 -> <窗口
"""
from __future__ import annotations

import json
import math
import os
from collections import Counter
from math import log, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LAM = 0.723607
OUT = {}


def lam1_bound(k, q):
    v = (k - 1) * (k - 2) * (1 - q)
    if v < 0:
        return float("nan")
    return ((k - 1) - sqrt(v)) / k


def qstar(k):
    num = 1 + (k - 1) * LAM ** 2
    corr = (k - 1) * (k - 2) * (LAM ** 2 / (1 - LAM))
    return (num - corr) / k


def peak_dim(q):
    best, bd = -1, None
    for D in range(2, 400):
        F = D * (D - 1) / 2 * q ** D
        if F > best:
            best, bd = F, D
    return bd


def check_bound(n=900000, seed=23):
    rng = np.random.default_rng(seed)
    worst = {}
    for _ in range(n):
        k = int(rng.integers(3, 16))
        w = np.sort(rng.random(k) ** rng.choice([0.4, 1.0, 2.5, 5.0]))[::-1]
        w = w / w.sum()
        q = float(np.sum(w ** 2))
        t3 = w[:3]
        l1 = float(t3[0] / t3.sum())
        b = lam1_bound(k, q)
        exc = l1 - b
        if exc > worst.get(k, (-9, None))[0]:
            worst[k] = (exc, {"q": round(q, 5), "lam1": round(l1, 5),
                              "bound": round(b, 5)})
    return {"per_k": {str(k): {"excess": round(v[0], 9), **v[1]}
                      for k, v in sorted(worst.items())},
            "max_excess": max(v[0] for v in worst.values())}


def check_tight(k, q):
    """等号构造：w1 由 λ*(q) 定，其余 k-1 块等权。"""
    b = lam1_bound(k, q)
    T = (k - 1) * sqrt((1 - q) / ((k - 1) * (k - 2))) if k > 2 else None
    if T is None:
        return None
    w1 = b * T
    rest = (1 - w1) / (k - 1)
    w = np.array([w1] + [rest] * (k - 1))
    if w.min() < -1e-12:
        return None
    l1 = w[0] / w[:3].sum()
    return {"k": k, "q": q, "bound": round(b, 6), "constructed": round(l1, 6),
            "q_ach": round(float(np.sum(w ** 2)), 6),
            "w": [round(float(x), 5) for x in w]}


# ---------------------------------------------------------------- 修剪测试
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


def reseed(hist):
    out = [set() for _ in hist]
    for site, book in enumerate(hist):
        for w in book:
            for step in (1, -1):
                s = w + (step,)
                if sum(s) != 0:
                    out[site].add(s)
    return out


def prune(hist, keep_top):
    """D222 口径：毁灭时每个 P_i 只保留"最高"若干亚层（按长度）。"""
    out = []
    for book in hist:
        if not book:
            out.append(set())
            continue
        lens = sorted({len(w) for w in book}, reverse=True)[:keep_top]
        out.append({w for w in book if len(w) in lens})
    return out


def run_pruned(gens, keep_top, period=3,
               seeds=("+", "-", "++", "--", "+++", "---")):
    active = [{wfs(t)} for t in seeds]
    hist = [set() for _ in seeds]
    trace = []
    for g in range(1, gens + 1):
        active = evolve(active, hist)
        if g % period == 0:
            active = reseed(hist)
            hist = prune(hist, keep_top)      # 只有在代际毁灭时才修剪
        cnt = Counter()
        for book in hist:
            for w in book:
                cnt[w] += 1
        tot = sum(cnt.values())
        if tot == 0:
            break
        ws = sorted(cnt.values(), reverse=True)
        q = sum((x / tot) ** 2 for x in ws)
        l1 = ws[0] / sum(ws[:3]) if len(ws) >= 3 else 1.0
        trace.append({"gen": g, "records": tot, "blocks": len(ws),
                      "q": round(q, 6), "lam1": round(l1, 5)})
    return trace


if __name__ == "__main__":
    print("(R51-1) 界检验 ...")
    OUT["bound_check"] = check_bound()
    print("  最大超出 = %.3e" % OUT["bound_check"]["max_excess"])
    for k, d in OUT["bound_check"]["per_k"].items():
        print("   k=%-3s 超出=%+.2e  q=%.4f λ1=%.5f 界=%.5f"
              % (k, d["excess"], d["q"], d["lam1"], d["bound"]))

    print("\n等号构造 ...")
    OUT["tight"] = [x for x in (check_tight(3, 0.6), check_tight(4, 0.7753),
                                check_tight(10, 0.5916), check_tight(6, 0.70))
                    if x]
    for r in OUT["tight"]:
        print("  k=%-3d q=%.4f 界=%.6f 构造值=%.6f (实际q=%.4f)"
              % (r["k"], r["q"], r["bound"], r["constructed"], r["q_ach"]))

    print("\n(R51-2) 门限表 ...")
    rows = []
    for k in range(3, 16):
        qs = qstar(k)
        rows.append({"k": k, "q_star": round(qs, 5),
                     "feasible": qs <= 1.0,
                     "peak_D": peak_dim(min(qs, 0.99999)) if qs <= 1 else None,
                     "peak_D_is_4": (peak_dim(min(qs, 0.99999)) == 4) if qs <= 1 else False})
    OUT["thresholds"] = rows
    for r in rows:
        print("  k=%-3d 语境性要求 q >= %-8.4f 可行=%-5s 账本峰 D=%s %s"
              % (r["k"], r["q_star"], r["feasible"], r["peak_D"],
                 "<== 峰在4" if r["peak_D_is_4"] else ""))
    best = [r for r in rows if r["feasible"] and r["peak_D_is_4"]]
    OUT["viable_k"] = [r["k"] for r in best]
    print("  ⟹ 能同时满足语境性 + 四维峰的块数:", OUT["viable_k"])

    print("\n(R51-2) 代际修剪测试 ...")
    OUT["pruning"] = {}
    for keep in [1, 2, 3, 5, 10, 20, None]:
        tr = run_pruned(24, keep if keep else 999)
        OUT["pruning"][str(keep)] = tr
        if tr:
            f = tr[-1]
            print("  keep_top=%-5s 末代: 记录=%-6d 块=%-6d q=%-9.6f λ1=%.5f 生存窗=%s"
                  % (keep, f["records"], f["blocks"], f["q"], f["lam1"],
                     0.5 < f["q"] < 0.6))
    with open(os.path.join(HERE, "R51_threshold_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("-> R51_threshold_results.json")
