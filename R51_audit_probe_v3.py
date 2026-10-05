#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R51 审计探针 v3：把 R44 的界推广到 k 块，并检验 R43 两层族

待验命题（本文推导，任务：独立复算 + 反例搜索）
  设 k 块归一谱 w（降序），q_k = Σ w_i²，
  顶三归一 λ1 = w1/T, T = w1+w2+w3（KCBS 3 维归约的第一分量）。
  R44 对 k=3 给出 λ1 ≤ (1+√(2q-1))/2。
  推广 (R51-1)：  λ1 ≤ (1 + √((k-1)(k q_k - 1))) / k ,  等号当
        w2 = w3 = ... = w_{k-1} = (T-w1)/(k-2),  w_k = 1-T.
  推论 (R51-2)：  语境性 (λ1 > λ* = 0.723607) ⇒ q_k ≥ q_k* = (1+(k-1)λ*²)/k
        k=3 -> 0.6 ;  k=4 -> 0.9518 ;  k=5 -> 1.0388 > 1 ;  k≥5 -> 空
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LAM_STAR = 0.723607
MU = np.array([math.sqrt(5), (5 - math.sqrt(5)) / 2, (5 - math.sqrt(5)) / 2])
OUT = {}


def top3_l1(w):
    s = np.sort(np.asarray(w, float))[::-1][:3]
    return float(s[0] / s.sum())


def qk(w):
    return float(np.sum(np.asarray(w, float) ** 2))


def bound(k, q):
    v = (k - 1) * (k * q - 1)
    return (1 + math.sqrt(v)) / k if v >= 0 else float("nan")


# ---------------------------------------------------------------- 1
def check_bound(n=1500000, seed=17):
    """随机搜索：所有 k 上 λ1 是否 ≤ bound(k,q_k)。"""
    rng = np.random.default_rng(seed)
    worst = {}
    for _ in range(n):
        k = int(rng.integers(3, 13))
        w = rng.random(k) ** rng.choice([0.3, 1.0, 3.0, 6.0])
        w = np.sort(w / w.sum())[::-1]
        q, l1, b = qk(w), top3_l1(w), bound(k, qk(w))
        exc = l1 - b
        if exc > worst.get(k, (-9, None))[0]:
            worst[k] = (exc, {"q": round(q, 6), "lam1": round(l1, 6),
                              "bound": round(b, 6), "w": [round(x, 5) for x in w]})
    OUT["1_bound_check"] = {
        "claim": "lam1 <= (1+sqrt((k-1)(k q -1)))/k",
        "per_k": {str(k): {"max_excess": round(v[0], 9), **v[1]}
                  for k, v in sorted(worst.items())},
        "global_max_excess": max(v[0] for v in worst.values()),
    }


# ---------------------------------------------------------------- 2
def check_tightness(k, q, seed=19):
    """在给定 (k,q) 上做局部优化，看能否达到等号。"""
    rng = np.random.default_rng(seed)
    # 用 R51-1 的参数化解作为起点
    lam = bound(k, q)
    if not (0 < lam < 1):
        return None
    T = math.sqrt((k - 1) / (k * q - 1)) if k * q > 1 else None
    if T is None:
        return None
    w1 = lam * T
    rest = (T - w1) / (k - 2) if k > 2 else 0.0
    w = np.array([w1] + [rest] * (k - 2) + [1 - T])
    if w.min() < -1e-12:
        return None
    w = np.sort(np.maximum(w, 0))[::-1]
    w = w / w.sum()
    got = top3_l1(w)
    # 局部细化
    for _ in range(20000):
        i, j = rng.integers(0, k, 2)
        if i == j:
            continue
        c = w.copy()
        s = rng.normal(0, 0.002)
        c[i] = max(c[i] + s, 0)
        c[j] = max(c[j] - s, 0)
        c = np.sort(c)[::-1]
        c = c / c.sum()
        if abs(qk(c) - q) > abs(qk(w) - q):
            continue
        if top3_l1(c) > top3_l1(w):
            w = c
    return {"k": k, "q_target": q, "analytic_bound": round(bound(k, q), 6),
            "optimized": round(top3_l1(w), 6), "q_achieved": round(qk(w), 6),
            "w": [round(float(x), 6) for x in w],
            "construct_w": [round(float(x), 6) for x in w]}


def threshold_table():
    rows = []
    for k in range(3, 10):
        qs = (1 + (k - 1) * LAM_STAR ** 2) / k
        rows.append({"k": k, "q_k_star": round(qs, 6),
                     "feasible": qs <= 1.0,
                     "ledger_peak_D": peak_dim(min(qs, 0.9999)) if qs <= 1 else None})
    return rows


def peak_dim(q):
    best, bd = -1, None
    for D in range(2, 200):
        F = D * (D - 1) / 2 * q ** D
        if F > best:
            best, bd = F, D
    return bd


# ---------------------------------------------------------------- 3
def two_layer(p, alpha, K):
    tail = np.array([a ** (-alpha) for a in range(1, K + 1)])
    w = np.concatenate([[p], (1 - p) * tail / tail.sum()])
    return np.sort(w)[::-1]


def check_two_layer():
    rows = []
    for p in [0.75, 0.80, 0.85, 0.90, 0.95]:
        for alpha in [1.0, 1.5, 2.0, 3.0]:
            for K in [8, 32, 128, 512]:
                w = two_layer(p, alpha, K)
                k = len(w)
                q, l1 = qk(w), top3_l1(w)
                b = bound(k, q)
                rows.append({"p": p, "alpha": alpha, "K": K, "k": k,
                             "q": round(q, 5), "lam1": round(l1, 5),
                             "bound": round(b, 5),
                             "contextual": l1 > LAM_STAR,
                             "S_max": round(float(np.dot(
                                 np.sort(w)[::-1][:3] /
                                 np.sort(w)[::-1][:3].sum(), MU)), 5),
                             "violates_bound": l1 > b + 1e-9})
    ctx = [r for r in rows if r["contextual"]]
    OUT["3_two_layer"] = {
        "n": len(rows), "n_contextual": len(ctx),
        "n_bound_violations": sum(r["violates_bound"] for r in rows),
        "min_q_among_contextual": min((r["q"] for r in ctx), default=None),
        "rows_sample": rows[:12],
    }


if __name__ == "__main__":
    print("1  界检验（随机搜索）..."); check_bound()
    print("2  等号紧性 ...")
    OUT["2_tightness"] = [x for x in (check_tightness(3, 0.6),
                                      check_tightness(4, 0.9518),
                                      check_tightness(5, 0.99)) if x]
    print("3  两层族 ..."); check_two_layer()
    OUT["2_thresholds"] = threshold_table()

    b = OUT["1_bound_check"]
    print("\n--- 1 界检验：最大超出 = %.3e ---" % b["global_max_excess"])
    for k, d in b["per_k"].items():
        print("  k=%-3s 超出=%+.2e  例: q=%.4f λ1=%.5f 界=%.5f"
              % (k, d["max_excess"], d["q"], d["lam1"], d["bound"]))
    print("--- 2 等号紧性 ---")
    for r in OUT["2_tightness"]:
        print("  k=%d q=%.4f: 解析界=%.6f 优化值=%.6f (q=%.4f)"
              % (r["k"], r["q_target"], r["analytic_bound"],
                 r["optimized"], r["q_achieved"]))
    print("--- 2 门限表 ---")
    for r in OUT["2_thresholds"]:
        print("  k=%d  语境性要求 q >= %.4f  可行=%s  账本峰 D=%s"
              % (r["k"], r["q_k_star"], r["feasible"], r["ledger_peak_D"]))
    d = OUT["3_two_layer"]
    print("--- 3 两层族 ---")
    print("  %d 组：语境 %d 组，违反界 %d 组，语境组最小 q=%s"
          % (d["n"], d["n_contextual"], d["n_bound_violations"],
             d["min_q_among_contextual"]))
    for r in d["rows_sample"][:6]:
        print("   p=%.2f a=%.1f K=%-4d k=%-4d q=%.4f λ1=%.4f 界=%.4f S=%.4f 语境=%s"
              % (r["p"], r["alpha"], r["K"], r["k"], r["q"], r["lam1"],
                 r["bound"], r["S_max"], r["contextual"]))
    with open(os.path.join(HERE, "R51_audit_v3_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("\n-> R51_audit_v3_results.json")
