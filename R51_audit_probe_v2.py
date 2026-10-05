#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R51 审计探针 v2：严格区分"全谱 q_d"与"顶三归一 q3 / λ1"

关键更正（v1 的错）：R44-2 的界是对**归一化**谱说的。
  全 k 块谱 w：q_k = Σ w_i²
  顶三归一谱 λ：λ_i = w_i / T，T = 前三块之和 ⇒ q3 = Σ λ_i² = (Σ_top3 w²)/T²
  显然 q3 ≠ q_k。KCBS 语境性用 λ1（顶三归一的第一分量）。
"""
from __future__ import annotations

import json
import os
from math import sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MU = np.array([sqrt(5), (5 - sqrt(5)) / 2, (5 - sqrt(5)) / 2])
OUT = {}


def smax(lam3):
    """3 维归一谱的最大 KCBS 值（降序配对，von Neumann 迹不等式）。"""
    return float(np.dot(np.sort(np.asarray(lam3))[::-1], MU))


def top3(w):
    s = sorted(w, reverse=True)[:3]
    t = sum(s)
    return np.array([x / t for x in s])


def lam1_bound(q):
    return (1 + sqrt(max(2 * q - 1, 0))) / 2


# ------------------------------------------------------------------ A
def audit_A(n=600000, seed=3):
    """(A1) 3 块谱：λ1 上确界是否 = (1+√(2q-1))/2，且 S_max(q) 闭式是否成立。
       (A2) 4/5/6 块谱：顶三归一 λ1 与**全谱** q_k 的关系（R44 未覆盖处）。"""
    rng = np.random.default_rng(seed)
    a1_worst, a2_worst = {}, {}
    for _ in range(n):
        k3 = 3
        w = rng.random(k3) ** rng.choice([0.4, 1.0, 2.5, 5.0])
        w = w / w.sum()
        q = float(np.sum(w ** 2))
        b = round(q, 3)
        a1_worst[b] = max(a1_worst.get(b, -1), float(np.max(w)))

        k = int(rng.integers(4, 9))
        w = rng.random(k) ** rng.choice([0.4, 1.0, 2.5, 5.0])
        w = w / w.sum()
        qk = float(np.sum(w ** 2))
        l1 = float(np.max(top3(w)))
        b = round(qk, 3)
        a2_worst[b] = max(a2_worst.get(b, -1), l1)

    rows1 = []
    for b in sorted(a1_worst):
        if b < 0.34:
            continue
        rows1.append({"q3": b, "lam1_num": round(a1_worst[b], 6),
                      "lam1_closed": round(lam1_bound(b), 6),
                      "excess": round(a1_worst[b] - lam1_bound(b), 9)})
    rows2 = []
    for b in sorted(a2_worst):
        if b < 0.34:
            continue
        rows2.append({"q_k": b, "lam1_top3_num": round(a2_worst[b], 6),
                      "lam1_bound_at_q_k": round(lam1_bound(b), 6),
                      "excess_over_bound": round(a2_worst[b] - lam1_bound(b), 6),
                      "contextual_at_q<0.6": (b < 0.6 and a2_worst[b] > 0.723607)})
    OUT["A1_three_block"] = {
        "claim": "R44-2: lam1_max(q)= (1+sqrt(2q-1))/2 ; S_max(3/5)=2",
        "rows": rows1,
        "max_excess": max(r["excess"] for r in rows1),
        "S_at_q_0.6": round(smax([lam1_bound(0.6),
                                  (1 - lam1_bound(0.6)) / 2,
                                  (1 - lam1_bound(0.6)) / 2]), 8),
    }
    OUT["A2_multiblock"] = {
        "claim": "R44 未覆盖：全谱 q_k<3/5 时顶三归一 λ1 能否 >0.7236？",
        "rows": rows2,
        "max_lam1_below_qk_0.6": max([r["lam1_top3_num"] for r in rows2
                                      if r["q_k"] < 0.6], default=None),
        "any_contextual_below_0.6": any(r["contextual_at_q<0.6"] for r in rows2),
    }


# ------------------------------------------------------------------ B
def audit_B(steps=400, seed=5):
    """直接优化：max λ1(顶三归一) s.t. 全 4 块谱 q4 = Q。"""
    rng = np.random.default_rng(seed)
    res = {}
    for Q in [0.45, 0.50, 0.55, 0.56, 0.58, 0.595, 0.60, 0.62, 0.65, 0.70]:
        best, arg = -1, None
        for _ in range(120000):
            # 参数化 4 块谱：w1>=w2>=w3>=w4
            a = rng.random(3)
            w = np.sort(np.concatenate([a, [1.0]]))[::-1]
            w = w / w.sum()
            q = float(np.sum(w ** 2))
            if abs(q - Q) > 0.004:
                continue
            l1 = float(np.max(top3(w)))
            if l1 > best:
                best, arg = l1, w.tolist()
        res[Q] = {"max_lam1": round(best, 6),
                  "bound_from_A1_at_same_q": round(lam1_bound(Q), 6),
                  "bound_from_R44_at_q3max": None,
                  "argmax_w": [round(x, 5) for x in arg] if arg else None,
                  "contextual": best > 0.723607}
    OUT["B_four_block_opt"] = {
        "note": "4 块全谱 q4 固定时，顶三归一 λ1 的上确界（随机 + 局部细化）",
        "rows": res,
    }
    # 局部细化：对最优附近做坐标上升
    refine = {}
    for Q, d in res.items():
        if d["argmax_w"] is None:
            continue
        w = np.array(d["argmax_w"])
        for _ in range(4000):
            i, j = rng.integers(0, 4, 2)
            if i == j:
                continue
            step = rng.normal(0, 0.004)
            cand = w.copy()
            cand[i] = max(cand[i] + step, 1e-6)
            cand[j] = max(cand[j] - step, 1e-6)
            cand = np.sort(cand)[::-1]
            cand = cand / cand.sum()
            if abs(np.sum(cand ** 2) - Q) > abs(np.sum(w ** 2) - Q) + 1e-9:
                continue
            if np.max(top3(cand)) > np.max(top3(w)):
                w = cand
        refine[Q] = {"max_lam1": round(float(np.max(top3(w))), 6),
                     "q": round(float(np.sum(w ** 2)), 6),
                     "w": [round(float(x), 5) for x in w]}
    OUT["B_four_block_refined"] = refine


# ------------------------------------------------------------------ C
def audit_C():
    """R49-3 的门限是充分条件；R44 的界把它变成必要条件——列出两者差异。"""
    lam1_star = 0.723607
    q_needed = (2 * lam1_star ** 2 - 2 * lam1_star + 1)   # q = λ1²+(1-λ1)² at equality
    rows = []
    for l1 in [0.70, 0.72, 0.723607, 0.75, 0.80]:
        qmin = l1 ** 2 + (1 - l1) ** 2          # worst case: rest on one block
        # 最小 q 使该 λ1 可能（其余集中在同一点）
        rows.append({"lam1": l1,
                     "lambda1_over_top3": round(l1, 6),
                     "min_q3_for_this_lam1": round(qmin, 6),
                     "S_max_if_rest_split_evenly": round(
                         smax([l1, (1 - l1) / 2, (1 - l1) / 2]), 6)})
    OUT["C_threshold"] = {
        "R49_sufficient_lambda1": lam1_star,
        "q3_at_which_S=2": 0.6,
        "lambda1_at_q3_0.6": round(lam1_bound(0.6), 6),
        "R44_necessary_condition": "contextual => q3 >= 3/5 (>=0.6)",
        "rows": rows,
    }


if __name__ == "__main__":
    print("A  3 块界与多块对照 ..."); audit_A()
    print("B  4 块优化 ...");        audit_B()
    print("C  门限 ...");            audit_C()
    a1 = OUT["A1_three_block"]
    print("\n--- A1 3 块：λ1 上确界 vs 闭式 ---")
    for r in a1["rows"]:
        print("  q3=%.3f  λ1数值=%.6f  闭式=%.6f  超出=%+.2e"
              % (r["q3"], r["lam1_num"], r["lam1_closed"], r["excess"]))
    print("  最大超出:", a1["max_excess"], " S_max(q=0.6)=", a1["S_at_q_0.6"])
    a2 = OUT["A2_multiblock"]
    print("--- A2 多块（4..8 块）---")
    for r in a2["rows"]:
        print("  q_k=%.3f  顶三归一λ1=%.6f  同q闭式界=%.6f  超出=%+.4f %s"
              % (r["q_k"], r["lam1_top3_num"], r["lam1_bound_at_q_k"],
                 r["excess_over_bound"],
                 "<-- 语境且 q<0.6" if r["contextual_at_q<0.6"] else ""))
    print("  q<0.6 时最大 λ1:", a2["max_lam1_below_qk_0.6"],
          " 是否存在语境:", a2["any_contextual_below_0.6"])
    print("--- B 4 块最优 ---")
    for Q, d in OUT["B_four_block_opt"].items():
        print("  q4=%.3f  maxλ1=%.6f  语境=%s" % (Q, d["max_lam1"], d["contextual"]))
    print("--- B 细化 ---")
    for Q, d in OUT["B_four_block_refined"].items():
        print("  q4=%.3f  maxλ1=%.6f  w=%s" % (d["q"], d["max_lam1"], d["w"]))
    with open(os.path.join(HERE, "R51_audit_v2_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("\n-> R51_audit_v2_results.json")
