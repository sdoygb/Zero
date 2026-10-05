#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R51 审计探针：三条判据的独立复算

A  R44 no-go 的普适性：闭式 S_max(q) 与"任意谱"数值上确界的一致性
B  R49 秩判据用**精确有理数**重算（原探针用浮点 SVD）
C  账本形式逃逸出口：把峰搬进 q>3/5 需要多重度多陡
D  两层族的 q 落点：与 R32 生存窗口 (1/2,3/5) 的冲突量
"""
from __future__ import annotations

import json
import os
from fractions import Fraction
from itertools import product
from math import log, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MU = [sqrt(5), (5 - sqrt(5)) / 2, (5 - sqrt(5)) / 2]
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
          53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103]
OUT = {}


def smax(lam):
    """S_max = <lam_desc, mu>，lam 自动降序。"""
    v = sorted(lam, reverse=True)[:3]
    if len(v) < 3:                      # 不足 3 块：补零（3 维归约）
        v = v + [0.0] * (3 - len(v))
    return sum(a * b for a, b in zip(v, MU))


def q_of(lam):
    return sum(x * x for x in lam)


# ---------------------------------------------------------------- A
def audit_A(n=400000, seed=7):
    """给定 q，S_max 的解析上确界是否被任意谱达到/超过。"""
    rng = np.random.default_rng(seed)
    worst = {}                                  # q-bucket -> max S
    for _ in range(n):
        k = int(rng.integers(3, 12))
        w = rng.random(k) ** rng.choice([0.5, 1.0, 2.0, 4.0])
        w = w / w.sum()
        q, s = q_of(w), smax(w)
        b = round(q, 2)
        if s > worst.get(b, -1):
            worst[b] = s
    rows = []
    for b in sorted(worst):
        if b <= 0.15:
            continue
        closed = MU[1] + (MU[0] - MU[1]) * (1 + sqrt(max(2 * b - 1, 0))) / 2
        rows.append({"q": b, "s_num": round(worst[b], 6),
                     "s_closed": round(closed, 6),
                     "excess": round(worst[b] - closed, 9)})
    OUT["A_universality"] = {
        "note": "s_num 不得超过 s_closed（R44-2 上界）；q<3/5 时 s_num 应 <2",
        "rows": rows,
        "max_excess": max(r["excess"] for r in rows),
        "max_S_below_3_5": max((r["s_num"] for r in rows if r["q"] < 0.595),
                               default=None),
    }


# ---------------------------------------------------------------- B
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
    """整数行向量的精确秩（Q 上 Gauss 消元，Fraction）。"""
    M = [[Fraction(x) for x in r] for r in rows]
    rank, ncol = 0, len(M[0]) if M else 0
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


def rank_of(w, exact=True):
    W = [Fraction(str(x)) for x in w]
    den = 1
    for x in W:
        den = den * x.denominator // np.gcd(den, x.denominator)
    nums = [int(x * den) for x in W]
    vecs = []
    for i in range(len(nums) - 1):
        a, b = exp_vec_int(nums[i]), exp_vec_int(nums[i + 1])
        if a is None or b is None:
            return None
        vecs.append([y - x for x, y in zip(a, b)])
    if exact:
        return exact_rank(vecs)
    return int(np.linalg.matrix_rank(np.array(vecs, float), tol=1e-9))


def audit_B():
    cases = [([2, 5, 20, 100], 2), ([144, 36, 16, 9], 2),
             ([1, 2, 4, 8], 1), ([16, 2, 3, 5], 3)]
    rows = []
    for w, want in cases:
        re_, rf = rank_of(w, True), rank_of(w, False)
        rows.append({"w": w, "claimed": want, "exact_rank": re_,
                     "float_rank": rf, "agree": re_ == rf,
                     "ok": re_ == want})
    # 随机整数权重下的精确 vs 浮点对照
    rng = np.random.default_rng(11)
    mism, tot = 0, 4000
    for _ in range(tot):
        k = 4
        w = [int(rng.integers(1, 400)) for _ in range(k)]
        a, b = rank_of(w, True), rank_of(w, False)
        if a != b:
            mism += 1
    OUT["B_exact_rank"] = {"cases": rows, "random_mismatch": mism,
                           "random_total": tot}


# ---------------------------------------------------------------- C
def m_pair(D):  return D * (D - 1) // 2
def m_dict(D):  return (D + 1) * D // 2


def audit_C():
    """峰在 D=4 当且仅当 q ∈ (m4/(m4+m3), m4/(m4+m5))。
    要覆盖 q=3/5 需上端 > 3/5 ⇔ m5/m4 < 5/3。"""
    rows = []
    for name, m in [("C(D,2)", m_pair), ("C(D+1,2)", m_dict),
                    ("D^3", lambda D: D ** 3), ("2^D", lambda D: 2 ** D),
                    ("D!", lambda D: float(np.math.factorial(D)) if hasattr(np, "math")
                     else __import__("math").factorial(D))]:
        m3, m4, m5 = m(3), m(4), m(5)
        lo, hi = m4 / (m4 + m3), m4 / (m4 + m5)
        rows.append({"mult": name, "ratio_m5_m4": round(m5 / m4, 4),
                     "window": [round(lo, 6), round(hi, 6)],
                     "covers_0.6": hi > 0.6,
                     "peak_at_4_for_q=0.6": lo < 0.6 < hi})
    need = [D for D in range(3, 9)
            if m_pair(4) / (m_pair(4) + m_pair(D)) > 0.6]
    OUT["C_ledger_escape"] = {
        "condition": "peak at D=4 covers q=3/5  <=>  m(5)/m(4) < 5/3",
        "rows": rows,
        "min_steepness_for_0.6": "m5/m4 < 1.6667",
        "standard_pair_carrier_m5_m4": 10 / 6,
    }


# ---------------------------------------------------------------- D
def two_layer(p, alpha, K):
    tail = [a ** (-alpha) for a in range(1, K + 1)]
    s = sum(tail)
    w = [p] + [(1 - p) * t / s for t in tail]
    return w


def audit_D():
    rows = []
    for p in [0.75, 0.80, 0.85, 0.90]:
        for alpha in [1.0, 2.0, 3.0]:
            for K in [8, 32, 256]:
                w = two_layer(p, alpha, K)
                q, s = q_of(w), smax(w)
                rows.append({"p": p, "alpha": alpha, "K": K,
                             "q": round(q, 5), "S_max": round(s, 5),
                             "contextual": s > 2,
                             "peak_dim": peak_dim(q),
                             "survives": 0.5 < q < 0.6})
    ok = [r for r in rows if r["contextual"] and r["survives"]]
    OUT["D_two_layer"] = {
        "rows": rows[:14],
        "n_contextual_and_surviving": len(ok),
        "n_total": len(rows),
        "q_range": [min(r["q"] for r in rows), max(r["q"] for r in rows)],
        "note": "R32 生存要求峰在 D=4 ⇔ q∈(0.5,0.6)；语境性要求 q>=0.6",
    }


def peak_dim(q):
    best, bd = -1, None
    for D in range(2, 60):
        F = m_pair(D) * q ** D
        if F > best:
            best, bd = F, D
    return bd


if __name__ == "__main__":
    print("A  R44 普适性 ..."); audit_A()
    print("B  R49 精确秩 ...");  audit_B()
    print("C  账本逃逸 ...");    audit_C()
    print("D  两层族 ...");      audit_D()
    print("\n--- A ---")
    print("  最大超出闭式上界:", OUT["A_universality"]["max_excess"])
    print("  q<0.595 的最大 S_max:", OUT["A_universality"]["max_S_below_3_5"])
    print("--- B ---")
    for r in OUT["B_exact_rank"]["cases"]:
        print("  w=%-18s 声称秩=%d 精确秩=%d 浮点秩=%d %s"
              % (r["w"], r["claimed"], r["exact_rank"], r["float_rank"],
                 "OK" if r["ok"] else "**不符**"))
    print("  随机 4 块整数权重精确/浮点不符:", OUT["B_exact_rank"]["random_mismatch"],
          "/", OUT["B_exact_rank"]["random_total"])
    print("--- C ---")
    for r in OUT["C_ledger_escape"]["rows"]:
        print("  m=%-9s m5/m4=%-7s 窗口=%-22s 覆盖0.6=%s"
              % (r["mult"], r["ratio_m5_m4"], r["window"], r["covers_0.6"]))
    print("--- D ---")
    d = OUT["D_two_layer"]
    print("  两层族 q 范围:", d["q_range"], " 同时语境+生存:",
          d["n_contextual_and_surviving"], "/", d["n_total"])
    for r in d["rows"][:6]:
        print("   p=%.2f a=%.1f K=%-4d q=%.5f S=%.4f 峰=D%d"
              % (r["p"], r["alpha"], r["K"], r["q"], r["S_max"], r["peak_dim"]))
    with open(os.path.join(HERE, "R51_audit_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("\n-> R51_audit_results.json")
