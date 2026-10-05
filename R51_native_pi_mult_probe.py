#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R51 探针 III：在原生 wipe_reseed 上测**全局闭合零层 Z\* 的多重度**

D211 第 3 步：Z* = N^{([C])}，"若同一模式闭合两次，则 n_α 增加"。
所以 π 的最强候选不是"distinct 词均权"也不是"长度类均权"，而是
        ω_c ∝ n_c = 闭合类 c 的**闭合记录多重度**（对所有站点求和）。
这既是 L1 层内的原生记账（历史层是记录，不是均匀先验），也是 R50 §1 说的 L1 对象。

本探针同时测四种读出候选（必须并列报，R50 §3 #5）：
  R1  每 (站点, 精确词) 一条记录      —— 记录多重度的最小口径
  R2  闭合**事件**多重度（含重闭合）  —— D211 "闭合两次则 n 增加"
  R3  长度类（把同一长度的类合并）
  R4  旋转类（[w] 商）
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
          53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]
OUT = {}
LIMIT = int(os.environ.get("LIMIT", "60000"))


def wfs(t): return tuple(1 if c == "+" else -1 for c in t)


def canon(w): return min(w[i:] + w[:i] for i in range(len(w)))


def evolve(active, hist, zero_mult, zero_set):
    """一步：每词派生 w±1；闭合者写入 L1 与 Z*（多重度累加）。"""
    nxt = [set() for _ in active]
    for site, paths in enumerate(active):
        for w in paths:
            for step in (1, -1):
                c = w + (step,)
                if sum(c) == 0:
                    hist[site].add(c)
                    k = canon(c)
                    zero_mult[k] += 1                 # R2：闭合事件多重度
                    zero_set.add(k)
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


def stats(counts):
    tot = sum(counts.values())
    ws = sorted(counts.values(), reverse=True)
    q = sum((x / tot) ** 2 for x in ws)
    l1 = ws[0] / sum(ws[:3])
    return tot, len(ws), q, l1


def gates(q, l1, k, counts):
    need = k - 1
    rank = None
    if k <= 300:
        nums = sorted(counts.values(), reverse=True)
        vecs, ok = [], True
        for i in range(len(nums) - 1):
            a, b = exp_vec(nums[i]), exp_vec(nums[i + 1])
            if a is None or b is None:
                ok = False
                break
            vecs.append([y - x for x, y in zip(a, b)])
        if ok:
            rank = exact_rank(vecs)
    return {"q_in_window": bool(0.5 < q < 0.6), "contextual": bool(l1 > LAM_STAR),
            "rank": rank, "needed_rank": need, "dense": rank == need if rank is not None else None}


def run(gens, period=3, seeds=("+", "-", "++", "--", "+++", "---")):
    active = [{wfs(t)} for t in seeds]
    hist = [set() for _ in seeds]
    zmult, zset = Counter(), set()
    trace = []
    for g in range(1, gens + 1):
        active = evolve(active, hist, zmult, zset)
        if g % period == 0:
            active = reseed(hist)
        # R1：每 (站点,词) 一条
        r1 = Counter()
        for book in hist:
            for w in book:
                r1[w] += 1
        cycR4 = Counter()
        for w, c in r1.items():
            cycR4[canon(w)] += c
        lenR3 = Counter()
        for w, c in r1.items():
            lenR3[len(w)] += c
        act = sum(len(s) for s in active)
        row = {"gen": g, "active": act, "records": sum(r1.values()),
               "z_mult_total": sum(zmult.values()), "z_classes": len(zset)}
        for nm, cnt in [("R1_record", r1), ("R2_eventmult", zmult),
                        ("R3_length", lenR3), ("R4_rotation", cycR4)]:
            t, k, q, l1 = stats(cnt)
            row[nm] = {"blocks": k, "q": round(q, 6), "lam1": round(l1, 5)}
        trace.append(row)
        if act > LIMIT:
            break
    return hist, zmult, zset, trace


if __name__ == "__main__":
    seeds = ("+", "-", "++", "--", "+++", "---")
    print("原生 wipe_reseed：多重度读出（上限 %d 活动词）..." % LIMIT)
    hist, zmult, zset, trace = run(40, 3, seeds)
    print("  gen  active  records  Z*多重 Z*类   q(R2事件) λ1(R2)  q(R1词) λ1(R1)")
    for t in trace:
        if t["gen"] % 1 == 0:
            print("  %3d %7d %8d %8d %6d   %.6f  %.5f   %.6f  %.5f"
                  % (t["gen"], t["active"], t["records"], t["z_mult_total"],
                     t["z_classes"], t["R2_eventmult"]["q"], t["R2_eventmult"]["lam1"],
                     t["R1_record"]["q"], t["R1_record"]["lam1"]))
    OUT["trace"] = trace
    print("\n--- 末代四种读出的两把尺子 ---")
    last = trace[-1]
    for nm, label in [("R1_record", "R1 每(站点,词)一条"),
                      ("R2_eventmult", "R2 闭合事件多重度"),
                      ("R3_length", "R3 长度类"),
                      ("R4_rotation", "R4 旋转类")]:
        d = last[nm]
        g = gates(d["q"], d["lam1"], d["blocks"],
                  zmult if nm == "R2_eventmult" else
                  Counter({w: 1 for w in zset}) if nm == "R4_rotation" else Counter())
        d["gates"] = g
        print("  %-22s 块=%-6d q=%-9.6f λ1=%-8.5f 生存窗=%-5s 语境=%-5s 秩=%s/%d"
              % (label, d["blocks"], d["q"], d["lam1"], g["q_in_window"],
                 g["contextual"], g["rank"], g["needed_rank"]))
    OUT["final_gen"] = last["gen"]
    # 幂律诊断：R2 的事件多重度分布
    ws = sorted(zmult.values(), reverse=True)
    if len(ws) > 5:
        xs = [log(i + 1) for i in range(len(ws))]
        ys = [log(max(v, 1)) for v in ws]
        n = len(xs)
        mx, my = sum(xs) / n, sum(ys) / n
        den = sum((a - mx) ** 2 for a in xs)
        slope = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / den if den else None
        print("\nR2 事件多重度排名-log 斜率（幂律指数）= %.3f" % slope)
        print("  最大/最小多重度 = %d / %d，类数 = %d" % (ws[0], ws[-1], len(ws)))
        OUT["powerlaw_slope_R2"] = round(slope, 4)
    with open(os.path.join(HERE, "R51_native_pi_mult_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("-> R51_native_pi_mult_results.json")
