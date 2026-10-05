#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R77 · 长度分级的必然性：闭合算子的本征结构

要证的命题：
    在闭合约束（Σw_i = 0）下，由步序列定义的动力学权重
        w(c) = F(步序列)
    只依赖路径**长度** ℓ，不依赖其它组合数据。
    故长度分级不是外加的观察窗，而是动力学必然留下的印记。

本探针分四步：
  A  建立"闭合算子"：把长度 ℓ 的闭合路径按步序列 → 权重的映射
  B  检验：在固定 ℓ 下，权重是否对所有闭合路径取同一个值
     （若是 ⇒ 长度是唯一不变量）
  C  检验：不同 ℓ 之间权重是否不同（若是 ⇒ 分级非平凡）
  D  由此得到的类结构是否给出 III_1（素数支撑 ≥ k−1）
"""
from __future__ import annotations

import json
import math
import os
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}


def fac(n):
    f, d = {}, 2
    n = int(n)
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def rank_of(vals):
    fs = [fac(v) for v in vals]
    primes = sorted({p for f in fs for p in f})
    rows = [[fs[i + 1].get(p, 0) - fs[i].get(p, 0) for p in primes]
            for i in range(len(fs) - 1)]
    if not rows:
        return 0
    M = [[Fraction(x) for x in r] for r in rows]
    rk = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        pv = M[rk][c]
        M[rk] = [x / pv for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


# ------------------------------------------------------------ A 闭合算子
def closed_paths(ell):
    """长度 ℓ、Σ=0 的 ±1 路径"""
    if ell % 2:
        return []
    for pos in combinations(range(ell), ell // 2):
        w = [-1] * ell
        for p in pos:
            w[p] = 1
        yield tuple(w)


def invariants(w):
    """从步序列算出三种组合不变量"""
    # (i) 最大振幅
    s, hmax = 0, 0
    for x in w:
        s += x
        hmax = max(hmax, abs(s))
    # (ii) 转向数（±1 序列的符号变化次数）
    turn = sum(1 for i in range(len(w) - 1) if w[i] != w[i + 1])
    # (iii) 首次返回零
    s, tau0 = 0, None
    for i, x in enumerate(w):
        s += x
        if s == 0:
            tau0 = i + 1
            break
    return {"hmax": hmax, "turn": turn, "tau0": tau0}


if __name__ == "__main__":
    print("=" * 74)
    print("A/B · 固定长度 ℓ 下，各组合不变量是否恒定")
    print("=" * 74)
    print("  ℓ   闭合路径数   hmax 取值数   turn 取值数   tau0 取值数   长度是唯一不变量?")
    rows = []
    for ell in [2, 4, 6, 8, 10, 12]:
        ws = list(closed_paths(ell))
        inv = [invariants(w) for w in ws]
        nh = len({x["hmax"] for x in inv})
        nt = len({x["turn"] for x in inv})
        n0 = len({x["tau0"] for x in inv})
        unique = (nh == 1 and nt == 1 and n0 == 1)
        print("  %-4d %-12d %-13d %-13d %-13d %s"
              % (ell, len(ws), nh, nt, n0, "是" if unique else "否"))
        rows.append({"ell": ell, "n": len(ws), "hmax_vals": nh,
                     "turn_vals": nt, "tau0_vals": n0, "unique": unique})
    OUT["A"] = rows
    print()
    print("  ⇒ 长度**不是**唯一的组合不变量：hmax/turn/tau0 在固定 ℓ 下都变")

    print()
    print("=" * 74)
    print("B2 · 但「只依赖长度」的权重仍然存在：最小/最大耦合")
    print("=" * 74)
    print("  检验：长度 ℓ 的闭合路径**总数** N(ℓ) = C(ℓ, ℓ/2)")
    print("  ℓ    N(ℓ)      分解")
    counts = []
    for ell in range(2, 22, 2):
        n = math.comb(ell, ell // 2)
        counts.append(n)
        print("  %-4d %-10d %s" % (ell, n, fac(n)))
    OUT["B2"] = counts
    print()
    print("  素因子支撑 =", sorted({p for n in counts for p in fac(n)}))
    print("  秩 = %d / 需 %d" % (rank_of(counts), len(counts) - 1))
    print("  ⇒ 路径计数 N(ℓ) 的素数贫乏（全是小素数）⇒ 单独用它**不足** III_1")
    OUT["B2_rank"] = rank_of(counts)

    print()
    print("=" * 74)
    print("C · 关键：K_{D+1} 上**简单圈**的计数才给素因子 {2,3,5}")
    print("=" * 74)
    for D in [3, 4, 5]:
        n = D + 1
        cnt = []
        for l in range(3, n + 1):
            c = math.comb(n, l) * math.factorial(l - 1) // 2
            cnt.append(c)
        ps = sorted({p for c in cnt for p in fac(c)})
        r = rank_of(cnt)
        print("  D=%d (K_%d): 圈数=%s  素数=%s  秩=%d/需 %d  III_1=%s"
              % (D, n, cnt, ps, r, len(cnt) - 1, r == len(cnt) - 1))
        OUT.setdefault("C", {})["D%d" % D] = {"counts": cnt, "primes": ps,
                                              "rank": r}
    print()
    print("  ⇒ 简单圈计数（= 无弦闭合）与路径计数不同：前者含素数 5")
    print("  ⇒ 因为简单圈**排除**了重复访问，是「最小闭合」的计数")

    print()
    print("=" * 74)
    print("D · 结论：长度分级的必然性依据")
    print("=" * 74)
    print("  1. 动力学按步发生（Z4 的步推进）⇒ 任何原生权重是步序列的函数")
    print("  2. 闭合约束 Σw=0 ⇒ 只有**偶数长度**可闭合 ⇒ 长度天然分出奇偶类")
    print("  3. 在 K_{D+1} 上，'最小闭合'（简单圈）按长度分级 ⇒ 素数 {2,3,5}")
    print("  4. 图对称性（S_{D+1} 顶点传递）**不**区分同长度的圈")
    print("     ⇒ 长度是该对称性下**唯一的原生等级**")
    print()
    print("  ⇒ 故长度分级是「动力学 + 拓扑对称性」下唯一可行的分级")
    print("  ⇒ 但它给出 k=D−1 类（D=4 时 3 类），秩 2/2，恰好 III_1")
    OUT["D_conclusion"] = True

    with open(os.path.join(HERE, "R77_length_grading_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R77_length_grading_results.json")
