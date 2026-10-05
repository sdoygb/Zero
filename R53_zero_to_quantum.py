#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R53 · 按"零→分层→成块"的设计，撤掉 L3 后逼出量子力学

设计（用户的图像，逐条落地）：
  (1) 零乱动         : 词 w = ±1 序列，从零出发（Sum w_i = 0 时闭合）
  (2) 分一层一层     : 高度 h = 部分和 S_k 的层级；一层 = 等高度面
  (3) 竖着看高度一致 : 同一 h 上的嵌套片段属于同一层
  (4) 横着看是一层东西: 层的成员 = 所有词的嵌套片段
  (5) 既然是集合     : 层内片段之间由局域补偿移动连通 -> 个体动力学
  (6) 聚合成块       : 密度结构决定"哪几层独立成块"，其余合并
  (7) 撤掉 L3        : π 不是输入，它 = 上述密度结构的函数

然后从 π 逼出量子力学：
  ω_c (推前测度) -> K = -log ω (模 Hamiltonian) -> 模流 -> KMS
  -> 非对易代数 M2 张量 C^k 上的态 -> GNS -> 语境性

判据（三门）：
  门1  q = Sum ω² ∈ (0.5, 0.6)       (R32 生存窗 -> 账本峰 D=4)
  门2  λ1 = ω_max / 前三和 > 0.723607 (KCBS 语境性)
  门3  秩 = k-1                       (R49 稠密判据 -> III_1)
"""
from __future__ import annotations

import json
import math
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
LAM = 0.723607
MU = (math.sqrt(5), (5 - math.sqrt(5)) / 2, (5 - math.sqrt(5)) / 2)
OUT = {}


# ------------------------------------------------------------ (1) 零乱动
def closed_words(L):
    """从零出发、回到零的 ±1 词（Σw=0）。"""
    for pos in combinations(range(L), L // 2):
        w = [-1] * L
        for p in pos:
            w[p] = 1
        yield tuple(w)


# ------------------------------------------------------------ (2)(3) 分层
def partial_sums(w):
    s, out = 0, []
    for x in w:
        s += x
        out.append(s)
    return out


def excursion_profile(w):
    """词在高度轴上的占据：层 h 上停留了几步（h = |S_k| 或 S_k 实现的高度）。"""
    ps = partial_sums(w)
    return Counter(abs(s) for s in ps)


def layer_density(words):
    """每一层 h 的"云密度" = **顶到该高度的云数**（= 最高浮出恰为 h 的词数）。
    这才是"停在某高度的一层云"，不是"穿过该高度的所有路径"。"""
    dens = defaultdict(int)
    for w in words:
        hmax = max(abs(s) for s in partial_sums(w))
        dens[hmax] += 1
    return dict(sorted(dens.items()))


# ------------------------------------------------------------ (6) 成块
def make_blocks(words, dens):
    """按原生密度规则成块：
    找密度极大点 h*；h > h* 的各层独立成块（内层，稀疏→向外递减支），
    h <= h* 的各层合并为一块（外层，递增支）。"""
    hs = sorted(dens)
    hstar = max(hs, key=lambda h: dens[h])
    cnt = Counter()
    for w in words:
        hmax = max(abs(s) for s in partial_sums(w))
        key = hmax if hmax > hstar else 0        # 0 = 合并的外层
        cnt[key] += 1
    return cnt, hstar


# ------------------------------------------------------------ 判据
def exp_vec(k, primes):
    v, kk = [], k
    for p in primes:
        e = 0
        while kk % p == 0:
            kk //= p
            e += 1
        v.append(e)
    return None if kk != 1 else v


def factorize(n):
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


def exact_rank(vals):
    fs = [factorize(v) for v in vals]
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


def ledger_peak(q):
    """F_D = C(D,2) q^D，相邻比 F_D/F_{D-1} = D/(D-2) * q。
    峰在 D  <=>  D/(D-2) * q > 1  且  (D+1)/(D-1) * q < 1。
    校验：q=5/9 时 (5/3)(5/9)=25/27<1 且 2q=10/9>1 => 峰 D=4（与 R32 §2 一致）。"""
    for D in range(2, 400):
        if D == 2:
            if 3 * q < 1:
                return 2
            continue
        if D * q > (D - 2) and (D + 1) * q < (D - 1):
            return D
    return None


def kcbs_smax(weights):
    """S_max = <降序权重, mu>（von Neumann 迹不等式，R37）。
    注意：必须先把顶三**归一**（KCBS 是 3 维归约上的判据）。"""
    v = sorted(weights, reverse=True)[:3]
    s = sum(v)
    if s <= 0:
        return 0.0
    v = [x / s for x in v]
    while len(v) < 3:
        v.append(0.0)
    return sum(a * b for a, b in zip(v, MU))


def analyse(cnt, label):
    tot = sum(cnt.values())
    ws = sorted((v / tot for v in cnt.values()), reverse=True)
    k = len(ws)
    q = sum(x * x for x in ws)
    t3 = ws[:3]
    lam1 = t3[0] / sum(t3) if len(t3) >= 3 else 1.0
    r = exact_rank(sorted(cnt.values(), reverse=True))
    smax = kcbs_smax(ws)
    return {"label": label, "blocks": k, "sizes": sorted(cnt.values(), reverse=True),
            "weights": ws, "q": round(q, 8),
            "lam1": round(lam1, 8), "rank": r, "needed_rank": k - 1,
            "S_max": round(smax, 6),
            "gate1_survival_window": bool(0.5 < q < 0.6),
            "gate2_contextuality": bool(smax > 2.0),
            "gate3_dense": bool(r == k - 1),
            "ledger_peak_D": ledger_peak(q)}


if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    words = list(closed_words(L))
    dens = layer_density(words)
    OUT["L"] = L
    OUT["n_words"] = len(words)
    OUT["layer_density"] = dens

    print("=" * 74)
    print("(1)(2)(3) 零乱动 -> 分层：L=%d，闭合词 %d 个" % (L, len(words)))
    print("=" * 74)
    print("  层高 h : 云密度（在该高度停留的词数占比）")
    tot_w = len(words)
    for h in sorted(dens):
        bar = "#" * int(50 * dens[h] / max(dens.values()))
        print("    %-4d : %6d  (%.4f)  %s" % (h, dens[h], dens[h] / tot_w, bar))

    cnt, hstar = make_blocks(words, dens)
    OUT["h_star"] = hstar
    print()
    print("(6) 成块：密度极大点在 h* = %d" % hstar)
    print("    规则：h > h* 的各层独立成块；h <= h* 的各层合并为一块")
    print("    块结构：", dict(sorted(cnt.items(), key=lambda kv: -kv[1])))

    a = analyse(cnt, "density-rule")
    OUT["result"] = a
    print()
    print("=" * 74)
    print("量子力学逼出：从 pi 到 omega 到 K 到模流")
    print("=" * 74)
    N = sum(cnt.values())
    ws = sorted((v / N for v in cnt.values()), reverse=True)
    print("  omega (推前测度) =", [round(x, 6) for x in ws])
    K = [-math.log(x) for x in ws]
    print("  K = -log omega    =", [round(x, 5) for x in K])
    print("  K 跨度            = %.5f  (非平凡: %s)" % (max(K) - min(K), max(K) - min(K) > 1e-6))
    print()
    print("  模谱（相邻对数比）=", [round(math.log(ws[i] / ws[i + 1]), 4)
                                for i in range(len(ws) - 1)])
    print()
    print("  门1 q        = %.6f   窗 (0.5, 0.6)          %s"
          % (a["q"], "通过" if a["gate1_survival_window"] else "未过"))
    print("  门2 S_max    = %.6f   > 2 (语境性)           %s"
          % (a["S_max"], "通过" if a["gate2_contextuality"] else "未过"))
    print("       lambda1  = %.6f   > 0.723607             %s"
          % (a["lam1"], "通过" if a["lam1"] > LAM else "未过"))
    print("  门3 秩       = %d / 需 %d                  %s"
          % (a["rank"], a["needed_rank"], "通过" if a["gate3_dense"] else "未过"))
    print("  账本峰 D     = %s" % a["ledger_peak_D"])
    allpass = (a["gate1_survival_window"] and a["gate2_contextuality"]
               and a["gate3_dense"])
    print()
    print("  >>> 三门%s" % ("全过：量子力学从 Z* 的嵌套结构中长出，L3 已被撤掉"
                          if allpass else "未全过"))
    OUT["all_pass"] = allpass

    with open(os.path.join(HERE, "R53_zero_to_quantum_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R53_zero_to_quantum_results.json")
