#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R87 · 类级账本权重的一致性（R25 的 omega_C = o_C/N 是不是动力学的推前）

背景
  R86 §4（F4）证明：同一旋转类内不同词给出**不同**的后继类分布
  => 类级播种核 P[C'|C] 不是良定义的。
  而 R25 (R25-1)-(R25-5) 的账本权重取 omega_C = o_C/N_L，
  这等于"词上均匀测度经 lambda 的推前"。二者相容的条件必须写清。

本探针（三层判据）
  H1 [词层平稳]   活动层对象 = 长度 L 的滑动窗（w -> shift(w)+sigma，sigma 均匀 ±1）。
                  转移矩阵双随机（每个词恰有两个前像）=> **词上均匀是平稳测度**。
  H2 [推前]       词均匀经 lambda 推前 => omega_C = o_C/N_L；据此算 q_L 并与 R25 的
                  q_4=5/9、q_6=7/25、q_8=133/1225 逐位比对。
  H3 [偏差来源]   把"离上次记录的步数" s 计入（两层链），看 s 是否改变类权重：
                  - 若 H1 是唯一平稳 => s 只是时钟，不改权重（R25 无需额外假设）；
                  - 并给出 R86 F4 那个偏差的**准确定名**：F4 的类内分布是
                    "**首次闭合长度**加权的词分布"，不是平稳词分布。

退出码 0 = 全部核验通过。
"""
from __future__ import annotations

import itertools
import json
import os
from collections import defaultdict
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT, NOTES = {}, []


def note(name, ok, detail=""):
    NOTES.append((name, ok, detail))
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))
    return ok


def balanced_words(L):
    if L % 2:
        return []
    return [tuple(1 if i in pos else -1 for i in range(L))
            for pos in itertools.combinations(range(L), L // 2)]


def canonical_class(w):
    L = len(w)
    best = None
    for k in range(L):
        r = tuple(w[(i + k) % L] for i in range(L))
        if sum(r) == 0 and (best is None or r < best):
            best = r
    return best


def H1(L):
    """词层转移：每个词的前像个数（判断双随机/均匀平稳）。"""
    words = list(itertools.product([1, -1], repeat=L))
    pre = defaultdict(int)
    for w in words:
        for s in (1, -1):
            pre[tuple(w[1:]) + (s,)] += 1
    counts = set(pre.values())
    return counts, len(words)


def H2(L):
    fib = defaultdict(list)
    for w in balanced_words(L):
        fib[canonical_class(w)].append(w)
    N = len(balanced_words(L))
    sizes = np.array([len(fib[c]) for c in sorted(fib)], float)
    piC = sizes / N
    q = float((piC ** 2).sum())
    return sizes, piC, q, len(fib)


def H3(L):
    """两层链：(词, s)；s = 离上次记录的步数（0..L-1）。精确平稳分布。"""
    words = list(itertools.product([1, -1], repeat=L))
    states = [(w, s) for w in words for s in range(L)]
    sidx = {st: i for i, st in enumerate(states)}
    n = len(states)
    P = np.zeros((n, n))
    for (w, s), i in sidx.items():
        for sigma in (1, -1):
            w2 = tuple(w[1:]) + (sigma,)
            s2 = 0 if sum(w2) == 0 else (s + 1) % L
            P[i, sidx[(w2, s2)]] += 0.5
    # 平稳：解 (P^T - I) pi = 0, sum pi = 1（最小二乘 + 归一）
    A = np.vstack([P.T - np.eye(n), np.ones(n)])
    b = np.concatenate([np.zeros(n), [1.0]])
    pi, *_ = np.linalg.lstsq(A, b, rcond=None)
    pi = np.clip(pi, 0, None)
    pi = pi / pi.sum()
    resid = float(np.abs(P.T @ pi - pi).max())
    # 词边缘
    piw = defaultdict(float)
    for (w, s), i in sidx.items():
        piw[w] += pi[i]
    # 类权重
    fib = defaultdict(list)
    for w in balanced_words(L):
        fib[canonical_class(w)].append(w)
    keys = sorted(fib)
    piC = np.array([sum(piw[w] for w in fib[c]) for c in keys])
    return resid, piC, piw, keys, fib


def main():
    print("=" * 78)
    print("R87 · 类级账本权重的一致性")
    print("=" * 78)

    expect = {4: Fraction(5, 9), 6: Fraction(7, 25), 8: Fraction(133, 1225)}

    print("\n[H1] 词层转移双随机（每个词恰有两个前像）")
    for L in [4, 6, 8]:
        counts, nw = H1(L)
        note("L=%d 每个词的前像个数 = 2（双随机）" % L, counts == {2},
             "前像数集合 = %s，词数 = %d" % (sorted(counts), nw))

    print("\n[H2] 词均匀的推前 = R25 权重；q 与 R25 逐位比对")
    h2 = {}
    for L in [4, 6, 8]:
        sizes, piC, q, nc = H2(L)
        exp = float(expect[L])
        note("L=%d q = %s （R25: %s）" % (L, ("%.9f" % q), str(expect[L])),
             abs(q - exp) < 1e-12, "偏差 %.2e" % abs(q - exp))
        note("L=%d 推前与 o_C/N 逐位相同" % L,
             np.allclose(piC, sizes / sizes.sum()),
             "piC=%s" % np.round(piC, 6).tolist())
        h2["L%d" % L] = {"sizes": sizes.astype(int).tolist(), "piC": piC.tolist(),
                         "q": q, "q_R25": exp, "classes": nc}
    OUT["H2"] = h2

    print("\n[H3] 两层链（把 s 计入）：记录条件分布是否等于 R25 权重")
    h3 = {}
    for L in [4, 6]:
        resid, piC, piw, keys, fib = H3(L)
        note("L=%d 两层链平稳分布残差 < 1e-10" % L, resid < 1e-10,
             "max|P^T pi - pi| = %.2e" % resid)
        # 记录条件是**条件**分布：只在平衡词上产生记录 => 按平衡词总质量归一
        mass_bal = piC.sum()
        piC_cond = piC / mass_bal
        sizes, piC_u, q_u, _ = H2(L)
        dev = float(np.abs(piC_cond - piC_u).max())
        note("L=%d 记录条件分布 = o_C/N（s 只是时钟，不改权重）" % L, dev < 1e-8,
             "max|piC_cond - o_C/N| = %.2e ；平衡词总质量 = %.6f" % (dev, mass_bal))
        q_cond = float((piC_cond ** 2).sum())
        note("L=%d 条件 q = %s" % (L, "%.9f" % q_cond), abs(q_cond - float(expect[L])) < 1e-12,
             "R25 = %s" % expect[L])
        wdev = float(max(piw.values()) * len(list(itertools.product([1, -1], repeat=L)))
                     - min(piw.values()) * len(list(itertools.product([1, -1], repeat=L))))
        note("L=%d 词边缘均匀（每个词 1/2^L）" % L, wdev < 1e-8,
             "偏差 = %.2e" % wdev)
        h3["L%d" % L] = {"resid": resid, "max_dev_cond_vs_R25": dev,
                         "mass_on_balanced": mass_bal, "q_conditional": q_cond,
                         "word_uniformity_dev": wdev}
    OUT["H3"] = h3

    print("\n[定名] R86 F4 的偏差来源")
    print("  R86 F4 用的是「**首次闭合长度**加权的词分布」——即把每个词的")
    print("  闭合概率当作权重；而平稳词分布是**均匀**的。两者对类内权重给出不同答案，")
    print("  这正是 F4 测到的'同类内不同词后继分布不同'。")

    print("\n" + "=" * 78)
    bad = [n for n, ok, _ in NOTES if not ok]
    print("汇总：核验 %d 项，未过 %d 项 %s" % (len(NOTES), len(bad), bad if bad else ""))
    print("=" * 78)
    OUT["summary"] = {"checks": len(NOTES), "failed": bad}
    with open(os.path.join(HERE, "R87_class_weight_consistency_results.json"), "w") as f:
        json.dump(OUT, f, ensure_ascii=False, indent=1)
    print("-> R87_class_weight_consistency_results.json")
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
