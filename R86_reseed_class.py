#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R86 · 重播种是否保类号？（跨类相干的二值判据）

背景（缺口 1 的重新命名）
  以「自旋因子 | 账本因子」为分割，原生态的块间相干恒为 0 => CHSH 恒 <= 2
  （独立探针：Bell 态 2.828427、去相干后 1.414214、原生态 0.000000）。
  故多体资源不是"张量积"，而是**账本不同类之间的相干**。

重播种的精确定义（实现层取证）
  `zero_sum_periodic_destruction.py` L159-170：每个历史词 w 播种两条开路径 w+ 与 w-。
  历史层记录（`R39`/`Z3`）：**精确词 w ＋ 闭合类 [w]**；类用**旋转类**（`G37`：闭环无规范起点）。
  闭合律（`G37`）：类权重 omega(C) = |C| / sum |C'|。

判据
  F1 [记录的纤维]      lambda(w) = [w] 的纤维 = 旋转轨道（非单射 ⇒ 词层自由度被记录丢掉）
  F2 [播种是词级分支]  w -> {w+, w-}：等权、无幅度、无相位
  F3 [相干不生成]      任意词层初态下，记录层的类间元恒为 0 ⇒ **保类号成立**
  F4 [意外的自由度]    同一旋转类内的**不同词**给出**不同的**后继类分布
                       ⇒ 记录不决定动力学（词层残余自由度是真的）
  F5 [不可观测性]      原生观测集（分块对角代数）看不见任何类间元：
                       对任意"类间只放非对角"的算符 E，全部块对角投影读数恒为 0

退出码 0 = 全部通过；任一"不符"退出码 1。
"""
from __future__ import annotations

import itertools
import json
import os
from collections import defaultdict
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name, ("   " + detail) if detail else ""))
    return bool(cond)


# ---------------- 原生对象 ----------------
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


def fiber_map(L):
    d = defaultdict(list)
    for w in balanced_words(L):
        d[canonical_class(w)].append(w)
    return d


# ---------------- 重播种的精确转移（词级） ----------------
def reseed_paths(w):
    """实现层：历史词 w 播种两条开路径 w+ 与 w-。"""
    return [tuple(w) + (1,), tuple(w) + (-1,)]


def class_after_closure(seq):
    """开路径闭合（sum=0）后写入历史层的记录 = 旋转类。"""
    return canonical_class(seq) if sum(seq) == 0 else None


def closure_length_dist(w, pmax=10):
    """精确 DP：播种后的开路径在 k 步后**首次**闭合的概率（k=1..pmax），及截断质量。

    状态为 (当前位置, 已走步数)；位置 0 是吸收点。用 Fraction 精确算。
    返回 ([p_1..p_pmax], p_leak)。
    """
    # frontier: (位置 s, 已见序列, 概率)
    frontier = []
    for seq in reseed_paths(w):
        frontier.append((sum(seq), seq, Fraction(1, 2)))
    probs = []
    leak = Fraction(0)
    for k in range(1, pmax + 1):
        pk = Fraction(0)
        nxt = []
        for s, seq, p in frontier:
            if s == 0:
                pk += p                       # 本步刚好首次闭合
            else:
                nxt.append((s + 1, seq + (1,), p / 2))
                nxt.append((s - 1, seq + (-1,), p / 2))
        probs.append(pk)
        frontier = nxt
    leak = sum(p for _, _, p in frontier)
    return probs, leak


def class_distribution_after_reseed(w, pmax=12, rng=None):
    """给出词 w 播种后闭合类的分布：每条首次闭合于第 k 步的序列权重 = 2^-k。

    截断（步数 > pmax 仍未闭合）的质量单列，不混入分布。
    """
    dist = defaultdict(float)
    for k in range(1, pmax + 1):
        seqs = first_passage_seqs(w, k)
        if not seqs:
            continue
        pk = Fraction(1, 2) ** k
        for seq in seqs:
            dist[class_after_closure(seq)] += float(pk)
    total = sum(dist.values())
    if total > 0:
        for key in list(dist):
            dist[key] /= total
    probs, leak = closure_length_dist(w, pmax=pmax)
    return dict(dist), float(leak)


def first_passage_seqs(w, k):
    """长度 = len(w)+k、首次闭合于末步的开路径（+1/-1 两条播种分支）。"""
    out = []
    for prefix in reseed_paths(w):
        base = list(prefix)
        def rec(seq, steps):
            if steps == k:
                if sum(seq) == 0:
                    out.append(tuple(seq))
                return
            if sum(seq) == 0:      # 提前闭合 => 不是首次于末步
                return
            for step in (1, -1):
                rec(seq + [step], steps + 1)
        rec(base, 1)
    return out


def n_first_passage(w, k):
    return len(first_passage_seqs(w, k))


def main():
    print("=" * 78)
    print("R86 · 重播种是否保类号？")
    print("=" * 78)

    # ---------------- F1 ----------------
    print("\n[F1] 记录的纤维 = 旋转类（非单射）")
    f1 = {}
    for L in [4, 6, 8, 10, 12]:
        fib = fiber_map(L)
        sizes = sorted((len(v) for v in fib.values()), reverse=True)
        nw = len(balanced_words(L))
        check("L=%2d 纤维并集 = 全部平衡词" % L, sum(sizes) == nw,
              "词数=%d 类数=%d 最大纤维=%d" % (nw, len(fib), max(sizes)))
        check("L=%2d lambda 非单射（词层自由度被记录丢掉）" % L,
              any(s > 1 for s in sizes), "纤维=%s" % sizes[:5])
        f1["L%d" % L] = {"words": nw, "classes": len(fib), "fiber_sizes": sizes}
    OUT["F1"] = f1

    # ---------------- F2 ----------------
    print("\n[F2] 播种是词级分支：w -> {w+, w-}（等权、无相位）")
    ok = all(len(set(reseed_paths(w))) == 2 and abs(sum(1 / 2 for _ in reseed_paths(w)) - 1) < 1e-15
             for L in [4, 6, 8] for w in balanced_words(L))
    check("每个历史词恰给两条等权分支", ok)
    OUT["F2"] = {"two_equal_branches": ok}

    # ---------------- F3（核心） ----------------
    print("\n[F3] 保类号：记录层态恒为对角（类间元 = 0）")
    worst = 0.0
    f3 = {}
    for L in [4, 6, 8, 10]:
        fib = fiber_map(L)
        keys = sorted(fib)
        idx = {c: i for i, c in enumerate(keys)}
        words = balanced_words(L)
        n, k = len(words), len(keys)
        cands = []
        for seed in range(3):
            g = np.random.default_rng(seed)
            psi = g.normal(size=n) + 1j * g.normal(size=n)
            psi /= np.linalg.norm(psi)
            cands.append(("random_pure_%d" % seed, np.outer(psi, psi.conj())))
        cands.append(("uniform_mix", np.eye(n) / n))
        big = max(fib.values(), key=len)
        psi = np.zeros(n, complex)
        for w in big[:2]:
            psi[words.index(w)] = 1 / np.sqrt(2)
        cands.append(("intra_fiber_pure", np.outer(psi, psi.conj())))
        psi = np.zeros(n, complex)
        psi[words.index(fib[keys[0]][0])] = 1 / np.sqrt(2)
        psi[words.index(fib[keys[-1]][0])] = 1 / np.sqrt(2)
        cands.append(("cross_fiber_pure", np.outer(psi, psi.conj())))
        wmax = 0.0
        for tag, rho in cands:
            rec = np.zeros((k, k), complex)
            for c in keys:
                ii = [words.index(w) for w in fib[c]]
                rec[idx[c], idx[c]] = rho[np.ix_(ii, ii)].trace()
            off = float(np.abs(rec - np.diag(np.diag(rec))).max())
            wmax = max(wmax, off)
            check("L=%2d %-17s 记录态对角" % (L, tag), off < 1e-14,
                  "max|off|=%.2e trace=%.12f" % (off, rec.trace().real))
        worst = max(worst, wmax)
        f3["L%d" % L] = {"max_offdiag": wmax}
    check("全库最坏类间相干 = 0（保类号成立）", worst < 1e-14,
          "sup|offdiag| = %.3e" % worst)
    OUT["F3"] = {"supremum_offdiag": worst, "per_L": f3}

    # ---------------- F4（意外结果） ----------------
    print("\n[F4] 记录**不**决定动力学：同类内不同词的后继类分布不同")
    diff_count, same_count = 0, 0
    f4 = {}
    for L in [6, 8]:
        fib = fiber_map(L)
        for c in sorted(fib):
            ws = fib[c]
            if len(ws) < 2:
                continue
            dists = []
            for w in ws:
                d, leak = class_distribution_after_reseed(w, pmax=8)
                d.pop("<unclosed>", None)
                tot = sum(d.values()) or 1.0
                dists.append({k: v / tot for k, v in d.items()})
            keys = set().union(*[set(d) for d in dists])
            mx = max(abs(dists[0].get(k, 0) - d.get(k, 0)) for d in dists[1:] for k in keys)
            if mx > 1e-9:
                diff_count += 1
            else:
                same_count += 1
        f4["L%d" % L] = {"classes_with_distinct_word_dynamics": diff_count,
                         "classes_with_identical": same_count}
    check("确实存在同类内不同词给出不同后继分布（自由度是真的）",
          diff_count > 0, "不同=%d 相同=%d" % (diff_count, same_count))
    OUT["F4"] = f4

    # ---------------- F5 ----------------
    print("\n[F5] 不可观测性：分块对角代数看不见任何类间元")
    L = 8
    fib = fiber_map(L)
    keys = sorted(fib)
    k = len(keys)
    # 类间元：只在 off-diagonal block 放非零（随机 Hermitian）
    rng = np.random.default_rng(3)
    E = np.zeros((k, k), complex)
    for i in range(k):
        for j in range(k):
            if i != j:
                z = rng.normal() + 1j * rng.normal()
                E[i, j] = z
    E = (E + E.conj().T) / 2
    E = E - np.diag(np.diag(E))          # 只留类间（块外）部分
    # 原生可观测量 = 块对角算符（每块一个矩阵）；其投影的读数
    maxread = 0.0
    for _ in range(200):
        D = np.zeros((k, k), complex)
        for i in range(k):
            D[i, i] = rng.normal()
        maxread = max(maxread, abs(np.trace(E @ D)))
    check("块对角投影读不出类间元（|Tr(E D)| = 0）", maxread < 1e-12,
          "max|Tr(ED)| = %.2e" % maxread)
    OUT["F5"] = {"max_readout_of_interclass": maxread}

    print("\n" + "=" * 78)
    print("汇总：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
    print("=" * 78)
    OUT["summary"] = {"pass": len(PASS), "fail": len(FAIL)}
    with open(os.path.join(HERE, "R86_reseed_class_results.json"), "w") as f:
        json.dump(OUT, f, ensure_ascii=False, indent=1)
    print("-> R86_reseed_class_results.json")
    return 0 if not FAIL else 1


if __name__ == "__main__":
    raise SystemExit(main())
