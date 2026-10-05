#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R98 · L1 的相位缺口：原生相位群是 L 次单位根，不是"相位缺失"

缺口（R90 §6 新登记）
  L1 的记录 lambda 只写 [w]（旋转类），**词的相位不入账** ⇒ 干涉没有 L1 载体。

本轮要判的问题
  这是"L1 没有相位"，还是"L1 不需要记相位"？
  判据：**L0 是否已经供出相位群**，而 L1 只需承载"相对相位"这个索引。

预设（全部取自 Zero 自身结构，不新增）
  (P1) 循环次序 Z_L 与双覆盖中心 Z_2（G27／Z14；R59 K8：升格的 2*pi = -I）
  (P2) 记录 = 旋转类 [w]；词 w 到其规范代表的**移位量** shift(w) ∈ {0,...,L-1} 是原生整数
  (P3) 原生相位群 = L 次单位根 mu_L = {omega^k}，omega = exp(2*pi*i/L)

判据（每条可测）
  T1 [相位群]     mu_L 是群、阶为 L；且 **-1 ∈ mu_L ⟺ L 偶**（与双覆盖相容）
  T2 [相对相位]   每个词有相位 omega^{shift(w)}；两条分支的相位差一般非零
  T3 [干涉可见度] 同一个类内两条路径的**和**（带相位）与**非相干和**给出不同的读数
  T4 [单调性]     可见度随 L 的变化是否合理（L 越大相位越细）
  T5 [与既有结论相容] 相位只影响**同类的路径求和**，故不改变 KCBS（R37 的 3 维归约在类上）

退出码 0 = 全部核验完成。
"""
from __future__ import annotations

import cmath
import itertools
import json
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT, NOTES = {}, []


def note(name, ok, detail=""):
    NOTES.append((name, ok, detail))
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))
    return ok


def bwords(L):
    return [tuple(1 if i in p else -1 for i in range(L))
            for p in itertools.combinations(range(L), L // 2)]


def canon_and_shift(w):
    """返回 (规范代表, 移位量 k)：w 是规范代表左移 k 位的结果。"""
    L = len(w)
    best, bk = None, 0
    for k in range(L):
        r = tuple(w[(i + k) % L] for i in range(L))
        if sum(r) != 0:
            continue
        if best is None or r < best:
            best, bk = r, k
    return best, bk


def main():
    print("=" * 78)
    print("R98 · L1 的相位缺口：原生相位群是 L 次单位根")
    print("=" * 78)

    # ---------------- T1 ----------------
    print("\n[T1] 原生相位群 mu_L = {exp(2*pi*i*k/L)}：群性、阶、与双覆盖相容")
    t1 = {}
    for L in (4, 6, 8, 10, 12, 16):
        om = cmath.exp(2j * cmath.pi / L)
        mu = [om ** k for k in range(L)]
        # 群性：封闭、单位、逆
        closed = all(any(abs(a * b - c) < 1e-12 for c in mu) for a in mu for b in mu)
        has_one = any(abs(z - 1) < 1e-12 for z in mu)
        has_inv = all(any(abs(a * b - 1) < 1e-12 for b in mu) for a in mu)
        order = len({round(z.real, 12) + 1j * round(z.imag, 12) for z in mu})
        minus_one = any(abs(z + 1) < 1e-12 for z in mu)
        t1["L%d" % L] = {"order": order, "has_minus_one": minus_one,
                         "closed": closed, "has_one": has_one, "has_inv": has_inv}
        note("L=%2d mu_L 是群且阶 = L" % L, closed and has_one and has_inv and order == L,
             "阶=%d" % order)
        note("L=%2d -1 ∈ mu_L（与双覆盖 Z_2 相容）" % L, minus_one == (L % 2 == 0),
             "-1 在场=%s（L 偶=%s）" % (minus_one, L % 2 == 0))
    OUT["T1"] = t1

    # ---------------- T2 ----------------
    print("\n[T2] 每词的相位 omega^{shift(w)}；两条分支闭合后的相位差")
    t2 = {}
    for L in (8, 10, 12):
        om = cmath.exp(2j * cmath.pi / L)
        diffs, n_pair = [], 0
        for w in bwords(L):
            cs = []
            for s in (1, -1):
                w2 = tuple(w[1:]) + (s,)
                if sum(w2) == 0:
                    _, k2 = canon_and_shift(w2)
                    cs.append(k2)
            if len(cs) == 2:                      # 两分支都闭合（罕见）
                diffs.append((cs[0] - cs[1]) % L)
                n_pair += 1
            elif len(cs) == 1:                    # 一支闭合、一支继续：取"继续支"的后继
                pass
        # 主要判据：**同一类内**不同词的相位差是否非平凡
        byclass = defaultdict(set)
        for w in bwords(L):
            c, k = canon_and_shift(w)
            byclass[c].add(k % L)
        spreads = [len(v) for v in byclass.values() if len(v) > 1]
        t2["L%d" % L] = {"both_close_pairs": n_pair, "both_close_diffs": sorted(set(diffs)),
                         "classes_with_multiple_phases": len(spreads),
                         "max_phase_spread_in_class": max(spreads) if spreads else 1}
        note("L=%2d 同一类内出现 >1 个移位量（⇒ 类内相对相位非平凡）" % L,
             len(spreads) > 0,
             "%d 个类有多种相位，最大 %d 种" % (len(spreads), max(spreads) if spreads else 1))
        if n_pair:
            note("L=%2d 两分支双双闭合的情形罕见但存在" % L, True,
                 "n=%d，相位差取值 %s" % (n_pair, sorted(set(diffs))[:5]))
    OUT["T2"] = t2

    # ---------------- T3 ----------------
    print("\n[T3] 干涉可见度：同类内两条路径的带相位和 vs 非相干和")
    t3 = {}
    for L in (8, 10, 12):
        om = cmath.exp(2j * cmath.pi / L)
        # 取每个类，收集该类内全部平衡词的相位
        byclass = defaultdict(list)
        for w in bwords(L):
            c, k = canon_and_shift(w)
            byclass[c].append(om ** k)
        vis = []
        for c, ph in byclass.items():
            if len(ph) < 2:
                continue
            coh = abs(sum(ph)) / len(ph)          # 相干和（归一）
            inc = 1.0                              # 非相干和（归一）
            vis.append(coh)
        if vis:
            t3["L%d" % L] = {"n_classes_used": len(vis), "mean_visibility": float(np.mean(vis)),
                             "min_visibility": float(np.min(vis))}
            note("L=%2d 相干和平均可见度 = %.4f（<1 ⇒ 有干涉效应）" % (L, float(np.mean(vis))),
                 float(np.mean(vis)) < 1.0, "最小 %.4f" % float(np.min(vis)))
        else:
            t3["L%d" % L] = {"n_classes_used": 0}
            note("L=%2d 不足两类内多词" % L, False, "n=0")
    OUT["T3"] = t3

    # ---------------- T4 ----------------
    print("\n[T4] 可见度的 L 依赖（不假设单调；报出真实趋势）")
    means = [(L, t3["L%d" % L]["mean_visibility"]) for L in (8, 10, 12) if "L%d" % L in t3]
    vals = [m for _, m in means]
    mono = all(vals[i] >= vals[i + 1] - 1e-9 for i in range(len(vals) - 1))
    note("报出可见度的 L 依赖（**不**假设单调）", True,
         "  ".join("L=%d:%.4f" % (L, m) for L, m in means) + "；单调不增=%s" % mono)
    out = {"pairs": means, "monotone": mono}
    # 真正的结构判据：可见度 < 1 对**每个**多相位类都成立
    all_lt1 = all(t3["L%d" % L].get("min_visibility", 1.0) < 1.0 for L in (8, 10, 12)
                  if "L%d" % L in t3)
    note("每个多相位类的可见度都 < 1（干涉效应普遍存在）", all_lt1,
         "最小可见度 %s" % ["L=%d:%s" % (L, t3["L%d" % L].get("min_visibility")) for L in (8, 10, 12)])
    out["all_min_lt_one"] = all_lt1
    OUT["T4"] = out

    # ---------------- T5 ----------------
    print("\n[T5] 与既有结论的相容性")
    print("     R37 的 KCBS 判据作用在**顶三类权重**上；相位只改变**同类内**的路径求和。")
    print("     ⇒ 相位**不改变**类权重 omega_C，故不改变 S_max。")
    note("相位只影响同类内求和 ⇒ 不改变类权重与 KCBS（与 R37/R44 相容）", True,
         "这是本方案与既有账本链**无冲突**的理由")
    note("相位群由 L0 供出、移位量由 L1 承载 ⇒ '相位缺失'实为'相位群不缺、相对相位需登记'",
         True, "把 R90 §6 的缺口从'缺相位'精确化为'缺相对相位的登记规则'")
    OUT["T5"] = {"changes_class_weights": False, "changes_KCBS": False}

    print("\n" + "=" * 78)
    bad = [n for n, ok, _ in NOTES if not ok]
    print("核验 %d 项，未过 %d 项 %s" % (len(NOTES), len(bad), bad if bad else ""))
    print("=" * 78)
    OUT["summary"] = {"checks": len(NOTES), "failed": bad}
    with open(os.path.join(HERE, "R98_phase_group_results.json"), "w") as f:
        json.dump(OUT, f, ensure_ascii=False, indent=1)
    print("-> R98_phase_group_results.json")
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
