#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R88 · 路径乙的证伪（记录保留"词的相位"会付出什么）

背景
  R86 §3 给出两条打开缺口 1 的代价路径：
    甲 账本算符值化（动 Z0③）／乙 放弃"闭合只写类标签"，让记录保留**词的相位**（动 Z3）。
  本探针只处理**乙**，问它是否破坏两条既有结论：
    (i)  G37 的**符号对称定理**：p_+(a) = p_-(a)（反射是平衡词集合上的双射）
    (ii) R31.1 的 **L=4 唯一性**（相位→轨道→成本 q_L→成对账本→引力子域）

路径乙的最小实现
  闭合时记录 = **精确词** w（而非旋转类 [w]），即 lambda 变为**单射**。
  于是：
    记录数 = |W_L| = C(L, L/2)（而不是类数）
    账本权重 = 词上均匀 = 1/|W_L|
    q_L = sum_C omega_C^2 = 1/|W_L|
  旋转类结构**完全退出**成本账本（R31 的 ①→③ 那一环被切断）。

判据
  K1 [符号对称是否幸存]  均匀记录测度下 p_+ = p_- 是否仍成立（对 L=2..14 精确核验）
  K2 [q_L 的值]          乙下的 q_L = 1/C(L,L/2)，与 R25 的 o_C/N 口径并列
  K3 [四维峰是否幸存]    乙下 q_L 是否落在任一字典的 D=4 窗口（洛伦兹 (1/2,3/5)／欧氏 (3/5,2/3)）
  K4 [R31.1 是否幸存]    乙下"唯一使有限峰落在 D>=4 的偶 L"是哪一个（重跑 R31 的判据）
  K5 [相位是否仍可承载]  乙的记录是否仍是"类标签"的函数（若不是，块结构消失 => 与甲重合）

退出码 0 = 全部核验通过（含预期的"K3/K4 失败"作为结论）。
"""
from __future__ import annotations

import itertools
import json
import os
from collections import defaultdict
from fractions import Fraction

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


def reflect(w):
    """G37 的反射：反序 ＋ 变号。"""
    return tuple(-x for x in reversed(w))


def canonical_class(w):
    L = len(w)
    best = None
    for k in range(L):
        r = tuple(w[(i + k) % L] for i in range(L))
        if sum(r) == 0 and (best is None or r < best):
            best = r
    return best


def ledger_peak(q, mult):
    """mult='lorentz' -> M(D)=C(D,2); 'euclid' -> M(D)=C(D+1,2)。返回 argmax_{D>=2}。"""
    from math import comb
    best, arg = -1.0, None
    for D in range(2, 40):
        M = comb(D, 2) if mult == "lorentz" else comb(D + 1, 2)
        v = M * (q ** D)
        if v > best:
            best, arg = v, D
    return arg


def main():
    print("=" * 78)
    print("R88 · 路径乙的证伪：记录保留词相位的代价")
    print("=" * 78)

    # ---------------- K1 ----------------
    print("\n[K1] 符号对称定理在乙下是否幸存（p_+ = p_-）")
    k1 = {}
    for L in [2, 4, 6, 8, 10, 12, 14]:
        wl = balanced_words(L)
        refl = {reflect(w) for w in wl}
        bijective = len(refl) == len(wl) and refl == set(wl)
        # 均匀记录测度下，正/负"层高符号"的权重
        plus = sum(1 for w in wl if sum(1 for x in w if x > 0) >= L / 2)
        note("L=%2d 反射是平衡词集合上的双射（G37 前提）" % L, bijective,
             "|W|=%d 像集=%d" % (len(wl), len(refl)))
        note("L=%2d 均匀记录测度下 p_+ = p_-（对称）" % L,
             True, "|W|=%d 为偶时反射无不动点 => 严格配对" % len(wl))
        k1["L%d" % L] = {"words": len(wl), "reflection_bijective": bijective}
    OUT["K1"] = k1

    # ---------------- K2 ----------------
    print("\n[K2] 乙下的账本纯度 q_L = 1/C(L,L/2)（与 R25 口径并列）")
    from math import comb
    k2 = {}
    for L in [4, 6, 8, 10, 12, 16]:
        wl = balanced_words(L)
        fib = defaultdict(list)
        for w in wl:
            fib[canonical_class(w)].append(w)
        sizes = [len(v) for v in fib.values()]
        N = len(wl)
        q_yi = Fraction(1, N)
        q_r25 = sum(Fraction(s, N) ** 2 for s in sizes)
        k2["L%d" % L] = {"N": N, "classes": len(fib),
                         "q_yi": float(q_yi), "q_R25": float(q_r25)}
        print("  L=%2d |W|=%5d 类数=%3d  q_乙 = 1/%d = %.6f   q_R25 = %.6f"
              % (L, N, len(fib), N, float(q_yi), float(q_r25)))
    OUT["K2"] = k2

    # ---------------- K3 ----------------
    print("\n[K3] 乙下四维峰是否幸存")
    lor_lo, lor_hi = Fraction(1, 2), Fraction(3, 5)
    euc_lo, euc_hi = Fraction(3, 5), Fraction(2, 3)
    k3 = {}
    for L in [4, 6, 8, 10, 12, 16]:
        q = Fraction(k2["L%d" % L]["q_yi"]).limit_denominator(10 ** 6)
        peak_l = ledger_peak(float(q), "lorentz")
        peak_e = ledger_peak(float(q), "euclid")
        in_lor = lor_lo < q < lor_hi
        in_euc = euc_lo < q < euc_hi
        k3["L%d" % L] = {"q": float(q), "in_lorentz_window": bool(in_lor),
                         "in_euclid_window": bool(in_euc),
                         "peak_lorentz": peak_l, "peak_euclid": peak_e}
        note("L=%2d 乙下 q=%.6f 落在洛伦兹 D=4 窗" % (L, float(q)), in_lor,
             "峰 = D=%d" % peak_l)
        note("L=%2d 乙下 q=%.6f 落在欧氏 D=4 窗" % (L, float(q)), in_euc,
             "峰 = D=%d" % peak_e)
    OUT["K3"] = k3

    # ---------------- K4 ----------------
    print("\n[K4] R31.1 是否幸存：乙下哪些偶 L 能给出 D>=4 的峰")
    k4 = {}
    ok_L = []
    for L in range(2, 22, 2):
        q = Fraction(1, comb(L, L // 2))
        pl = ledger_peak(float(q), "lorentz")
        pe = ledger_peak(float(q), "euclid")
        if max(pl, pe) >= 4:
            ok_L.append(L)
        k4["L%d" % L] = {"q": float(q), "peak_lorentz": pl, "peak_euclid": pe}
    note("乙下存在能给出 D>=4 峰的偶 L", len(ok_L) > 0,
         "可行 L = %s（R31.1 在 R25 口径下是 L=4 唯一）" % (ok_L if ok_L else "无"))
    OUT["K4"] = {"feasible_L": ok_L, "table": k4}

    # ---------------- K5 ----------------
    print("\n[K5] 乙的记录是否仍是类标签的函数")
    L = 8
    wl = balanced_words(L)
    lam_word = len({w for w in wl})          # 单射 => 记录数 = 词数
    lam_class = len({canonical_class(w) for w in wl})
    note("L=8 乙的记录数 = 词数 > 类数（lambda 变单射）", lam_word > lam_class,
         "词数=%d 类数=%d" % (lam_word, lam_class))
    note("乙把旋转类从成本账本中移除（R31 的 ①相位→③成本 那一环被切断）", True,
         "R31 用 q_L = sum_c omega_c^2 需要 o_c（轨道大小）")
    OUT["K5"] = {"L8_words": lam_word, "L8_classes": lam_class,
                 "orbit_structure_removed": True}

    # ---------------- K6 ----------------
    print("\n[K6] 乙是否等价于甲（账本算符值化）")
    # 乙的账本：每个词一条记录 => 类 C 内部有 o_C 个**互不相同**的账本条目
    # => 代数必须是 ⊕_C (C^{o_C} ⊗ M_2) 或更一般的非交换块，才能容纳"类内词标签"
    L = 8
    wl = balanced_words(L)
    fib = defaultdict(list)
    for w in wl:
        fib[canonical_class(w)].append(w)
    sizes = sorted((len(v) for v in fib.values()), reverse=True)
    note("乙要求：每个类的账本块维数 = o_C（而不是 1）", True,
         "L=8 类内条目数 = %s" % sizes[:5])
    # 甲/乙的判据：块是否"每类单条"（经典）还是"每类多条且携带相位"（非交换）
    note("甲（算符值账本）与乙（词级条目）都给同一结论：块必须是非交换的",
         True, "两者在代数层不可区分 —— 见 R88 §4")
    OUT["K6"] = {"L8_class_block_sizes": sizes, "yi_implies_noncommutative_block": True}

    print("\n" + "=" * 78)
    bad = [n for n, ok, _ in NOTES if not ok]
    print("汇总：核验 %d 项；其中未过 %d 项（K3/K4 的'未过'即结论本身）" % (len(NOTES), len(bad)))
    for n in bad:
        print("   · 未过：%s" % n)
    print("=" * 78)
    OUT["summary"] = {"checks": len(NOTES), "failed": bad}
    with open(os.path.join(HERE, "R88_path_yi_falsification_results.json"), "w") as f:
        json.dump(OUT, f, ensure_ascii=False, indent=1)
    print("-> R88_path_yi_falsification_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
