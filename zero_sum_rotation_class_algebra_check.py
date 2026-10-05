#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_sum_rotation_class_algebra_check.py —— Zero 文章①【旋转类代数】的核验
=========================================================================
独立实断言（全部**从零重算**，不 import 任何 Zero 程序）：
  F1  规范代表 canon(w) 是类不变量（同类的词给同一个 canon）
  F2  项链计数闭式（含 d | gcd(L, L/2) 限制）与枚举逐项一致
  F3  内部归零数 r(w) 的分布与文章表格逐项一致，且各行之和 = K(L)
  F4  三套闭合代数 M = 1 / 1+r / 2^r 的总重数与文章表格一致
  F5  文章写明了"三套只是代数、不是物理后代数"（与②的 no-go 呼应）
"""
import io
import itertools
import os
import sys
from collections import Counter
from math import comb, gcd

HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")
FAIL = 0
N = 0


def check(name, cond, detail="", level="ind"):
    global FAIL, N
    N += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    tag = ("v" if ok else "x") if (level == "ind" or LEDGER_MODE == "A" or not ok) else "i"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))


def head(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


def rd(p):
    return io.open(os.path.join(HERE, p), encoding="utf-8").read()


DOC = rd("zero_sum_rotation_class_algebra.md")


def canon(w):
    return min(w[i:] + w[:i] for i in range(len(w)))


def r_count(w):
    total = returns = 0
    for v in w[:-1]:
        total += v
        if total == 0:
            returns += 1
    return returns


def phi(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def K_closed(L):
    g = gcd(L, L // 2)
    return sum(phi(d) * comb(L // d, L // (2 * d)) for d in range(1, g + 1) if g % d == 0) // L


def words(L):
    return [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]


# ---------------------------------------------------------------- F1
head("F1  canon(w) 是旋转类不变量")
bad = tot = 0
for L in (4, 6, 8):
    for w in words(L):
        for k in range(L):
            rot = w[k:] + w[:k]
            tot += 1
            if canon(rot) != canon(w):
                bad += 1
check("全部旋转与原名给同一个 canon", bad == 0, "测 %d 个旋转，反例 %d" % (tot, bad))
check("canon 的结果是原词的某个旋转", all(
    canon(w) in [w[k:] + w[:k] for k in range(len(w))] for w in words(6)))

# ---------------------------------------------------------------- F2
head("F2  项链计数闭式 vs 枚举（文章 §2 的表）")
REF_K = {2: 1, 4: 2, 6: 4, 8: 10, 10: 26, 12: 80}
ok_all = True
for L, ref in sorted(REF_K.items()):
    enum = len({canon(w) for w in words(L)})
    clo = K_closed(L)
    ok = (enum == ref == clo)
    ok_all &= ok
    check("L=%-3d 枚举 %-3d = 闭式 %-3d = 文章 %-3d" % (L, enum, clo, ref), ok)
_L10 = 10
_g = gcd(_L10, _L10 // 2)
check("闭式带 d | gcd(L,L/2) 限制（L=10：正确 26；误取 d | L 给 27）",
      sum(phi(d) * comb(_L10 // d, _L10 // (2 * d))
          for d in range(1, _g + 1) if _g % d == 0) // _L10 == 26)
check("文章写明该限制", "d\\mid\\gcd" in DOC or "d\\,\\mid\\,\\gcd" in DOC)

# ---------------------------------------------------------------- F3
head("F3  r(w) 分布（文章 §3 的表）")
REF_R = {2: {0: 1}, 4: {0: 1, 1: 1}, 6: {0: 2, 1: 1, 2: 1},
         8: {0: 5, 1: 3, 2: 1, 3: 1}, 12: {0: 37, 1: 21, 2: 13, 3: 7, 4: 1, 5: 1}}
for L, ref in sorted(REF_R.items()):
    dist = dict(sorted(Counter(r_count(c) for c in {canon(w) for w in words(L)}).items()))
    check("L=%-3d r 分布 = %s" % (L, ref), dist == ref, "实算 %s" % dist)
    check("L=%-3d 各行之和 = K(L) = %d" % (L, REF_K[L]), sum(dist.values()) == REF_K[L])

# ---------------------------------------------------------------- F4
head("F4  三套闭合代数总重数（文章 §4 的表）")
REF_M = {4: (2, 3, 3), 8: (10, 18, 23), 12: (80, 157, 235)}
for L, (m1, m2, m3) in sorted(REF_M.items()):
    rs = [r_count(c) for c in {canon(w) for w in words(L)}]
    got = (len(rs), sum(1 + x for x in rs), sum(2 ** x for x in rs))
    check("L=%-3d (sterile, single_cut, all_cuts) = %s" % (L, (m1, m2, m3)),
          got == (m1, m2, m3), "实算 %s" % (got,))
check("文章给出三条公式 M=1 / 1+r / 2^r",
      "M_{\\rm sterile}=1" in DOC and "1+r(w)" in DOC and "2^{\\,r(w)}" in DOC)

# ---------------------------------------------------------------- F5
head("F5  文章与②的 no-go 呼应")
check("文章指出三套不是物理后代数", "不是物理后代数" in DOC or "分支程序计数" in DOC)
check("文章指向繁殖转移定理一文",
      "zero_sum_reproduction_transition_theorems.md" in DOC)

# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
