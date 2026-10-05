#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z12_check.py —— 【E1 收口：有限值标号；类规模因子律；均匀化】核验
==================================================================
独立实断言：
  F1  因子律：类规模取值集 = {2} ∪ {T 的偶因子 ≥ 4}（T=2..16）
  F2  L=4 值集 = {2,4}，归一化 {0.75,1.5}
  F3  对比度上界 2^k（k=4 → 16）
  F4  均匀标号 ⟹ 精确平坦（Z0③ 回归）
  F5  随机二值纹理 ⟹ 块间散度 = 1/√块（CLT 率）
  F6  周期二值纹理 ⟹ 粗粒化精确均匀
  F7  二值场：w 取值数受限（交替标号给 2 个值）
  F8  文档结论与引文在位
"""
import io
import itertools
import os
import sys
from collections import Counter
from math import comb

import numpy as np
from scipy.sparse import csr_matrix

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


DOC = io.open(os.path.join(HERE, "Z12_geometric_input_closed_binary_labeling.md"),
              encoding="utf-8").read()


# ---------------------------------------------------------------- 工具
def balanced(T):
    return [w for w in itertools.product((1, -1), repeat=T) if sum(w) == 0]


def rotclass(w):
    T = len(w)
    return min(tuple(w[i:] + w[:i]) for i in range(T))


def size_set(T):
    cls = {}
    for w in balanced(T):
        cls.setdefault(rotclass(w), 0)
        cls[rotclass(w)] += 1
    return sorted(set(cls.values()))


def predicted(T):
    return sorted({2} | {d for d in range(4, T + 1, 2) if T % d == 0})


def ring(c):
    n = len(c)
    r, cc, v = [], [], []
    for i in range(n):
        j = (i + 1) % n
        r += [i, j]
        cc += [j, i]
        v += [c[i], c[i]]
    return csr_matrix((v, (r, cc)), shape=(n, n))


def W_sum(A, k):
    P = A.copy()
    W = A.multiply(P) * 2
    for m in range(3, k + 1):
        P = (P @ A).tocsr()
        W = W + m * A.multiply(P)
    return W.tocsr()


def wedges(c, k):
    n = len(c)
    W = W_sum(ring(c), k)
    return np.array([W[i, (i + 1) % n] for i in range(n)])


# ---------------------------------------------------------------- F1
head("F1  因子律：类规模取值集 = {2} ∪ {T 的偶因子 ≥ 4}")
for T in (2, 4, 6, 8, 10, 12, 14, 16):
    got, pred = size_set(T), predicted(T)
    check("T=%2d：实测 %s = 预测 %s" % (T, got, pred), got == pred)
check("2 恒在取值集中（交替类）", all(2 in size_set(T) for T in (2, 4, 6, 8, 10, 12, 14, 16)))

# ---------------------------------------------------------------- F2
head("F2  L=4 下的值集 = {2,4}")
avail = []
for T in (2, 4):
    cls = {}
    for w in balanced(T):
        cls.setdefault(rotclass(w), 0)
        cls[rotclass(w)] += 1
    avail += list(cls.values())
vals = sorted(set(avail))
check("可用词长 T≤L=4 的类规模集合 = {2,4}", vals == [2, 4], "%s" % avail)
mean = sum(avail) / len(avail)
norm = sorted(set(round(v / mean, 6) for v in avail))
check("归一化取值 = {0.75, 1.5}（均值 8/3）", abs(mean - 8 / 3) < 1e-12 and norm == [0.75, 1.5],
      "均值 %.6f  值 %s" % (mean, norm))
check("取值数 = 2（二值）", len(norm) == 2)
check("范围 = 2.00", abs(vals[-1] / vals[0] - 2.0) < 1e-12)

# ---------------------------------------------------------------- F3
head("F3  对比度上界 (max/min)^k = 2^k")
bounds = {k: (vals[-1] / vals[0]) ** k for k in (2, 4, 8)}
check("k=2 → 4", abs(bounds[2] - 4) < 1e-9, "%.1f" % bounds[2])
check("k=4 → 16（L=4 的上界）", abs(bounds[4] - 16) < 1e-9, "%.1f" % bounds[4])
check("k=8 → 256", abs(bounds[8] - 256) < 1e-9, "%.1f" % bounds[8])

# ---------------------------------------------------------------- F4
head("F4  均匀标号 ⟹ 精确平坦（Z0③ 回归）")
w = wedges(np.ones(64), 4)
check("均匀标号：边权相对差 = 0", float((w.max() - w.min()) / w.mean()) == 0.0,
      "%.1e" % float((w.max() - w.min()) / w.mean()))

# ---------------------------------------------------------------- F5
head("F5  随机二值纹理 ⟹ 块间散度 = 1/√块")
Nn = 4096
rng = np.random.default_rng(3)
c = np.where(rng.random(Nn) < 0.5, 0.75, 1.5)
w = wedges(c, 4)
rates = []
for B in (4, 8, 16, 32, 64, 128):
    nb = Nn // B
    blocks = np.array([w[i * B:(i + 1) * B].mean() for i in range(nb)])
    s = float(blocks.std() / blocks.mean())
    rates.append((B, s, 1.0 / np.sqrt(B)))
    print("      块=%3d  散度=%.4f  1/√块=%.4f" % (B, s, 1.0 / np.sqrt(B)))
check("散度随块减小", all(rates[i][1] > rates[i + 1][1] for i in range(len(rates) - 1)))
check("散度 ≈ 1/√块（相对偏差 < 20%）",
      all(abs(s / p - 1) < 0.2 for _, s, p in rates),
      "最大偏差 %.1f%%" % (100 * max(abs(s / p - 1) for _, s, p in rates)))

# ---------------------------------------------------------------- F6
head("F6  周期二值纹理 ⟹ 粗粒化精确均匀")
for p in (2, 4, 8):
    cc = np.array([0.75 if (i % p) < p / 2 else 1.5 for i in range(1024)], float)
    ww = wedges(cc, 4)
    B = p * 4
    nb = 1024 // B
    blocks = np.array([ww[i * B:(i + 1) * B].mean() for i in range(nb)])
    spread = float((blocks.max() - blocks.min()) / blocks.mean())
    check("周期 p=%d：块平均极差/均值 = 0（精确均匀）" % p, spread == 0.0, "%.1e" % spread)
    check("周期 p=%d：局部对比度 > 1（局部非均匀）" % p, float(ww.max() / ww.min()) > 1.0,
          "%.3f" % float(ww.max() / ww.min()))

# ---------------------------------------------------------------- F7
head("F7  二值场：w 的取值数受限")
c_alt = np.array([0.75 if i % 2 == 0 else 1.5 for i in range(64)], float)
w_alt = wedges(c_alt, 4)
check("交替标号：w 取值数 = 2", len(set(np.round(w_alt, 9))) == 2,
      "%d" % len(set(np.round(w_alt, 9))))
check("交替标号：对比度 < 上界 16", float(w_alt.max() / w_alt.min()) < 16.0,
      "%.3f" % float(w_alt.max() / w_alt.min()))
c_rnd = np.where(rng.random(64) < 0.5, 0.75, 1.5)
check("随机标号：对比度仍 < 上界 16", float(wedges(c_rnd, 4).max() / wedges(c_rnd, 4).min()) < 16.0,
      "%.3f" % float(wedges(c_rnd, 4).max() / wedges(c_rnd, 4).min()))

# ---------------------------------------------------------------- F8
head("F8  文档结论与引文在位")
check("结论盒：值集 ＝ {2}∪偶因子／二值", "偶因子" in DOC and "二值" in DOC and "均匀化" in DOC)
check("写明 L=4 ⟹ 二值", "二值" in DOC and "0.75" in DOC)
check("写明 E1 收口（自由场 → 有限值标号）", "收口" in DOC and "有限值标号" in DOC)
check("写明均匀化三例（均匀／周期／随机）", "均匀化" in DOC and "周期" in DOC and "随机" in DOC)
check("登记带价条款 Z-E5", "Z-E5" in DOC and "买" in DOC and "价" in DOC)
check("给出预言（共形因子只取两值）", "预言" in DOC and "两个值" in DOC)
check("诚实边界在位（归纳律／L=4 等级／对齐）", "归纳" in DOC and "对齐" in DOC)
for fn, keys in [
    ("G61_locking_the_five_integers.md", ["L=4"]),
    ("G46_k_is_the_lifetime.md", ["顶点传递"]),
    ("zero_sum_rotation_class_algebra.md", ["旋转类"]),
    ("Z11_correction_k_is_fixed_and_age_is_not_site.md", ["与寿命无关"]),
    ("Z9_pi_filter_and_lifetime_fork.md", ["旋转类"]),
    ("Z8_native_scale_field_candidate.md", ["推前重数"]),
    ("Z7_embedding_input_explicit_dictionary.md", ["闭式字典"]),
    ("G58_I2a_resolved_as_embedding_input.md", ["Γ-收敛"]),
]:
    p = os.path.join(HERE, fn)
    txt = io.open(p, encoding="utf-8").read() if os.path.exists(p) else ""
    miss = [k for k in keys if k not in txt]
    check("引文在位：%s" % fn, os.path.exists(p) and not miss, "缺 %s" % miss if miss else "")

# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
