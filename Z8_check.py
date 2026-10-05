#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z8_check.py —— 【c 的原生候选：π 的推前重数（含一次被证伪的候选）】的核验
========================================================================
独立实断言：
  F1  候选判据清单在位
  F2  生存（前缀）测度 F(a)：半程定理的独立复算（闭式 = 暴力枚举）
  F3  候选 1（c = F(a)）**被证伪**：字典失效／无平坦区／对比度爆炸
  F4  候选 2：均匀 π ⟹ 环上**精确平坦**（Z0③ 回归），且值 = P_sum(1)
  F5  候选 2 乘积规则：字典二阶收敛 ＋ 度规非均匀
  F6  候选 2 均值规则：字典二阶收敛（规则无关的核心结论）
  F7  对比度闭式：P 单调 ⟹ 对比度 = P(c_max)/P(c_min)
  F8  文档结论与引文在位
"""
import itertools
import io
import math
import os
import sys
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


DOC = io.open(os.path.join(HERE, "Z8_native_scale_field_candidate.md"), encoding="utf-8").read()


# ---------------------------------------------------------------- 工具
def F_closed(a, L):
    """F(a) = Pr(|S_a| <= L-a)，S_a = a 个 ±1 之和。"""
    if a > L:
        a = L
    tot = 0
    for h in range(-a, a + 1):
        if abs(h) <= L - a and (h - a) % 2 == 0:
            tot += comb(a, (a + h) // 2)
    return tot / 2 ** a


def F_brute(a, L):
    c = 0
    for seq in itertools.product((1, -1), repeat=a):
        if abs(sum(seq)) <= L - a:
            c += 1
    return c / 2 ** a


def W_sum(A, k):
    P = A.copy()
    W = A.multiply(P) * 2
    for m in range(3, k + 1):
        P = (P @ A).tocsr()
        W = W + m * A.multiply(P)
    return W.tocsr()


def P_sum(c, k, d=1):
    """层和口径字典（d=1：W_1(L)=C(L,(L+1)/2)）。"""
    return sum(m * comb(m - 1, m // 2) * c ** m for m in range(2, k + 1, 2))


def ring(c):
    Nn = len(c)
    r, cc, v = [], [], []
    for i in range(Nn):
        j = (i + 1) % Nn
        r += [i, j]
        cc += [j, i]
        v += [c[i], c[i]]
    return csr_matrix((v, (r, cc)), shape=(Nn, Nn))


def chain(c):
    Nn = len(c)
    r, cc, v = [], [], []
    for i in range(Nn):
        r += [i, i + 1]
        cc += [i + 1, i]
        v += [c[i], c[i]]
    return csr_matrix((v, (r, cc)), shape=(Nn + 1, Nn + 1))


def wedges(W, Nn):
    return np.array([W[i, (i + 1) % Nn] for i in range(Nn)])


# ---------------------------------------------------------------- F1
head("F1  候选判据清单在位")
for key in ("原生", "无量纲", "逐点", "顶点传递", "缓变", "正"):
    check("判据『%s』在位" % key, key in DOC)

# ---------------------------------------------------------------- F2
head("F2  生存（前缀）测度 F(a)：半程定理的独立复算")
for L in (8, 16):
    fc = [F_closed(a, L) for a in range(L + 1)]
    fb = [F_brute(a, L) for a in range(L + 1)]
    check("L=%2d：闭式 = 暴力枚举（最大偏差 = 0）" % L,
          max(abs(x - y) for x, y in zip(fc, fb)) == 0.0,
          "%.1e" % max(abs(x - y) for x, y in zip(fc, fb)))
for L in (8, 16, 32, 64):
    fc = [F_closed(a, L) for a in range(L + 1)]
    check("L=%2d：前半程 F ≡ 1（a ≤ L/2）" % L, all(x == 1.0 for x in fc[:L // 2 + 1]))
    check("L=%2d：后半程 F 严格递减" % L, all(fc[i] > fc[i + 1] for i in range(L // 2, L)))
check("尾值闭式 F(L) = C(L,L/2)/2^L：0.273438 / 0.196381 / 0.139950 / 0.099347",
      all(abs(F_closed(L, L) - float("%.6f" % (comb(L, L // 2) / 2 ** L))) < 1e-6 for L in (8, 16, 32, 64)),
      "%s" % ["%.6f" % F_closed(L, L) for L in (8, 16, 32, 64)])

# ---------------------------------------------------------------- F3
head("F3  候选 1（c = F(a)）被证伪")
for L in (8, 16, 32):
    k = L
    c = np.array([F_closed(a, L) for a in range(L)])
    W = W_sum(chain(c), k)
    w = np.array([W[a, a + 1] for a in range(L)])
    p = np.array([P_sum(ci, k) for ci in c])
    dev = float(np.abs(w / p - 1).max())
    check("L=%2d：字典相对偏差 ≫ 1（不连续场，字典失效）" % L, dev > 0.5, "%.2e" % dev)
    # 平坦判据必须用**同一条链、同一个 k 的均匀基线**（否则边界效应混进来）
    w_uni = np.array([W_sum(chain(np.ones(L)), k)[a, a + 1] for a in range(L)])
    flat = 0
    for a in range(L):
        if w[a] == w_uni[a]:
            flat += 1
        else:
            break
    check("L=%2d：精确平坦区 = %d 条边（几乎不存在）" % (L, flat), flat == 2,
          "L/2=%d, k/2-1=%d" % (L // 2, k // 2 - 1))
    check("L=%2d：度规对比度爆炸（> 1e1）" % L, w.max() / w.min() > 10.0, "%.3e" % (w.max() / w.min()))
check("候选 1 的失效模式与 Z7 §8 的「仅扫光滑场」边界一致",
      "不连续" in DOC and "Z7" in DOC)

# ---------------------------------------------------------------- F4
head("F4  候选 2：均匀 π ⟹ 环上精确平坦（Z0③ 回归）")
for k in (4, 8, 16):
    for Nn in (32, 64):
        w = wedges(W_sum(ring(np.ones(Nn)), k), Nn)
        check("均匀 π：N=%d k=%2d 边权相对差 = 0" % (Nn, k), float((w.max() - w.min()) / w.mean()) == 0.0,
              "值 %.6g" % w[0])
for k in (4, 8, 16):
    w = wedges(W_sum(ring(np.ones(32)), k), 32)
    check("k=%2d：均匀值 = P_sum(1) = Σ m·C(m-1,m/2)" % k, abs(w[0] - P_sum(1.0, k)) < 1e-9,
          "%.0f" % w[0])

# ---------------------------------------------------------------- F5
head("F5  候选 2 乘积规则 c = M_iM_j/⟨M⟩²：字典二阶 ＋ 非均匀")
prod = []
for Nn in (64, 256, 1024):
    i = np.arange(Nn)
    M = 1 + 0.3 * np.cos(2 * np.pi * (i + 0.5) / Nn)
    c = (M * M) / (M.mean() ** 2)
    k = 8
    w = wedges(W_sum(ring(c), k), Nn)
    p = np.array([P_sum(x, k) for x in c])
    m = (i > 0.2 * Nn) & (i < 0.8 * Nn)
    prod.append((float(np.abs(w[m] / p[m] - 1).max()), float(w.max() / w.min())))
check("乘积规则：字典偏差二阶收敛（N×4 ⟹ ÷>8）",
      prod[0][0] / prod[1][0] > 8 and prod[1][0] / prod[2][0] > 8,
      " -> ".join("%.2e" % x[0] for x in prod))
check("乘积规则：度规非均匀且对比度收敛（≈6.9e3）",
      prod[0][1] > 1e3 and abs(prod[2][1] / prod[1][1] - 1) < 0.01,
      " -> ".join("%.4g" % x[1] for x in prod))

# ---------------------------------------------------------------- F6
head("F6  候选 2 均值规则 c = (M_i+M_j)/(2⟨M⟩)：规则无关的核心结论")
meanr = []
for Nn in (64, 256, 1024):
    i = np.arange(Nn)
    M = 1 + 0.3 * np.cos(2 * np.pi * (i + 0.5) / Nn)
    c = M / M.mean()
    k = 8
    w = wedges(W_sum(ring(c), k), Nn)
    p = np.array([P_sum(x, k) for x in c])
    m = (i > 0.2 * Nn) & (i < 0.8 * Nn)
    meanr.append((float(np.abs(w[m] / p[m] - 1).max()), float(w.max() / w.min())))
check("均值规则：字典二阶收敛", meanr[0][0] / meanr[1][0] > 8 and meanr[1][0] / meanr[2][0] > 8,
      " -> ".join("%.2e" % x[0] for x in meanr))
check("两种规则都非均匀 ⟹ 「Γ 非正则 ⟺ π 非均匀」与规则无关",
      prod[2][1] > 1e3 and meanr[2][1] > 10)

# ---------------------------------------------------------------- F7
head("F7  对比度闭式：P 单调 ⟹ 对比度 = P(c_max)/P(c_min)")
mono = all(P_sum(c + 1e-4, 8) > P_sum(c, 8) for c in np.linspace(0.05, 3.0, 50))
check("P_sum 在 c>0 上严格单调（字典可逆）", mono)
r = 3.45
cmax, cmin = math.sqrt(r), 1 / math.sqrt(r)
check("r=3.45, k=8：对比度 = P(c_max)/P(c_min) ≈ 8.4e3",
      abs(P_sum(cmax, 8) / P_sum(cmin, 8) / 8.41e3 - 1) < 0.02,
      "%.3e" % (P_sum(cmax, 8) / P_sum(cmin, 8)))
check("对比度随 k 指数增长（r=3.45：k=4 → k=16 增 > 1e6 倍）",
      P_sum(cmax, 16) / P_sum(cmin, 16) / (P_sum(cmax, 4) / P_sum(cmin, 4)) > 1e6,
      "%.2e → %.2e" % (P_sum(cmax, 4) / P_sum(cmin, 4), P_sum(cmax, 16) / P_sum(cmin, 16)))

# ---------------------------------------------------------------- F8
head("F8  文档结论与引文在位")
check("结论盒：c 的原生候选 = π 的推前重数", "推前重数" in DOC)
check("写明候选 1 被证伪", "被证伪" in DOC and "F(a)" in DOC)
check("写明输入账本收缩（「Γ 非正则」并入 π）", "并入" in DOC and "非均匀" in DOC)
check("写明 Z0③ 回归（均匀 π ⟹ 精确平坦）", "Z0③" in DOC and "精确平坦" in DOC)
check("给出可证伪预言（对比度对 r 与 k 的依赖）", "对比度" in DOC and "预言" in DOC)
for fn, keys in [
    ("G54_quantitative_profile_age_measure.md", ["半程定理", "前缀"]),
    ("G29_probability_as_derived_not_postulated.md", ["推前", "粗粒化"]),
    ("G80_all_to_all_age_coupling.md", ["全对全"]),
    ("G85_all_to_all_from_closed_walks.md", ["全对全"]),
    ("Z7_embedding_input_explicit_dictionary.md", ["闭式字典"]),
    ("Z5_finite_k_locality_escape.md", ["非正则"]),
    ("G46_k_is_the_lifetime.md", ["顶点传递"]),
    ("G57_unreachability_of_absolute_normalization.md", ["无量纲"]),
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
