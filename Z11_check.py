#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z11_check.py —— 【更正 Z10 §3：k=L 是固定步数；年龄 ≠ 站点】核验
==================================================================
独立实断言：
  F1  k 固定（=4,8,16）＋空间格 N：大范围场的字典偏差与对比度（多项式、非退化）
  F2  k = N（＝"年龄＝站点"场景）：退化（复现 Z10 §3 的数）
  F3  对比度 = range^k 的增长律（多项式于 range、指数于 k）
  F4  "年龄＝站点 ⟹ k ∝ N" 的环节核验（含 G46/G61/G58 的引文）
  F5  更正后的 C7 表述与账本更新在位
  F6  文档结论与引文在位
"""
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


DOC = io.open(os.path.join(HERE, "Z11_correction_k_is_fixed_and_age_is_not_site.md"),
              encoding="utf-8").read()
Z10 = io.open(os.path.join(HERE, "Z10_position_field_and_amplitude_criterion.md"),
              encoding="utf-8").read()


# ---------------------------------------------------------------- 工具
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


def Pk(c, k):
    return sum(m * comb(m - 1, m // 2) * c ** m for m in range(2, k + 1, 2))


def field_of_range(Nn, rng):
    eps = (rng - 1) / (rng + 1)
    i = np.arange(Nn)
    return 1 + eps * np.cos(2 * np.pi * (i + 0.5) / Nn)


def measure(c, k):
    Nn = len(c)
    W = W_sum(ring(c), k)
    w = np.array([W[j, (j + 1) % Nn] for j in range(Nn)])
    p = np.array([Pk(x, k) for x in c])
    core = [j for j in range(Nn) if 0.3 * Nn < j < 0.7 * Nn]
    return float(np.abs(w[core] / p[core] - 1).max()), float(w.max() / w.min())


def Pcut(i, L):
    if i <= 0 or i >= L:
        return 1.0
    if i % 2:
        return 0.0
    return comb(i, i // 2) * comb(L - i, (L - i) // 2) / comb(L, L // 2)


# ---------------------------------------------------------------- F1
head("F1  k 固定：大范围场非退化（字典偏差小、对比度多项式）")
NS = 512
fixed = {}
for k in (4, 8, 16):
    for rng in (3, 10, 50):
        dev, ctr = measure(field_of_range(NS, rng), k)
        fixed[(k, rng)] = (dev, ctr)
        print("      k=%2d 范围=%2d  dev=%.1e  ctr=%.1e" % (k, rng, dev, ctr))
check("k=4：范围 50 的场仍非退化（dev < 1e-3，ctr < 1e5）",
      fixed[(4, 50)][0] < 1e-3 and fixed[(4, 50)][1] < 1e5,
      "dev=%.1e ctr=%.1e" % fixed[(4, 50)])
check("k=8：范围 50 的场 dev < 1e-2", fixed[(8, 50)][0] < 1e-2, "%.1e" % fixed[(8, 50)][0])
check("所有固定 k 的偏差都 < 5e-2（字典有效）",
      all(v[0] < 5e-2 for v in fixed.values()),
      "max dev = %.1e" % max(v[0] for v in fixed.values()))
check("k=4 时对比度随 range 单调（多项式放大）",
      fixed[(4, 3)][1] < fixed[(4, 10)][1] < fixed[(4, 50)][1],
      "%s" % ["%.1e" % fixed[(4, r)][1] for r in (3, 10, 50)])

# ---------------------------------------------------------------- F2
head("F2  k = N（＝年龄＝站点）：退化（复现 Z10 §3）")
deg = []
for L in (64, 128, 256):
    Nn = L // 2
    c = np.array([Pcut(2 * j, L) for j in range(Nn)])
    c = c / c.mean()
    dev, ctr = measure(c, Nn)
    deg.append((L, Nn, dev, ctr))
    print("      L=%3d N=%3d  dev=%.2e  ctr=%.2e" % (L, Nn, dev, ctr))
check("k=N：偏差随 L 暴涨（>1）", deg[0][2] > 1 and deg[-1][2] > deg[0][2] * 10,
      " -> ".join("%.2e" % d[2] for d in deg))
check("k=N：对比度爆炸（>1e15）", all(d[3] > 1e15 for d in deg),
      " -> ".join("%.2e" % d[3] for d in deg))
check("两场景差别是读法而非场（固定 k 的 dev 全部 ≪ 1，k=N 的全部 ≫ 1）",
      max(v[0] for v in fixed.values()) < 0.05 and min(d[2] for d in deg) > 1.0)

# ---------------------------------------------------------------- F3
head("F3  对比度增长律：多项式于 range、指数于 k")
rows = []
for rng in (3, 10, 50):
    row = []
    for k in (2, 4, 8, 16, 32):
        row.append(measure(field_of_range(NS, rng), k)[1])
    rows.append((rng, row))
    print("      范围=%2d  " % rng + "  ".join("k=%2d:%.1e" % (k, v) for k, v in zip((2, 4, 8, 16, 32), row)))
check("固定范围：对比度随 k 单调增（指数于 k）",
      all(all(r[1][i] < r[1][i + 1] for i in range(4)) for r in rows))
check("固定 k：对比度随 range 单调增（多项式于 range）",
      all(rows[0][1][i] < rows[1][1][i] < rows[2][1][i] for i in range(5)))
sl = float(np.polyfit([2, 4, 8, 16, 32], np.log(rows[0][1]), 1)[0])
check("范围 3 时 ln(对比度) 对 k 线性，斜率 ≈ ln3=1.0986",
      abs(sl - math.log(3)) < 0.15, "斜率 %.4f（ln3=%.4f）" % (sl, math.log(3)))

# ---------------------------------------------------------------- F4
head("F4  『年龄＝站点 ⟹ k ∝ N』：环节核验")
G46 = io.open(os.path.join(HERE, "G46_k_is_the_lifetime.md"), encoding="utf-8").read()
G58 = io.open(os.path.join(HERE, "G58_I2a_resolved_as_embedding_input.md"), encoding="utf-8").read()
G61 = io.open(os.path.join(HERE, "G61_locking_the_five_integers.md"), encoding="utf-8").read()
Z0 = io.open(os.path.join(HERE, "Z0_zero_never_rests_single_axiom.md"), encoding="utf-8").read()
check("环节①：Z0 说年龄＝词长", "年龄 = 词长" in Z0 or "年龄＝词长" in Z0)
check("环节②：G46 说 k 被逼成 L", "k=L" in G46 or "$k=L$" in G46)
check("环节③：G58 说公理内读法是『L 固定步数』", "固定步数" in G58 and "量纲常数" in G58)
check("环节④：G61 锁 L=4（【条件】）", "$L=4$" in G61 or "L=4" in G61)
check("⇒ 若年龄＝站点，则细化空间 ⟹ L∝N ⟹ k∝N ⟹ §§F2 的退化",
      "年龄 ≠ 站点" in DOC and "k=N" in DOC)

# ---------------------------------------------------------------- F5
head("F5  更正后的 C7 与账本更新在位")
check("C7 新表述：对比度 = range^k 有界", "range}(c)^{\\,k}" in DOC or "range}{}^{k}" in DOC or "\\mathrm{range}(c)^{\\,k}" in DOC)
check("写明 Z10 §3 的数值正确、错的是场景", "读法更正" in DOC or "场景" in DOC)
check("写明 Z8 候选 2 存活", "存活" in DOC and "候选 2" in DOC)
check("写明 Z9 的 τ_i 分叉结掉（(a) 存活／(b) 出局）", "分叉" in DOC and "(a) 存活" in DOC)
check("写明对 I5 的新约束（站点与寿命无关）", "与寿命无关" in DOC)
check("Z10 已加更正横幅", "Z11" in Z10 and "更正" in Z10)

# ---------------------------------------------------------------- F6
head("F6  文档结论与引文在位")
check("结论盒：更正 + 新约束", "不是公理内的读法" in DOC and "该识别被排除" in DOC)
check("三种读法表在位", "固定步数" in DOC and "propto N" in DOC and "量纲常数" in DOC)
check("诚实边界在位（读法更正／$L=4$ 的等级／$k$ 未知）", "读法更正" in DOC and "余量" in DOC)
for fn, keys in [
    ("G58_I2a_resolved_as_embedding_input.md", ["固定步数", "量纲常数"]),
    ("G61_locking_the_five_integers.md", ["L=4"]),
    ("G46_k_is_the_lifetime.md", ["顶点传递"]),
    ("D_arc/D214_local_zero_sum_transport.md", ["额外识别"]),
    ("Z9_pi_filter_and_lifetime_fork.md", ["分叉"]),
    ("Z8_native_scale_field_candidate.md", ["推前重数"]),
    ("Z7_embedding_input_explicit_dictionary.md", ["闭式字典"]),
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
