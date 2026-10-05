#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z10_check.py —— 【位置型原生场：D214 的环与归零概率；新判据 C7】核验
====================================================================
独立实断言：
  F1  归零概率闭式 = 暴力枚举（L=4..14）
  F2  宇称：奇数位置精确为 0（宇称倍增被迫）
  F3  批量缓变：体内相对步长 ~ O(1/L)（拟合指数 ≈ -1）
  F4  端点幂律边界层 ~ 1/(2(j+1))
  F5  场的范围随 L 按 √L 增长
  F6  固定 k：字典二阶收敛（每加倍 ÷4）
  F7  k = N（= L 步）：退化（偏差暴涨、对比度爆炸）
  F8  r(w) 分布对账（类层面 = 语料；词层面均值 = Σ_i P_i）
  F9  文档结论与引文在位
"""
import itertools
import io
import math
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


DOC = io.open(os.path.join(HERE, "Z10_position_field_and_amplitude_criterion.md"),
              encoding="utf-8").read()


# ---------------------------------------------------------------- 工具
def Pcut(i, L):
    """归零（切割）概率：Pr(S_i = 0 | 长度 L 的平衡词)。"""
    if i <= 0 or i >= L:
        return 1.0
    if i % 2:
        return 0.0
    return comb(i, i // 2) * comb(L - i, (L - i) // 2) / comb(L, L // 2)


def balanced(L):
    return [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]


def rotclass(w):
    L = len(w)
    return min(tuple(w[i:] + w[:i]) for i in range(L))


def rcount(w):
    s = 0
    c = 0
    for i, x in enumerate(w, 1):
        s += x
        if s == 0 and i < len(w):
            c += 1
    return c


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


def field(L):
    """宇称倍增环上的归一化归零概率场。"""
    Nn = L // 2
    c = np.array([Pcut(2 * j, L) for j in range(Nn)])
    return c / c.mean(), Nn


# ---------------------------------------------------------------- F1
head("F1  归零概率闭式 = 暴力枚举")
for L in (4, 6, 8, 10, 12, 14):
    W = balanced(L)
    cnt = Counter()
    for w in W:
        s = 0
        for i, x in enumerate(w, 1):
            s += x
            if s == 0:
                cnt[i] += 1
    dev = max(abs(cnt.get(i, 0) / len(W) - Pcut(i, L)) for i in range(1, L + 1))
    check("L=%2d：max|暴力-闭式| = 0" % L, dev == 0.0, "%.1e" % dev)

# ---------------------------------------------------------------- F2
head("F2  宇称：奇数位置精确为 0 ⟹ 宇称倍增被迫")
for L in (8, 64, 256):
    odd = max(Pcut(i, L) for i in range(1, L, 2))
    check("L=%3d：奇数位置 P 恒为 0" % L, odd == 0.0, "%.1e" % odd)
check("偶数位置非零（有效图 = 宇称倍增环）", all(Pcut(2 * j, 64) > 0 for j in range(1, 31)))

# ---------------------------------------------------------------- F3
head("F3  批量缓变：体内相对步长 ~ O(1/L)")
res = []
for L in (64, 128, 256, 512, 1024):
    P = [Pcut(2 * j, L) for j in range(L // 2 + 1)]
    bulk = [j for j in range(L // 2) if 0.2 * L / 2 < j < 0.8 * L / 2]
    step = max(abs(P[j + 1] - P[j]) / P[j] for j in bulk)
    res.append((L, step))
    print("      L=%4d 体内最大相对步长 = %.5f" % (L, step))
slope = float(np.polyfit(np.log([r[0] for r in res]), np.log([r[1] for r in res]), 1)[0])
check("拟合指数 ≈ -1（O(1/L) 缓变）", abs(slope + 1) < 0.15, "斜率 %.3f" % slope)
check("步长单调递减", all(res[i][1] > res[i + 1][1] for i in range(len(res) - 1)))

# ---------------------------------------------------------------- F4
head("F4  端点幂律边界层 ~ 1/(2(j+1))")
for L in (256, 1024):
    P = [Pcut(2 * j, L) for j in range(L // 2 + 1)]
    rel = [abs(P[j + 1] - P[j]) / P[j] for j in range(6)]
    pred = [1.0 / (2 * (j + 1)) for j in range(6)]
    check("L=%4d：前 6 个相对步长 ≈ 1/(2(j+1))（偏差 < 0.02）" % L,
          max(abs(a - b) for a, b in zip(rel, pred)) < 0.02,
          "%s" % ["%.3f" % x for x in rel])

# ---------------------------------------------------------------- F5
head("F5  场的范围随 L 按 √L 增长")
rngs = []
for L in (64, 128, 256, 512):
    c, Nn = field(L)
    rngs.append((L, float(c.max() / c.min())))
    print("      L=%4d 范围 = %.2f" % (L, rngs[-1][1]))
check("范围单调增长", all(rngs[i][1] < rngs[i + 1][1] for i in range(len(rngs) - 1)),
      "%s" % ["%.2f" % r[1] for r in rngs])
check("增长慢于 L（≈√L）", rngs[-1][1] / rngs[0][1] < math.sqrt(rngs[-1][0] / rngs[0][0]) * 1.5,
      "×%.2f（√(L比)=%.2f）" % (rngs[-1][1] / rngs[0][1], math.sqrt(rngs[-1][0] / rngs[0][0])))

# ---------------------------------------------------------------- F6
head("F6  固定 k：字典二阶收敛")
fixed = []
for L in (64, 128, 256, 512):
    c, Nn = field(L)
    k = 4
    W = W_sum(ring(c), k)
    w = np.array([W[j, (j + 1) % Nn] for j in range(Nn)])
    p = np.array([Pk(x, k) for x in c])
    core = [j for j in range(Nn) if 0.4 * Nn < j < 0.6 * Nn]
    fixed.append(float(np.abs(w[core] / p[core] - 1).max()))
print("      k=4 体内字典偏差 = %s" % ["%.2e" % x for x in fixed])
check("固定 k=4：字典偏差随 L 二阶下降（每加倍 ÷>3）",
      all(fixed[i] / fixed[i + 1] > 3 for i in range(len(fixed) - 1)),
      " -> ".join("%.2e" % x for x in fixed))

# ---------------------------------------------------------------- F7
head("F7  k = N（= L 步）：退化")
deg = []
for L in (64, 128, 256):
    c, Nn = field(L)
    k = Nn
    W = W_sum(ring(c), k)
    w = np.array([W[j, (j + 1) % Nn] for j in range(Nn)])
    p = np.array([Pk(x, k) for x in c])
    core = [j for j in range(Nn) if 0.4 * Nn < j < 0.6 * Nn]
    dev = float(np.abs(w[core] / p[core] - 1).max())
    ctr = float(w.max() / w.min())
    deg.append((L, dev, ctr))
    print("      L=%4d k=%3d 偏差=%.2e 对比度=%.2e" % (L, k, dev, ctr))
check("k=N：字典偏差随 L 暴涨（> 1）", deg[0][1] > 1 and deg[-1][1] > deg[0][1] * 10,
      " -> ".join("%.2e" % d[1] for d in deg))
check("k=N：度规对比度爆炸（> 1e15）", all(d[2] > 1e15 for d in deg),
      " -> ".join("%.2e" % d[2] for d in deg))
check("⇒ 新判据 C7（有界振幅）", "C7" in DOC and "有界振幅" in DOC)

# ---------------------------------------------------------------- F8
head("F8  r(w) 分布对账")
EXPECT = {8: {0: 5, 1: 3, 2: 1, 3: 1}, 12: {0: 37, 1: 21, 2: 13, 3: 7, 4: 1, 5: 1}}
for L in (8, 12):
    W = balanced(L)
    cls = {}
    wd = Counter()
    for w in W:
        r = rcount(w)
        cls[rotclass(w)] = r
        wd[r] += 1
    cd = Counter(cls.values())
    check("L=%2d：类层面 r 分布 = 语料" % L, cd == Counter(EXPECT[L]), "%s" % dict(sorted(cd.items())))
    mean_w = sum(k * v for k, v in wd.items()) / sum(wd.values())
    span = sum(Pcut(i, L) for i in range(1, L))
    check("L=%2d：词层面均值 E[r] = Σ_i P_i" % L, abs(mean_w - span) < 1e-12,
          "%.6f vs %.6f" % (mean_w, span))

# ---------------------------------------------------------------- F9
head("F9  文档结论与引文在位")
check("结论盒：位置型场存在但死于 k=L", "位置型原生场存在" in DOC and "退化度规" in DOC)
check("写明站点图 = D214 的词位环（识别部分 = I5）", "环图" in DOC and "I5" in DOC)
check("写明逃出引理 79 的理由（位置条件化）", "位置条件化" in DOC)
check("写明宇称倍增是被迫的", "宇称倍增" in DOC)
check("写明 C7 与判据清单更新到七条", "C7" in DOC and "七条" in DOC)
check("诚实边界在位（替身／识别）", "替身" in DOC and "识别" in DOC)
for fn, keys in [
    ("D214_local_zero_sum_transport.md", ["环图", "额外识别"]),
    ("G54_quantitative_profile_age_measure.md", ["引理 79"]),
    ("G59_I7_settled_native_cone_and_its_residue.md", ["宇称"]),
    ("zero_sum_rotation_class_algebra.md", ["内部归零数"]),
    ("G46_k_is_the_lifetime.md", ["k=L", "顶点传递"]),
    ("Z9_pi_filter_and_lifetime_fork.md", ["年龄筛"]),
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
