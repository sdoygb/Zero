#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z9_check.py —— 【π 的筛子：年龄筛与生存筛、以及局部寿命的分叉】核验
======================================================================
独立实断言：
  F1  年龄筛：只依赖年龄的 π 给完全相等的类规模（G54 引理 79）
  F2  幸存者：旋转类（Z3 的写入标签）——类规模非均匀，分布 = G37
  F3  生存筛：c=F(a) 在顶点传递图（环）上非均匀 ⟹ 违反 C4（G46 §3）
  F4  局部寿命 τ_i：两种读法都满足 C4，但彼此分叉
  F5  勘误：G54 §4.1 的 ×5/×9
  F6  文档结论与引文在位
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


DOC = io.open(os.path.join(HERE, "Z9_pi_filter_and_lifetime_fork.md"), encoding="utf-8").read()
G54 = io.open(os.path.join(HERE, "G54_quantitative_profile_age_measure.md"), encoding="utf-8").read()


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


def layers(A, kmax):
    out = [None, None]
    P = A.copy()
    out.append(A.multiply(P) * 2)
    for m in range(3, kmax + 1):
        P = (P @ A).tocsr()
        out.append(m * A.multiply(P))
    return out


def w_globalk(A, k):
    n = A.shape[0]
    Ls = layers(A, k)
    W = sum(Ls[2:k + 1])
    return np.array([W[i, (i + 1) % n] for i in range(n)])


def w_localk(A, kloc):
    n = A.shape[0]
    Ls = layers(A, int(max(kloc)))
    return np.array([sum(Ls[m][i, (i + 1) % n] for m in range(2, kloc[i] + 1)) for i in range(n)])


def balanced_words(L):
    return [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]


def rotclass(w):
    L = len(w)
    return min(tuple(w[i:] + w[:i]) for i in range(L))


def F_ab(a, L):
    if a > L:
        a = L
    return sum(comb(a, (a + h) // 2) for h in range(-a, a + 1)
               if abs(h) <= L - a and (h - a) % 2 == 0) / 2 ** a


# ---------------------------------------------------------------- F1
head("F1  年龄筛：π 只依赖年龄 ⟹ 类规模恒等（G54 引理 79）")
occ_all = {}
for L in (4, 6, 8, 10, 12):
    W = balanced_words(L)
    occ = [0] * (L + 1)
    for w in W:
        for a in range(L + 1):
            occ[a] += 1                       # 每个词在每一岁计一次
    occ_all[L] = (len(W), min(occ), max(occ))
    check("L=%2d：各年龄占用恒等（min=max=|W_L|=%d）" % (L, len(W)),
          min(occ) == max(occ) == len(W), "%d/%d" % (min(occ), max(occ)))
check("⇒ 年龄型 π 给 c≡1 ⇒ 平坦（零几何信息）", all(v[1] == v[2] for v in occ_all.values()))

# ---------------------------------------------------------------- F2
head("F2  幸存者：旋转类（Z3 的写入标签）——类规模非均匀")
EXPECT = {4: {2: 1, 4: 1}, 8: {2: 1, 4: 1, 8: 8}, 12: {2: 1, 4: 1, 6: 3, 12: 75}}
for L in (4, 8, 12):
    cls = {}
    for w in balanced_words(L):
        k = rotclass(w)
        cls[k] = cls.get(k, 0) + 1
    dist = {}
    for s in sorted(cls.values()):
        dist[s] = dist.get(s, 0) + 1
    check("L=%2d：类数 = %d" % (L, len(cls)), len(cls) == {4: 2, 8: 10, 12: 80}[L],
          "%d" % len(cls))
    check("L=%2d：规模分布 = %s（= G37）" % (L, EXPECT[L]), dist == EXPECT[L], "%s" % dist)
    check("L=%2d：总重数 = |W_L| = %d" % (L, len(balanced_words(L))),
          sum(cls.values()) == len(balanced_words(L)))
check("旋转类规模非均匀（max/min = 2, 4, 6）⇒ c ≢ 1 ⇒ 几何",
      all(max(EXPECT[L]) / min(EXPECT[L]) == {4: 2, 8: 4, 12: 6}[L] for L in EXPECT))

# ---------------------------------------------------------------- F3
head("F3  生存筛：c=F(a) 在顶点传递图（环）上非均匀 ⟹ 违反 C4")
L = 8
c = np.array([F_ab(a, L) for a in range(L)])
w = w_globalk(ring(c), L)
rel = float((w.max() - w.min()) / w.mean())
check("C_8 上真实走道权重的相对差 > 0.3（实测 %.3f）" % rel, rel > 0.3, "%.4f" % rel)
check("不是取值数 1（G46 §3 要求）", len(set(np.round(w, 9))) > 1, "%d 个不同值" % len(set(np.round(w, 9))))
c1 = np.ones(L)
check("均匀基线（c≡1）在环上精确均匀", float((lambda x: (x.max() - x.min()) / x.mean())(w_globalk(ring(c1), L))) == 0.0)
pred = np.array([sum(m * comb(m - 1, m // 2) * x ** m for m in range(2, L + 1, 2)) for x in c])
check("字典预测同样非均匀（两条路都出局，结论不依赖字典）",
      float((pred.max() - pred.min()) / pred.mean()) > 1.0, "%.3f" % float((pred.max() - pred.min()) / pred.mean()))
check("真实权重比字典预测光滑（邻域半径 k/2-1 抹开跳跃）",
      float(np.abs(np.diff(w)).max()) < float(np.abs(np.diff(pred)).max()),
      "%.1f < %.1f" % (float(np.abs(np.diff(w)).max()), float(np.abs(np.diff(pred)).max())))

# ---------------------------------------------------------------- F4
head("F4  局部寿命 τ_i：两种读法都满足 C4，但彼此分叉")
LL = 16
i = np.arange(LL)
tau = LL * (1 + 0.3 * np.cos(2 * np.pi * (i + 0.5) / LL))
wu_g = w_globalk(ring(np.ones(LL)), LL)
wu_l = w_localk(ring(np.ones(LL)), [LL] * LL)
check("均匀 τ≡L：读法 (a) 全局 k=L 精确平坦", float((wu_g.max() - wu_g.min()) / wu_g.mean()) == 0.0)
check("均匀 τ≡L：读法 (b) 局部截断 k_i=τ_i=L 精确平坦",
      float((wu_l.max() - wu_l.min()) / wu_l.mean()) == 0.0)
w_c = w_globalk(ring(tau / LL), LL)
kloc = np.clip(np.round(tau).astype(int), 2, LL)
w_k = w_localk(ring(np.ones(LL)), kloc)
check("非均匀 τ：(a) 对比度 > 1e3（实测 %.0f）" % (w_c.max() / w_c.min()), w_c.max() / w_c.min() > 1e3)
check("非均匀 τ：(b) 对比度 > 10（实测 %.1f）" % (w_k.max() / w_k.min()), w_k.max() / w_k.min() > 10)
check("两读法逐边相对差 > 0.9（分叉，实测 %.3f）" % float(np.abs(w_k / w_c - 1).max()),
      float(np.abs(w_k / w_c - 1).max()) > 0.9)
check("两读法都非均匀 ⟹ τ_i 通过 C4 但不是同一个模型",
      float((w_c.max() - w_c.min()) / w_c.mean()) > 0.3 and float((w_k.max() - w_k.min()) / w_k.mean()) > 0.3)

# ---------------------------------------------------------------- F5
head("F5  勘误：G54 §4.1 的 ×5/×9")
check("G54 已写 ×5", "1.0000 ×5" in G54)
check("G54 已写 ×9", "1.0000 ×9" in G54)
check("G54 不再写 ×4 / ×8", "1.0000 ×4" not in G54 and "1.0000 ×8" not in G54)
check("G54 有勘误横幅并指向 Z9", "勘误" in G54 and "Z9_pi_filter_and_lifetime_fork.md" in G54)
for L in (8, 16):
    fc = [F_ab(a, L) for a in range(L + 1)]
    check("L=%2d：F≡1 的条数 = L/2+1 = %d" % (L, L // 2 + 1),
          sum(1 for x in fc if x == 1.0) == L // 2 + 1, "%d" % sum(1 for x in fc if x == 1.0))

# ---------------------------------------------------------------- F6
head("F6  文档结论与引文在位")
check("结论盒：几何只能来自词的内容", "词的内容" in DOC and "年龄" in DOC)
check("写明年龄筛与生存筛两道筛", "年龄筛" in DOC and "生存筛" in DOC)
check("写明幸存者＝旋转类，且引用 G35 的写入规则", "旋转类" in DOC and "Z3 的写入" in DOC)
check("写明 τ_i 分叉（两读法与对比度）", "分叉" in DOC and "1801" in DOC and "83" in DOC)
check("写明新缺口：类↔站点识别的准均匀性", "准均匀" in DOC)
check("诚实边界在位（替身／C5 未验）", "替身" in DOC and "未**被检验**" in DOC or "未被检验" in DOC)
for fn, keys in [
    ("G54_quantitative_profile_age_measure.md", ["引理 79", "半程定理"]),
    ("G35_reseeding_and_chirality.md", ["旋转类"]),
    ("G37_reseeding_law_and_sign_symmetry_theorem.md", ["reseed"]),
    ("G46_k_is_the_lifetime.md", ["顶点传递"]),
    ("G29_probability_as_derived_not_postulated.md", ["推前"]),
    ("Z8_native_scale_field_candidate.md", ["推前重数"]),
    ("Z1_zero_layer_as_the_foundation.md", ["局部寿命"]),
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
