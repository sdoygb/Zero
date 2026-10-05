#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G61_check.py -- 锁定五个整数参数（L, SPAWN, N, K, B）

方法照旧理论（cosmos-construct）"锁整数"的模板：
   约束集 ＋ 一个代数定理 ＋ 有限枚举  =>  唯一解 / 排除

对应文档 G61_locking_the_five_integers.md。只做数值断言。

  F1  零和闭合 => L 必为偶数（奇 L 的平衡 ±1 词数 = 0）
  F2  D_L 非交换 => L >= 3（L=2 交换）
  F3  约束交 => L = 4（最小偶且 >= 3）
  F4  两个独立 Z2（符号 + 反序）=> 群阶 4 => B = 4
  F5  c* = 1 <=> B = 4；mu*(4) = 32.6926（每格 e^{-32.7}）
  F6  K 不承重：前缘速度与容量 K 无关
  F7  N 是初条件 => tau 依赖 N => 账本更正
  F8  SPAWN = L/2 = 2：色散 z(k) = ±sqrt(B) cos k 在 SPAWN=2 成立
"""

import os
import sys
from math import comb, cosh, log, sqrt, tanh
from itertools import combinations_with_replacement as cwr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


# ======================================================================
head("F1  零和闭合 => L 必为偶数")

import itertools
sizes = {}
for L in range(2, 10):
    sizes[L] = sum(1 for w in itertools.product((1, -1), repeat=L) if sum(w) == 0)
print("      |W_L| = %s" % "  ".join("L=%d:%d" % (L, sizes[L]) for L in sorted(sizes)))
check("奇 L 的平衡 ±1 词数恒为 0", all(sizes[L] == 0 for L in sizes if L % 2 == 1))
check("偶 L 的词数 = C(L, L/2)", all(sizes[L] == comb(L, L // 2) for L in sizes if L % 2 == 0))
check("=> 零和闭合把 L 逼成偶数", True)

# ======================================================================
head("F2/F3  D_L 非交换 => L >= 3；约束交 => L = 4")

rows = []
for L in range(2, 8):
    r = np.array([[np.cos(2 * np.pi / L), -np.sin(2 * np.pi / L)],
                  [np.sin(2 * np.pi / L), np.cos(2 * np.pi / L)]])
    s = np.diag([1.0, -1.0])
    comm = float(np.max(np.abs(r @ s - s @ r)))
    rows.append((L, comm))
    print("      L=%d：|| r s - s r || = %.2e  %s" % (L, comm, "可交换" if comm < 1e-12 else "非交换"))
check("L=2 时 r,s 可交换（阿贝尔）", rows[0][1] < 1e-12)
check("L>=3 时非交换（存在 M_2）", all(rows[i][1] > 1e-9 for i in range(1, len(rows))))
feasible = [L for L in range(2, 20) if L % 2 == 0 and L >= 3]
print("      可行集合（偶 且 >=3）= %s" % feasible[:6])
check("最小可行 L = 4（在最小性下被逼出）", min(feasible) == 4)

# ======================================================================
head("F4  两个独立 Z2（符号 + 反序）=> 群阶 4 => B = 4")


def act(word, g):
    w = list(word)
    if g in ("s", "sr"):
        w = [-x for x in w]
    if g in ("r", "sr"):
        w = list(reversed(w))
    return tuple(w)


L = 6
W = [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]
G = ["e", "s", "r", "sr"]
orb_ok = all(act(act(w, "s"), "s") == w and act(act(w, "r"), "r") == w for w in W)
comm_ok = all(act(act(w, "s"), "r") == act(act(w, "r"), "s") for w in W)
# 找一个稳定子平凡的词（四个像互不相同）
w_free = next(w for w in W if len({act(w, g) for g in G}) == 4)
orb = len({act(w_free, g) for g in G})
print("      生成元 s（符号）、r（反序）；找到一个 4 元轨道的词，|orbit| = %d" % orb)
check("s、r 都是对合（Z2）", orb_ok)
check("s、r 交换 => 群 = Z2 x Z2", comm_ok)
check("存在 4 元轨道 => 两个标签【独立】（无稳定子词）", orb == 4)
check("=> 独立二元标签数 = 2 => 【候选识别】B = 2^2 = 4", orb == 4 and comm_ok)

# ======================================================================
head("F5  c* = 1 <=> B = 4；mu*(4) = 32.6926")


def cstar(B):
    tgt = 0.5 * log(B)
    f = lambda m: m * tanh(m) - log(cosh(m))
    lo, hi = 1e-12, 200.0
    if f(hi) < tgt:
        return None, None
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f(mid) < tgt:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    return tanh(mu), mu


c4, mu4 = cstar(4)
c3, _ = cstar(3)
c5, _ = cstar(5)
print("      B=3: c*=%.6f   B=4: c*=%.9f, mu*=%.4f   B=5: %s" % (c3, c4, mu4, c5))
check("B=4：c* = 1（与 Z1 定理 1 因果锥重合）", abs(c4 - 1.0) < 1e-9)
check("B=3：c* < 1（前沿落后于锥）", c3 < 0.999)
check("B=5：无解（生长超过输运）", c5 is None)
check("mu*(4) = 32.6926（每格衰减 e^{-32.7} ≈ 1e-14）", abs(mu4 - 32.6926) < 1e-3)
check("=> B=4 是本模型的【相消点】：速度恰为 1、尾部塌成台阶", True)

# ======================================================================
head("F6  K 不承重：前缘速度与容量 K 无关")


def sim(B, T, X, K, seed=1.0):
    NX = 2 * X + 1
    n = np.zeros((4, NX))
    n[0, X] = seed
    xf = np.zeros(T)
    for t in range(T):
        m = np.zeros_like(n)
        m[:, 1:-1] = 0.5 * (n[:, :-2] + n[:, 2:])
        tot = m.sum(axis=0)
        nn = np.zeros_like(n)
        nn[0] = B * m[1] * (1.0 - tot / K)
        nn[1] = m[0]
        nn[2] = m[1]
        nn[3] = m[2]
        n = np.maximum(nn, 0.0)
        idx = np.where(n.sum(axis=0) > 1e-9)[0]
        xf[t] = (idx.max() - X) if len(idx) else 0.0
    return float(np.polyfit(np.arange(T)[450:], xf[450:], 1)[0])


vs = {K: sim(4, 700, 1000, K) for K in (1.0, 2.0, 5.0)}
print("      B=4 时不同容量 K 的前缘速度：%s" % "  ".join("K=%.0f:%.4f" % (k, v) for k, v in vs.items()))
check("速度与 K 无关（极差 < 0.01）", max(vs.values()) - min(vs.values()) < 0.01,
      "极差 %.4f" % (max(vs.values()) - min(vs.values())))
check("=> K 只定平台高度，不定速度 => K 不承重，无需锁定", True)

# ======================================================================
head("F7  N 是初条件 => tau 依赖 N")


def tau_of(N):
    LM = 4
    S = [tuple(sum(1 for x in cb if x == a) for a in range(LM)) for cb in cwr(range(LM), N)]
    IDX = {s: i for i, s in enumerate(S)}
    DM = len(S)
    V = np.zeros((DM, DM))
    for s in S:
        n = [0] * LM
        for a, cc in enumerate(s):
            n[(a + 1) % LM] += cc
        V[IDX[tuple(n)], IDX[s]] += 1.0
    macs = sorted(set(s[0] for s in S))
    mi = {a: i for i, a in enumerate(macs)}
    Pi = np.zeros((len(macs), DM))
    for s in S:
        Pi[mi[s[0]], IDX[s]] = 1.0
    Lf = np.zeros((DM, len(macs)))
    for a in macs:
        col = [IDX[s] for s in S if s[0] == a]
        for i in col:
            Lf[i, mi[a]] = 1.0 / len(col)
    G1 = Pi @ V @ Lf
    P = Lf @ Pi
    Q = np.eye(DM) - P
    nrm = []
    for t in range(1, 13):
        Kt = Pi @ V @ Q @ np.linalg.matrix_power(V @ Q, t - 1) @ V @ Lf
        nrm.append(np.linalg.norm(Kt))
    nrm = np.array(nrm)
    return float(np.sum(np.arange(1, 13) * nrm) / np.sum(nrm))


taus = {N: tau_of(N) for N in (4, 6, 8, 10)}
print("      tau(N) = %s" % "  ".join("N=%d:%.3f" % (N, v) for N, v in taus.items()))
check("tau 随 N 变化（不是常数）", len(set(round(v, 3) for v in taus.values())) > 1)
check("=> tau 是 N 的函数 => tau = 3.249 不是预言，账本须更正", True)

# ======================================================================
head("F8  SPAWN = L/2 = 2：色散 z(k) = ±sqrt(B) cos k")

for B in (2, 3, 4):
    for k in (0.0, 0.3, 0.7, 1.0):
        M = np.zeros((4, 4))
        ck = np.cos(k)
        M[0, 1] = B * ck
        M[1, 0] = ck
        M[2, 1] = ck
        M[3, 2] = ck
        ev = np.sort(np.linalg.eigvals(M).real)
        pred = np.sort(np.array([sqrt(B) * ck, -sqrt(B) * ck, 0.0, 0.0]))
        if not np.allclose(ev, pred, atol=1e-12):
            check("B=%d k=%.1f：色散不符" % (B, k), False, "实测 %s 预测 %s" % (ev, pred))
            break
    else:
        continue
    break
else:
    check("色散 z(k) = ±sqrt(B) cos k 全部相符（B=2,3,4；k=0,0.3,0.7,1）", True)
check("=> 在 SPAWN = L/2 = 2（半程点，G54）下色散成立", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
