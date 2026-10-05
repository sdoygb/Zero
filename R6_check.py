#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R6_check.py -- H5 字典误差界（定理 R6.3）的独立核验。

对应文档 R6_h5_dictionary_error_bound.md。只读文件，不写回任何既有文件。

  F1  精确游走数 W_d 与 Z7 表一致
  F2  反射—反转配对：一阶项精确为零（暴力枚举游程）
  F3  周期环面 1D/2D/3D：O(a^2) 与常数稳定性
  F4  尖锐速率：C^inf -> 2、C^{0,1} -> 1、C^{0,alpha} -> alpha、间断 -> 0
  F5  常数随 k 的经验增长落在 [2.5, 3.3]
  F6  文档诚实边界：把 H5 写成定理并写死三条前提

数值是旁证，不替代文档中的引理 R6.1–R6.2 与定理 R6.3 的证明。
"""

import io
import os
import sys
from itertools import product

import numpy as np
from scipy.sparse import csr_matrix

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


DOC = io.open(
    os.path.join(HERE, "R6_h5_dictionary_error_bound.md"), encoding="utf-8"
).read()


def anchor(*tokens):
    return all(tok in DOC for tok in tokens)


# ======================================================================
# 精确游走数 W_d(L)：Z^d 上从 0 到 e_1 的 L 步游走数
def walk_count(d, L):
    if L < 0:
        return 0
    cur = {(0,) * d: 1}
    for _ in range(L):
        nxt = {}
        for pos, val in cur.items():
            for ax in range(d):
                for s in (-1, 1):
                    q = list(pos)
                    q[ax] += s
                    q = tuple(q)
                    nxt[q] = nxt.get(q, 0) + val
        cur = nxt
    return cur.get((1,) + (0,) * (d - 1), 0)


def prefactor(d, k):
    return {m: m * walk_count(d, m - 1) for m in range(2, k + 1)}


# ----------------------------------------------------------------------
# 周期环面上的加权邻接矩阵与闭环权
def build_periodic(dim, N, c, k):
    n = N ** dim

    def idx(p):
        return sum((p[t] % N) * (N ** t) for t in range(dim))

    rows, cols, vals = [], [], []
    for p in product(range(N), repeat=dim):
        a = idx(p)
        for ax in range(dim):
            q = list(p)
            q[ax] = (p[ax] + 1) % N
            q = tuple(q)
            b = idx(q)
            if a == b:
                continue
            mid = [p[t] / N for t in range(dim)]
            mid[ax] = ((p[ax] + 0.5) % N) / N
            rows.append(a)
            cols.append(b)
            vals.append(c(mid[0]) if dim == 1 else c(tuple(mid)))
    A = csr_matrix((vals, (rows, cols)), shape=(n, n))
    A = A + A.T  # 每条边只登记一次，转置补回对称项，不再除 2
    P = A.copy()
    W = A.multiply(P) * 2
    for m in range(3, k + 1):
        P = (P @ A).tocsr()
        W = W + m * A.multiply(P)
    return W, idx


def rel_error_1d(N, c, k):
    """周期 1D：逐边返回 |W_e/P(c_e)-1|。"""
    pf = prefactor(1, k)
    W, idx = build_periodic(1, N, c, k)
    xe = ((np.arange(N) + 0.5) % N) / N
    ce = c(xe)
    P = sum(pf[m] * ce ** m for m in range(2, k + 1))
    diag1 = np.array([W[i, (i + 1) % N] for i in range(N)])
    return np.abs(diag1 / P - 1)


def rel_error_2d(N, c, k):
    pf = prefactor(2, k)
    W, idx = build_periodic(2, N, c, k)
    errs = []
    for i in range(N - 1):
        for j in range(N):
            a = idx((i, j))
            ce = c(((i + 0.5) / N, j / N))
            P = sum(pf[m] * ce ** m for m in range(2, k + 1))
            errs.append(abs(W[a, idx((i + 1, j))] / P - 1))
    return np.array(errs)


# ======================================================================
head("F1  精确游走数 W_d 与 Z7 表一致")

W1 = [walk_count(1, L) for L in (1, 3, 5, 7)]
W2 = [walk_count(2, L) for L in (1, 3, 5, 7)]
W3 = [walk_count(3, L) for L in (1, 3, 5, 7)]
print("      W1=%s  W2=%s  W3=%s" % (W1, W2, W3))
check("W1=[1,3,10,35]", W1 == [1, 3, 10, 35])
check("W2=[1,9,100,1225]", W2 == [1, 9, 100, 1225])
check("W3=[1,15,310,7455]", W3 == [1, 15, 310, 7455])
check("偶数步游走数为 0（二部性）", all(walk_count(d, 2) == 0 for d in (1, 2, 3)))


# ======================================================================
head("F2  反射—反转配对：一阶项精确为零（暴力枚举）")


def walks_from_to(d, steps, target):
    """枚举 Z^d 中从 0 到 target 的 steps 步游走（位移序列）。"""
    out = []
    for seq in product(range(d), repeat=steps):
        for signs in product((-1, 1), repeat=steps):
            pos = [0] * d
            ok = True
            for ax, s in zip(seq, signs):
                pos[ax] += s
            if tuple(pos) == target:
                out.append(list(zip(seq, signs)))
    return out


def midpoint_sum_zero(d, steps):
    """Σ_walks Σ_j (mid_j - x_e)，x_e = e_1/2；返回精确有理数 = 0。"""
    target = (1,) + (0,) * (d - 1)
    xe = [0.5] + [0.0] * (d - 1)
    total = np.zeros(d)
    for w in walks_from_to(d, steps, target):
        pos = [0.0] * d
        for ax, s in w:
            edge_mid = list(pos)
            edge_mid[ax] += 0.5 * s
            for k2 in range(d):
                total[k2] += edge_mid[k2] - xe[k2]
            pos[ax] += s
    return total


ok_zero = True
for d in (1, 2, 3):
    for steps in (1, 3, 5):
        t = midpoint_sum_zero(d, steps)
        if np.max(np.abs(t)) > 1e-12:
            ok_zero = False
        print("      d=%d steps=%d  Σ(mid-x_e)=%s  |w|=%d" %
              (d, steps, np.round(t, 12).tolist(), walk_count(d, steps)))
check("一阶项对所有 (d,steps) 精确为零（引理 R6.2）", ok_zero)


# ======================================================================
head("F3  周期环面 1D/2D/3D：O(a^2) 与常数稳定性")

c1 = lambda x: 1 + 0.3 * np.cos(2 * np.pi * x)
k = 8
Ns = (64, 128, 256, 512, 1024, 2048)
rels = [float(rel_error_1d(N, c1, k).max()) for N in Ns]
C1 = [r * N ** 2 for r, N in zip(rels, Ns)]
print("      1D max rel = %s" % ["%.3e" % r for r in rels])
print("      1D rel*N^2 = %s" % ["%.3f" % c for c in C1])
check("1D 误差按 O(a^2) 缩小（相邻比值约 4）",
      all(3.6 < rels[i] / rels[i + 1] < 4.4 for i in range(len(rels) - 1)))
check("1D 常数 rel*N^2 收敛（首尾相对差 < 1%）",
      abs(C1[0] - C1[-1]) / C1[-1] < 0.01)

c2 = lambda p: 1 + 0.2 * np.cos(2 * np.pi * p[0]) + 0.15 * np.cos(2 * np.pi * p[1])
Ns2 = (16, 32, 64, 128)
rels2 = [float(rel_error_2d(N, c2, 6).max()) for N in Ns2]
C2 = [r * N ** 2 for r, N in zip(rels2, Ns2)]
print("      2D max rel = %s" % ["%.3e" % r for r in rels2])
print("      2D rel*N^2 = %s" % ["%.3f" % c for c in C2])
check("2D 误差按 O(a^2) 缩小", all(rels2[i] / rels2[i + 1] > 3.6 for i in range(len(rels2) - 1)))
check("2D 常数 rel*N^2 收敛（首尾相对差 < 3%）",
      abs(C2[0] - C2[-1]) / C2[-1] < 0.03)


# ======================================================================
head("F4  尖锐速率：C^inf / C^{0,1} / C^{0,alpha} / 间断")


def circdist(x, a):
    d = np.abs(x - a)
    return np.minimum(d, 1 - d)


def fit_exponent(seq, Ns):
    return [np.log(seq[i] / seq[i + 1]) / np.log(Ns[i + 1] / Ns[i])
            for i in range(len(seq) - 1)]


Nsh = (512, 1024, 2048, 4096, 8192)

smooth = lambda x: 1 + 0.3 * np.cos(2 * np.pi * x)
rs = [float(rel_error_1d(N, smooth, k).max()) for N in Nsh]
exp_s = fit_exponent(rs, Nsh)[-1]
print("      smooth exponent=%.3f" % exp_s)
check("C^inf 指数 ≈ 2", abs(exp_s - 2.0) < 0.05)

kink = lambda x: 1 + 0.3 * np.cos(2 * np.pi * x) + 0.2 * circdist(x, 0.5)
rk = [float(rel_error_1d(N, kink, k).max()) for N in Nsh]
exp_k = fit_exponent(rk, Nsh)[-1]
print("      C^{0,1} kink exponent=%.3f" % exp_k)
check("C^{0,1} 指数 ≈ 1", abs(exp_k - 1.0) < 0.12)

for al, tol in ((0.5, 0.08), (0.25, 0.06)):
    ca = lambda x, a=al: 1 + 0.3 * circdist(x, 0.5) ** a
    ra = [float(rel_error_1d(N, ca, k).max()) for N in Nsh]
    ea = fit_exponent(ra, Nsh)[-1]
    print("      Holder alpha=%.2f exponent=%.3f" % (al, ea))
    check("C^{0,alpha=%.2f} 指数 ≈ alpha" % al, abs(ea - al) < tol)


def step_field(x, jump=0.37):
    return 1 + 0.3 * np.cos(2 * np.pi * x) + 0.3 * (x > jump)


rd = [float(rel_error_1d(N, step_field, k).max()) for N in (512, 1024, 2048, 4096)]
print("      discontinuous rel = %s" % ["%.3e" % r for r in rd])
check("间断场 O(1)：误差不随 N 下降", all(r > 0.5 for r in rd)
      and abs(rd[-1] - rd[0]) / rd[0] < 0.05)


# ======================================================================
head("F5  常数随 k 的经验增长落在 [2.5,3.3]")

Cs = {}
for kk in (4, 8, 12):
    r = float(rel_error_1d(1024, c1, kk).max())
    Cs[kk] = r * 1024 ** 2
print("      C(k) = %s" % {kk: round(v, 2) for kk, v in Cs.items()})
exp1 = np.log(Cs[8] / Cs[4]) / np.log(8 / 4)
exp2 = np.log(Cs[12] / Cs[8]) / np.log(12 / 8)
print("      empirical exponents = %.2f, %.2f" % (exp1, exp2))
check("常数经验指数在 [2.5,3.3]", 2.5 <= exp1 <= 3.3 and 2.5 <= exp2 <= 3.3)


# ======================================================================
head("F6  文档诚实边界")

check("H5 写成定理而非猜测",
      anchor("定理 R6.3", "反射—反转配对", "精确消去", "可以证明，不再是数值猜测"))
check("速率尖锐写明", anchor("速率是", "尖锐", "C^{1,1}", "C^{0,\\alpha}"))
check("三条前提写死", anchor("必须固定", "周期", "连续", "c_{\\min}"))
check("间断反例在位", anchor("间断", "O(1)", "失败"))
check("数值与证明分离", anchor("数值是旁证不是证明", "不替代"))
check("文档存在且非空", len(DOC) > 2000)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
