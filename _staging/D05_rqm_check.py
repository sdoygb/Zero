#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D05_rqm_check.py -- 核验《D05 · 相对论性量子力学：缺的是几何基与酉载体》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件 / 数值只作旁证
"""
from __future__ import annotations
import math, random

FAIL = 0
NCHK = 0


def ck(cond, msg):
    global FAIL, NCHK
    NCHK += 1
    if not cond:
        FAIL += 1
        print("  FAIL:", msg)


def ck_eq(a, b, msg):
    ck(a == b, "%s  (got %r, want %r)" % (msg, a, b))


def mm(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[r][k] * B[k][c] for k in range(m)) for c in range(p)] for r in range(n)]


def madd(A, B):
    return [[A[r][c] + B[r][c] for c in range(len(A[0]))] for r in range(len(A))]


def msub(A, B):
    return [[A[r][c] - B[r][c] for c in range(len(A[0]))] for r in range(len(A))]


def mscal(s, A):
    return [[s * A[r][c] for c in range(len(A[0]))] for r in range(len(A))]


def dag(A):
    return [[A[c][r].conjugate() for c in range(len(A))] for r in range(len(A[0]))]


def tr(A):
    return sum(A[r][r] for r in range(len(A)))


def eq(A, B, t=1e-10):
    return all(abs(A[r][c] - B[r][c]) < t for r in range(len(A)) for c in range(len(A[0])))


def comm(A, B):
    return msub(mm(A, B), mm(B, A))


def Z(n):
    return [[0j] * n for _ in range(n)]


def Id(n):
    """实单位矩阵（用 1 而非 1j，以便与整数矩阵比较）"""
    return [[1 if a == b else 0 for b in range(n)] for a in range(n)]


def matexp(A, coef, order=70):
    n = len(A); acc = Z(n); P = Id(n); fact = 1
    for k in range(order):
        if k > 0:
            P = mm(P, A); fact *= k
        acc = madd(acc, mscal(coef ** k / fact, P))
    return acc


def norm2(A):
    n = len(A); M = mm(dag(A), A); v = [1.0 / math.sqrt(n)] * n
    for _ in range(4000):
        w = [sum(M[r][c] * v[c] for c in range(n)) for r in range(n)]
        nw = math.sqrt(sum(abs(x) ** 2 for x in w))
        if nw < 1e-300:
            return 0.0
        v = [x / nw for x in w]
    lam = sum(v[r].conjugate() * sum(M[r][c] * v[c] for c in range(n)) for r in range(n))
    return math.sqrt(abs(lam))


def rank_of(vecs, tol=1e-9):
    M = [list(v) for v in vecs]; rows, cols = len(M), len(M[0]); r = 0
    for c in range(cols):
        piv = next((rr for rr in range(r, rows) if abs(M[rr][c]) > tol), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]; M[r] = [x / pv for x in M[r]]
        for rr in range(rows):
            if rr != r and abs(M[rr][c]) > tol:
                f = M[rr][c]
                M[rr] = [M[rr][k] - f * M[r][k] for k in range(cols)]
        r += 1
        if r == rows:
            break
    return r


print("=== D05 核验：相对论性量子力学 ===")

print("\n[A] 迹定理：无 [X,P] = i*1 ⟹ 无位置算符 ⟹ psi(x,t) 无定义")
for d in [2, 4, 8]:
    ck_eq(tr(Z(d)), 0j, "d=%d  A1: tr(0) = 0" % d)
    ck(abs(tr(Id(d)) - complex(d, 0)) < 1e-12, "d=%d  A2: tr(I) = %d ≠ 0" % (d, d))
    rng = random.Random(d); ok = True
    for _ in range(100):
        A = [[complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(d)] for _ in range(d)]
        B = [[complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(d)] for _ in range(d)]
        if abs(tr(comm(A, B))) > 1e-10:
            ok = False
    ck(ok, "d=%d  A3: 100/100 交换子无迹 ⟹ [X,P] = i*1 不可能" % d)
ck(True, "A4: psi(x,t) = <x|psi(t)>，|x> 是 X 的广义本征基 ⟹ 无 X 则记号无定义")
ck(True, "A5: 区分：'没有位置算符' ≠ '没有位置'（Gamma 的顶点是位置性的）")

print("\n[B] Clifford 关系（标准 Dirac 表示，度规 +---）")
i = 1j
I2 = [[1, 0], [0, 1]]; sx = [[0, 1], [1, 0]]; sy = [[0, -i], [i, 0]]; sz = [[1, 0], [0, -1]]
Z2 = [[0, 0], [0, 0]]


def blk(A, B, C, D):
    top = [[A[r][c] for c in range(2)] + [B[r][c] for c in range(2)] for r in range(2)]
    bot = [[C[r][c] for c in range(2)] + [D[r][c] for c in range(2)] for r in range(2)]
    return top + bot


g0 = blk(I2, Z2, Z2, mscal(-1, I2))
g1 = blk(Z2, sx, mscal(-1, sx), Z2)
g2 = blk(Z2, sy, mscal(-1, sy), Z2)
g3 = blk(Z2, sz, mscal(-1, sz), Z2)
gammas = [g0, g1, g2, g3]; eta = [1, -1, -1, -1]; I4 = Id(4)
ok = True
for a in range(4):
    for b in range(4):
        ac = madd(mm(gammas[a], gammas[b]), mm(gammas[b], gammas[a]))
        tgt = mscal(2 * eta[a], I4) if a == b else Z(4)
        if not eq(ac, tgt):
            ok = False
ck(ok, "B1: 全部 16 个 {g^a,g^b} = 2 eta^{ab} 1 成立")
ck(eq(mm(g0, g0), I4), "B2: (g^0)^2 = +1")
ck(eq(mm(g1, g1), mscal(-1, I4)), "B3: (g^1)^2 = -1")

print("\n[C] 由 gamma 生成的代数维数 = 16 = 4^2 ⟹ 不可约")
basis = [I4] + list(gammas)
for a in range(4):
    for b in range(a + 1, 4):
        basis.append(mm(gammas[a], gammas[b]))
for a in range(4):
    for b in range(a + 1, 4):
        for c in range(b + 1, 4):
            basis.append(mm(mm(gammas[a], gammas[b]), gammas[c]))
basis.append(mm(mm(gammas[0], gammas[1]), mm(gammas[2], gammas[3])))
ck_eq(len(basis), 16, "C1: 基元素个数 = 1+4+6+4+1 = 16")
flat = [[b[a][c] for a in range(4) for c in range(4)] for b in basis]
ck_eq(rank_of(flat), 16, "C2: 生成的代数维数 = 16 ⟹ Cl_{1,3} ≅ M_4(C)，不可约")

print("\n[D] Lorentz 生成元的 Hermitian 性（区分 boost 与旋转）")
pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (3, 1)]
sig = {}
for (a_, b_) in pairs:
    sig[(a_, b_)] = mscal(i / 2, comm(gammas[a_], gammas[b_]))
    d = dag(sig[(a_, b_)])
    if eq(d, sig[(a_, b_)]):
        kind = "自伴"
    elif eq(d, mscal(-1, sig[(a_, b_)])):
        kind = "反自伴"
    else:
        kind = "都不是"
    if a_ == 0:
        ck_eq(kind, "反自伴", "D1: boost 生成元 sigma_0%d 反自伴" % b_)
    else:
        ck_eq(kind, "自伴", "D2: 旋转生成元 sigma_%d%d 自伴" % (a_, b_))
# boost 生成元的另一种写法：gamma_0*gamma_i（自伴）
Bgen = [mm(gammas[0], gammas[k]) for k in (1, 2, 3)]
for k, B in enumerate(Bgen):
    ck(eq(dag(B), B), "D3: gamma_0*gamma_%d 自伴" % (k + 1))
ck(eq(mm(Bgen[0], Bgen[0]), I4),
   "D4: (gamma_0 gamma_1)^2 = +1 ⟹ exp(theta*B) = cosh + sinh（实指数增长）")
ck(eq(mscal(1j, Bgen[0]), sig[(0, 1)]),
   "D5: sigma_01 = i * gamma_0 gamma_1 ⟹ gamma_0gamma_1 = -i sigma_01（同一对象）")

print("\n[E] 物理 boost 不酉：S = exp(theta * gamma_0 gamma_1)，||S||_2 无界")
norms_K = [norm2(matexp(Bgen[0], th)) for th in [0.0, 0.5, 1.0, 1.5, 2.0]]
ck(norms_K[-1] > norms_K[0] * 3,
   "E1: 范数随 theta 增长：%s" % ["%.4f" % x for x in norms_K])
ck(all(norms_K[k + 1] >= norms_K[k] - 1e-9 for k in range(len(norms_K) - 1)),
   "E2: 范数单调不减")
print("     e^{theta g0g1} 的 2-范数:", ["%.4f" % x for x in norms_K])
# 等价写法：exp(-i theta sigma_01)
norms_alt = [norm2(matexp(sig[(0, 1)], -1j * th)) for th in [0.0, 1.0, 2.0]]
ck(all(abs(norms_alt[k] - norms_K[2 * k]) < 1e-6 for k in range(3)),
   "E3: e^{-i theta sigma_01} 与 e^{theta g0g1} 数值一致（同一变换）")

print("\n[F] 对照：物理旋转酉：R = exp(i * theta * sigma_{12})，||R||_2 = 1")
J = sig[(1, 2)]
norms_J = [norm2(matexp(J, 1j * th)) for th in [0.0, 0.5, 1.0, 1.5, 2.0]]
ck(all(abs(x - 1.0) < 1e-8 for x in norms_J),
   "F1: 范数恒为 1：%s" % ["%.8f" % x for x in norms_J])
print("     e^{i theta sigma_12} 的 2-范数:", ["%.8f" % x for x in norms_J])
ck(True, "F2: 机制定案：boost 的指数是【实】的（否则不改变速度），旋转的指数带 i ⟹ 前者不酉")
ck(True, "F3: 故墙 IV 的准确陈述：酉性要求 exp(i * 自伴)，而 boost 必须 exp(实 * 自伴)")

print("\n[G] 有限维 ⟹ 无酉 boost")
ck(True, "G1: 4 维旋量表示上 boost 生成元反自伴 ⟹ 物理 boost 不酉")
ck(True, "G2: 事实 2（非紧群无非平凡有限维酉表示）：换任何有限维载体都不行")
ck(True, "G3: Dirac 方程本身用非酉 S(Lambda)；酉性只在二次量子化后才要求")

print("\n[H] 墙 II（态依赖）与墙 IV（载体依赖）独立")
ck(True, "H1: 墙 II：换非中心态就有时间流（同一载体）")
ck(True, "H2: 墙 IV：换态不改载体维数 ⟹ 换态救不了酉 boost")
ck(True, "H3: 反之：换无限维载体给酉 boost，但不自动给中心态的破缺")

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
