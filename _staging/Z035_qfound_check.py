#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D12_qfound_check.py -- 核验《D12 · 量子基础检查》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件 / 数值只作旁证

术语（标准物理）：Bell 不等式；CHSH 量；Tsirelson 界；乘积态 product state；
  单态 singlet state；中心态 central state；块标量观测量 block-scalar observable；
  交换子代数 commutative subalgebra；Gelfand 对偶；Kochen–Specker 定理；语境性 contextuality；
  联合分布 joint distribution。

分组：
  A  CHSH 的经典界（乘积态给 0；确定性模型给 |S|<=2）
  B  单态达到 Tsirelson 界 2*sqrt2（最优角，精确）
  C  随机抽样不超 Tsirelson 界
  D  中心态的结构与"朴素因子化"的【否证】
  E  块标量观测量构成交换子代数（与一切对易）
  F  非交换性 => 语境性可能；交换性 => 无语境性（Gelfand）
"""
from __future__ import annotations

import math
import random

FAIL = 0
NCHK = 0


def ck(cond, msg):
    global FAIL, NCHK
    NCHK += 1
    if not cond:
        FAIL += 1
        print("  FAIL:", msg)


# ---------------------------------------------------------------- 线性代数
def mm(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]


def dag(A):
    n, m = len(A), len(A[0])
    return [[A[j][i].conjugate() for j in range(n)] for i in range(m)]


def kron(A, B):
    n, m = len(A), len(A[0])
    p, q = len(B), len(B[0])
    return [[A[i // p][j // q] * B[i % p][j % q] for j in range(m * q)] for i in range(n * p)]


def tr(A):
    return sum(A[i][i] for i in range(len(A)))


def eye(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def comm_norm(A, B):
    X, Y = mm(A, B), mm(B, A)
    return max(abs(X[i][j] - Y[i][j]) for i in range(len(A)) for j in range(len(A)))


sx = [[0, 1], [1, 0]]
sy = [[0, -1j], [1j, 0]]
sz = [[1, 0], [0, -1]]


def sigma(n):
    return [[n[0] * sx[i][j] + n[1] * sy[i][j] + n[2] * sz[i][j] for j in range(2)] for i in range(2)]


def direction(theta, phi):
    return (math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta))


def rho_of(psi):
    return [[psi[i] * psi[j].conjugate() for j in range(len(psi))] for i in range(len(psi))]


def corr(rho, a, b):
    """E(a,b) = <psi| sigma(a) ⊗ sigma(b) |psi>"""
    A = kron(sigma(a), eye(2))
    B = kron(eye(2), sigma(b))
    return tr(mm(mm(rho, A), B)).real


print("=== D12 核验：量子基础检查 ===")

# ---------------------------------------------------------------- A
print("\n[A] CHSH 的经典界")
prod = rho_of([1, 0, 0, 0])          # |00>，乘积态
a = direction(math.pi / 2, 0)
ap = direction(math.pi / 2, math.pi / 2)
b = direction(math.pi / 2, math.pi / 4)
bp = direction(math.pi / 2, 3 * math.pi / 4)


def chsh(rho, a, ap, b, bp):
    return corr(rho, a, b) - corr(rho, a, bp) + corr(rho, ap, b) + corr(rho, ap, bp)


S_prod = chsh(prod, a, ap, b, bp)
ck(abs(S_prod) < 1e-12, "A1: 乘积态 |00> 的 CHSH 量 = 0（远低于经典界 2）")
# 确定性局域模型：A_i, B_j ∈ {±1}，S 的上界
best = 0
for A1 in (1, -1):
    for A2 in (1, -1):
        for B1 in (1, -1):
            for B2 in (1, -1):
                s = A1 * B1 - A1 * B2 + A2 * B1 + A2 * B2
                best = max(best, abs(s))
ck(best == 2, "A2: 确定性局域模型的 max|S| = 2（经典的 CHSH 界，实得 %d）" % best)
print("     乘积态 S = %.1f ; 确定性局域模型 max|S| = %d" % (S_prod, best))

# ---------------------------------------------------------------- B

# ---------------------------------------------------------------- A2'
print("\n[A2'] Tsirelson 的算子恒等式 S_op^2 = 4I + [A,A']⊗[B,B']")
def addm(*Ms):
    n, m = len(Ms[0]), len(Ms[0][0])
    return [[sum(M[i][j] for M in Ms) for j in range(m)] for i in range(n)]
def subm(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def smulm(c, A):
    return [[c * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def Sop_of(A, Ap, B, Bp):
    return addm(kron(A, B), smulm(-1, kron(A, Bp)), kron(Ap, B), kron(Ap, Bp))


def rand_dir():
    return direction(math.acos(random.uniform(-1, 1)), random.uniform(0, 2 * math.pi))


def comm_mat(A, B):
    return subm(mm(A, B), mm(B, A))


def maxeig(M):
    n = len(M)
    v = [1.0] * n
    for _ in range(600):
        w = [sum(M[i][j] * v[j] for j in range(n)) for i in range(n)]
        nrm = math.sqrt(sum(x * x for x in w))
        if nrm < 1e-15:
            return 0.0
        v = [x / nrm for x in w]
    w = [sum(M[i][j] * v[j] for j in range(n)) for i in range(n)]
    return sum(v[i] * w[i] for i in range(n))


Sop = Sop_of(sigma(a), sigma(ap), sigma(b), sigma(bp))
lhs = mm(Sop, Sop)
rhs = addm(smulm(4, eye(4)), kron(comm_mat(sigma(a), sigma(ap)), comm_mat(sigma(b), sigma(bp))))
ck(all(abs(lhs[i][j] - rhs[i][j]) < 1e-8 for i in range(4) for j in range(4)),
   "A2'a: 最优角上 S_op^2 = 4I + [A,A']⊗[B,B']（精确）")


worst = 0.0
for _ in range(60):
    Am = sigma(rand_dir()); Apm = sigma(rand_dir())
    Bm = sigma(rand_dir()); Bpm = sigma(rand_dir())
    sq = mm(Sop_of(Am, Apm, Bm, Bpm), Sop_of(Am, Apm, Bm, Bpm))
    re = [[sq[i][j].real for j in range(4)] for i in range(4)]
    worst = max(worst, maxeig(re))
ck(worst <= 8 + 1e-6, "A2'b: 60 组随机方向上 ||S_op||^2 = %.6f ≤ 8" % worst)
print("     恒等式精确成立；60 组随机方向 max ||S_op||^2 = %.6f ≤ 8" % worst)

print("\n[B] 单态达到 Tsirelson 界 2√2")
sing = rho_of([0, 1 / math.sqrt(2), -1 / math.sqrt(2), 0])
ck(abs(corr(sing, a, b) + sum(a[i] * b[i] for i in range(3))) < 1e-9,
   "B1: 单态上 E(a,b) = -a·b（标准结果，本文数值确认）")
S_sing = chsh(sing, a, ap, b, bp)
ck(abs(abs(S_sing) - 2 * math.sqrt(2)) < 1e-9,
   "B2: 单态的 |S| = 2√2 = %.10f（实得 %.10f）" % (2 * math.sqrt(2), abs(S_sing)))
ck(abs(S_sing) > 2, "B3: 单态违反经典界（|S| = %.6f > 2）" % abs(S_sing))
print("     E(a,b)=%.6f E(a,b')=%.6f E(a',b)=%.6f E(a',b')=%.6f" % (
    corr(sing, a, b), corr(sing, a, bp), corr(sing, ap, b), corr(sing, ap, bp)))
print("     S = %.10f  (|S| = 2√2 = %.10f)" % (S_sing, 2 * math.sqrt(2)))

# ---------------------------------------------------------------- C
print("\n[C] 随机抽样不超 Tsirelson 界")
random.seed(12)


bestS = 0.0
for _ in range(4000):
    s = abs(chsh(sing, rand_dir(), rand_dir(), rand_dir(), rand_dir()))
    bestS = max(bestS, s)
ck(bestS <= 2 * math.sqrt(2) + 1e-9,
   "C1: 4000 次抽样的最大 |S| = %.6f ≤ 2√2" % bestS)
print("     抽样最大 |S| = %.6f（界 %.6f）" % (bestS, 2 * math.sqrt(2)))

# ---------------------------------------------------------------- D
print("\n[D] 中心态：朴素因子化的【否证】")


def rho_central(w1, w2):
    return [[w1 / 2, 0, 0, 0], [0, w1 / 2, 0, 0], [0, 0, w2 / 2, 0], [0, 0, 0, w2 / 2]]


def rnd4():
    return [[complex(random.uniform(-1, 1), random.uniform(-1, 1)) for _ in range(4)] for _ in range(4)]


rho_c = rho_central(0.3, 0.7)
worst = 0.0
for _ in range(200):
    A, B = rnd4(), rnd4()
    AB = tr(mm(mm(rho_c, A), B))
    pr = tr(mm(rho_c, A)) * tr(mm(rho_c, B))
    worst = max(worst, abs(AB - pr))
# 朴素主张 "<AB> = <A><B>" 被否证
ck(worst > 0.1, "D1: 中心态上一般 A,B 不因子化（最大偏差 %.4f > 0.1）" % worst)
print("     中心态上 max|<AB> - <A><B>| = %.4f ⟹ 朴素因子化【假】" % worst)
# 但块标量观测量与一切对易
def blk(a1, a2):
    return [[a1, 0, 0, 0], [0, a1, 0, 0], [0, 0, a2, 0], [0, 0, 0, a2]]


# 块标量算符与【分块对角】算符对易
Tdiag = [[1, 0.2, 0, 0], [0.3, 2, 0, 0], [0, 0, 5, 0.1], [0, 0, 0.4, -1]]
ck(comm_norm(blk(1, -1), Tdiag) < 1e-12,
   "D2a: 块标量算符与分块对角算符对易")
# 但【不】与含块外相干的算符对易（=> 中心只是分块标量，不含一般算符）
Toff = [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]
ck(comm_norm(blk(1, -1), Toff) > 1e-9,
   "D2b: 块标量算符【不】与含块间相干的算符对易（故非与一切对易）")
ck(comm_norm(blk(1, -1), rho_central(0.3, 0.7)) < 1e-12,
   "D2c: 块标量算符与账本态对易（账本态是中心元）")
print("     块标量：与分块对角对易、与账本态对易；但不与块间相干对易")
# 中心态的 <AB> 与 <A><B> 在【块标量】上的关系：仍不因子化
AB = tr(mm(mm(rho_c, blk(1, -1)), blk(1, -1)))
pr = tr(mm(rho_c, blk(1, -1))) * tr(mm(rho_c, blk(1, -1)))
ck(abs(AB - pr) > 0.1, "D3: 即使限于块标量，也不因子化（差 %.4f）" % abs(AB - pr))
print("     块标量上 <AB>=%.4f, <A><B>=%.4f ⟹ 是经典混合，不是独立" % (AB, pr))

# ---------------------------------------------------------------- E
print("\n[E] 块标量观测量构成交换子代数（两两对易）")
for (a1, a2) in [(1.0, -1.0), (0.2, 0.7), (3.0, 0.0)]:
    for (b1, b2) in [(1.0, -1.0), (0.3, 0.9), (0.0, 2.0)]:
        ck(comm_norm(blk(a1, a2), blk(b1, b2)) < 1e-12,
           "E1: [blk%s, blk%s] = 0" % ((a1, a2), (b1, b2)))
print("     所有块标量（= 中心）算符两两对易 ⟹ 交换子代数")

# ---------------------------------------------------------------- F
print("\n[F] 非交换性 ⇒ 语境性可能；交换性 ⇒ 无语境性")
for (u, v) in [((1, 0, 0), (0, 1, 0)), ((1, 0, 0), (0, 0, 1)), ((0, 1, 0), (0, 0, 1))]:
    ck(comm_norm(sigma(u), sigma(v)) > 1e-9,
       "F1: 同一块内不同方向的赝自旋不可对易 [σ%s, σ%s]" % (u, v))
ck(comm_norm(sigma((1, 0, 0)), sigma((1, 0, 0))) < 1e-12, "F2: 同方向对易")
print("     同一块内：不同方向不可对易 ⟹ Kochen–Specker 型矛盾【可能】")
print("     块标量（中心）：全对易 ⟹ Gelfand 对偶 C(X) ⟹ 联合分布存在 ⟹ 无语境性")

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
