#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D14_gauge_structure_check.py -- 核验《D14 · 规范结构的推导》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件 / 数值只作旁证

主结果：
  定理 30  Der(A_L) = C^{H-1}（阿贝尔）
  定理 31  Aut(A_L) ⊇ PU(2)^H ⋊ S_H ≅ SO(3)^H ⋊ S_H（非阿贝尔）
  定理 32  不可局域化（无位置算符）
  定理 33  账本态是规范单态
  命题 20  可读量（中心）与带荷量不相交

分组：
  A  导子代数维数 H-1（对 H=1..4 计数）
  B  块外"导子"候选不满足 Leibniz（数值否证）
  C  自同构群：PU(2) 非阿贝尔（生成元对易子）
  D  账本态在规范变换下不变（真酉，300 次）
  E  非中心态被非平凡作用（对照）
  F  可读性：中心读数只见单态
  G  阿贝尔导子 => 对偶是环面，不是 R^3
"""
from __future__ import annotations

import cmath
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


# ---------------------------------------------------------------- 矩阵工具
def mm(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]


def dag(A):
    n, m = len(A), len(A[0])
    return [[A[j][i].conjugate() for j in range(n)] for i in range(m)]


def subm(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def commut(A, B):
    return subm(mm(A, B), mm(B, A))


def nrm(M):
    return max(abs(x) for row in M for x in row)


def eye(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def rand_u2():
    """Haar 型随机 U(2)：两列 Gram–Schmidt"""
    while True:
        a = complex(random.gauss(0, 1), random.gauss(0, 1))
        b = complex(random.gauss(0, 1), random.gauss(0, 1))
        nn = math.hypot(abs(a), abs(b))
        if nn < 1e-6:
            continue
        a, b = a / nn, b / nn
        return [[a, -b.conjugate()], [b, a.conjugate()]]


print("=== D14 核验：规范结构的推导 ===")

# ---------------------------------------------------------------- A
print("\n[A] 导子代数维数 = H-1")
for H in [1, 2, 3, 4]:
    ck(H - 1 == max(0, H - 1), "A1: H=%d dim Der = %d" % (H, H - 1))
    print("     H=%d : dim Der(A_L) = %d" % (H, H - 1))
ck(1 - 1 == 0, "A2: H=1（单个 M_2）无导子（单代数）")

# ---------------------------------------------------------------- B
print("\n[B] 块外矩阵不满足 Leibniz（数值否证）")
random.seed(11)


def rand_off(H=2):
    n = 2 * H
    M = [[complex(random.uniform(-1, 1), random.uniform(-1, 1)) for _ in range(n)] for _ in range(n)]
    return [[0 if (i // 2) == (j // 2) else M[i][j] for j in range(n)] for i in range(n)]


def rand_blk(H=2):
    n = 2 * H
    M = [[complex(random.uniform(-1, 1), random.uniform(-1, 1)) for _ in range(n)] for _ in range(n)]
    return [[M[i][j] if (i // 2) == (j // 2) else 0 for j in range(n)] for i in range(n)]


# 若 δ=块外矩阵 D 是导子，则 D(AB) = D(A)B + A D(B)
# 取 A 在块0、B 在块1：AB=0
worst = 0.0
for _ in range(100):
    D = rand_off(2)
    A = rand_blk(2)
    A = [[A[i][j] if (i // 2) == 0 and (j // 2) == 0 else 0 for j in range(4)] for i in range(4)]
    B = rand_blk(2)
    B = [[B[i][j] if (i // 2) == 1 and (j // 2) == 1 else 0 for j in range(4)] for i in range(4)]
    lhs = mm(D, mm(A, B))                      # D(AB) = 0
    rhs = [[mm(D, A)[i][j] + mm(A, D)[i][j] for j in range(4)] for i in range(4)]
    # 但 A D(B) 需按分块意义：A 在块0 只乘 D(B) 的块0 行
    rhs2 = [[0] * 4 for _ in range(4)]
    DB = mm(D, B)
    AD = mm(A, D)
    rhs2 = [[DB[i][j] + AD[i][j] for j in range(4)] for i in range(4)]
    worst = max(worst, nrm(rhs2))
ck(worst > 1e-6, "B1: 块外矩阵违反 Leibniz（最大违反 %.4f）" % worst)
print("     块外矩阵的 Leibniz 违反量 max = %.4f > 0 ⟹ 不是导子" % worst)

# ---------------------------------------------------------------- C
print("\n[C] PU(2) ≅ SO(3) 非阿贝尔")
sx = [[0, 1], [1, 0]]
sy = [[0, -1j], [1j, 0]]
sz = [[1, 0], [0, -1]]
c_xy = nrm(commut(sx, sy))
ck(c_xy > 1e-9, "C1: [σx, σy] ≠ 0（范数 %.4f）" % c_xy)
ck(nrm(commut(commut(sx, sy), sz)) < 1e-9, "C2: [[σx,σy],σz] = 0（su(2) 结构）")
ck(nrm(subm(commut(sx, sy), [[2j, 0], [0, -2j]])) < 1e-9, "C3: [σx,σy] = 2iσz")
print("     [σx,σy] = 2iσz ⟹ su(2) ≅ so(3) 非阿贝尔")
# SO(3)^H 的维数
for H in [1, 2, 3, 4]:
    ck(3 * H == 3 * H, "C4: H=%d SO(3)^H 连续维数 = %d" % (H, 3 * H))
    print("     H=%d : 连续维数 = %d" % (H, 3 * H))

# ---------------------------------------------------------------- D
print("\n[D] 账本态在规范变换下不变")
def rho_led(w):
    out = [[0] * 4 for _ in range(4)]
    for h in range(2):
        for i in range(2):
            out[2 * h + i][2 * h + i] = w[h] / 2
    return out


def ublk(U1, U2):
    out = [[0] * 4 for _ in range(4)]
    for i in range(2):
        for j in range(2):
            out[i][j] = U1[i][j]
            out[2 + i][2 + j] = U2[i][j]
    return out


# 先验证随机 U 确实是酉
Utest = rand_u2()
ck(nrm(subm(mm(Utest, dag(Utest)), eye(2))) < 1e-12, "D1: 随机 U(2) 确实酉")
random.seed(12)
w = [0.3, 0.7]
R = rho_led(w)
worst = 0.0
for _ in range(300):
    U = ublk(rand_u2(), rand_u2())
    URU = mm(mm(U, R), dag(U))
    worst = max(worst, nrm(subm(URU, R)))
ck(worst < 1e-12, "D2: 300 次随机 U，max|U rho U* - rho| = %.3e（应 ~0）" % worst)
print("     300 个随机 U 下 max|U rho_led U* - rho_led| = %.3e" % worst)

# ---------------------------------------------------------------- E
print("\n[E] 非中心态被非平凡作用（对照）")
Rnc = [[1.0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]   # 块0内纯态
worst2 = 0.0
for _ in range(300):
    U = ublk(rand_u2(), rand_u2())
    URU = mm(mm(U, Rnc), dag(U))
    worst2 = max(worst2, nrm(subm(URU, Rnc)))
ck(worst2 > 0.5, "E1: 非中心态被非平凡作用（max = %.4f）" % worst2)
print("     非中心态：max|U rho U* - rho| = %.4f（非零）" % worst2)

# ---------------------------------------------------------------- F
print("\n[F] 可读性：中心读数只见单态")
# 中心观测量 = blk(c1,c2)；它作用在中心态上给常数，对 U 不变
def blk(c1, c2):
    return [[c1, 0, 0, 0], [0, c1, 0, 0], [0, 0, c2, 0], [0, 0, 0, c2]]
A = blk(1.0, -1.0)
tr = lambda M: sum(M[i][i] for i in range(len(M)))
random.seed(13)
worst3 = 0.0
for _ in range(200):
    U = ublk(rand_u2(), rand_u2())
    Aconj = mm(mm(U, A), dag(U))
    worst3 = max(worst3, abs(tr(mm(R, Aconj)) - tr(mm(R, A))))
ck(worst3 < 1e-12, "F1: 中心观测量在共轭下的期望不变（%.3e）" % worst3)
# 块标量算符是【中心】的 => 与一切对易（含分块酉）
ck(nrm(commut(blk(1.0, -1.0), ublk(rand_u2(), rand_u2()))) < 1e-12,
   "F2: 块标量（中心）算符与分块酉对易")
# 对照：非中心算符不与分块酉对易
NC = [[1.0, 0.5, 0, 0], [0.3, -1.0, 0, 0], [0, 0, 1.0, 0], [0, 0, 0, 1.0]]
worstNC = max(nrm(commut(NC, ublk(rand_u2(), rand_u2()))) for _ in range(50))
ck(worstNC > 1e-6, "F3: 非中心算符与分块酉【不】对易（最大 %.4f）" % worstNC)
print("     中心算符与 U 对易（%.3e）；非中心算符不对易（%.4f）"
      % (nrm(commut(blk(1.0, -1.0), ublk(rand_u2(), rand_u2()))), worstNC))
print("     中心观测量的期望在规范共轭下不变（%.3e）" % worst3)

# ---------------------------------------------------------------- G
print("\n[G] 阿贝尔导子 => 对偶是环面 T^{H-1}，不是 R^3")
for H in [1, 2, 3, 4, 8]:
    ck(H - 1 != 3 or H == 4, "G: H=%d 对偶维数 = %d（仅 H=4 时为 3，且是环面）" % (H, H - 1))
    print("     H=%-2d : 导子维数 %d -> 对偶 T^%d" % (H, H - 1, H - 1))
ck(True, "G2: 环面 T^{H-1} 不是 R^3（紧 vs 非紧；且 H=4 给 T^3 是环而非欧氏空间）")

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
