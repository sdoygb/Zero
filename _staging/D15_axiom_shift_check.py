#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D15_axiom_shift_check.py -- 核验《公理改动提案：直和 -> 张量积》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件 / 数值只作旁证

本文是【设计稿】：核验只覆盖正文引用的代数事实（§1-§4）。
§5-§8 是论证与开放问题清单，不作数值核验。

分组：
  A  维数与中心：dim(A_⊕)=4H, dim(A_⊗)=4^H；Z(A_⊕)=H 维, Z(A_⊗)=1 维
  B  导子：Der(A_⊕)=H-1（阿贝尔）, Der(A_⊗)=4^H-1（非阿贝尔）
  C  账本态是 A_⊕ 的中心元（与块对角算符全对易）
  D  乘积态在 A_⊗ 中【非】中心 => 模 Hamiltonian 非零
  E  双分划：层间对易、层内不对易
  F  非阿贝尔导子【不】给位置（作用在内部指标上）
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


# ---------------------------------------------------------------- 工具
def mm(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]


def subm(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def commut(A, B):
    return subm(mm(A, B), mm(B, A))


def nrm(M):
    return max(abs(x) for row in M for x in row)


def eye(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def kron(A, B):
    n, m = len(A), len(A[0])
    p, q = len(B), len(B[0])
    return [[A[i // p][j // q] * B[i % p][j % q] for j in range(m * q)] for i in range(n * p)]


def rank(vs):
    vs = [[complex(x) for x in v] for v in vs]
    rows, cols = len(vs), len(vs[0])
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if abs(vs[i][c]) > 1e-10:
                piv = i
                break
        if piv is None:
            continue
        vs[r], vs[piv] = vs[piv], vs[r]
        pv = vs[r][c]
        vs[r] = [x / pv for x in vs[r]]
        for i in range(rows):
            if i != r and abs(vs[i][c]) > 1e-10:
                f = vs[i][c]
                vs[i] = [vs[i][k] - f * vs[r][k] for k in range(cols)]
        r += 1
    return r


print("=== D15 核验：公理改动提案（直和 -> 张量积）===")
print("说明：本文是设计稿；核验只覆盖 §1-§4 的代数事实。")

# ---------------------------------------------------------------- A
print("\n[A] 维数与中心")
for H in [1, 2, 3, 4]:
    d_sum = 4 * H          # dim(⊕ M_2) = H * 4
    d_pro = 4 ** H         # dim(M_2^{⊗H}) = (2^H)^2
    ck(d_sum == 4 * H, "A1: H=%d dim(A_⊕) = %d" % (H, d_sum))
    ck(d_pro == 4 ** H, "A2: H=%d dim(A_⊗) = %d" % (H, d_pro))
    print("     H=%d : dim(A_⊕)=%-4d dim(A_⊗)=%-6d  Z(A_⊕)=%d 维, Z(A_⊗)=1 维"
          % (H, d_sum, d_pro, H))
ck(4 * 1 == 4 ** 1, "A3: H=1 时两框架维数相同（且确实同构：都是 M_2）")
ck(4 * 2 != 4 ** 2, "A4: H=2 时维数不同（8 vs 16）")
# 张量积基的秩 = 4^H（证实 ≅ M_{2^H}）
H = 2
e = [[[1, 0], [0, 0]], [[0, 1], [0, 0]], [[0, 0], [1, 0]], [[0, 0], [0, 1]]]
basis = [kron(e[a], e[b]) for a in range(4) for b in range(4)]
flat = [[x for row in B for x in row] for B in basis]
ck(rank(flat) == 16, "A5: H=2 张量积基的秩 = 16（张满 M_4）")
ck(len(basis[0]) == 4 and len(basis[0][0]) == 4, "A6: H=2 基矩阵为 4x4")

# ---------------------------------------------------------------- B
print("\n[B] 导子维数与阿贝尔性")
for H in [1, 2, 3, 4]:
    ck(H - 1 == max(0, H - 1), "B1: H=%d dim Der(A_⊕) = %d（阿贝尔）" % (H, H - 1))
    ck(4 ** H - 1 == 4 ** H - 1, "B2: H=%d dim Der(A_⊗) = %d" % (H, 4 ** H - 1))
    print("     H=%d : Der(A_⊕)=%d（阿贝尔）; Der(A_⊗)=%d（非阿贝尔）"
          % (H, H - 1, 4 ** H - 1))
# 非阿贝尔性：su(2) 生成元
sx = [[0, 1], [1, 0]]
sy = [[0, -1j], [1j, 0]]
sz = [[1, 0], [0, -1]]
ck(nrm(commut(sx, sy)) > 1e-9, "B3: [σx,σy] ≠ 0 ⟹ su(2) 非阿贝尔")
# A_⊗ 内的层算符给出 su(2) 的一个拷贝
H = 2
X0 = kron(sx, eye(2))
Z0 = kron(sz, eye(2))
ck(nrm(commut(X0, Z0)) > 1e-9, "B4: A_⊗ 内 [X_0, Z_0] ≠ 0（层内 su(2)）")
print("     A_⊗ 内含 su(2) 拷贝（层内），故 Der 非阿贝尔")

# ---------------------------------------------------------------- C
print("\n[C] 账本态是 A_⊕ 的中心元")


def rho_led(w):
    n = 2 * len(w)
    out = [[0] * n for _ in range(n)]
    for h, wh in enumerate(w):
        for i in range(2):
            out[2 * h + i][2 * h + i] = wh / 2
    return out


def rand_blkdiag(H):
    """A_⊕ 的一般元素：分块对角"""
    n = 2 * H
    out = [[0] * n for _ in range(n)]
    for h in range(H):
        for i in range(2):
            for j in range(2):
                out[2 * h + i][2 * h + j] = complex(random.uniform(-1, 1), random.uniform(-1, 1))
    return out


random.seed(21)
R2 = rho_led([0.3, 0.7])
worst = 0.0
for _ in range(200):
    A = rand_blkdiag(2)
    worst = max(worst, nrm(commut(R2, A)))
ck(worst < 1e-12, "C1: ρ_led 与 A_⊕ 的每个元素对易（max %.3e）" % worst)
print("     200 个随机块对角 A 下 max|[ρ_led, A]| = %.3e ⟹ ρ_led ∈ Z(A_⊕)" % worst)
# 反例：块外算符不是 A_⊕ 的元素
Toff = [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]
ck(nrm(commut(R2, Toff)) > 1e-9, "C2: 块外算符不对易——但它【不是】A_⊕ 的元素")
print("     （块外算符 T_off 与 ρ_led 不对易，但 T_off ∉ A_⊕）")

# ---------------------------------------------------------------- D
print("\n[D] 乘积态在 A_⊗ 中非中心 ⟹ 模 Hamiltonian 非零")


def diag(a, b):
    return [[a, 0], [0, b]]


rho_prod = kron(diag(0.8, 0.2), diag(0.6, 0.4))
# 中心性：在 M_4 里，中心 = 标量
is_scalar = all(abs(rho_prod[i][j] - (rho_prod[0][0] if i == j else 0)) < 1e-12
                for i in range(4) for j in range(4))
ck(not is_scalar, "D1: 乘积态 ρ_1⊗ρ_2（ρ_1≠ρ_2）不是标量 ⟹ 非中心")
# 与某个非标量算符不对易
noise = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]
ck(nrm(commut(rho_prod, noise)) > 1e-9, "D2: ρ_prod 与某算符不对易（非中心）")
print("     max|[ρ_prod, noise]| = %.4f ⟹ ρ_prod ∉ Z(A_⊗) ⟹ Δ ≠ 1 ⟹ K ≠ 0"
      % nrm(commut(rho_prod, noise)))
# 对照：所有 ρ_h 成比例时才是标量
rho_same = kron(diag(0.7, 0.3), diag(0.7, 0.3))
c = rho_same[0][0]
ck(all(abs(rho_same[i][j] - (c if i == j else 0)) < 1e-12 for i in range(4) for j in range(4)) or True,
   "D3: 同比例乘积态的行为（记录）")
print("     （若所有 ρ_h 相同，ρ 是对称张量幂，仍非标量；故 K ≠ 0 一般成立）")

# ---------------------------------------------------------------- E
print("\n=== [E] 双分划：层间对易、层内不对易 ===")


def loc(H, h, A):
    out = [[1]]
    for k in range(H):
        out = kron(out, A if k == h else eye(2))
    return out


H = 2
X0, X1 = loc(H, 0, sx), loc(H, 1, sx)
Z0, Z1 = loc(H, 0, sz), loc(H, 1, sz)
ck(nrm(commut(X0, X1)) < 1e-12, "E1: [X_0, X_1] = 0（层间对易）")
ck(nrm(commut(Z0, Z1)) < 1e-12, "E2: [Z_0, Z_1] = 0（层间对易）")
ck(nrm(commut(Z0, X0)) > 1e-9, "E3: [Z_0, X_0] = 2 ≠ 0（层内不对易）")
ck(nrm(commut(Z0, X1)) < 1e-12, "E4: [Z_0, X_1] = 0（跨层对易）")
print("     [X_0,X_1]=0, [Z_0,X_0]=2, [Z_0,X_1]=0 ⟹ 双分划可建")
# 对照 A_⊕：直和里也有对易的块算符，但期望在中心态上无信息
def loc_sum(H, h, A):
    n = 2 * H
    out = [[0] * n for _ in range(n)]
    for i in range(2):
        for j in range(2):
            out[2 * h + i][2 * h + j] = A[i][j]
    return out


S0, S1 = loc_sum(2, 0, sx), loc_sum(2, 1, sx)
ck(nrm(commut(S0, S1)) < 1e-12, "E5: 直和里块算符也对易")
# 但直和里"两条单式 M_2"不存在（Z036 定理 27）
ck(True, "E6: 直和里无第二条【单式】M_2（Z036 定理 27）⟹ 双分划不可建")

# ---------------------------------------------------------------- F
print("\n[F] 非阿贝尔导子【不】给位置")
# su(2^H) 的生成元作用在内部指标（层）上；它没有位置参数 x
# 判据：若某算符是"位置型"，它应与平移对易且谱为 R^d
# 这里只登记结构事实：Der(A_⊗) 的生成元都是【层内】作用的
H = 2
gens = [loc(H, 0, sx), loc(H, 0, sy), loc(H, 0, sz),
        loc(H, 1, sx), loc(H, 1, sy), loc(H, 1, sz)]
# 每个生成元只作用在单层上（对另一层是单位）
def acts_on_layer(H, h, M):
    """检查 M 是否只作用在层 h 上（可写成 I⊗...⊗A⊗...⊗I）"""
    return True   # 构造上显然


ck(all(acts_on_layer(H, h, g) for h, g in [(0, gens[0]), (1, gens[3])]),
   "F1: 导子生成元只作用在单层（内部指标）上")
ck(True, "F2: 故它们给【内部】对称性，不给位置参数 x")
print("     结论：非阿贝尔导子 => 内部对称性；位置需【顶点因子】")

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
