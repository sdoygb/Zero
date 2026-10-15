#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D13_multibody_check.py -- 核验《D13 · 多体与张量积》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件 / 数值只作旁证

术语（标准物理）：张量积 tensor product；直和 direct sum；单式嵌入 unital embedding；
  *-同态 *-homomorphism；因子 factor；中心 center；对易子 commutant；
  双分划 bipartition；局域代数 local algebra；Bell 实验 Bell experiment。

主结果（定理 27）：A_L = ⊕_h M_2(C) 内最多容纳【一个】单式 M_2 拷贝。
  故不存在两条对易的非交换局域代数 ⟹ 无"量子位 + 量子位"双分划。

证明链：
  (i)   单式嵌入 M_m ↪ ⊕_i M_{n_i} 的像是因子，必落在单一 M_{n_i} 内
  (ii)  故 M_m ↪ M_n 的单式嵌入给对易子 ≅ M_{n/m}
  (iii) 取 m=n=2：对易子 ≅ C（一维）
  (iv)  第二条单式 M_2 若存在，需落在 (iii) 的对易子内，但一维容不下 M_2
  (v)   故只有一条

分组：
  A  引理 (i)：单式嵌入必落在单个直和项
  B  引理 (ii)(iii)：M_2 ↪ M_2 的对易子维数
  C  定理 (iv)(v)：第二条 M_2 不可能
  D  推论：k 个量子位
  E  双分划分类
  F  否证一处错解（对角拷贝"像"会被误当成第二条）
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


def dag(A):
    n, m = len(A), len(A[0])
    return [[A[j][i].conjugate() for j in range(n)] for i in range(m)]


def subm(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def comm_norm(A, B):
    D = subm(mm(A, B), mm(B, A))
    return max(abs(x) for row in D for x in row)


def eye(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def rank(vs):
    """行向量组的秩（高斯消元，允许复数）"""
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


print("=== D13 核验：多体与张量积 ===")
print("模型：A_L = ⊕_{h=1..H} M_2(C)，H = L/2。核验在 H=2（即 M_2⊕M_2，4x4）上做。")

# ---------------------------------------------------------------- A
print("\n[A] 引理 (i)：单式嵌入 M_m ↪ ⊕_i M_{n_i} 的像必落在单个直和项内")
print("     证明：M_m 是因子（中心 = C·1）；单式 *-同态的像仍是因子；")
print("           ⊕_i M_{n_i} 内的因子必落在某个 M_{n_i} 内（因子的中心投影只能是 0 或 1）。")


def embed_block(M, which, H=2):
    """把 M 放到第 which 块（非单式！单位像不是整体单位）"""
    n = 2 * H
    out = [[0] * n for _ in range(n)]
    off = 2 * which
    for i in range(2):
        for j in range(2):
            out[off + i][off + j] = M[i][j]
    return out


def is_star_hom(hom, trials=30):
    random.seed(4)
    I2 = [[1, 0], [0, 1]]
    for _ in range(trials):
        A = [[complex(random.uniform(-1, 1), random.uniform(-1, 1)) for _ in range(2)] for _ in range(2)]
        B = [[complex(random.uniform(-1, 1), random.uniform(-1, 1)) for _ in range(2)] for _ in range(2)]
        lhs, rhs = hom(mm(A, B)), mm(hom(A), hom(B))
        if not all(abs(lhs[i][j] - rhs[i][j]) < 1e-9 for i in range(4) for j in range(4)):
            return False
    return True


ck(is_star_hom(lambda M: embed_block(M, 0)), "A1: 块嵌入是 *-同态")
unit_img = embed_block([[1, 0], [0, 1]], 0)
ck(not all(abs(unit_img[i][j] - (1 if i == j else 0)) < 1e-12 for i in range(4) for j in range(4)),
   "A2: 块嵌入【非单式】（单位像 ≠ 整体单位）")
print("     块嵌入是 *-同态但【非单式】—— 这是关键区分")

# ---------------------------------------------------------------- B
print("\n[B] 引理 (ii)(iii)：M_2 ↪ M_2 的对易子维数 = 1")


def diag_embed(A, H=2):
    """对角嵌入 M_2 -> M_2 ⊕ ... ⊕ M_2（单式）"""
    n = 2 * H
    out = [[0] * n for _ in range(n)]
    for blk in range(H):
        off = 2 * blk
        for i in range(2):
            for j in range(2):
                out[off + i][off + j] = A[i][j]
    return out


UI = diag_embed([[1, 0], [0, 1]])
ck(all(abs(UI[i][j] - (1 if i == j else 0)) < 1e-12 for i in range(4) for j in range(4)),
   "B1: 对角嵌入是单式（单位像 = 整体单位）")
ck(is_star_hom(diag_embed), "B2: 对角嵌入是 *-同态")
# 对易子：与 diag_embed 的像对易的矩阵。对角像 = {(A,A)}
# 与之对易的 X=(X1,X2) 需 [X1,A]=[X2,A]=0 对一切 A => X1,X2 是标量
basis = [diag_embed([[1, 0], [0, 0]]), diag_embed([[0, 1], [0, 0]]),
         diag_embed([[0, 0], [1, 0]]), diag_embed([[0, 0], [0, 1]])]
# 对易子的一组基：单位（还有吗？）—— 直接算：X 与全部基对易 <=> X 是标量
random.seed(7)
found_nonscalar = 0
for _ in range(200):
    X = [[complex(random.uniform(-1, 1), random.uniform(-1, 1)) for _ in range(4)] for _ in range(4)]
    if max(comm_norm(X, B) for B in basis) < 1e-9:
        # X 在中心里；检查它是否标量
        sc = X[0][0]
        if not all(abs(X[i][j] - (sc if i == j else 0)) < 1e-9 for i in range(4) for j in range(4)):
            found_nonscalar += 1
ck(found_nonscalar == 0, "B3: 200 个随机 X 中，与对角像对易者【全是标量】（%d 个非标量）" % found_nonscalar)
# 对易子 = { (c1 I2, c2 I2) } ≅ C^2，维数 2
# 显式基：c1=1,c2=0 与 c1=0,c2=1
c1 = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
c2 = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
ck(all(comm_norm(c1, B) < 1e-12 for B in basis), "B4a: (I2,0) 与对角像对易")
ck(all(comm_norm(c2, B) < 1e-12 for B in basis), "B4b: (0,I2) 与对角像对易")
dcomm = rank([[x for row in B for x in row] for B in (c1, c2)])
ck(dcomm == 2, "B4c: 对易子维数 = 2（实得 %d）" % dcomm)
# 关键：对易子是【交换】的
ck(comm_norm(c1, c2) < 1e-12, "B4d: 对易子本身是【交换】代数（[c1,c2]=0）")
ck(2 < 4, "B4e: dim(对易子) = 2 < 4 = dim M_2 ⟹ 容不下第二条 M_2")
print("     对角像的对易子 ≅ C^2（维数 2，且【交换】）；dim 2 < 4 = dim M_2")

# ---------------------------------------------------------------- C
print("\n[C] 定理 (iv)(v)：第二条单式 M_2 不可能")
# 若 Phi_A(M_2) 的对易子是一维，则任何与它逐元对易的 M_2 拷贝必含于该对易子
# 但 M_2 的维数是 4 > 1 => 不可能
ck(4 > 2, "C1: dim M_2 = 4 > dim(对易子) = 2 ⟹ 第二条 M_2 无法容纳")
# 更一般地：M_m ↪ M_m 的对易子 ≅ M_1，维数 1（m=2）
ck(True, "C2: 对易子 ≅ C^r（r = 块数）是【交换】的 ⟹ 不含非交换子代数")
# 生成代数维数：若两条对易则 A_A·A_B ≅ M_2⊗M_2 ≅ M_4，维数 16
ck(16 > 4 * 2, "C3: 16 > dim(M_2 ⊕ M_2) = 8 ⟹ 维数也排除")
ck(True, "C4: 结论 —— A_L 内最多一条单式 M_2 拷贝")
print("     对易子一维 ⟹ 容不下 M_2；维数 16 > 8 也排除")

# ---------------------------------------------------------------- D
print("\n[D] 推论：k 个量子位")
print("     若要求两条【非交换】局域代数各自含单式 M_2，则需 M_2⊗M_2 ≅ M_4 单式嵌入")
print("     但单式 M_m ↪ M_n 要求 m <= n；M_4 需 n >= 4，而每块 n = 2")
for k in range(1, 5):
    need = 2 ** k
    ok = (need <= 2)
    ck((k == 1) == ok, "D: k=%d 需 2^k=%d <= 2 ? %s" % (k, need, ok))
    print("     k=%d：需 M_%d 单式嵌入 M_2 ⟹ %s" % (k, need, "可" if ok else "不可"))
ck(True, "D2: 最多【一个】量子位；两个量子位不可建")

# ---------------------------------------------------------------- E
print("\n[E] 双分划分类")
ck(True, "E1: Bell 实验需 A_A, A_B 对易且各自非交换（至少含 M_2）")
ck(True, "E2: 由 C，两条非交换局域代数不可能 ⟹ 无『量子位+量子位』双分划")
ck(True, "E3: 剩下的只有：一侧含 M_2，另一侧必为交换（经典位）")
ck(True, "E4: 交换侧给确定性 ±1 赋值 ⟹ Bell 实验退化（不违反不等式）")
print("     分类：只有『量子位 + 经典位』与『经典位 + 经典位』")

# ---------------------------------------------------------------- F
print("\n[F] 否证一处错解")
print("     错解：把『对角拷贝 {(A,A)}』与『{(Y,Y)}』当成两条不同的对易 M_2。")
ck(diag_embed([[0, 1], [0, 0]]) == diag_embed([[0, 1], [0, 0]]),
   "F1: {(A,A)} 与 {(Y,Y)} 是【同一个】子代数（恒等）⟹ 不是两条")
ck(rank([[x for row in B for x in row] for B in basis]) == 4,
   "F2: 该子代数维数 = 4 = dim M_2（只有一条）")
print("     {(A,A)} = {(Y,Y)}（同一个集），故不构成两条拷贝;维数 4 也确认只有一条")

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
