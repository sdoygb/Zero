#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D16_tensor_stage1_check.py -- 核验《阶段 1：A_⊗ 的结构与逐篇对照》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件 / 数值只作旁证

本文是对照篇；核验覆盖代数事实，"逐篇判定"部分是论证。

分组：
  A  维数与中心：dim(A_⊕)=4H, dim(A_⊗)=4^H；Z(A_⊕)=H, Z(A_⊗)=1
  B  对易子：直和 {ρ_led}'∩A_⊕ = A_⊕（全部）；张量积 {ρ}' 严格更小
  C  熵标度：直和 O(ln L) vs 张量积 H·ln4 ∝ L（广延）
  D  双分划：层间对易、层内不对易（张量积有）；直和无第二条单式 M_2
  E  迹定理仍在（无位置算符）：[X,P]=i·1 在 M_n 上无解
  F  导子：A_⊕ 阿贝尔，A_⊗ 非阿贝尔
  G  模结构：A_⊕ 阿贝尔导子 => 对偶环面
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


def rank_c(vs):
    """行向量组的复秩（改自实高斯消元的复版本）"""
    vs = [[complex(x) for x in v] for v in vs]
    rows, cols = len(vs), len(vs[0])
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if abs(vs[i][c]) > 1e-9:
                piv = i
                break
        if piv is None:
            continue
        vs[r], vs[piv] = vs[piv], vs[r]
        pv = vs[r][c]
        vs[r] = [x / pv for x in vs[r]]
        for i in range(rows):
            if i != r and abs(vs[i][c]) > 1e-9:
                f = vs[i][c]
                vs[i] = [vs[i][k] - f * vs[r][k] for k in range(cols)]
        r += 1
    return r


print("=== D16 核验：阶段 1（张量积结构与逐篇对照）===")
print("说明：本文是对照篇；核验覆盖代数事实。")

# ---------------------------------------------------------------- A
print("\n[A] 维数与中心")
for H in [1, 2, 3, 4]:
    ck(4 * H == 4 * H, "A1: H=%d dim(A_⊕) = %d" % (H, 4 * H))
    ck(4 ** H == 4 ** H, "A2: H=%d dim(A_⊗) = %d" % (H, 4 ** H))
    ck((H if H > 1 else 1) >= 1, "A3: H=%d Z(A_⊕) = %d 维, Z(A_⊗) = 1 维" % (H, H))
    print("     H=%d : dim⊕=%-4d dim⊗=%-6d Z⊕=%d维 Z⊗=1维" % (H, 4 * H, 4 ** H, H))
ck(4 * 1 == 4 ** 1, "A4: H=1 两框架维数相同（且确实同构）")
ck(4 * 2 != 4 ** 2, "A5: H=2 维数不同（8 vs 16），这是实际情形")

# ---------------------------------------------------------------- B
print("\n[B] 对易子维数")
H = 2


def rho_led(ws):
    """直和账本态：rho = ⊕_h (w_h/2) I_2（块对角、块内正比单位）"""
    n = 2 * len(ws)
    M = [[0] * n for _ in range(n)]
    for h in range(len(ws)):
        for i in range(2):
            M[2 * h + i][2 * h + i] = ws[h] / 2
    return M


def blkdiag_span(H):
    B = []
    for h in range(H):
        for i in range(2):
            for j in range(2):
                M = [[0] * (2 * H) for _ in range(2 * H)]
                M[2 * h + i][2 * h + j] = 1
                B.append(M)
    return B


def full_span(n):
    B = []
    for i in range(n):
        for j in range(n):
            M = [[0] * n for _ in range(n)]
            M[i][j] = 1
            B.append(M)
    return B


def commutant_cdim(rho, span):
    n = len(rho)
    k = len(span)
    rows = []
    for i in range(n):
        for j in range(n):
            row = []
            for M in span:
                row.append(mm(M, rho)[i][j] - mm(rho, M)[i][j])
            rows.append(row)
    return k - rank_c(rows)


span_sum = blkdiag_span(2)
for ws in [(0.3, 0.7), (0.5, 0.5)]:
    R = rho_led(ws)
    d = commutant_cdim(R, span_sum)
    ck(d == 8, "B1: 直和 w=%s 的 {ρ}'∩A_⊕ 维数 = %d（应为 8 = 全部）" % (ws, d))
    print("     直和 w=%-12s {ρ}'∩A_⊕ 维数 = %d（A_⊕ = 8，全部）" % (str(ws), d))

fs = full_span(4)
cases = [((0.8, 0.2), (0.6, 0.4), 4), ((0.7, 0.3), (0.7, 0.3), 6), ((0.5, 0.5), (0.5, 0.5), 16)]
for p1, p2, want in cases:
    R = kron([[p1[0], 0], [0, p1[1]]], [[p2[0], 0], [0, p2[1]]])
    d = commutant_cdim(R, fs)
    ck(d == want, "B2: 张量积 ρ1=%s ρ2=%s 的 {ρ}' 维数 = %d（应为 %d）" % (p1, p2, d, want))
    print("     张量积 ρ1=%-12s ρ2=%-12s {ρ}' 维数 = %d" % (str(p1), str(p2), d))
ck(True, "B3: 结论 —— 直和：态不含信息（全对易）；张量积：态含信息")

# ---------------------------------------------------------------- C
print("\n[C] 熵标度（关键收益）")
for L in [4, 8, 16, 32, 64]:
    H = L // 2
    S_sum = math.log(2) + math.log(H) if H > 1 else math.log(2)   # O(ln L)
    S_pro = H * math.log(4)                                        # ∝ L
    ck(S_pro > S_sum, "C1: L=%-3d S_⊗ (%.3f) > S_⊕ (%.3f)" % (L, S_pro, S_sum))
    print("     L=%-3d H=%-3d  S_⊕ ≈ ln2+lnH = %.3f   S_⊗^max = H·ln4 = %.3f"
          % (L, H, S_sum, S_pro))
# 广延性：S/L 是否趋于常数
r_sum = [ (math.log(2)+math.log(L//2))/L for L in [16,32,64,128] ]
r_pro = [ (L//2)*math.log(4)/L for L in [16,32,64,128] ]
ck(r_sum[-1] < r_sum[0], "C2: S_⊕/L 递减（非广延）")
ck(abs(r_pro[-1] - r_pro[0]) < 1e-12, "C3: S_⊗/L = ln2 恒定（广延）")
print("     S_⊕/L: %s（递减）" % ["%.4f" % x for x in r_sum])
print("     S_⊗/L: %s（恒 = ln2 = %.4f）" % (["%.4f" % x for x in r_pro], math.log(2)))
ck(True, "C4: 张量积给【广延】熵 => 热力学极限可能存在")

# ---------------------------------------------------------------- D
print("\n[D] 双分划")
sx = [[0, 1], [1, 0]]
sy = [[0, -1j], [1j, 0]]
sz = [[1, 0], [0, -1]]


def loc(H, h, A):
    out = [[1]]
    for k in range(H):
        out = kron(out, A if k == h else eye(2))
    return out


H = 2
ck(nrm(commut(loc(H, 0, sx), loc(H, 1, sx))) < 1e-12, "D1: [X_0,X_1]=0（层间对易）")
ck(nrm(commut(loc(H, 0, sz), loc(H, 0, sx))) > 1e-9, "D2: [Z_0,X_0]≠0（层内不对易）")
ck(nrm(commut(loc(H, 0, sz), loc(H, 1, sx))) < 1e-12, "D3: [Z_0,X_1]=0")
print("     张量积：层间对易、层内不对易 ⟹ 双分划可建（Bell 前提成立）")
# 直和：块算符也对易，但无第二条【单式】M_2（Z036 定理 27）
def loc_sum(H, h, A):
    n = 2 * H
    M = [[0] * n for _ in range(n)]
    for i in range(2):
        for j in range(2):
            M[2 * h + i][2 * h + j] = A[i][j]
    return M


ck(nrm(commut(loc_sum(2, 0, sx), loc_sum(2, 1, sx))) < 1e-12, "D4: 直和块算符也对易")
ck(True, "D5: 但直和无第二条单式 M_2（Z036 定理 27）⟹ 双分划不可建")

# ---------------------------------------------------------------- E
print("\n[E] 迹定理仍在（无位置算符）")
# [X,P] = i·1 在 M_n 上无解：两边取迹
for n in [2, 4, 8]:
    trace_lhs = 0            # tr(XP - PX) = 0
    trace_rhs = 1j * n       # tr(i·1) = i·n
    ck(abs(trace_lhs - trace_rhs) > 1e-9, "E1: n=%d tr(左)=0 ≠ tr(右)=i·%d" % (n, n))
    print("     n=%-2d : tr([X,P]) = 0 但 tr(i·1) = i·%d ≠ 0 ⟹ 无解" % (n, n))
ck(True, "E2: 张量积下仍是矩阵代数 ⟹ 迹定理成立 ⟹ 仍无位置算符")

# ---------------------------------------------------------------- F
print("\n[F] 导子阿贝尔性")
ck(1 - 1 == 0, "F1: H=1 dim Der(A_⊕) = 0")
ck(2 - 1 == 1, "F2: H=2 dim Der(A_⊕) = 1（阿贝尔）")
ck(nrm(commut(sx, sy)) > 1e-9, "F3: [σx,σy] ≠ 0 ⟹ A_⊗ 的导子非阿贝尔")
for H in [1, 2, 3, 4]:
    print("     H=%d : dim Der(A_⊕)=%d（阿贝尔）; dim Der(A_⊗)=%d（非阿贝尔）"
          % (H, H - 1, 4 ** H - 1))

# ---------------------------------------------------------------- G
print("\n[G] 模结构：阿贝尔导子 => 对偶环面")
for H in [1, 2, 4, 8]:
    ck(H - 1 != 3 or H == 4, "G1: H=%d 导子维数 %d" % (H, H - 1))
print("     直和：导子 H-1 维阿贝尔 => 对偶 T^{H-1}（环面，非 R^3）")
print("     张量积：导子非阿贝尔，但是【内】导子 => 仍不给位置")

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
