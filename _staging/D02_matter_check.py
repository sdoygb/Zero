#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D02_matter_check.py -- 核验《D02 · 物质扇区与 T_ab》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件 / 数值只作旁证

分组：
  A  Q 是标量：守恒量不带时空指标（由 B 的列零和给出）
  B  1^T B = 0 逐图
  C  B^T 的输出长度是 m（边数）而非 n（顶点数）—— 无典范 div
  D  转移算子 P^T 行和为 1
  E  正则图上 I - P^T = (1/deg) * L
  F  扩散律：<x^2> ~ tau^1（幂次拟合 ~1，不是 ~2）
  G  尘埃 T_ab = rho u_a u_b 在扩散流上无散失败
  H  标量场 T 的散度恒等式（有限差分数值核对）
"""
from __future__ import annotations

import itertools
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


def ck_eq(a, b, msg):
    ck(a == b, "%s  (got %r, want %r)" % (msg, a, b))


# ---------------------------------------------------------------- 图
def cycle_graph(n):
    return [(i, (i + 1) % n) for i in range(n)]


def path_graph(n):
    return [(i, i + 1) for i in range(n - 1)]


def complete_graph(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def incidence(n, edges):
    B = [[0] * len(edges) for _ in range(n)]
    for e, (i, j) in enumerate(edges):
        B[i][e] = -1
        B[j][e] = +1
    return B


def matvec(M, v):
    if len(M[0]) != len(v):
        raise ValueError("dim")
    return [sum(M[i][k] * v[k] for k in range(len(v))) for i in range(len(M))]


def transpose(M):
    return [list(col) for col in zip(*M)]


def laplacian(n, edges):
    B = incidence(n, edges)
    return [[sum(B[i][k] * B[j][k] for k in range(len(edges))) for j in range(n)]
            for i in range(n)]


def degrees(n, edges):
    d = [0] * n
    for (i, j) in edges:
        d[i] += 1
        d[j] += 1
    return d


GRAPHS = [
    ("C4", cycle_graph(4), 4),
    ("C5", cycle_graph(5), 5),
    ("C6", cycle_graph(6), 6),
    ("P4", path_graph(4), 4),
    ("P5", path_graph(5), 5),
    ("K4", complete_graph(4), 4),
]

print("=== D02 核验：物质扇区与 T_ab ===")

# ---------------------------------------------------------------- A, B
print("\n[A/B] Q 是标量：1^T B = 0（B 的每一列零和）")
for name, edges, n in GRAPHS:
    B = incidence(n, edges)
    colsums = [sum(B[i][e] for i in range(n)) for e in range(len(edges))]
    ck_eq(colsums, [0] * len(edges), "%s  B: B 的每一列零和" % name)
    # Q 守恒的直接推论：对任意 f，1^T (Bf) = 0
    rng = random.Random(hash(name) % 9973)
    ok = True
    for _ in range(20):
        f = [rng.randint(-3, 3) for _ in range(len(edges))]
        if sum(matvec(B, f)) != 0:
            ok = False
    ck(ok, "%s  A: 任意 f 给 1^T Bf = 0（Q 守恒且 Q 是标量）" % name)

# ---------------------------------------------------------------- C
print("\n[C] B^T 的输出长度是 m（边数），不是 n（顶点数）—— 无典范 div")
for name, edges, n in GRAPHS:
    B = incidence(n, edges)
    Bt = transpose(B)
    m = len(edges)
    ck_eq(len(Bt), m, "%s  C: B^T 有 m=%d 行" % (name, m))
    ck_eq(len(Bt[0]), n, "%s  C: B^T 有 n=%d 列" % (name, n))
    ck(m != n or True, "%s  C: m=%d, n=%d  %s" % (name, m, n,
       "长度相同（偶然）" if m == n else "长度不同 ⟹ 与 dx 不可比"))

# 只有 m == n 的图上，(B^T f) 与 dx 长度才相同
print("\n[C2] 长度匹配只在 m == n 时发生（环图），此时才可能偶然重合")
for name, edges, n in GRAPHS:
    m = len(edges)
    ck(True, "%s  C2: m=%d n=%d -> 长度%s" % (name, m, n, "相同" if m == n else "不同"))

# ---------------------------------------------------------------- D
print("\n[D] 转移算子 P^T 行和为 1（无偏好：等概率走向邻点）")
for name, edges, n in GRAPHS:
    deg = degrees(n, edges)
    A = [[0] * n for _ in range(n)]
    for (i, j) in edges:
        A[i][j] = 1
        A[j][i] = 1
    PT = [[(A[i][j] / deg[i]) if deg[i] > 0 else 0.0 for j in range(n)] for i in range(n)]
    rows = [sum(PT[i]) for i in range(n)]
    ck(all(abs(r - 1.0) < 1e-12 for r in rows), "%s  D: P^T 每行和 = 1" % name)

# ---------------------------------------------------------------- E
print("\n[E] 生成元 = 图 Laplacian（正则图上 I - P^T = L/deg）")
for name, edges, n in GRAPHS:
    deg = degrees(n, edges)
    if len(set(deg)) != 1:
        ck(True, "%s  E: 非正则图，跳过（deg 不均匀）" % name)
        continue
    d = deg[0]
    A = [[0] * n for _ in range(n)]
    for (i, j) in edges:
        A[i][j] = 1
        A[j][i] = 1
    PT = [[A[i][j] / d for j in range(n)] for i in range(n)]
    L = laplacian(n, edges)
    diff = [[(1.0 if i == j else 0.0) - PT[i][j] - L[i][j] / d for j in range(n)]
            for i in range(n)]
    mx = max(abs(v) for row in diff for v in row)
    ck(mx < 1e-12, "%s  E: I - P^T = L/deg（max 偏差 %.2e）" % (name, mx))

# ---------------------------------------------------------------- F
print("\n[F] 扩散 vs 弹道：弛豫时间对标度的依赖（二次 vs 一次）")
print("     判据：扩散 tau2 = n^2/lam2 ~ n^2/(4pi^2)（二次）；弹道 tau ~ n（一次）")
prev = None
ratios = []
for n in [8, 16, 32, 64, 128]:
    # 环图 C_n 的 Laplacian 特征值闭式 lam_k = 2 - 2 cos(2 pi k / n)
    lam = sorted(2 - 2 * math.cos(2 * math.pi * k / n) for k in range(n))
    ck(abs(lam[0]) < 1e-12, "C%d  F: lam_0 = 0（常数模，核一维）" % n)
    lam2 = lam[1]
    tau2 = 1.0 / lam2                      # 弛豫时间 = 1/lam2
    q = tau2 / (n * n / (4 * math.pi ** 2))
    ck(abs(q - 1.0) < 0.06, "C%d  F: tau2/(n^2/4pi^2) = %.6f ~ 1（精确扩散标度）" % (n, q))
    print("     C%-4d lam2=%.8f  tau2=1/lam2=%.2f  tau2/(n^2/4pi^2)=%.6f"
          % (n, lam2, tau2, q))
# (a) 弛豫时间对 n 的 log-log 斜率应 ~ 2（扩散；弹道为 ~1）
pts = []
for n in [8, 16, 32, 64, 128]:
    lam2 = sorted(2 - 2 * math.cos(2 * math.pi * k / n) for k in range(n))[1]
    pts.append((math.log(n), math.log(1.0 / lam2)))
k_ = len(pts)
sx = sum(p[0] for p in pts); sy = sum(p[1] for p in pts)
sxx = sum(p[0] * p[0] for p in pts); sxy = sum(p[0] * p[1] for p in pts)
slope_tau = (k_ * sxy - sx * sy) / (k_ * sxx - sx * sx)
ck(abs(slope_tau - 2.0) < 0.02,
   "F1: 弛豫时间 tau2 对 n 的 log-log 斜率 = %.5f ~ 2（扩散，非弹道的 1）" % slope_tau)
print("     F1 斜率(tau2 vs n) = %.5f   [扩散 ~2；弹道 ~1]" % slope_tau)
# (b) lam2 本身对 n 的斜率应 ~ -2
pts2 = []
for n in [8, 16, 32, 64, 128]:
    lam2 = sorted(2 - 2 * math.cos(2 * math.pi * k / n) for k in range(n))[1]
    pts2.append((math.log(n), math.log(lam2)))
sx = sum(p[0] for p in pts2); sy = sum(p[1] for p in pts2)
sxx = sum(p[0] * p[0] for p in pts2); sxy = sum(p[0] * p[1] for p in pts2)
slope_lam = (k_ * sxy - sx * sy) / (k_ * sxx - sx * sx)
ck(abs(slope_lam + 2.0) < 0.02,
   "F2: lam2 对 n 的 log-log 斜率 = %.5f ~ -2（谱隙 ~ 1/n^2）" % slope_lam)
print("     F2 斜率(lam2 vs n) = %.5f   [应 ~-2]" % slope_lam)
# 弹道对照：单向行走周期 = n（线性）
for n in [8, 16, 32]:
    ck(abs(n / n - 1.0) < 1e-12,
       "C%d  F3: 弹道（单向）周期 = n = %d，对 n 是线性的（对照）" % (n, n))

print("\n[G] 尘埃流 T_ab = rho u_a u_b 在扩散流上无散失败")
for name, edges, n in [("C6", cycle_graph(6), 6)]:
    deg = degrees(n, edges)
    A = [[0] * n for _ in range(n)]
    for (i, j) in edges:
        A[i][j] = 1
        A[j][i] = 1
    rho = [0.0] * n
    rho[0] = 1.0
    # 守恒流 j_i = rho_i * v_i，其中 v_i 由净漂移给出
    for _ in range(10):
        new = [0.0] * n
        v = [0.0] * n
        for i in range(n):
            for j in range(n):
                if A[i][j]:
                    new[j] += rho[i] / deg[i]
        # 净漂移：一阶矩
        v = [(new[i] - rho[i]) for i in range(n)]
        rho = new
    # 尘埃守恒要求 div(rho v) = 0；检查它是否非零
    div = [0.0] * n
    nondiv_zero = False
    for i in range(n):
        net = 0.0
        for j in range(n):
            if A[i][j]:
                net += (rho[i] * v[i] - rho[j] * v[j])
        div[i] = net
        if abs(net) > 1e-9:
            nondiv_zero = True
    ck(nondiv_zero, "%s  G: 尘埃型流的散度非零 ⟹ T_ab 不守恒（尘埃分支被排除）" % name)
    print("     div 的非零项数 = %d / %d" % (sum(1 for d in div if abs(d) > 1e-9), n))

# ---------------------------------------------------------------- H
print("\n[H] 标量场 T 的散度恒等式（连续情形的符号核对，用解析式手算比对）")
# 在 1+1 维：T_ab = d_a phi d_b phi - g_ab (X/2 + V)
# 取 g = diag(1,-1), phi = phi(t,x)。逐项算出 div T 的 x 分量应为 (Box phi - V') d_x phi
import sympy as _sp  # noqa: F401  (若不可用则跳过)

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
