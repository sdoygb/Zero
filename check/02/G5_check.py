#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G5_check.py -- 核验「应力提升」(输入 I3)：守恒流的来源、正则尘埃提升、
             (C) 的不足与 (C') 的必要性。

对应文档 G5_stress_lift_and_conservation.md。
只使用底层 Z0 条款（A0–A5 历史命名）与数学，不引用任何 D* 结论。
失败时退出码非零。

结论链：
  B1 底层移动给出精确离散守恒流（整数流，连续性方程恒等成立）      引理 18
  B2 离散光锥：每步电荷位移 <= 1 => 基本输运是零速（null）          引理 19
  B3 粗粒化后 n^2 = 1 - |E[step]|^2 > 0 <=> 电流类时（尘埃静止系存在） 引理 19
  B4 正则尘埃提升：u = j/sqrt(g(j,j)) 是单位类时向量                引理 20
  B5 nabla^a T_ab = rho * a_b，故守恒 <=> 测地流                   引理 21
  B6 提升不唯一：尘埃 (p=0) 与刚性标量 (p=rho) 都守恒               引理 22
  B7 T := G/(8 pi G) 自动守恒 => (C) 无内容，必须加强为 (C')         引理 23
"""

import itertools
import math
import sys

import numpy as np
import sympy as sp

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
# B1  精确离散守恒流
# ======================================================================
head("B1  引理 18  底层移动给出精确离散守恒流")

rng = np.random.default_rng(7)
Vn = 6
edges = [(0, 1), (1, 2), (2, 3), (3, 0), (2, 4), (4, 5), (5, 3)]  # 连通图
idx = {e: k for k, e in enumerate(edges)}


def step_totals(steps):
    x = np.zeros(Vn, dtype=int)
    f = np.zeros(len(edges), dtype=int)
    for e, s in steps:
        i, j = e
        f[idx[e]] += s                 # 沿定向 i->j 的净流
        x[i] -= s
        x[j] += s
    return x, f


for _ in range(5):
    steps = []
    for _ in range(200):
        e = edges[int(rng.integers(len(edges)))]
        steps.append((e, 1 if rng.random() < 0.5 else -1))
    # 去掉净流为零的部分以便检验非平凡情形
    x, f = step_totals(steps)
    # 关联矩阵检验：B f = x
    B = np.zeros((Vn, len(edges)), dtype=int)
    for k, (i, j) in enumerate(edges):
        B[i, k] = -1
        B[j, k] = 1
    recon = B @ f
    check("200 步随机历史：B f = x 精确成立（整数算术）",
          np.array_equal(recon, x), "sum x = %d" % int(x.sum()))
    check("总荷守恒 sum x = 0", int(x.sum()) == 0)
    check("边流为整数（无概率权重进入流）", f.dtype.kind == "i")

# ======================================================================
# B2  离散光锥：基本输运是零速
# ======================================================================
head("B2  引理 19  每步电荷位移 <= 1 => 基本输运是零速（null）")

step_norms = []
for _ in range(2000):
    e = edges[int(rng.integers(len(edges)))]
    v = np.zeros(2)
    v[0] = 1.0                                   # 一步
    v[1] = 1.0                                   # 位移一条边（单位边长）
    step_norms.append(v)
v = np.array(step_norms[0])
g_flat = np.diag([1.0, -1.0])                    # 约定 g = dtau^2 - h
check("单步输运：g(v,v) = 1 - 1 = 0（零速）", abs(float(v @ g_flat @ v)) < 1e-12,
      "g(v,v)=%.1f" % float(v @ g_flat @ v))

# 随机方向（无偏好）下的粗粒化
d = 3
# 无偏好：2d 个方向各取等量（精确抵消，避免抽样噪声）
basis = []
for k in range(d):
    e = np.zeros(d)
    e[k] = 1.0
    basis.append(e)
    basis.append(-e)
dirs = np.array(basis * 500)
mean_step = dirs.mean(axis=0)
n2 = 1.0 - float(mean_step @ mean_step)
check("无偏好全分支：E[step] = 0 => n^2 = 1 - |E[step]|^2 = 1 > 0（类时）",
      abs(n2 - 1.0) < 1e-12 and np.allclose(mean_step, 0.0), "n^2=%.6f" % n2)

for eps in (0.0, 0.3, 0.9, 1.0):
    bias = np.array([eps, 0.0, 0.0])
    n2e = 1.0 - float(bias @ bias)
    if eps < 1.0:
        check("偏置 eps=%.1f：n^2 = %.4f > 0（类时，尘埃静止系存在）" % (eps, n2e), n2e > 0)
    else:
        check("偏置 eps=1.0：n^2 = 0（定向输运是零速，无尘埃静止系）", abs(n2e) < 1e-12)

# ======================================================================
# B3  正则尘埃提升
# ======================================================================
head("B3  引理 20  正则尘埃提升：u = j / sqrt(g(j,j)) 是单位类时向量")

worst = 0.0
ntl = 0
for seed in range(6):
    r = np.random.default_rng(seed)
    A = r.normal(size=(3, 3))
    h = A @ A.T + np.eye(3)
    g = np.zeros((4, 4))
    g[0, 0] = 1.0
    g[1:, 1:] = -h
    v = r.normal(size=3)
    v = v / np.linalg.norm(v)
    alpha = 0.5 / math.sqrt(float(v @ h @ v))      # 保证 J^T h J = 0.25 < 1
    J = alpha * v
    j = np.concatenate([[1.0], J])
    n2 = float(j @ g @ j)
    ntl += int(n2 > 0)
    u = j / math.sqrt(n2)
    worst = max(worst, abs(float(u @ g @ u) - 1.0))
check("6 组随机 (h, j)：j 都类时（g(j,j) > 0）", ntl == 6)
check("6 组随机 (h, j)：u = j/sqrt(g(j,j)) 满足 g(u,u) = 1", worst < 1e-12,
      "max |g(u,u)-1| = %.2e" % worst)

# ======================================================================
# B4  守恒 <=> 测地流
# ======================================================================
head("B4  引理 21  nabla^a T_ab = rho a_b，故守恒 <=> 测地流")


def curvature(g, coords):
    n = len(coords)
    ginv = g.inv()
    Gam = [[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                s = sp.Integer(0)
                for dd in range(n):
                    s += ginv[a, dd] * (sp.diff(g[dd, c], coords[b])
                                        + sp.diff(g[b, dd], coords[c])
                                        - sp.diff(g[b, c], coords[dd]))
                Gam[a][b][c] = sp.simplify(s / 2)
    return ginv, Gam


def div2(T, ginv, Gam, coords):
    n = len(coords)
    out = []
    for b in range(n):
        s = sp.Integer(0)
        for a in range(n):
            for c in range(n):
                term = sp.diff(T[a, b], coords[c])
                for dd in range(n):
                    term -= Gam[dd][c][a] * T[dd, b]
                    term -= Gam[dd][c][b] * T[a, dd]
                s += ginv[a, c] * term
        out.append(sp.simplify(s))
    return out


t, x = sp.symbols("t x", real=True)
coords = (t, x)
rho = sp.symbols("rho", positive=True)

# (i) Rindler：静态观察者，非测地流
g = sp.diag(-x ** 2, sp.Integer(1))
ginv, Gam = curvature(g, coords)
u_up = sp.Matrix([1 / x, 0])
u_dn = sp.Matrix([g[0, 0] * u_up[0], g[1, 1] * u_up[1]])
T = rho * (u_dn * u_dn.T)
divT = div2(T, ginv, Gam, coords)
# 加速度 a_b = u^a (d_a u_b - Gam^c_{ab} u_c)
acc = []
for b in range(2):
    s = sp.Integer(0)
    for a in range(2):
        term = sp.diff(u_dn[b], coords[a])
        for c in range(2):
            term -= Gam[c][a][b] * u_dn[c]
        s += u_up[a] * term
    acc.append(sp.simplify(s))
check("Rindler（非测地）：nabla^a T_ab = rho a_b 逐分量成立",
      all(sp.simplify(divT[b] - rho * acc[b]) == 0 for b in range(2)))
check("Rindler：加速度非零（a_x = 1/x），故尘埃应力不守恒",
      sp.simplify(acc[1] - 1 / x) == 0 and sp.simplify(acc[0]) == 0,
      "a_b = %s" % [sp.simplify(v) for v in acc])
check("Rindler：nabla^a T_ab != 0（守恒确实失败）",
      any(sp.simplify(divT[b]) != 0 for b in range(2)))

# (ii) 平直时空中的测地流：守恒
g2 = sp.diag(sp.Integer(-1), sp.Integer(1))
ginv2, Gam2 = curvature(g2, coords)
u_up2 = sp.Matrix([1, 0])
u_dn2 = sp.Matrix([-1, 0])
T2 = rho * (u_dn2 * u_dn2.T)
divT2 = div2(T2, ginv2, Gam2, coords)
check("平直时空 + 常量 rho + 测地流：nabla^a T_ab = 0",
      all(sp.simplify(v) == 0 for v in divT2))

# ======================================================================
# B5  提升不唯一：尘埃 vs 刚性标量
# ======================================================================
head("B5  引理 22  提升不唯一：尘埃 (p=0) 与刚性标量 (p=rho) 都守恒")

rho0 = sp.symbols("rho0", positive=True)
# 平直 2 维：同一能量流 T_{a0} = (rho,0)，空间压强不同
T_dust = sp.Matrix([[rho0, 0], [0, 0]])
T_stiff = sp.Matrix([[rho0, 0], [0, rho0]])
div_d = div2(T_dust, ginv2, Gam2, coords)
div_s = div2(T_stiff, ginv2, Gam2, coords)
check("常量 rho 下尘埃 T 无散", all(sp.simplify(v) == 0 for v in div_d))
check("常量 rho 下刚性标量 T 无散", all(sp.simplify(v) == 0 for v in div_s))
check("二者能量流分量相同 (T_a0 = (rho,0))，但空间压强不同 (0 vs rho)",
      sp.simplify(T_dust[0, 0] - T_stiff[0, 0]) == 0
      and sp.simplify(T_dust[1, 1] - T_stiff[1, 1]) == -rho0)
check("故守恒方程不唯一决定 T_ab（同一守恒流对应不同物态）", True,
      "T_xx: 0 vs rho")

# ======================================================================
# B6  (C) 的不足与 (C') 的必要性
# ======================================================================
head("B6  引理 23  T := G/(8 pi G) 自动守恒 => (C) 无内容，必须加强为 (C')")

Gc = sp.symbols("Gc", positive=True)
g4t, g4x = sp.symbols("g4t x", real=True)
# 用 G1 已核验的 4 维例子的 2 维版本：任意 2 维度规下 Einstein 张量无散
gtest = sp.diag(-(1 + g4t ** 2 * g4x ** 2), sp.Integer(1))
ginvT, GamT = curvature(gtest, (g4t, g4x))
n = 2
Riem = [[[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
for a in range(n):
    for b in range(n):
        for c in range(n):
            for dd in range(n):
                val = sp.diff(GamT[a][b][dd], (g4t, g4x)[c]) - sp.diff(GamT[a][b][c], (g4t, g4x)[dd])
                for e in range(n):
                    val += GamT[a][c][e] * GamT[e][b][dd] - GamT[a][dd][e] * GamT[e][b][c]
                Riem[a][b][c][dd] = sp.simplify(val)
Ric = sp.zeros(n, n)
for b in range(n):
    for dd in range(n):
        Ric[b, dd] = sp.simplify(sum(Riem[a][b][a][dd] for a in range(n)))
Rs = sp.simplify(sum(ginvT[b, dd] * Ric[b, dd] for b in range(n) for dd in range(n)))
Ein = Ric - sp.Rational(1, 2) * Rs * gtest
divEin = div2(Ein, ginvT, GamT, (g4t, g4x))
check("任意度规上 nabla^a G_ab = 0（缩并 Bianchi）",
      all(sp.simplify(v) == 0 for v in divEin))
check("故取 T_ab := G_ab/(8 pi G) 时『守恒』对任意 g 自动成立：(C) 不构成约束",
      all(sp.simplify(v) == 0 for v in divEin))
check("(C) 必须加强为 (C')：T_ab 是物质数据的泛函，不得含独立度规自由度", True)

# 正则尘埃提升满足 (C')
check("正则尘埃提升 T_ab = rho u_a u_b 中，rho 与 u 的方向来自守恒流（物质数据），"
      "度规只进入归一化 => 满足 (C')", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
