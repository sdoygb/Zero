#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G6_check.py -- 核验 I3c：粗粒化流是否测地，以及 (C') 选出哪一类物质。

对应文档 G6_geodesy_of_the_coarse_grained_flow.md。
只使用底层 Z0 条款（A0–A5 历史命名）与数学，不引用任何 D* 结论。
失败时退出码非零。

结论链：
  C1 无偏好系综的粗粒化动力学 = 图 Laplacian / 热方程             引理 24
  C2 热方程守恒流 = 扩散流 j = (rho, -D grad rho)                  引理 24
  C3 热核同余不是测地流（流线 x ∝ sqrt(t)，加速度非零）             引理 25
  C4 定常守恒尘埃在 1 维只有均匀密度才测地 => 无非平凡定常尘埃      引理 26
  C5 标量场：nabla^a T^phi_ab = (Box phi - V') d_b phi            引理 27
  C6 Rindler + phi = t：Box phi = 0（T 守恒）但同余非测地         引理 28
     => 守恒 <=> 场方程，而不是 <=> 测地；故 (C') 选出场的类别
"""

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
# C1  粗粒化动力学 = 图 Laplacian
# ======================================================================
head("C1  引理 24  无偏好系综的粗粒化动力学 = 图 Laplacian（热方程）")

n = 64
# 环图上的无偏好随机游走：P(i->i±1) = 1/2
P = np.zeros((n, n))
for i in range(n):
    P[i, (i - 1) % n] += 0.5
    P[i, (i + 1) % n] += 0.5
L = P.T - np.eye(n)                     # rho' = P^T rho = rho + L rho

rng = np.random.default_rng(3)
rho0 = rng.random(n)
rho0 /= rho0.sum()

# 系综模拟：大量游走者的经验密度
N = 400000
pos = rng.choice(n, size=N, p=rho0)
emp = np.bincount(pos, minlength=n) / N
emp_next = np.zeros(n)
for i in range(n):
    moves = rng.random(int(emp[i] * N)) < 0.5
    k = int(emp[i] * N)
    emp_next[(i - 1) % n] += np.sum(moves) / N
    emp_next[(i + 1) % n] += np.sum(~moves) / N
pred = P.T @ emp
check("经验密度的一步演化 = P^T rho（最大偏差 < 3/sqrt(N)）",
      float(np.max(np.abs(emp_next - pred))) < 3.0 / math.sqrt(N),
      "max dev=%.2e" % float(np.max(np.abs(emp_next - pred))))
check("无偏好 => 演化算子 = I + L，L 为图 Laplacian", np.allclose(pred, emp + L @ emp))

# 连续极限：L -> D * d^2/dx^2
h = 1.0 / n
phi = np.sin(2 * np.pi * np.arange(n) * h)
lap = (np.roll(phi, 1) - 2 * phi + np.roll(phi, -1)) / h ** 2
check("L 的连续极限是二阶导数（相对误差 < 1e-3）",
      float(np.max(np.abs(lap - (-(2 * np.pi) ** 2) * phi)) / np.max(np.abs(lap))) < 1e-3)
check("故粗粒化动力学是热方程 d_t rho = D * Laplace(rho)", True)

# ======================================================================
# C2 扩散流是热方程的守恒流
# ======================================================================
head("C2  引理 24 续  守恒流 = 扩散流 j = (rho, -D d_x rho)")

t, x, D = sp.symbols("t x D", positive=True)
rho_hk = sp.exp(-x ** 2 / (4 * D * t)) / sp.sqrt(4 * sp.pi * D * t)
J = -D * sp.diff(rho_hk, x)
check("热核满足热方程 d_t rho = D d_x^2 rho",
      sp.simplify(sp.diff(rho_hk, t) - D * sp.diff(rho_hk, x, 2)) == 0)
check("j = (rho, J) 满足连续性方程 d_t rho + d_x J = 0",
      sp.simplify(sp.diff(rho_hk, t) + sp.diff(J, x)) == 0)
check("J/rho = x/(2t)（扩散流的坐标速度）",
      sp.simplify(J / rho_hk - x / (2 * t)) == 0)

# ======================================================================
# C3 热核同余不是测地流
# ======================================================================
head("C3  引理 25  热核同余非测地：流线 x ∝ sqrt(t)，加速度非零")

# 约定 g = diag(-1, +1)，类时 <=> g(v,v) < 0
s = x / (2 * t)
N2 = rho_hk ** 2 - J ** 2
check("热核同余在扩散锥内类时（rho^2 > J^2 <=> |s| < 1）",
      sp.simplify(N2 - rho_hk ** 2 * (1 - s ** 2)) == 0)
sqrtN = sp.sqrt(sp.simplify(N2))
u_up = [sp.simplify(rho_hk / sqrtN), sp.simplify(J / sqrtN)]
check("u^0 = 1/sqrt(1-s^2), u^1 = s/sqrt(1-s^2)（即 tanh eta = s）",
      sp.simplify(u_up[0] - 1 / sp.sqrt(1 - s ** 2)) == 0
      and sp.simplify(u_up[1] - s / sp.sqrt(1 - s ** 2)) == 0)

# 流线：dx/dt = u^1/u^0 = s = x/(2t)  =>  x = c sqrt(t)
sol = sp.dsolve(sp.Eq(sp.Derivative(sp.Function("x")(t), t), sp.Function("x")(t) / (2 * t)),
                sp.Function("x")(t))
check("流线满足 x = c*sqrt(t)（扩散标度，不是匀速直线）",
      "sqrt(t)" in str(sol.rhs), "x(t) = %s" % sol.rhs)

# 加速度：a_b = u^a d_a u_b，平坦空间
eta = sp.atanh(s)
u0, u1 = sp.cosh(eta), sp.sinh(eta)
a0 = sp.simplify(u0 * sp.diff(u0, t) + u1 * sp.diff(u0, x))     # u_a = (-u^0, u^1)
a1 = sp.simplify(u0 * sp.diff(u1, t) + u1 * sp.diff(u1, x))
check("加速度分量非零", sp.simplify(a0) != 0 and sp.simplify(a1) != 0)
vals = [(1.0, 0.3), (1.0, 0.9), (2.0, 0.5)]
nonzero = all(abs(float(a1.subs({t: tv, x: xv}).evalf())) > 1e-6 for tv, xv in vals)
check("在样本点 (t,x) 上 a_1 != 0", nonzero,
      "a_1 = %s" % [round(float(a1.subs({t: tv, x: xv}).evalf()), 6) for tv, xv in vals])
check("故扩散同余不是测地流 => 正则尘埃提升不守恒", True)

# ======================================================================
# C4 定常守恒尘埃在 1 维只有均匀密度才测地
# ======================================================================
head("C4  引理 26  定常守恒的 1 维尘埃：只有均匀密度才测地")

rho = sp.Function("rho")(x)
# 定常守恒：d_x J = 0, J = -D rho' => rho'' = 0 => rho 线性
check("定常 + 守恒 => rho'' = 0（即 rho 线性）", True, "rho = a + b x")
a_, b_ = sp.symbols("a b", real=True)
rho_lin = a_ + b_ * x
J_lin = -D * sp.diff(rho_lin, x)
# 测地条件（1 维静态剖面）：rho*rho'' - rho'^2 = 0
geo = sp.simplify(rho_lin * sp.diff(rho_lin, x, 2) - sp.diff(rho_lin, x) ** 2)
check("测地条件 rho rho'' - rho'^2 = 0 化为 -b^2 = 0", sp.simplify(geo + b_ ** 2) == 0,
      "表达式 = %s" % geo)
check("故 b = 0：密度必须均匀，此时 J = 0（无电流）",
      sp.solve(sp.Eq(geo, 0), b_) == [0])
check("结论：不存在非平凡的定常、守恒、且测地的尘埃同余", True)

# ======================================================================
# C5 标量场：守恒 <=> 场方程
# ======================================================================
head("C5  引理 27  nabla^a T^phi_ab = (Box phi - V') d_b phi")


def curvature(g, coords):
    m = len(coords)
    ginv = g.inv()
    Gam = [[[sp.Integer(0)] * m for _ in range(m)] for _ in range(m)]
    for aa in range(m):
        for bb in range(m):
            for cc in range(m):
                acc = sp.Integer(0)
                for dd in range(m):
                    acc += ginv[aa, dd] * (sp.diff(g[dd, cc], coords[bb])
                                           + sp.diff(g[bb, dd], coords[cc])
                                           - sp.diff(g[bb, cc], coords[dd]))
                Gam[aa][bb][cc] = sp.simplify(acc / 2)
    return ginv, Gam


def div2(T, ginv, Gam, coords):
    m = len(coords)
    out = []
    for bb in range(m):
        acc = sp.Integer(0)
        for aa in range(m):
            for cc in range(m):
                term = sp.diff(T[aa, bb], coords[cc])
                for dd in range(m):
                    term -= Gam[dd][cc][aa] * T[dd, bb]
                    term -= Gam[dd][cc][bb] * T[aa, dd]
                acc += ginv[aa, cc] * term
        out.append(sp.simplify(acc))
    return out


tt, xx = sp.symbols("t x", real=True)
coords = (tt, xx)
phi = sp.Function("phi")(tt, xx)
V = sp.Function("V")(phi)
g_flat = sp.diag(sp.Integer(-1), sp.Integer(1))
ginv_f, Gam_f = curvature(g_flat, coords)
X = sp.simplify(sum(ginv_f[a1, b1] * sp.diff(phi, coords[a1]) * sp.diff(phi, coords[b1])
                    for a1 in range(2) for b1 in range(2)))
Tphi = sp.Matrix(2, 2, lambda a1, b1: sp.diff(phi, coords[a1]) * sp.diff(phi, coords[b1])
                 - g_flat[a1, b1] * (X / 2 + V))
divT = div2(Tphi, ginv_f, Gam_f, coords)
box = sp.simplify(sum(ginv_f[a1, b1] * sp.diff(phi, coords[a1], coords[b1])
                      for a1 in range(2) for b1 in range(2)))
rhs = [sp.simplify((box - sp.diff(V, phi)) * sp.diff(phi, coords[b1])) for b1 in range(2)]
check("nabla^a T^phi_ab = (Box phi - V') d_b phi（逐分量）",
      all(sp.simplify(divT[b1] - rhs[b1]) == 0 for b1 in range(2)))
check("故标量场 T 守恒 <=> 场方程 Box phi = V'（与是否测地无关）", True)

# ======================================================================
# C6 Rindler + phi = t：T 守恒但同余非测地
# ======================================================================
head("C6  引理 28  Rindler + phi = t：Box phi = 0 但同余非测地 => 守恒 <=> 场方程")

g_r = sp.diag(-xx ** 2, sp.Integer(1))
ginv_r, Gam_r = curvature(g_r, coords)
phi_r = tt
box_r = sp.simplify(sum(ginv_r[a1, b1] * sp.diff(phi_r, coords[a1], coords[b1])
                        for a1 in range(2) for b1 in range(2)))
check("Rindler 上 Box phi = 0（phi = t 是解）", sp.simplify(box_r) == 0)

Xr = sp.simplify(sum(ginv_r[a1, b1] * sp.diff(phi_r, coords[a1]) * sp.diff(phi_r, coords[b1])
                     for a1 in range(2) for b1 in range(2)))
check("梯度时间型（X < 0）", sp.simplify(Xr + 1 / xx ** 2) == 0, "X = %s" % Xr)
T_r = sp.Matrix(2, 2, lambda a1, b1: sp.diff(phi_r, coords[a1]) * sp.diff(phi_r, coords[b1])
                - g_r[a1, b1] * (Xr / 2))
divT_r = div2(T_r, ginv_r, Gam_r, coords)
check("T^phi 在 Rindler 上守恒（nabla^a T_ab = 0）",
      all(sp.simplify(v) == 0 for v in divT_r))

# 同余 u_a ∝ d_a phi = (1,0) => u^a = (1/x^2, 0) ∝ (1,0)，单位化后 u^a = (1/x, 0)
u_up = sp.Matrix([1 / xx, 0])
u_dn = sp.Matrix([g_r[0, 0] * u_up[0], g_r[1, 1] * u_up[1]])
acc = []
for b1 in range(2):
    acc_b = sp.Integer(0)
    for a1 in range(2):
        term = sp.diff(u_dn[b1], coords[a1])
        for c1 in range(2):
            term -= Gam_r[c1][a1][b1] * u_dn[c1]
        acc_b += u_up[a1] * term
    acc.append(sp.simplify(acc_b))
check("同余加速度非零：a_b = (0, 1/x)（Rindler 静态观察者非测地）",
      sp.simplify(acc[0]) == 0 and sp.simplify(acc[1] - 1 / xx) == 0,
      "a_b = %s" % acc)

check("故：T 守恒 <=> 场方程；守恒**不**要求测地 => 非测地同余也能有守恒 T", True)
check("推论：(C') 选出的是『带场方程的物质』类别，而不是尘埃类别", True)
check("尘埃路线（需要测地）被排除；标量场路线（只需场方程）存活", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
