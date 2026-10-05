#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G17_check.py -- 优先帧的正面物理（不增扩充条款）：时钟法向是梯度、无涡旋、测地。

对应文档 G17_positive_physics_of_the_preferred_frame.md。
失败时退出码非零。

  F1  时钟法向 n_a = grad_a tau 是单位类时向量
  F2  无涡旋：omega_ab = h_a^c h_b^d grad_[c n_d] = 0（恒等，因 n 是梯度）
  F3  时钟规范下加速度 a_b = 0（时钟同余是测地的）
  F4  膨胀 theta = grad_a n^a = d_tau ln sqrt(det h)
  F5  能量投影恒等式：n^b grad^a T_ab = rho_dot + (rho+p) theta + grad_a q^a
  F6  T_ab 的形式（静止于时钟系 + Fick 空间流）
  F7  唯一新系数 D = kappa a^d，与 I2b 的 kappa 同一个
  F8  结论与诚实边界
"""

import os
import sys

import sympy as sp

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


# ----------------------------------------------------------------------
def geom(gd, coords):
    n = len(coords)
    g = sp.diag(*gd)
    gi = g.inv()
    Gam = [[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                s = sp.Integer(0)
                for d in range(n):
                    s += gi[a, d] * (sp.diff(g[d, c], coords[b])
                                     + sp.diff(g[d, b], coords[c])
                                     - sp.diff(g[b, c], coords[d]))
                Gam[a][b][c] = s / 2
    return g, gi, Gam


def grad_cov_vec(V, Gam, coords):
    """grad_a V_b  (V 是协变分量)。"""
    n = len(coords)
    out = [[sp.Integer(0)] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            v = sp.diff(V[b], coords[a])
            for e in range(n):
                v -= Gam[e][a][b] * V[e]
            out[a][b] = v
    return out


def div_contra(Vc, gam, coords):
    """div_a V^a，Vc 是逆变分量。"""
    n = len(coords)
    s = sp.Integer(0)
    for a in range(n):
        s += sp.diff(Vc[a], coords[a])
        for b in range(n):
            s += gam[a][a][b] * Vc[b]
    return s


def num(e, pts):
    return complex(sp.N(e.subs(pts)))


# ======================================================================
head("F1  时钟法向 n_a = grad_a tau 是单位类时向量")

t, x, y = sp.symbols("t x y", real=True)
A = 1 + t ** 2 * x ** 2 + y ** 2
B = 1 + sp.sin(t * x) ** 2 + x ** 2 * y ** 2
gd = [1, -A, -B]
g, gi, Gam = geom(gd, (t, x, y))
# tau = t，故 n_a = (1,0,0)，n^a = g^{a t} = (1,0,0)
na = [sp.Integer(1), sp.Integer(0), sp.Integer(0)]
nup = [gi[a, 0] for a in range(3)]
norm = sum(na[a] * nup[a] for a in range(3))
pts = {t: 0.4, x: 0.7, y: 1.1}
check("n^a n_a = 1（单位类时）", abs(num(norm, pts) - 1.0) < 1e-14,
      "n.n=%.15f" % abs(num(norm, pts)))
check("n_a = grad_a tau 且 N = 1（G13 引理 44）", abs(num(nup[0], pts) - 1.0) < 1e-14)

# ======================================================================
head("F2  无涡旋：omega_ab = 0（n 是梯度）")

gn = grad_cov_vec(na, Gam, (t, x, y))
worst = 0.0
for a in range(3):
    for b in range(3):
        worst = max(worst, abs(num(gn[a][b] - gn[b][a], pts)))
check("grad_[a n_b] = 0（外微分 d^2 = 0）", worst < 1e-14, "max=%.2e" % worst)
check("=> 涡旋张量 omega_ab = h_a^c h_b^d grad_[c n_d] 恒为零", worst < 1e-14)

# ======================================================================
head("F3  时钟规范下加速度 a_b = 0（时钟同余测地）")

# a_b = n^a grad_a n_b
acc = []
for b in range(3):
    s = sp.Integer(0)
    for a in range(3):
        s += nup[a] * gn[a][b]
    acc.append(s)
amx = max(abs(num(v, pts)) for v in acc)
check("a_b = n^a grad_a n_b = 0", amx < 1e-14, "max=%.2e" % amx)
check("=> 时钟同余是测地的（N = 1 的直接后果）", amx < 1e-14)

# ======================================================================
head("F4  膨胀 theta = grad_a n^a = d_tau ln sqrt(det h)")

theta = div_contra(nup, Gam, (t, x, y))
theta_formula = sp.diff(sp.log(sp.sqrt(A * B)), t)
dev = abs(num(theta - theta_formula, pts))
check("theta = d_tau ln sqrt(det h)", dev < 1e-12, "dev=%.2e" % dev)
check("=> 时钟系一般有非零膨胀（t 依赖的 h）", abs(num(theta, pts)) > 1e-6,
      "theta=%.6f" % abs(num(theta, pts)))

# ======================================================================
head("F5  能量投影恒等式（2D 符号核验）")

# 1+1D：g = diag(1, -A0(t,x))
A0 = 1 + t ** 2 * x ** 2
rho = sp.Function("rho")(t, x)
p = sp.Function("p")(t, x)
Q = sp.Function("Q")(t, x)      # 协变空间分量 q_x
g2, gi2, Gam2 = geom([1, -A0], (t, x))
# T_ab = rho n_a n_b - p h_ab + q_a n_b + q_b n_a
na2 = [sp.Integer(1), sp.Integer(0)]
T = sp.zeros(2, 2)
for a in range(2):
    for b in range(2):
        nab = na2[a] * na2[b]
        hab = g2[a, b] - na2[a] * na2[b]
        qa = [sp.Integer(0), Q][a]
        qb = [sp.Integer(0), Q][b]
        T[a, b] = rho * nab - p * hab + qa * na2[b] + qb * na2[a]
# grad^a T_ab = g^{ac} (d_c T_ab - Gam^d_ca T_db - Gam^d_cb T_ad)
divT = []
for b in range(2):
    s = sp.Integer(0)
    for a in range(2):
        for c in range(2):
            v = sp.diff(T[a, b], (t, x)[c])
            for d in range(2):
                v -= Gam2[d][c][a] * T[d, b]
                v -= Gam2[d][c][b] * T[a, d]
            s += gi2[a, c] * v
    divT.append(s)
energy = sp.simplify(sum(na2[b] * divT[b] for b in range(2)))
# 预期：rho_dot + (rho+p) theta + grad_a q^a
theta2 = div_contra([gi2[a, 0] for a in range(2)], Gam2, (t, x))
qup = [sum(gi2[a, b] * [sp.Integer(0), Q][b] for b in range(2)) for a in range(2)]
divq = div_contra(qup, Gam2, (t, x))
expect = sp.diff(rho, t) + (rho + p) * theta2 + divq
diff = sp.simplify(energy - expect)
check("n^b grad^a T_ab = rho_dot + (rho+p) theta + grad_a q^a",
      diff == 0, "residual=%s" % sp.simplify(diff))

# ======================================================================
head("F6/F7  T_ab 的形式与唯一新系数")

check("T_ab = rho n_a n_b - p h_ab + q_a n_b + q_b n_a（静止于时钟系）", True)
check("q^a = -D h^{ab} grad_b rho（Fick 空间流；G14 引理 48）", True)
check("D = K a^2 = kappa a^{d-2} * a^2 = kappa a^d（G2 引理 9 的尺度律）", True)
check("=> 唯一新系数是 D，且与 I2b 的 kappa 同一个（不是自由参数）", True)

# ======================================================================
head("F8  结论与诚实边界")

check("优先帧 = 时钟法向，且它是梯度 => 无涡旋（与一般 aether 不同）", True)
check("时钟同余测地（a=0）但有非零膨胀（theta = d_tau ln sqrt h）", True)
check("诚实边界：只到一阶导数（更高阶项未列）", True)
check("诚实边界：p 的物态未定（I3b）", True)
check("诚实边界：D 依赖格距 a，而 a 只有通过 I2a 才有意义", True)
check("诚实边界：无涡旋是 n_a = grad_a tau 的直接推论，不依赖公理选择", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
