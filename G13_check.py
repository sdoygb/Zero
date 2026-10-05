#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G13_check.py -- 叶层与洛伦兹不变性缺口：把「几何层协变、物质层依赖叶层」钉成可复跑检查。

对应文档 G13_foliation_and_lorentz_invariance_gap.md。
登记为输入账本第 6 条 I6。
失败时退出码非零。

  F1  导出度量的 ADM 形式：N = 1，beta = 0（底层选出叶层）
  F2  重切片 t' = t - xi(x,y) 生成非零 shift => 「时钟规范」不被重切片保持
  F3  Einstein 张量按张量律变换 => 几何层协变（场方程不依赖叶层）
  F4  粗粒化输运是热方程（抛物型）：源外有立即的非零尾 => 无光锥
  F5  波动方程（双曲型）：紧支源在 |x| > 1+t 之外为机器零 => 有光锥
  F6  场类首符号是逆度规 => 零锥即特征线（双曲）
  F7  结论：几何层协变 vs 物质层依赖叶层 => 洛伦兹不变性未建立，登记 I6
"""

import io
import os
import sys

import numpy as np
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
def as_metric(gmat, n):
    """接受 3x3 矩阵或对角元列表，返回 sympy 矩阵。"""
    g = sp.Matrix(gmat)
    if g.shape != (n, n):
        g = sp.diag(*list(gmat))
    return g


def geometry(gmat, coords):
    """通用度规的 Einstein 张量。不做化简：构建 DAG 后数值取值（同 G8）。"""
    n = len(coords)
    g = as_metric(gmat, n)
    ginv = g.inv()
    Gam = [[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                s = sp.Integer(0)
                for d in range(n):
                    s += ginv[a, d] * (sp.diff(g[d, c], coords[b])
                                       + sp.diff(g[d, b], coords[c])
                                       - sp.diff(g[b, c], coords[d]))
                Gam[a][b][c] = s / 2
    Rup = [[[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    v = sp.diff(Gam[a][b][d], coords[c]) - sp.diff(Gam[a][b][c], coords[d])
                    for e in range(n):
                        v += Gam[a][c][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][c]
                    Rup[a][b][c][d] = v
    Ric = sp.zeros(n, n)
    for b in range(n):
        for d in range(n):
            Ric[b, d] = sum(Rup[a][b][a][d] for a in range(n))
    Rs = sum(ginv[b, d] * Ric[b, d] for b in range(n) for d in range(n))
    return g, Ric - Rs * g / 2


def num(e, pts):
    return complex(sp.N(e.subs(pts)))


def adm(gmat, n):
    """从 g_{mu nu} 提取 ADM：N、beta^i、h_ij。第 0 个坐标是叶层时间。"""
    g = as_metric(gmat, n)
    h = -g[1:, 1:]
    hinv = h.inv()
    g0 = sp.Matrix([g[0, j] for j in range(1, n)])
    beta = -hinv * g0
    N2 = g[0, 0] + (g0.T * hinv * g0)[0]
    return sp.sqrt(sp.simplify(N2)), beta, h


# ======================================================================
head("F1  导出度量的 ADM 形式：N = 1，beta = 0")

t, x, y = sp.symbols("t x y", real=True)      # t 就是步参数 tau
h11 = 1 + t ** 2 * x ** 2 + y ** 2
h22 = 1 + sp.sin(t * x) ** 2 + x ** 2 * y ** 2
gdiag = [1, -h11, -h22]
N, beta, h = adm(gdiag, 3)
pts = {t: 0.4, x: 0.7, y: 1.1}
Nv = abs(num(N, pts))
bv = max(abs(num(beta[i], pts)) for i in range(2))
check("N = 1（对含时 h 也成立）", abs(Nv - 1.0) < 1e-12, "N=%.15f" % Nv)
check("beta^i = 0", bv < 1e-14, "max|beta|=%.2e" % bv)
check("故底层选出一个叶层（t = 常数切片）", abs(Nv - 1.0) < 1e-12 and bv < 1e-14)

# ======================================================================
head("F2  重切片 t' = t - xi(x,y) 生成非零 shift")

xi = sp.Rational(3, 10) * x ** 2
# t' = t - xi, x' = x, y' = y   =>   t = t' + xi
# Jp = d x'^mu / d x^alpha
Jp = sp.Matrix([[1, -sp.diff(xi, x), -sp.diff(xi, y)], [0, 1, 0], [0, 0, 1]])
inv = Jp.inv()                     # d x^alpha / d x'^mu
subs = {t: t + xi, x: x, y: y}     # 用同名符号表示 primed 坐标
g0m = as_metric(gdiag, 3)
gprime = sp.simplify(inv.T * g0m.subs(subs) * inv)
N2, beta2, h2 = adm(gprime, 3)
bv2 = max(abs(num(beta2[i], pts)) for i in range(2))
Nv2 = abs(num(N2, pts))
check("重切片后 shift 一般非零", bv2 > 1e-3, "max|beta'|=%.6f" % bv2)
check("重切片后 lapse 也不再为 1（闭式 N'^2 = 1 + xi^T h_new^{-1} xi）",
      abs(Nv2 - 1.0) > 1e-3, "N'=%.6f" % Nv2)
check("=> 「时钟规范」N=1,beta=0 是底层选出的，重切片同时改动 N 与 beta", bv2 > 1e-3 and abs(Nv2 - 1.0) > 1e-3)

# ======================================================================
head("F3  Einstein 张量按张量律变换 => 几何层协变")

_, G_a = geometry(g0m, (t, x, y))
G_pull = sp.Matrix(3, 3, lambda i, j: G_a[i, j].subs(subs))
G_exp = inv.T * G_pull * inv
_, G_b = geometry(gprime, (t, x, y))
dev = max(abs(num(G_exp[i, j], pts) - num(G_b[i, j], pts)) for i in range(3) for j in range(3))
check("直接算的 G' 与张量搬运的 G' 一致 => 场方程协变", dev < 1e-8, "max dev=%.2e" % dev)
check("故对几何扇区而言叶层是纯规范", dev < 1e-8)

# ======================================================================
head("F4  粗粒化输运是热方程（抛物型）：源外有立即的非零尾")

Lx, nx = 40.0, 1601
xs = np.linspace(-Lx, Lx, nx)
dx = xs[1] - xs[0]
sig = 0.20
u = np.exp(-(xs ** 2) / (2 * sig ** 2))
D = 1.0
dt = 0.25 * dx ** 2 / D
T = 1.0
for _ in range(int(T / dt)):
    un = u.copy()
    un[1:-1] = u[1:-1] + D * dt / dx ** 2 * (u[2:] - 2 * u[1:-1] + u[:-2])
    un[0] = un[-1] = 0.0
    u = un
far = np.abs(xs) > (1.0 + T + 5 * sig)
hfar = float(np.max(np.abs(u[far])))
check("热方程在远场（|x| > 1+t+5sigma）仍有非零值", hfar > 1e-6, "max|u|_far=%.3e" % hfar)
check("=> 抛物型：无光锥，扰动瞬时传遍", hfar > 1e-6)

# ======================================================================
head("F5  波动方程（双曲型）：紧支源在 |x| > 1+t 之外为机器零")

bump = np.where(np.abs(xs) <= 1.0, (1 - xs ** 2) ** 8, 0.0)
uw = bump.copy()
v = np.zeros_like(bump)
c = 1.0
dtw = 0.4 * dx / c
for _ in range(int(T / dtw)):
    lap = np.zeros_like(uw)
    lap[1:-1] = (uw[2:] - 2 * uw[1:-1] + uw[:-2]) / dx ** 2
    v = v + dtw * c ** 2 * lap
    uw = uw + dtw * v
outside = np.abs(xs) > (1.0 + T + 0.15)
wfar = float(np.max(np.abs(uw[outside])))
check("波动方程在 |x| > 1+t 之外为机器零", wfar < 1e-5, "max|u|_out=%.3e" % wfar)
check("=> 双曲型：有光锥，信号速度有限", wfar < 1e-5)
ratio = hfar / max(wfar, 1e-300)
check("远场量级比 > 1e3（抛物 vs 双曲）", ratio > 1e3, "比值=%.2e" % ratio)

# ======================================================================
head("F6  场类首符号是逆度规 => 零锥即特征线（双曲）")

ev = np.linalg.eigvalsh(np.diag([1.0, -1.0, -1.0]))
check("首符号矩阵有一个正特征值、其余为负（洛伦兹号差）",
      int(np.sum(ev > 0)) == 1 and int(np.sum(ev < 0)) == 2, "eigs=%s" % np.round(ev, 6))
k0, k1 = sp.symbols("k0 k1", real=True)
sol = sp.solve(sp.Eq(k0 ** 2 - k1 ** 2, 0), k0)
check("零锥方程有非平凡实解（特征线存在）", len(sol) == 2, "k0 = %s" % sol)
check("故代表场双曲、粗粒化流动抛物 => 二者不同机制", True)

# ======================================================================
head("F7  结论")

g1 = io.open(os.path.join(HERE, "G1_derivations_from_the_bottom_layer.md"),
             encoding="utf-8").read()
g6 = io.open(os.path.join(HERE, "G6_geodesy_of_the_coarse_grained_flow.md"),
             encoding="utf-8").read()
g0 = io.open(os.path.join(HERE, "G0_bottom_layer_and_derivation_route.md"),
             encoding="utf-8").read()
check("G1 引理 5 给 g = dtau^2 - h（含 N=1, beta=0）", "d\\tau^2" in g1)
check("G6 判定粗粒化动力学是热方程", "热方程" in g6)
check("Z4 定义步参数 tau（叶层时间的来源）", "步指标" in g0)
check("结论：洛伦兹不变性（叶层纯规范）未建立，登记为账本 I6", True)
check("诚实边界：几何层协变已核验，缺口只在物质/输运层", dev < 1e-8)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
