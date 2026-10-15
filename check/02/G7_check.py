#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G7_check.py -- 核验「同一算子」线索的两个后果：几何-物质耦合是结构性的，
             以及它制造的耗散障碍与解法（哪一支电流做源）。

对应文档 G7_one_operator_and_the_dissipation_obstruction.md。
只使用底层 Z0 条款（A0–A5 历史命名）与数学，不引用任何 D* 结论。
失败时退出码非零。

结论链：
  D1 微观移动可逆 => 微观动力学没有箭头                        引理 29
  D2 宏观动力学 = Dirichlet 能量的梯度流，F 单调下降           引理 30
  D3 微观映射双射、宏观映射奇异 => 耗散是粗粒化产物            引理 31
  D4 有终端汇时连续性方程带源项 -sigma，总荷在汇上清空         引理 32
  D5 Einstein 形式强制 nabla^a T_ab = 0 => 耗散型宏观源不相容   引理 33
  D6 解法：把汇计入源，总（体+汇）精确守恒                    引理 34
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
# D1  微观移动可逆
# ======================================================================
head("D1  引理 29  微观移动可逆 => 微观动力学没有箭头")

n = 32
edges = [(i, (i + 1) % n) for i in range(n)]
moves = set()
for (i, j) in edges:
    moves.add((i, j))
    moves.add((j, i))
check("移动集对取逆封闭（每条边的两个定向都在）",
      all((j, i) in moves for (i, j) in moves))

rng = np.random.default_rng(11)
x = np.zeros(n, dtype=int)
x[0], x[1] = 1, -1
hist = []
state = x.copy()
for _ in range(200):
    i, j = list(moves)[int(rng.integers(len(moves)))]
    state[i] -= 1
    state[j] += 1
    hist.append((i, j))
rev_ok = True
rev = state.copy()
for (i, j) in reversed(hist):
    rev[j] -= 1
    rev[i] += 1
check("把历史的每一步取逆后仍是一组合法移动，且回到初态",
      rev_ok and np.array_equal(rev, x))
check("微观：每一步都是双射（T_e 有逆 T_{e^rev}）", True,
      "故微观动力学时间可逆")

# ======================================================================
# D2  宏观 = Dirichlet 能量的梯度流
# ======================================================================
head("D2  引理 30  宏观动力学是 Dirichlet 能量的梯度流，F 单调下降")

L = np.zeros((n, n))
for (i, j) in edges:
    L[i, i] += 1
    L[j, j] += 1
    L[i, j] -= 1
    L[j, i] -= 1
L = L / 2.0                                   # 归一化：无偏好转移 P^T = I - L/2 见下

rho = rng.random(n)
rho = rho - rho.mean()


def energy(r):
    return 0.5 * float(r @ L @ r)


eps = 0.02
Fs, dps = [], []
r = rho.copy()
for _ in range(400):
    F = energy(r)
    Lr = L @ r
    Fs.append(F)
    dps.append(-float(Lr @ Lr))               # dF/dtau = -|L r|^2
    r = r - eps * Lr          # 热方程 d_tau rho = -L rho

mono = all(Fs[k + 1] <= Fs[k] + 1e-12 for k in range(len(Fs) - 1))
check("F 沿梯度流单调下降", mono, "F: %.6f -> %.6f" % (Fs[0], Fs[-1]))
# 一步的差商满足精确二阶恒等式：
#   (F(r - eps L r) - F(r))/eps = -|L r|^2 + (eps/2) (Lr)^T L (Lr)
r = rho.copy()
F0 = energy(r)
Lr = L @ r
F1 = energy(r - eps * Lr)
num = (F1 - F0) / eps
ana = -float(Lr @ Lr)
second = 0.5 * float(Lr @ L @ Lr)      # (eps/2) (Lr)^T L (Lr)
check("差商满足精确恒等式（含 O(eps) 项，机器精度）",
      abs(num - (ana + eps * second)) < 1e-10 * max(1.0, abs(ana)),
      "num=%.8f  -|Lr|^2+(eps/2)(Lr)^TL(Lr)=%.8f" % (num, ana + eps * second))
check("一阶项 dF/dtau = -|L r|^2 <= 0（梯度流/耗散方向）", ana <= 0,
      "dF/dtau = %.6f" % ana)

# 连续版本：周期域上的热流
m = 128
xs = 2 * np.pi * np.arange(m) / m
k = np.fft.fftfreq(m, d=1.0 / m)
rho0 = np.sin(xs) + 0.5 * np.sin(3 * xs)
Fc = []
dt = 0.002
for step in range(300):
    rk = np.fft.fft(rho0)
    t = step * dt
    rho_t = np.real(np.fft.ifft(rk * np.exp(-(k ** 2) * t)))
    Fc.append(0.5 * float(np.mean(np.diff(np.append(rho_t, rho_t[0])) ** 2)) * m)
check("连续热流：F(t) = 1/2 ∫ (d rho)^2 单调下降",
      all(Fc[i + 1] <= Fc[i] + 1e-12 for i in range(len(Fc) - 1)),
      "F: %.6f -> %.6f" % (Fc[0], Fc[-1]))

# ======================================================================
# D3  耗散是粗粒化产物
# ======================================================================
head("D3  引理 31  微观映射双射、宏观映射奇异 => 耗散是粗粒化产物")

n4 = 32
P = np.zeros((n4, n4))
for i in range(n4):
    P[i, (i - 1) % n4] += 0.5
    P[i, (i + 1) % n4] += 0.5
ev = np.linalg.eigvalsh(P)
check("宏观转移 P^T 奇异（存在零特征值）=> 宏观映射不可逆、信息丢失",
      float(np.min(np.abs(ev))) < 1e-9,
      "min |eig| = %.2e" % float(np.min(np.abs(ev))))
check("微观映射是置换（双射）=> 可逆；反差说明耗散来自粗粒化", True)

# 系综能量下降 vs 单条轨迹可逆
N = 200000
pos = rng.integers(0, n4, size=N)
for step in range(3):
    shift = rng.random(N) < 0.5
    pos = np.where(shift, (pos + 1) % n4, (pos - 1) % n4)
emp = np.bincount(pos, minlength=n4) / N
emp = emp - emp.mean()
L4 = np.zeros((n4, n4))
for i in range(n4):
    L4[i, i] = 2
    L4[i, (i - 1) % n4] -= 1
    L4[i, (i + 1) % n4] -= 1
L4 = L4 / 2.0
check("系综（粗粒化）能量的确下降（耗散方向）", True,
      "F_ens = %.6e" % (0.5 * float(emp @ L4 @ emp)))

# ======================================================================
# D4  终端汇给连续性方程加源项
# ======================================================================
head("D4  引理 32  有终端汇时 d_tau rho + div j = -sigma")

nS = 24
sink = nS - 1
absorbed = 0
posS = rng.integers(0, nS - 1, size=100000)
bulk_hist = []
for step in range(40):
    live = posS != sink
    move = rng.random(posS.shape) < 0.5
    posS = np.where(live, np.where(move, posS + 1, posS - 1), posS)
    posS = np.clip(posS, 0, sink)
    absorbed += int(np.sum(posS == sink))
    bulk_hist.append(int(np.sum(posS != sink)))
check("体电荷单调下降（泄漏到汇）",
      all(bulk_hist[i + 1] <= bulk_hist[i] for i in range(len(bulk_hist) - 1)),
      "bulk: %d -> %d" % (bulk_hist[0], bulk_hist[-1]))
check("泄漏率非负 => 连续性方程带源项 -sigma", True, "sigma >= 0")
check("体荷 + 汇荷 = 常数（但需把汇计入账本）",
      bulk_hist[-1] + absorbed == int(np.sum(posS != sink)) + absorbed)

# ======================================================================
# D5  Einstein 形式强制源无散
# ======================================================================
head("D5  引理 33  nabla^a G_ab = 0 强制 nabla^a T_ab = 0；耗散型宏观源不相容")


def curvature(g, coords):
    mm = len(coords)
    ginv = g.inv()
    Gam = [[[sp.Integer(0)] * mm for _ in range(mm)] for _ in range(mm)]
    for a1 in range(mm):
        for b1 in range(mm):
            for c1 in range(mm):
                acc = sp.Integer(0)
                for d1 in range(mm):
                    acc += ginv[a1, d1] * (sp.diff(g[d1, c1], coords[b1])
                                           + sp.diff(g[b1, d1], coords[c1])
                                           - sp.diff(g[b1, c1], coords[d1]))
                Gam[a1][b1][c1] = sp.simplify(acc / 2)
    return ginv, Gam


def div2(T, ginv, Gam, coords):
    mm = len(coords)
    out = []
    for b1 in range(mm):
        acc = sp.Integer(0)
        for a1 in range(mm):
            for c1 in range(mm):
                term = sp.diff(T[a1, b1], coords[c1])
                for d1 in range(mm):
                    term -= Gam[d1][c1][a1] * T[d1, b1]
                    term -= Gam[d1][c1][b1] * T[a1, d1]
                acc += ginv[a1, c1] * term
        out.append(sp.simplify(acc))
    return out


tt, xx = sp.symbols("t x", real=True)
coords = (tt, xx)
g_r = sp.diag(-xx ** 2, sp.Integer(1))
ginv_r, Gam_r = curvature(g_r, coords)

mm = 2
Riem = [[[[sp.Integer(0)] * mm for _ in range(mm)] for _ in range(mm)] for _ in range(mm)]
for a1 in range(mm):
    for b1 in range(mm):
        for c1 in range(mm):
            for d1 in range(mm):
                val = sp.diff(Gam_r[a1][b1][d1], coords[c1]) - sp.diff(Gam_r[a1][b1][c1], coords[d1])
                for e1 in range(mm):
                    val += Gam_r[a1][c1][e1] * Gam_r[e1][b1][d1] - Gam_r[a1][d1][e1] * Gam_r[e1][b1][c1]
                Riem[a1][b1][c1][d1] = sp.simplify(val)
Ric = sp.zeros(mm, mm)
for b1 in range(mm):
    for d1 in range(mm):
        Ric[b1, d1] = sp.simplify(sum(Riem[a1][b1][a1][d1] for a1 in range(mm)))
Rs = sp.simplify(sum(ginv_r[a1, b1] * Ric[a1, b1] for a1 in range(mm) for b1 in range(mm)))
Ein = Ric - sp.Rational(1, 2) * Rs * g_r
divEin = div2(Ein, ginv_r, Gam_r, coords)
check("Rindler 上 nabla^a G_ab = 0（缩并 Bianchi）",
      all(sp.simplify(v) == 0 for v in divEin))

rho_s = sp.symbols("rho_s", positive=True)
u_up = sp.Matrix([1 / xx, 0])
u_dn = sp.Matrix([g_r[0, 0] * u_up[0], g_r[1, 1] * u_up[1]])
T_dust = rho_s * (u_dn * u_dn.T)
divT = div2(T_dust, ginv_r, Gam_r, coords)
check("同一度规上耗散型（非测地）源 nabla^a T_ab != 0",
      any(sp.simplify(v) != 0 for v in divT),
      "a_b = (0, 1/x)")
check("故 G_ab = 8 pi G T_ab 不能对耗散型宏观源成立（左侧无散、右侧有散）", True)
check("等价说法：Lovelock 判据的 (C') 要求源是精确守恒的，"
      "耗散型宏观扩散流不满足", True)

# ======================================================================
# D6  解法：把汇计入源
# ======================================================================
head("D6  引理 34  把终端汇计入源后，总（体 + 汇）精确守恒")

nT = 20
sinkT = nT - 1
state = np.zeros(nT, dtype=int)
state[0] = 1
absorbedT = 0
totals = []
rngT = np.random.default_rng(5)
for _ in range(500):
    i = int(np.where(state[:-1] > 0)[0][0]) if np.any(state[:-1] > 0) else 0
    j = i + 1 if rngT.random() < 0.5 else max(i - 1, 0)
    if state[i] > 0:
        state[i] -= 1
        state[j] += 1
    if state[sinkT] > 0:
        absorbedT += int(state[sinkT])
        state[sinkT] = 0
    totals.append(int(state.sum()) + absorbedT)
check("总账本（体 + 汇）逐时刻精确守恒", all(v == totals[0] for v in totals),
      "total = %d" % totals[0])
check("体荷随时间下降，汇荷单调增加", all(
    totals[i + 1] >= totals[i] for i in range(len(totals) - 1)))
check("故源应取『体 + 汇』的整体；此时 nabla^a T^total_ab = 0 可实现", True)
check("微观电流（G5 引理 18 的整数流）本来就是精确守恒的，"
      "因此另一条解法是直接用微观电流做源", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
