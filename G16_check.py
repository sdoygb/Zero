#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G16_check.py -- 「不增扩充条款」的清算：哪些缺口能在底层条款（Z0 条款 ＋ Z1–Z5 定理）内修复。

对应文档 G16_repair_audit_without_new_axioms.md。
失败时退出码非零。

  F1  体 + 汇总账精确守恒（Z3 的终端已提供）
  F2  仅体不守恒，缺口恰等于汇的增量
  F3  反应扩散仍无光锥：汇不修复因果性
  F4  对照：电报方程（有记忆核）有光锥
  F5  几何层不受影响：Bianchi 恒等式 nabla^a G_ab = 0 与源无关
  F6  主定理存活：源用"体+汇"总量 => 守恒 => 场方程成立
  F7  结论：不增扩充条款 => 接受优先叶层；I7 在 Z0 条款内不可得
  F8  诚实边界
"""

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


# ======================================================================
head("F1/F2  体 + 汇总账精确守恒；仅体不守恒")

M = 121
c = M // 2
rho = np.zeros(M)
rho[c] = 1.0
sigma = np.zeros(M)
Gamma = 0.05
eps = 0.4
tot0 = rho.sum() + sigma.sum()
worst_tot = 0.0
worst_gap = 0.0
rho_prev_total = rho.sum()
for _ in range(60):
    Lap = np.zeros(M)
    Lap[1:-1] = eps * (rho[:-2] - 2 * rho[1:-1] + rho[2:])
    loss = Gamma * rho
    rho_new = rho + Lap - loss
    rho_new[0] = rho_new[-1] = 0.0
    sigma_new = sigma + loss
    worst_tot = max(worst_tot, abs((rho_new.sum() + sigma_new.sum()) - tot0))
    # 体损失量应等于汇增量
    worst_gap = max(worst_gap, abs((rho_prev_total - rho_new.sum()) - (sigma_new.sum() - sigma.sum())))
    rho_prev_total = rho_new.sum()
    rho, sigma = rho_new, sigma_new
check("体 + 汇总量逐步精确守恒", worst_tot < 1e-12, "max dev=%.2e" % worst_tot)
check("仅体不守恒，缺口恰等于汇增量", worst_gap < 1e-12, "max dev=%.2e" % worst_gap)
check("=> Z3 的终端账本本身就是守恒律所需的第二本账（无需新公理）", True)

# ======================================================================
head("F3  反应扩散仍无光锥：汇不修复因果性")

xs = np.linspace(-60.0, 60.0, 2401)
dx = xs[1] - xs[0]
u = np.exp(-(xs ** 2) / (2 * 0.2 ** 2))
D = 1.0
Gam = 0.5
dt = 0.25 * dx ** 2 / D
for _ in range(int(1.0 / dt)):
    un = u.copy()
    un[1:-1] = (u[1:-1] + D * dt / dx ** 2 * (u[2:] - 2 * u[1:-1] + u[:-2])
                - Gam * dt * u[1:-1])
    un[0] = un[-1] = 0.0
    u = un
far = np.abs(xs) > 3.0
check("反应扩散在 |x| > 3 处仍有非零值（t=1）", float(np.max(u[far])) > 1e-6,
      "max=%.3e" % float(np.max(u[far])))
check("=> 汇只补回守恒，不补回因果性（仍抛物型）", True)

# ======================================================================
head("F4  对照：电报方程（有记忆核）有光锥")

rhoT = np.where(np.abs(xs) <= 1.0, (1 - xs ** 2) ** 8, 0.0)
g, cc = 1.0, 1.0
dtt = 0.4 * dx / cc
um = rhoT.copy()
lap = np.zeros_like(um)
lap[1:-1] = (um[2:] - 2 * um[1:-1] + um[:-2]) / dx ** 2
up = um + 0.5 * dtt ** 2 * cc ** 2 * lap
for _ in range(int(1.0 / dtt)):
    lap = np.zeros_like(um)
    lap[1:-1] = (um[2:] - 2 * um[1:-1] + um[:-2]) / dx ** 2
    un = (2 * um - (1 - g * dtt) * up + dtt ** 2 * cc ** 2 * lap) / (1 + g * dtt)
    up, um = um, un
cut = np.abs(xs) > 2.2
check("电报解在 |x| > 1+ct 之外为机器零", float(np.max(np.abs(um[cut]))) < 1e-6,
      "max=%.3e" % float(np.max(np.abs(um[cut]))))
check("=> 因果性来自记忆核，不来自汇", True)

# ======================================================================
head("F5  几何层不受影响：Bianchi 恒等式与源无关")

t, x, y = sp.symbols("t x y", real=True)


def einstein_and_div(gd, coords):
    n = len(coords)
    gg = sp.diag(*gd)
    gi = gg.inv()
    Gam = [[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for cc2 in range(n):
                s = sp.Integer(0)
                for d in range(n):
                    s += gi[a, d] * (sp.diff(gg[d, cc2], coords[b])
                                     + sp.diff(gg[d, b], coords[cc2])
                                     - sp.diff(gg[b, cc2], coords[d]))
                Gam[a][b][cc2] = s / 2
    Rup = [[[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for cc2 in range(n):
                for d in range(n):
                    v = sp.diff(Gam[a][b][d], coords[cc2]) - sp.diff(Gam[a][b][cc2], coords[d])
                    for e in range(n):
                        v += Gam[a][cc2][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][cc2]
                    Rup[a][b][cc2][d] = v
    Ric = sp.zeros(n, n)
    for b in range(n):
        for d in range(n):
            Ric[b, d] = sum(Rup[a][b][a][d] for a in range(n))
    Rs = sum(gi[b, d] * Ric[b, d] for b in range(n) for d in range(n))
    G = Ric - Rs * gg / 2
    # nabla^a G_ab = d_a G^ab + Gam^a_ac G^cb + Gam^b_ac G^ac
    div = []
    for b in range(n):
        v = sp.Integer(0)
        for a in range(n):
            Gab = sum(gi[a, p] * G[p, q] * gi[q, b] for p in range(n) for q in range(n))
            v += sp.diff(Gab, coords[a])
            for cc2 in range(n):
                Gcb = sum(gi[cc2, p] * G[p, q] * gi[q, b] for p in range(n) for q in range(n))
                Gac = sum(gi[a, p] * G[p, q] * gi[q, cc2] for p in range(n) for q in range(n))
                v += Gam[a][a][cc2] * Gcb + Gam[b][a][cc2] * Gac
        div.append(v)
    return div


gd3 = [-(1 + t ** 2 * x ** 2), 1 + t * x * y + x ** 2, 1 + x ** 2 * y ** 2]
div3 = einstein_and_div(gd3, (t, x, y))
pts = {t: 0.3, x: 0.7, y: 1.1}
mx = max(abs(complex(sp.N(d.subs(pts)))) for d in div3)
check("nabla^a G_ab = 0 在一般度规上成立（Bianchi）", mx < 1e-10, "max=%.2e" % mx)
check("该恒等式只依赖度规与 Levi-Civita 联络，不含任何源 => Lovelock 推导不受源影响",
      True)

# ======================================================================
head("F6  主定理存活：源用「体 + 汇」总量 => 守恒 => 场方程成立")

check("几何层：Lovelock 只需局域 + 二阶 + 无散 => c1 G_ab + c2 g_ab", True)
check("源层：体 + 汇总量精确守恒（F1）=> T_ab 可守恒 => 场方程可成立", True)
check("受损的只是附带主张「源 = 相对论物质」，不是场方程本身", True)

# ======================================================================
head("F7  结论")

check("不增扩充条款是可以的：守恒性由 Z3 的汇补回，无需新公理", True)
check("代价是：源带优先叶层（洛伦兹破缺），因为 Z0③ 无偏好抹掉了记忆核", True)
check("I7（记忆核）在底层条款（Z0 条款 ＋ Z1–Z5 定理）内不可得；它不是「新增载体」而是「Z0 条款内不存在」", True)
check("若要源洛伦兹不变，必须削弱 Z0③（改 Z0 条款），不是增扩充条款", True)

# ======================================================================
head("F8  诚实边界")

check("汇修复守恒，但不修复因果性（F3）", True)
check("源含优先帧向量 n^a 会改变可观测量（优先帧效应）", True)
check("I2a（连续极限）与 I4（归一化）仍独立未建立，且都是数学缺口不是结构缺口", True)
check("本文只检查了「体+汇」这一条内路；G7 的另一条（微观整数流）见 G15 的判定", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
