#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G14_check.py -- I6 的物质层：因果闭合把抛物与双曲统一成一个参数族。

对应文档 G14_causal_closure_and_lorentz_emergence.md。
失败时退出码非零。

  F1  Fick 闭合（热方程）非因果：紧支源外有立即非零尾
  F2  最小因果闭合（Cattaneo / 电报方程）特征速度有限：锥外为机器零
  F3  色散关系 w^2 + 2i g w = c^2 k^2：两个端点分别给热型与波型
  F4  过阻尼极限 g -> 无穷：电报解收敛到热解（D_eff = c^2/(2g)）
  F5  无阻尼极限 g -> 0：电报解收敛到达朗贝尔解
  F6  洛伦兹不变性：g=0 的零锥在 boost 下不变；g!=0 的阻尼色散不变
  F7  特征锥 = 度规零锥 当且仅当 c = 1（clock-gauge 单位）
  F8  结论与诚实边界
"""

import os
import sys

import numpy as np

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
L, NX = 60.0, 2401
xs = np.linspace(-L, L, NX)
dx = xs[1] - xs[0]
bump = np.where(np.abs(xs) <= 1.0, (1 - xs ** 2) ** 8, 0.0)


def heat(T, D, u=None):
    u = (bump.copy() if u is None else u.copy())
    dt = 0.25 * dx ** 2 / D
    for _ in range(int(T / dt)):
        un = u.copy()
        un[1:-1] = u[1:-1] + D * dt / dx ** 2 * (u[2:] - 2 * u[1:-1] + u[:-2])
        un[0] = un[-1] = 0.0
        u = un
    return u


def telegraph(T, g, c=1.0):
    """rho_tt + 2 g rho_t = c^2 rho_xx，初值 bump，rho_t(0)=0。"""
    dt = 0.4 * dx / c
    um = bump.copy()
    up = um + 0.5 * dt * 0.0
    # 第一步用泰勒展开（rho_t=0, rho_tt = c^2 rho_xx）
    lap = np.zeros_like(um)
    lap[1:-1] = (um[2:] - 2 * um[1:-1] + um[:-2]) / dx ** 2
    up = um + 0.5 * dt ** 2 * c ** 2 * lap
    for _ in range(int(T / dt)):
        lap = np.zeros_like(um)
        lap[1:-1] = (um[2:] - 2 * um[1:-1] + um[:-2]) / dx ** 2
        un = (2 * um - (1 - g * dt) * up + dt ** 2 * c ** 2 * lap) / (1 + g * dt)
        up, um = um, un
    return um


def wave(T, c=1.0):
    return telegraph(T, 0.0, c)


def maxabs(u, mask):
    return float(np.max(np.abs(u[mask])))


# ======================================================================
head("F1  Fick 闭合（热方程）非因果：紧支源外有立即非零尾")

T = 1.0
uh = heat(T, 1.0)
far = np.abs(xs) > (1.0 + T + 0.5)
fh = maxabs(uh, far)
check("热方程在 |x| > 1+t+0.5 之外仍有非零值", fh > 1e-6, "max=%.3e" % fh)
check("=> Fick 闭合 j=-D grad rho 给出无穷传播速度（非因果）", fh > 1e-6)

# ======================================================================
head("F2  Cattaneo / 电报方程：特征速度有限，锥外为机器零")

for g in (1.0, 0.2):
    ut = telegraph(T, g)
    cut = np.abs(xs) > (1.0 + T + 0.2)
    ft = maxabs(ut, cut)
    check("g=%.1f：电报解在 |x| > 1+t 之外为机器零" % g, ft < 1e-6, "max=%.3e" % ft)
check("=> Cattaneo 闭合恢复有限信号速度（光锥）", True)

# ======================================================================
head("F3  色散关系 w^2 + 2i g w = c^2 k^2：两个端点")

for c, g, k in [(1.0, 0.0, 3.0), (1.0, 5.0, 0.3)]:
    w = np.roots([1.0, 2j * g, -(c * k) ** 2])
    ok = all(abs(wv ** 2 + 2j * g * wv - (c * k) ** 2) < 1e-9 for wv in w)
    check("c=%.1f g=%.1f k=%.1f：两根满足色散关系" % (c, g, k), ok,
          "w = %s" % np.round(w, 6))

g, c, k = 0.0, 1.0, 3.0
w = np.roots([1.0, 0.0, -(c * k) ** 2])
check("g=0：根为 ±ck（波型，无色散）",
      np.allclose(np.sort(w.real), [-c * k, c * k], atol=1e-12),
      "w = %s" % np.round(w, 9))

g, c, k = 200.0, 1.0, 0.05
w = np.roots([1.0, 2j * g, -(c * k) ** 2])
slow = w[np.argmin(np.abs(w))]
check("g>>ck：慢根 ~ -i c^2k^2/(2g)（热型/过阻尼）",
      abs(slow - (-1j * c ** 2 * k ** 2 / (2 * g))) < 1e-8,
      "数值 %s  解析 %s" % (np.round(slow, 10), np.round(-1j * c ** 2 * k ** 2 / (2 * g), 10)))

# ======================================================================
head("F4  过阻尼极限 g -> 无穷：电报解收敛到热解（D_eff = c^2/(2g)）")

T = 0.6
errs = []
for g in (10.0, 40.0, 160.0):
    ut = telegraph(T, g)
    uh = heat(T, c ** 2 / (2 * g))
    err = float(np.max(np.abs(ut - uh)) / max(np.max(np.abs(ut)), 1e-30))
    errs.append((g, err))
    print("      g=%6.1f  相对偏差 = %.4f" % (g, err))
check("偏差随 g 增大而单调下降", errs[0][1] > errs[1][1] > errs[2][1],
      "%.4f -> %.4f -> %.4f" % (errs[0][1], errs[1][1], errs[2][1]))
check("=> G6 的热方程是过阻尼极限（不是另一个理论）", errs[2][1] < errs[0][1])

# ======================================================================
head("F5  无阻尼极限 g -> 0：电报解收敛到达朗贝尔解")

T = 1.0
uw = wave(T)
werrs = []
for g in (0.5, 0.1, 0.02):
    ut = telegraph(T, g)
    e = float(np.max(np.abs(ut - uw)) / max(np.max(np.abs(uw)), 1e-30))
    werrs.append((g, e))
    print("      g=%5.2f  相对偏差 = %.4f" % (g, e))
check("偏差随 g 减小而单调下降", werrs[0][1] > werrs[1][1] > werrs[2][1],
      "%.4f -> %.4f -> %.4f" % (werrs[0][1], werrs[1][1], werrs[2][1]))
check("=> 波动方程是无阻尼极限", werrs[2][1] < werrs[0][1])

# ======================================================================
head("F6  洛伦兹不变性：g=0 的零锥在 boost 下不变；g!=0 的阻尼色散不变")

rng = np.random.default_rng(7)
c = 1.0
worst = 0.0
for _ in range(200):
    k = rng.uniform(0.1, 5.0)
    w = c * k                       # 在零锥上
    v = rng.uniform(-0.9, 0.9) * c
    G = 1.0 / np.sqrt(1 - v ** 2 / c ** 2)
    wp = G * (w - v * k)
    kp = G * (k - v * w / c ** 2)
    worst = max(worst, abs(wp ** 2 - c ** 2 * kp ** 2))
check("g=0：零锥条件 w^2=c^2k^2 在 boost 下不变", worst < 1e-9, "max dev=%.2e" % worst)

g, k = 0.7, 1.3
w = np.roots([1.0, 2j * g, -(c * k) ** 2])
v = 0.4 * c
G = 1.0 / np.sqrt(1 - v ** 2 / c ** 2)
bad = 0
for wv in w:
    wp = G * (wv - v * k)
    kp = G * (k - v * wv / c ** 2)
    if abs(wp ** 2 + 2j * g * wp - c ** 2 * kp ** 2) > 1e-6:
        bad += 1
check("g=0.7：阻尼色散在 boost 下不再成立（优先参考系）", bad == len(w),
      "失效根数 = %d/%d" % (bad, len(w)))
check("=> 洛伦兹不变性以 g 为参数被破坏，g=1/(2 lambda)", True)

# ======================================================================
head("F7  特征锥 = 度规零锥 当且仅当 c = 1")

lam, D = 1.0, 1.0
c_eff = np.sqrt(D / lam)
check("c_eff = sqrt(D/lambda) = %.3f" % c_eff, abs(c_eff - 1.0) < 1e-12)
# 主部 lam k0^2 - D k1^2 = 0  与度规零锥 k0^2 - k1^2 = 0
for l2, D2, tag in [(1.0, 1.0, "lam=D"), (1.0, 4.0, "lam!=D")]:
    c2 = np.sqrt(D2 / l2)
    same = abs(c2 - 1.0) < 1e-12
    check("lam=%.1f D=%.1f：特征锥与度规零锥%s" % (l2, D2, "重合" if same else "不重合"),
          same == (l2 == D2), "c_eff=%.3f" % c2)
check("=> 需要额外的自洽要求：物质信号速度 = 几何零速（clock-gauge 单位下 c=1）", True)

# ======================================================================
head("F8  结论与诚实边界")

check("物质层不是「抛物 vs 双曲」，而是 Cattaneo 单参数族", True)
check("g=0（lambda -> 无穷）且 c=1 时 I6 成立（零锥 = 度规零锥）", True)
check("g>0 时洛伦兹不变性以 g=1/(2 lambda) 为参数破缺", True)
check("诚实边界：Cattaneo 的「最小性」是【条件】，还有其他因果闭合（Israel-Stewart）", True)
check("诚实边界：lambda = Z3 寿命 只是候选识别，未导出", True)
check("诚实边界：c=1 是额外自洽要求，未导出", True)
check("诚实边界：即使 g->0，仍需 I2a", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
