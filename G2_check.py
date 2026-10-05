#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G2_check.py -- 核验「局域连续极限」(输入 I2) 的分解与两个障碍。

对应文档 G2_local_continuum_limit.md。
只使用底层 Z0 条款（A0–A5 历史命名）与数学，不引用任何 D* 结论。
失败时退出码非零。

两组实验：
  网格差分型能量  -> 尺度律与极限值（引理 9）
  环面有效电阻    -> 距离型几何量的尺度障碍（引理 10）
环面用 Fourier 精确公式，无边界效应。
"""

import math
import sys

import numpy as np

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
# 一、网格差分型能量（引理 9）
# ======================================================================

def test_function(n, d):
    """f(x) = prod_k sin(pi x_k)，节点 x_k = (i+1)/(n+1)，边界取 0。"""
    h = 1.0 / (n + 1)
    axes = [(np.arange(n) + 1) * h for _ in range(d)]
    grids = np.meshgrid(*axes, indexing="ij")
    f = np.ones([n] * d)
    for g in grids:
        f = f * np.sin(np.pi * g)
    return f


def dirichlet_energy(n, d, K):
    """E = sum_{edges} K (f_i - f_j)^2，含与 f=0 边界的连线。"""
    f = test_function(n, d)
    total = 0.0
    for k in range(d):
        total += float(np.sum(np.diff(f, axis=k) ** 2))
        sl = [slice(None)] * d
        sl[k] = 0
        total += float(np.sum(f[tuple(sl)] ** 2))
        sl[k] = -1
        total += float(np.sum(f[tuple(sl)] ** 2))
    return K * total


def continuum_energy(d, samples=300):
    """int |grad f|^2 的数值积分与解析值 d pi^2 2^{-d}。"""
    xs = (np.arange(samples) + 0.5) / samples
    w = 1.0 / samples
    grids = np.meshgrid(*([xs] * d), indexing="ij")
    total = 0.0
    for k in range(d):
        term = np.pi ** 2 * np.cos(np.pi * grids[k]) ** 2
        for j in range(d):
            if j != k:
                term = term * np.sin(np.pi * grids[j]) ** 2
        total += float(np.sum(term)) * (w ** d)
    return total, d * math.pi ** 2 * (0.5 ** d)


def slope_sequence(d, alpha, ns):
    """E(a) ~ a^{alpha+2-d}(1 + c a)，故 s(a) := log2(E(2a)/E(a)) = m + k a。

    用三点 Richardson 外推 m = 2 s(a/2) - s(a) 消去 O(a) 修正。
    """
    a = [1.0 / (n + 1) for n in ns]
    E = [dirichlet_energy(n, d, a[i] ** alpha) for i, n in enumerate(ns)]
    s1 = math.log2(E[1] / E[0])
    s2 = math.log2(E[2] / E[1])
    return 2.0 * s2 - s1


head("D1  引理 9  尺度律：K=1 时 E 按 a^{d-2} 漂移；K=a^{d-2} 时收敛")

for d in (1, 2, 3):
    ns = (32, 64, 128)
    m_fixed = slope_sequence(d, 0.0, ns)
    m_scaled = slope_sequence(d, float(d - 2), ns)
    check(
        "d=%d  K=1：外推斜率 = d-2 = %d" % (d, d - 2),
        abs(m_fixed - (d - 2)) < 0.02,
        "extrapolated=%.5f" % m_fixed,
    )
    check(
        "d=%d  K=a^{d-2}：外推斜率 = 0（收敛）" % d,
        abs(m_scaled) < 0.02,
        "extrapolated=%.5f" % m_scaled,
    )

head("D2  引理 9 续  外推斜率 = d-2-alpha，唯一零点 alpha = d-2")

for d in (1, 2, 3):
    alphas = [-1.0, 0.0, 0.5, 1.0, 2.0] if d == 1 else (
        [-1.0, 0.0, 1.0, 2.0, 3.0] if d == 2 else [0.0, 0.5, 1.0, 1.5, 2.0]
    )
    rows = []
    for alpha in alphas:
        m = slope_sequence(d, alpha, (32, 64, 128))
        pred = d - 2 - alpha
        rows.append((alpha, m))
        check(
            "d=%d  alpha=%.1f：外推斜率 = d-2-alpha = %.3f" % (d, alpha, pred),
            abs(m - pred) < 0.02,
            "measured=%.5f" % m,
        )
    best = min(rows, key=lambda r: abs(r[1]))
    check(
        "d=%d  唯一使斜率消失的 alpha = d-2 = %d" % (d, d - 2),
        abs(best[0] - (d - 2)) < 1e-9 and abs(best[1]) < 0.02,
        "best alpha=%.1f" % best[0],
    )

head("D3  引理 9 续  正确标度下极限值等于局部型 int |grad f|^2")

for d in (1, 2, 3):
    ref, exact = continuum_energy(d)
    ns = (48, 96, 192) if d < 3 else (24, 48, 96)
    a = [1.0 / (n + 1) for n in ns]
    ratios = [dirichlet_energy(n, d, a[i] ** (d - 2)) / exact for i, n in enumerate(ns)]
    coef = np.polyfit(a, ratios, 1)
    intercept = float(coef[1])
    check(
        "d=%d  E(K=a^{d-2}) / int|grad f|^2 外推到 a->0 得 1" % d,
        abs(intercept - 1.0) < 0.01,
        "ratios=%s  intercept=%.5f" % ([round(r, 5) for r in ratios], intercept),
    )
    check(
        "d=%d  连续参照的解析值与数值积分一致" % d,
        abs(ref - exact) / exact < 1e-3,
        "quad=%.6f exact=%.6f" % (ref, exact),
    )

# ======================================================================
# 二、环面有效电阻（引理 10）
# ======================================================================

def torus_resistance(d, n, rvecs):
    """d 维 n^d 环面上的有效电阻，Fourier 精确公式。

    lambda_k = sum_i 2(1 - cos(2 pi k_i / n))，k != 0；
    R(r) = 2( S0/N - Re(ifftn(H))[r] )，H = 1/lambda（k=0 处置 0）。
    """
    lam = np.zeros([n] * d)
    for i in range(d):
        shape = [1] * d
        shape[i] = n
        ki = np.arange(n).reshape(shape)
        lam = lam + 2.0 * (1.0 - np.cos(2.0 * np.pi * ki / n))
    H = np.zeros([n] * d)
    mask = lam > 1e-12
    H[mask] = 1.0 / lam[mask]
    S0 = float(H.sum())
    N = float(n ** d)
    G = np.real(np.fft.ifftn(H))
    return [2.0 * (S0 / N - float(G[tuple(rv)])) for rv in rvecs]


head("D0  校验  环面上相邻节点的有效电阻收敛到无限格点值 R_adj = 1/d")

for d, n in ((1, 256), (2, 128), (3, 64)):
    rv = [0] * d
    rv[0] = 1
    R = torus_resistance(d, n, [tuple(rv)])[0]
    check("d=%d  R_adj -> 1/d = %.4f（有限尺寸修正 O(1/n)）" % (d, 1.0 / d),
          abs(R - 1.0 / d) < 5e-3, "measured=%.5f" % R)

head("D4a  引理 10  d=1：R(r) = r(n-r)/n 精确（线性，唯一的度量型维度）")

n = 256
for r in (1, 10, 40):
    got = torus_resistance(1, n, [(r,)])[0]
    pred = r * (n - r) / n
    check("d=1  n=%d r=%d：R = r(n-r)/n" % (n, r),
          abs(got - pred) < 1e-9, "R=%.10f pred=%.10f" % (got, pred))

head("D4b  引理 10  d=2：R(r) ~ (1/pi) ln r（对数，sqrt(R) 次线性）")

n2 = 1024
rs = [16, 32, 64, 128, 256]
vals = torus_resistance(2, n2, [(r, 0) for r in rs])
slope01 = (vals[1] - vals[0]) / math.log(2.0)
check("d=2  相邻倍 r 的 ln 系数 = 1/pi = %.5f" % (1.0 / math.pi),
      abs(slope01 - 1.0 / math.pi) < 2e-3,
      "measured=%.5f" % slope01)
check("d=2  R(r)/r -> 0（次线性增长，不是距离）",
      vals[-1] / rs[-1] < 0.02,
      "R(256)/256=%.5f" % (vals[-1] / rs[-1]))
check("d=2  R(r) 随 r 单调增长且无上界迹象",
      all(vals[i + 1] > vals[i] for i in range(len(vals) - 1)),
      "R=%s" % [round(v, 4) for v in vals])

head("D4c  引理 10  d=3：R(r) 饱和到有限常数（sqrt(R) 有界）")

n3 = 64
rs3 = [1, 2, 4, 8, 16]
vals3 = torus_resistance(3, n3, [(r, 0, 0) for r in rs3])
diffs3 = [vals3[i + 1] - vals3[i] for i in range(len(vals3) - 1)]
check("d=3  R(r) 单调递增但增量递减（饱和）",
      all(vals3[i + 1] > vals3[i] for i in range(len(vals3) - 1))
      and all(diffs3[i + 1] < diffs3[i] for i in range(len(diffs3) - 1)),
      "R=%s  diffs=%s" % ([round(v, 5) for v in vals3], [round(x, 5) for x in diffs3]))
check("d=3  R(r) 在 r<=16 上有界且增量递减 => 收敛到有限常数",
      vals3[-1] < 0.5 and all(diffs3[i + 1] < diffs3[i] for i in range(len(diffs3) - 1)),
      "R(16)=%.5f  diff(8->16)=%.5f" % (vals3[-1], diffs3[-1]))
n_indep = [torus_resistance(3, nn, [(4, 0, 0)])[0] for nn in (16, 32, 64)]
check("d=3  R(4) 与系统尺寸无关（收敛到无限格点值）",
      max(n_indep) - min(n_indep) < 5e-3,
      "R(4) at n=16,32,64: %s" % [round(v, 5) for v in n_indep])

check("三种尺度的定性差别：d=1 线性、d=2 对数、d=3 饱和",
      abs(torus_resistance(1, n, [(32,)])[0] / 32.0 - (n - 32) / n) < 1e-9
      and vals[-1] / rs[-1] < 0.02
      and vals3[-1] < 0.5,
      "d1 近似线性；d2 R/r=%.5f；d3 R(16)=%.4f" % (vals[-1] / rs[-1], vals3[-1]))

head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
