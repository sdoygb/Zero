#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
a2_bridge_check.py —— 热核系数 $a_2$ 的 $R$ 系数 $=1/6$ 的独立核验（Z142 §5.3 的作用量桥）

【核的是什么】
Z142 §5.3 的作用量桥用的是标准的 Seeley–DeWitt／Gilkey 展开：对 Laplace 型算子

    theta(t) = tr e^{-tP} ~ (4 pi t)^{-n/2} ( a_0 + a_1 t + a_2 t^2 + ... ),
    a_0 = Vol = ∫ sqrt(g),        a_1 = (1/6)∫ sqrt(g) R − ∫ sqrt(g) E .

本脚本**不用**该公式的文本，而是直接从**谱**算 heat trace，拟合出 $a_1/a_0$，再与 $R/6$ 比较——
故它核的是**常数 $1/6$ 本身**，不是重复抄公式。

【两个对象】
  A. 半径 $a$ 的球面 $S^2$：本征值 $l(l+1)/a^2$、重数 $2l+1$。$R = 2/a^2$，$Vol = 4\pi a^2$，
     $\int \sqrt g R = 8\pi$（Gauss–Bonnet，与 $a$ 无关）。
     ⟹ 预测 $\dfrac{a_1}{a_0} = \dfrac{(1/6)\cdot 8\pi}{4\pi a^2} = \dfrac{1}{3a^2} = \dfrac{R}{6}$。
  B. 平坦环面 $T^2$（边长 $L$）：$R=0$ ⟹ 预测 $a_1/a_0 = 0$（数值上按指数小）。
     这一条是**对照**：说明拟合出的斜率确实在测 $R$，而不是某个普适常数。

【口径】退出码 0 ＝ 全部判据成立（PASS）；1 ＝ 有判据不成立。本脚本**没有**"FAIL＝要报的判定"一类语义。
"""

import sys

import numpy as np

TOL = 2e-3          # 斜率比的相对容差（拟合截断误差量级）

results = []


def ck(cond, label, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {label}" + (f" —— {detail}" if detail else ""))
    results.append((bool(cond), label, detail))
    return cond


def theta_sphere(t, a, lmax=4000):
    """S^2(a) 的 heat trace：sum_l (2l+1) exp(-l(l+1)t/a^2)。"""
    l = np.arange(lmax + 1)
    x = l * (l + 1) * t / (a * a)
    w = np.where(x < 700.0, (2 * l + 1) * np.exp(-np.clip(x, 0, 700.0)), 0.0)
    return float(w.sum())


def theta_torus(t, L, mmax=400):
    """平坦环面 T^2（边长 L）的 heat trace：sum_{m,n} exp(-(2 pi/L)^2 (m^2+n^2) t)。"""
    m = np.arange(-mmax, mmax + 1)
    X, Y = np.meshgrid(m, m, indexing="ij")
    x = (2 * np.pi / L) ** 2 * (X ** 2 + Y ** 2) * t
    return float(np.where(x < 700.0, np.exp(-np.clip(x, 0, 700.0)), 0.0).sum())


def fit_a1_over_a0(trace_and_scale, ts, deg=3):
    """theta(t) -> F(t) = (4 pi t)^{n/2} theta(t) = a_0 + a_1 t + ...；返回 a_1/a_0。"""
    F = np.array([(4 * np.pi * t) ** (trace_and_scale[2] / 2.0) * trace_and_scale[0](t)
                  for t in ts])
    co = np.polyfit(ts, F, deg)[::-1]          # co[0] + co[1] t + ...
    return co[0], co[1], co[1] / co[0]


def head(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


# ────────────────────────────────────────────────────────── A：球面
head("A｜$S^2(a)$：$a_1/a_0$ 是否等于 $R/6=1/(3a^2)$")

TS = np.linspace(1e-3, 6e-3, 26)
for a in (1.0, 2.0, 3.0):
    a0, a1, ratio = fit_a1_over_a0((lambda t, a=a: theta_sphere(t, a), None, 2), TS)
    pred = 1.0 / (3 * a * a)
    rel = abs(ratio - pred) / pred
    line = (f"a={a:g}：拟合 a0={a0:.8f}（理论 4πa²={4*np.pi*a*a:.8f}）、"
            f"a1={a1:.8f}（理论 (1/6)∫R={4*np.pi/3:.8f}）\n"
            f"        a1/a0={ratio:.8f}   R/6={pred:.8f}   相对差={rel:.2e}")
    print("  " + line)
    ck(rel < TOL, f"A 球面 a={a:g}：a1/a0 = R/6", f"相对差 {rel:.2e}")

# ────────────────────────────────────────────────────────── B：平坦环面（对照）
head("B｜平坦 $T^2$（$R=0$，对照）：$a_1/a_0$ 是否趋于 0")

for L in (4.0, 8.0):
    a0, a1, ratio = fit_a1_over_a0((lambda t, L=L: theta_torus(t, L), None, 2), TS)
    rel = abs(ratio - 0.0)
    print(f"  L={L:g}：拟合 a0={a0:.8f}（理论 L²={L*L:.8f}）、a1={a1:.3e}、"
          f"a1/a0={ratio:.3e}")
    ck(rel < 1e-4, f"B 平坦环面 L={L:g}：R=0 ⟹ a1/a0≈0", f"|a1/a0|={rel:.2e}")

# ────────────────────────────────────────────────────────── 汇总
npass = sum(1 for ok, _, _ in results if ok)
print()
print("=" * 78)
print(f"合计 {len(results)} 条：PASS {npass} / FAIL {len(results) - npass}")
print("=" * 78)
print("结论：球面上 a1/a0 与 R/6 相符、平坦环面上为 0 ⟹ 热核第二系数中的 R 系数确是 1/6，")
print("      且它不是普适常数（随 R 变）。这就是 Z142 §5.3 作用量桥所引用的标准常数。")

sys.exit(0 if npass == len(results) else 1)
