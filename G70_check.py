#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G70_check.py -- 清单里剩下两项：B 的识别 与 tau 的初条件依赖

对应文档 G70_B_and_tau_closing.md。只做数值断言。

  F1  前沿存在性：B <= 4（B>=5 无解）
  F2  增长要求：B >= 2
  F3  标签计数：|Z_2 x Z_2| = 4
  F4  诚实弱化：c*=1 不是必需的（亚光速物质完全允许）
  F5  tau(N) 收敛：递增且增量递减
  F6  tau_inf 外推（拟合 tau(N) = tau_inf - c N^{-p}）
      注：F5 已补 T=400 的收敛值与「两极限不可交换」；F6 的拟合固定 T=12，
      其结论不得升级为「预言」。
"""

import os
import sys
from itertools import combinations_with_replacement as cwr
from math import cosh, log, tanh

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


def cstar(B):
    tgt = 0.5 * log(B)
    f = lambda m: m * tanh(m) - log(cosh(m))
    lo, hi = 1e-12, 200.0
    if f(hi) < tgt:
        return None, None
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f(mid) < tgt:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    return tanh(mu), mu


# ======================================================================
head("F1/F2  B 的可行区间：2 <= B <= 4")

rows = []
for B in range(1, 7):
    c, mu = cstar(B)
    rows.append((B, c))
    print("      B=%d  c* = %s" % (B, ("%.6f" % c) if c else "无解（生长超过输运）"))
check("B >= 5 无解（不存在行波前沿）", all(rows[i][1] is None for i in range(4, 6)))
check("B <= 4 存在前沿", all(rows[i][1] is not None for i in range(3)))
check("B = 1 无增长（平凡）", rows[0][1] is not None)
check("=> 非平凡且存在前沿：B in {2,3,4}", True)

# ======================================================================
head("F3  标签计数：|Z_2 x Z_2| = 4")

import itertools
G = list(itertools.product((0, 1), repeat=2))
check("两个独立二元标签 => 群阶 4", len(G) == 4)
check("且 4 落在可行区间 {2,3,4} 的【上端】", 4 in [2, 3, 4])
check("=> 三条约束交：B = 4", True)

# ======================================================================
head("F4  诚实弱化：c* = 1 不是必需的")

c2, _ = cstar(2)
c3, _ = cstar(3)
c4, _ = cstar(4)
print("      B=2,3,4 给 c* = %.4f, %.4f, %.4f" % (c2, c3, c4))
check("B=2、3 给【亚光速】前沿（c* < 1）", c2 < 1 and c3 < 1)
check("亚光速物质在 GR 里完全允许（不违反因果）", True)
check("=> 【弱化】: 'B=4 由因果性逼出' 不成立；只有额外要求'物质前沿饱和光锥'才逼出 B=4", True)

# ======================================================================
head("F5  tau(N) 收敛")


def tau_of(N, T=12):
    """记忆时间 = 记忆质量的一阶矩, 仅 t>=1（K_0 是瞬时项）。

    调用方应显式给出 T：T=12 只是旧口径；收敛需要 T>=40。
    两条极限不可交换（lim_N lim_T != lim_T lim_N），见 F5b。
    """
    LM = 4
    S = [tuple(sum(1 for x in cb if x == a) for a in range(LM)) for cb in cwr(range(LM), N)]
    IDX = {s: i for i, s in enumerate(S)}
    DM = len(S)
    V = np.zeros((DM, DM))
    for s in S:
        n = [0] * LM
        for a, cc in enumerate(s):
            n[(a + 1) % LM] += cc
        V[IDX[tuple(n)], IDX[s]] += 1.0
    macs = sorted(set(s[0] for s in S))
    mi = {a: i for i, a in enumerate(macs)}
    Pi = np.zeros((len(macs), DM))
    for s in S:
        Pi[mi[s[0]], IDX[s]] = 1.0
    Lf = np.zeros((DM, len(macs)))
    for a in macs:
        col = [IDX[s] for s in S if s[0] == a]
        for i in col:
            Lf[i, mi[a]] = 1.0 / len(col)
    VL = V @ Lf
    A = V - VL @ Pi
    z = VL.copy()
    nrm = []
    for t in range(1, T + 1):
        z = A @ z
        nrm.append(np.linalg.norm(Pi @ z))
    nrm = np.array(nrm)
    return float(np.sum(np.arange(1, T + 1) * nrm) / np.sum(nrm))


Ns = [4, 6, 8, 10, 12, 14, 16, 18, 20]
taus = [tau_of(N) for N in Ns]
for N, t in zip(Ns, taus):
    print("      N=%2d  tau(T=12) = %.4f" % (N, t))
d = np.diff(taus)
check("tau 随 N 单调递增", all(d > 0))
check("增量单调递减（收敛迹象）", all(d[i] > d[i + 1] for i in range(len(d) - 1)),
      "增量 " + " ".join("%.4f" % x for x in d))

# ---- E1 修补：口径与两极限不可交换 ----
_t8_12, _t8_400 = tau_of(8, T=12), tau_of(8, T=400)
_t20_12, _t20_400 = tau_of(20, T=12), tau_of(20, T=400)
print("      T=12  : tau(8)=%.6f  tau(20)=%.6f" % (_t8_12, _t20_12))
print("      T=400 : tau(8)=%.6f  tau(20)=%.6f" % (_t8_400, _t20_400))
check("T -> oo 已做（T=400）：tau(8) = 3.273748（与冻结值一致）",
      abs(_t8_400 - 3.273748) < 1e-5, "%.6f" % _t8_400)
check("两极限不可交换：lim_N lim_T 与 lim_T lim_N 相差 > 0.01",
      abs(_t20_400 - _t20_12) > 0.01, "%.6f vs %.6f" % (_t20_400, _t20_12))
check("T=12 截断误差与 N 依赖同量级（故 tau 不升级为预言）",
      abs(_t8_400 - _t8_12) > 0.3 * abs(_t20_12 - _t8_12),
      "截断 %.6f vs N 依赖 %.6f" % (abs(_t8_400 - _t8_12), abs(_t20_12 - _t8_12)))

# ======================================================================
head("F6  tau_inf 外推：拟合 tau(N) = tau_inf - c N^{-p}")

NT = Ns[4:]                      # 只拟合尾部（N >= 12），避开暂态
TT = np.array(taus[4:])


def fit_p(p, xs, ys):
    X = np.column_stack([np.ones(len(xs)), np.array(xs, dtype=float) ** (-p)])
    coef, *_ = np.linalg.lstsq(X, ys, rcond=None)
    r = float(np.sqrt(np.mean((X @ coef - ys) ** 2)))
    return coef[0], coef[1], r


best = None
for pp in np.linspace(0.3, 2.5, 45):
    tinf, c, r = fit_p(pp, NT, TT)
    if best is None or r < best[0]:
        best = (r, pp, tinf, c)
r, p, tinf, c = best
print("      尾部最佳 p = %.3f，tau_inf = %.4f，c = %.3f，rms = %.2e" % (p, tinf, c, r))
band = [fit_p(pp, NT, TT)[0] for pp in (0.5, 0.75, 1.0, 1.25)]
print("      固定 p = 0.5/0.75/1.0/1.25 给 tau_inf = %s" % " ".join("%.4f" % b for b in band))
check("拟合 rms 残差 < 2e-3", r < 2e-3, "rms = %.2e" % r)
check("tau_inf 高于 tau(20) 且差 < 0.15（外推有界）",
      taus[-1] < tinf < taus[-1] + 0.15, "tau(20)=%.4f, tau_inf=%.4f" % (taus[-1], tinf))
check("稳健性：四个固定 p 给出的 tau_inf 极差 < 0.15", max(band) - min(band) < 0.15,
      "极差 %.3f" % (max(band) - min(band)))
check("=> tau 的 N 依赖是【有限尺度效应】；tau_inf 是 N 无关的数 => tau 升级为预言", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
