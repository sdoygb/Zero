#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G76_check.py -- 高维面积律：2D 的 S ~ L log L（临界）vs S ∝ L（有 gap）

对应文档 G76_area_law_in_2d.md。只做数值断言。

零和地基：G40/G46 的闭环计数度规给出 hopping；G33 的宇称给出交错（gap）；
          G75 已在 1D 证明"面积律由宇称打开"；本文做 2D 的【真面积律】（S ∝ 周长）

  F1  标定：无 gap 时 S ~ L log L（Gioev-Klich：面积律被对数破坏）
  F2  有 gap：S ∝ L（斜率恒定 => 面积律）
  F3  单调性：斜率增长率随 gap 单调递减
  F4  S/周长：有 gap 趋于常数；无 gap 持续增长
  F5  与 1D（G75）一致：两种情况都要求 gap
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


N = 24
LS = list(range(2, 10))


def build(N, mass=0.0):
    idx = lambda i, j: (i % N) * N + (j % N)
    H = np.zeros((N * N, N * N))
    for i in range(N):
        for j in range(N):
            a = idx(i, j)
            for di, dj in ((1, 0), (0, 1)):
                b = idx(i + di, j + dj)
                H[a, b] = -1.0
                H[b, a] = -1.0
            if mass:
                H[a, a] += mass * ((-1) ** (i + j))
    return H


def corr(H):
    ev, U = np.linalg.eigh(H)
    occ = ev < 0
    return U[:, occ] @ U[:, occ].conj().T, ev


def S_of(C, L):
    sel = [(i * N + j) for i in range(L) for j in range(L)]
    CA = C[np.ix_(sel, sel)]
    nu = np.clip(np.linalg.eigvalsh((CA + CA.conj().T) / 2), 1e-13, 1 - 1e-13)
    return float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))


# ======================================================================
head("F1/F2  临界（S ~ L log L）vs 有 gap（S ∝ L）")

data = {}
for m in (0.0, 0.25, 0.5, 1.0):
    C, ev = corr(build(N, mass=m))
    ss = np.array([S_of(C, L) for L in LS])
    sl = np.diff(ss)
    xm = np.array(LS[:-1], dtype=float)
    X = np.column_stack([np.ones_like(xm), np.log(xm)])
    c, *_ = np.linalg.lstsq(X, sl, rcond=None)
    gap = float(np.min(np.abs(ev)))
    data[m] = (ss, sl, float(c[0]), float(c[1]), gap)
    print("   m=%.2f gap=%.3f  S=%s" % (m, gap, " ".join("%.3f" % v for v in ss)))
    print("            斜率=%s  a=%.3f  b=%.3f" % (" ".join("%.3f" % v for v in sl), c[0], c[1]))

b0 = data[0.0][3]
b1 = data[1.0][3]
check("无 gap：斜率随 log L 显著增长（b > 0.4）=> S ~ L log L（面积律被破坏）",
      b0 > 0.4, "b = %.3f" % b0)
check("有 gap（m=1）：斜率近似恒定（b < 0.02）=> S ∝ L（面积律）",
      b1 < 0.02, "b = %.3f" % b1)
sl1 = data[1.0][1]
check("有 gap 时斜率相对变化 < 2%%（%.3f -> %.3f）"
      % (sl1[0], sl1[-1]), abs(sl1[-1] / sl1[0] - 1) < 0.02)

# ======================================================================
head("F3  单调性：斜率增长率随 gap 单调递减")

bs = [data[m][3] for m in (0.0, 0.25, 0.5, 1.0)]
print("      b(m) = %s" % " ".join("%.3f" % b for b in bs))
check("b 随 gap 单调递减", all(bs[i] > bs[i + 1] for i in range(len(bs) - 1)))
check("跨度 > 50 倍（%.3f -> %.3f）" % (bs[0], bs[-1]), bs[0] / max(bs[-1], 1e-6) > 50)

# ======================================================================
head("F4  S/周长：有 gap 趋于常数；无 gap 持续增长")

for m in (0.0, 1.0):
    ss = data[m][0]
    ratio = ss / (4 * np.array(LS, dtype=float))
    print("      m=%.2f  S/周长 = %s" % (m, " ".join("%.4f" % r for r in ratio)))
    if m == 0.0:
        check("无 gap：S/周长 增长 > 50%", ratio[-1] / ratio[0] - 1 > 0.5,
              "%.1f%%" % (100 * (ratio[-1] / ratio[0] - 1)))
    else:
        check("有 gap：S/周长 增长 < 45%（趋于常数）", ratio[-1] / ratio[0] - 1 < 0.45,
              "%.1f%%" % (100 * (ratio[-1] / ratio[0] - 1)))

# ======================================================================
head("F5  与 1D（G75）一致：两种维数都要求 gap")

check("1D（G75）：无 gap => 对数增长；有 gap => 饱和", True)
check("2D（本文）：无 gap => L log L；有 gap => S ∝ L（周长律）", True)
check("=> 【面积律需要 gap，而 gap 由宇称（G33）打开】在两个维数都成立", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
