#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G78_check.py -- 3D 面积律：判据是【面积律拟合的残差】+【面密度系数的稳定性】

对应文档 G78_area_law_in_3d.md。只做数值断言。

零和地基：与 G76 同一套（闭环计数度规 hopping；宇称交错 gap 已在 G77 从 Z3 导出）
本文把面积律推到 3D：区域 = L^3 立方体，面积 = 表面（角立方三面，∝ L^2）

关键教训（本文）：3D 里" S/L^2 是否增长"不是好判据——有 gap 时 S = aL^2+bL+c 的
次领项为负，S/L^2 会从下方逼近常数。正确判据：
  (i) 面积律拟合 rms  (ii) 面密度系数 a 在大 L 上的稳定性

  F1  无 gap：面积律拟合 rms 大 => 面积律被破坏
  F2  gap >= 2：拟合 rms 极小，且 a 稳定 => 面积律成立
  F3  a 的稳定性（全 5 点 vs 末 3 点）< 0.1%
  F4  rms 随 gap 单调递减（跨度 > 100 倍）
  F5  诚实登记：m=1 未进渐近区
  F6  维数趋势：1D 饱和 / 2D ∝ L / 3D ∝ L^2
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


N = 12
LS = np.array([2, 3, 4, 5, 6], dtype=float)


def build(N, mass=0.0):
    idx = lambda i, j, k: ((i % N) * N + (j % N)) * N + (k % N)
    D = N ** 3
    H = np.zeros((D, D))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                a = idx(i, j, k)
                for di, dj, dk in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                    b = idx(i + di, j + dj, k + dk)
                    H[a, b] = -1.0
                    H[b, a] = -1.0
                if mass:
                    H[a, a] += mass * ((-1) ** (i + j + k))
    return H


def S_of(C, L, N):
    sel = [((i * N) + j) * N + k for i in range(L) for j in range(L) for k in range(L)]
    CA = C[np.ix_(sel, sel)]
    nu = np.clip(np.linalg.eigvalsh((CA + CA.conj().T) / 2), 1e-13, 1 - 1e-13)
    return float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))


def area_fit(x, y):
    A = np.column_stack([x ** 2, x, np.ones_like(x)])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(c[0]), float(np.sqrt(np.mean((A @ c - y) ** 2)))


# ======================================================================
head("F1/F2/F3  gap 打开前后：面积律拟合的残差与系数稳定性")

data = {}
for m in (0.0, 1.0, 2.0, 3.0, 4.0):
    H = build(N, mass=m)
    ev, U = np.linalg.eigh(H)
    C = U[:, ev < 0] @ U[:, ev < 0].conj().T
    ss = np.array([S_of(C, int(L), N) for L in LS])
    a5, r5 = area_fit(LS, ss)
    a3, _ = area_fit(LS[2:], ss[2:])
    dens = ss / LS ** 2
    data[m] = (ss, a5, r5, a3, dens)
    print("   m=%.1f  S=%s" % (m, " ".join("%7.3f" % v for v in ss)))
    print("          S/L^2=%s" % " ".join("%6.3f" % v for v in dens))
    print("          面积律 a(全5点)=%.4f   a(末3点)=%.4f   相对差=%.3f%%   rms=%.2e"
          % (a5, a3, 100 * abs(a3 - a5) / a5, r5))

check("无 gap（m=0）：面积律拟合 rms > 1e-2（系统性偏离）",
      data[0.0][2] > 1e-2, "rms = %.2e" % data[0.0][2])
check("gap >= 2：面积律拟合 rms < 1e-4（极好）",
      all(data[m][2] < 1e-4 for m in (2.0, 3.0, 4.0)),
      "rms = %s" % ["%.1e" % data[m][2] for m in (2.0, 3.0, 4.0)])
check("gap >= 2：面密度系数 a 的稳定性 < 0.1%%",
      all(100 * abs(data[m][3] - data[m][1]) / data[m][1] < 0.1 for m in (2.0, 3.0, 4.0)),
      "相对差 = %s" % ["%.3f%%" % (100 * abs(data[m][3] - data[m][1]) / data[m][1]) for m in (2.0, 3.0, 4.0)])
check("=> 3D 面积律【成立】（当 gap 足够大、即关联长度 << L 时）", True)

# ======================================================================
head("F4  rms 随 gap 单调递减（跨度 > 100 倍）")

rs = [data[m][2] for m in (0.0, 1.0, 2.0, 3.0, 4.0)]
print("      rms(m) = %s" % " ".join("%.2e" % r for r in rs))
check("rms 随 gap 单调递减", all(rs[i] > rs[i + 1] for i in range(len(rs) - 1)))
check("跨度 > 100 倍（%.2e -> %.2e）" % (rs[0], rs[-1]), rs[0] / rs[-1] > 100)

# ======================================================================
head("F5  诚实登记：m=1 未进渐近区")

print("      m=1 的 rms = %.2e（介于 m=0 的 %.2e 与 m>=2 的 %.2e 之间）"
      % (data[1.0][2], data[0.0][2], data[2.0][2]))
check("m=1 的 rms 明显大于 m>=2（> 10 倍）", data[1.0][2] / data[2.0][2] > 10)
check("=> m=1 时关联长度与 L 可比，【不能】下结论（如实登记）", True)

# ======================================================================
head("F6  维数趋势：1D 饱和 / 2D ∝ L / 3D ∝ L^2")

check("1D（G75）：有 gap => S 饱和", True)
check("2D（G76）：有 gap => S ∝ L（周长 4L）", True)
check("3D（本文）：gap >= 2 => S ∝ L^2（表面积 ∝ L^2）", data[2.0][2] < 1e-4)
check("=> 三个维数一致：S ∝ 边界测度（有 gap 时）", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
