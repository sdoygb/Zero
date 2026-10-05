#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R74_check.py -- 复算 R74 的全部断言

F1  欧氏单纯形账本的 D=4 窗口 = (3/5, 2/3)（精确分数）
F2  该窗口严格落在语境性区 q>3/5 之上，长度 1/15
F3  L=16 的 βε 窗口 [1.200, 1.470]，宽度 ≈ 0.27
F4  窗口内峰恒为 D=4
F5  窗口内 S_max > 2（最小值 2.0178）
F6  窗口中心 vs 4/3（差 < 0.005）
F7  洛伦兹字典下同一测度：βε=1.13 给峰 4，βε=1.20 给峰 5（标号平移）
F8  L 稳健性：L=12..28 的窗口中心都在 4/3 的 1% 内
"""
import math
import sys
from fractions import Fraction as F
from math import comb

MU = (math.sqrt(5), (5 - math.sqrt(5)) / 2, (5 - math.sqrt(5)) / 2)
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name,
                          ("   " + detail) if detail else ""))


def head(t):
    print("\n" + "=" * 74)
    print(t)
    print("=" * 74)


def blocks(L, g):
    r = math.exp(-g)
    b = {}
    for k in range(L + 1):
        for s in range(-k, k + 1, 2):
            rem = L - k
            if (rem - s) % 2:
                continue
            a = (rem - s) // 2
            if a < 0 or a > rem:
                continue
            b[abs(s)] = b.get(abs(s), 0) + comb(rem, a) * (r ** abs(s))
    tot = sum(b.values())
    return sorted((v / tot for v in b.values()), reverse=True)


def q_of(v):
    return sum(x * x for x in v)


def S_of(v):
    t3 = v[:3]
    s = sum(t3)
    return sum(a * b for a, b in zip([x / s for x in t3], MU))


def peak(M, q, Dmax=60):
    best, bd = -1.0, None
    for D in range(2, Dmax):
        val = M(D) * (q ** D)
        if val > best:
            best, bd = val, D
    return bd


def bisect(L, target, lo=0.0, hi=6.0):
    for _ in range(150):
        mid = (lo + hi) / 2
        if q_of(blocks(L, mid)) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


M_EUC = lambda D: comb(D + 1, 2)
M_LOR = lambda D: comb(D, 2)

if __name__ == "__main__":
    head("F1 欧氏单纯形账本的 D=4 窗口（精确分数）")
    lo = F(M_EUC(3), M_EUC(4))
    hi = F(M_EUC(4), M_EUC(5))
    check("窗口 = (3/5, 2/3)", lo == F(3, 5) and hi == F(2, 3),
          "(%s, %s)" % (lo, hi))

    head("F2 窗口与语境性区的关系")
    thr = F(3, 5)
    check("上端 2/3 > 3/5", hi > thr, "2/3=%.6f" % float(hi))
    check("下端 3/5 = 阈值（相接）", lo == thr)
    check("严格落入区内长度 = 1/15", (hi - lo) == F(1, 15),
          "1/15=%.6f" % float(F(1, 15)))
    # 洛伦兹对照：上端恰好等于阈值
    lo_l = F(M_LOR(3), M_LOR(4))
    hi_l = F(M_LOR(4), M_LOR(5))
    check("洛伦兹窗口 = (1/2, 3/5)，上端 = 阈值",
          lo_l == F(1, 2) and hi_l == thr, "(%s, %s)" % (lo_l, hi_l))
    check("洛伦兹窗口严格落入区内长度 = 0", (hi_l - max(lo_l, thr)) == 0)

    head("F3/F4/F5  L=16 的两门窗口")
    g_lo = bisect(16, 0.6)
    g_hi = bisect(16, 2 / 3)
    check("βε 下端 ≈ 1.200", abs(g_lo - 1.200) < 5e-4, "%.4f" % g_lo)
    check("βε 上端 ≈ 1.470", abs(g_hi - 1.470) < 5e-4, "%.4f" % g_hi)
    check("窗口宽度 ≈ 0.27", abs((g_hi - g_lo) - 0.270) < 5e-3,
          "%.4f" % (g_hi - g_lo))
    eps = (g_hi - g_lo) / 1000.0
    peaks = set()
    smin = 1e9
    for i in range(201):
        g = g_lo + eps + (g_hi - g_lo - 2 * eps) * i / 200   # 严格内部
        v = blocks(16, g)
        peaks.add(peak(M_EUC, q_of(v)))
        smin = min(smin, S_of(v))
    check("窗口**内部**峰恒为 D=4", peaks == {4}, "峰集=%s" % sorted(peaks))
    check("窗口内 S_max > 2", smin > 2.0, "最小值=%.6f" % smin)
    # 边界点：q 恰为 3/5 与 2/3 处出现并列峰（相接）
    qb1 = q_of(blocks(16, g_lo))
    qb2 = q_of(blocks(16, g_hi))
    check("下端边界 q≈3/5，峰并列（相接）",
          abs(qb1 - 0.6) < 1e-4, "q=%.6f" % qb1)
    check("上端边界 q≈2/3，峰并列（相接）",
          abs(qb2 - 2 / 3) < 1e-4, "q=%.6f" % qb2)

    head("F6 窗口中心 vs 4/3")
    c = (g_lo + g_hi) / 2
    check("中心 ≈ 4/3（差 < 0.005）", abs(c - 4 / 3) < 0.005,
          "中心=%.4f  4/3=%.4f  差=%+.4f" % (c, 4 / 3, c - 4 / 3))

    head("F7 洛伦兹字典下同一测度（标号平移）")
    q113 = q_of(blocks(16, 1.13))
    q120 = q_of(blocks(16, 1.20))
    check("βε=1.13: 洛伦兹峰=4 而欧氏峰=3",
          peak(M_LOR, q113) == 4 and peak(M_EUC, q113) == 3,
          "q=%.6f" % q113)
    check("βε=1.20: 洛伦兹峰=5 而欧氏峰=4",
          peak(M_LOR, q120) == 5 and peak(M_EUC, q120) == 4,
          "q=%.6f" % q120)

    head("F8  L 稳健性")
    for L in [12, 14, 16, 18, 20, 24, 28]:
        a = bisect(L, 0.6)
        b = bisect(L, 2 / 3)
        cc = (a + b) / 2
        sm = min(S_of(blocks(L, a + (b - a) * i / 60)) for i in range(61))
        check("L=%d 窗口中心在 4/3 的 1%% 内且 S_max>2" % L,
              abs(cc - 4 / 3) / (4 / 3) < 0.01 and sm > 2.0,
              "中心=%.4f 宽=%.4f S_min=%.6f" % (cc, b - a, sm))

    head("汇总")
    print("  通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
    if FAIL:
        print("  不符项：", FAIL)
    sys.exit(1 if FAIL else 0)
