#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R72_layer_audit.py -- 复算 R72 的三条选维链与层归属断言

F1  R23 演化层链：S_D^evo = a_L^D，a_L = C(L,L/2)/2^L < 1
F2  该链在 D>=4 上的最大者恒为 D=4（对 L=4..24）
F3  账本链：F_D = C(D,2) q^D，q in (0.5,0.6) => 峰 D=4
F4  窗口宽 Delta_D = 2/(D(D+1)) 在 D>=4 上最大者为 D=4
F5  号差：(1,m-1) 由随机正定 h 核验（m=2..5）
F6  G89 命题 2：原生对合给 (1,1) 当且仅当 m=3
F7  因果区间多项式的 d 阶差分恒为常数（d_s=1,2）
"""
import sys
import math
from itertools import product
from math import comb

import numpy as np

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name,
                          ("   " + detail) if detail else ""))


def head(t):
    print("\n" + "=" * 74)
    print(t)
    print("=" * 74)


# F1/F2 ---------------------------------------------------------------
def evo_survival(L, D):
    a = comb(L, L // 2) / 2 ** L
    return a ** D


if __name__ == "__main__":
    head("F1  R23 演化层链：a_L < 1")
    for L in [4, 6, 8, 16, 20, 24]:
        a = comb(L, L // 2) / 2 ** L
        check("a_L < 1 (L=%d)" % L, 0 < a < 1, "a=%.6f" % a)

    head("F2  该链在 D>=4 上的最大者恒为 D=4")
    for L in [4, 6, 8, 16, 20, 24]:
        vals = {D: evo_survival(L, D) for D in range(4, 9)}
        mx = max(vals, key=lambda k: vals[k])
        check("argmax_{D>=4} S_D^evo = 4 (L=%d)" % L, mx == 4,
              "argmax=%d  S_4=%.6e" % (mx, vals[4]))

    head("F3  账本链：q in (0.5,0.6) => 峰 D=4")
    def ledger_peak(q):
        for D in range(2, 300):
            if D == 2:
                if 3 * q < 1:
                    return 2
                continue
            if D * q > (D - 2) and (D + 1) * q < (D - 1):
                return D
        return None
    check("q=5/9 峰=4", ledger_peak(5 / 9) == 4, "峰=%s" % ledger_peak(5 / 9))
    check("q=0.5578 峰=4", ledger_peak(0.5578) == 4, "峰=%s" % ledger_peak(0.5578))
    check("q=0.45 峰=3", ledger_peak(0.45) == 3, "峰=%s" % ledger_peak(0.45))
    check("q=0.65 峰=5", ledger_peak(0.65) == 5, "峰=%s" % ledger_peak(0.65))

    head("F4  窗口宽 Delta_D = 2/(D(D+1))，D>=4 上最大者为 D=4")
    widths = {D: 2.0 / (D * (D + 1)) for D in range(4, 15)}
    mx = max(widths, key=lambda k: widths[k])
    check("argmax Delta_D = 4", mx == 4, "argmax=%d  Delta_4=%.6f" % (mx, widths[4]))
    check("Delta_D 随 D 严格递减",
          all(widths[D] > widths[D + 1] for D in range(4, 14)))

    head("F5  号差 (1,m-1) 由随机正定 h 核验")
    rng = np.random.default_rng(3)
    ok_all = True
    for m in [2, 3, 4, 5]:
        ok = 0
        for _ in range(200):
            A = rng.normal(size=(m, m))
            h = A @ A.T + 0.1 * np.eye(m)
            g = np.zeros((m + 1, m + 1))
            g[0, 0] = 1.0
            g[1:, 1:] = -h
            ev = np.linalg.eigvalsh(g)
            if int(np.sum(ev > 1e-10)) == 1 and int(np.sum(ev < -1e-10)) == m:
                ok += 1
        check("m=%d 号差 (1,%d) 200/200" % (m, m), ok == 200, "%d/200" % ok)
        ok_all = ok_all and ok == 200

    head("F6  G89 命题 2：原生对合给 (1,1) 当且仅当 m=3")
    def cycles(sig):
        m = len(sig)
        seen = [False] * m
        out = []
        for i in range(m):
            if seen[i]:
                continue
            c, j = 0, i
            while not seen[j]:
                seen[j] = True
                j = sig[j]
                c += 1
            out.append(c)
        return out
    from itertools import permutations
    for m in range(2, 8):
        n11 = 0
        for sig in permutations(range(m)):
            if any(sig[sig[i]] != i for i in range(m)):
                continue
            cyc = cycles(sig)
            for eps in (1, -1):
                dp = len(cyc) - 1 if eps == 1 else sum(1 for c in cyc if c % 2 == 0)
                dm = sum(1 for c in cyc if c % 2 == 0) if eps == 1 else len(cyc) - 1
                if (dp, dm) == (1, 1):
                    n11 += 1
        if m == 3:
            check("m=3 给出 (1,1) 的对合数 > 0", n11 > 0, "%d" % n11)
        else:
            check("m=%d 给出 (1,1) 的对合数 = 0" % m, n11 == 0, "%d" % n11)

    head("F7  因果区间多项式的 d 阶差分恒为常数")
    def cone(T, dims):
        return [set(p for p in product(range(-t, t + 1), repeat=dims)
                    if sum(abs(x) for x in p) <= t
                    and (t - sum(abs(x) for x in p)) % 2 == 0)
                for t in range(T + 1)]

    def interval(T, dtau, dx, dims):
        c = cone(T, dims)
        tgt = tuple((dx[i] if i < len(dx) else 0) for i in range(dims))
        return sum(1 for tau in range(dtau + 1) for p in c[tau]
                   if tuple(tgt[i] - p[i] for i in range(dims)) in c[dtau - tau])

    for ds in [1, 2]:
        T = 12
        s = [interval(T, dt, (0,) * ds, ds) for dt in range(0, T + 1, 2)]
        v = s[:]
        for _ in range(ds + 1):
            v = [v[i + 1] - v[i] for i in range(len(v) - 1)]
        check("d_s=%d 的 %d 阶差分为常数" % (ds, ds + 1),
              all(abs(x - v[0]) < 1e-9 for x in v), "差分值=%d" % v[0])
        expect = 2 ** ds
        check("d_s=%d 的差分值 = 2^d_s = %d" % (ds, expect), v[0] == expect,
              "实得 %d" % v[0])

    head("汇总")
    print("  通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
    if FAIL:
        print("  不符项：", FAIL)
    sys.exit(1 if FAIL else 0)
