#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R20_A1_probe.py —— A1（几何一致性）在 1+1D 自由费米子支线上的数值判定
=====================================================================
离线自足；不访问网络，不修改任何项目文件。

判定对象（R19 §5 的 A1）：
    模流与闭合旋转流是否满足 sigma_{theta/2 pi} = alpha_theta，
    等价地：模 Hamiltonian h = log((1-C)/C) 是否在连续极限里正比于
    几何 boost 的生成元 B（权函数 xi(x) = x(l-x)/l 取在键中点），
    且比例常数取到 A1 要求的 2 pi。

约定：站点 x_j = j（j=0..N-1），区间长 l = N-1，键 b 的键中点在 x=b+1/2，
      故 beta_b = (b+1/2)(l-b-1/2)/l（对 b -> l-1-b 对称）。

两个检验层：
  (a) 近邻权重 w_b = |h_{b,b+1}| 与 beta_b 的逐点比例（局部形状）；
  (b) 二次型尺度 lambda_f = q_h(f)/q_B(f)（连续极限的正确比较层）。

这是有限 N 数值判定，不是定理。
"""
import sys

import mpmath as mp

mp.mp.dps = 60

NS = (16, 24, 32, 48, 64, 80)
TRIM = 1

ASSERTIONS = 0
FAILURES = 0


def check(name, condition, detail=""):
    global ASSERTIONS, FAILURES
    ASSERTIONS += 1
    ok = bool(condition)
    if not ok:
        FAILURES += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))
    return ok


def fmt(x, digits=10):
    return mp.nstr(x, digits)


def build_h(N):
    """C_ij = sin(pi(i-j)/2)/(pi(i-j))（i!=j），C_ii=1/2；h=log((1-C)/C)。"""
    C = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            d = i - j
            C[i, j] = mp.mpf(1) / 2 if d == 0 else mp.sin(mp.pi * d / 2) / (mp.pi * d)
    ev, evec = mp.eigsy(C)
    h = mp.matrix(N, N)
    for k in range(N):
        val = mp.log((1 - ev[k]) / ev[k])
        for i in range(N):
            for j in range(N):
                h[i, j] += val * evec[i, k] * evec[j, k]
    return h


def beta_bond(N):
    """键中点的 xi：l = N-1，x_b = b+1/2。"""
    length = mp.mpf(N - 1)
    return [
        (mp.mpf(b) + mp.mpf(1) / 2) * (length - mp.mpf(b) - mp.mpf(1) / 2) / length
        for b in range(N - 1)
    ]


def test_functions(N):
    """若干光滑测试函数：f_i = phi(i/l)。"""
    l = mp.mpf(N - 1)
    us = [mp.mpf(i) / l for i in range(N)]
    return {
        "sin(pi u)": [mp.sin(mp.pi * u) for u in us],
        "u(1-u)": [u * (1 - u) for u in us],
        "u(1-u)^2": [(u * (1 - u)) ** 2 for u in us],
        "sin(2 pi u)": [mp.sin(2 * mp.pi * u) for u in us],
    }


def quad(M, vec, n):
    return mp.fsum(vec[i] * M[i, j] * vec[j] for i in range(n) for j in range(n))


def boost_matrix(N):
    """B：键权 beta_b 的 hopping 矩阵（无对角）。"""
    beta = beta_bond(N)
    B = mp.matrix(N, N)
    for b in range(N - 1):
        B[b, b + 1] = beta[b]
        B[b + 1, b] = beta[b]
    return B


def fit_linear(xs, ys):
    n = len(xs)
    sx, sy = mp.fsum(xs), mp.fsum(ys)
    sxx = mp.fsum(x * x for x in xs)
    sxy = mp.fsum(xs[i] * ys[i] for i in range(n))
    det = n * sxx - sx * sx
    b = (n * sxy - sx * sy) / det
    return (sy - b * sx) / n, b


def main():
    print("=" * 72)
    print("R20 A1 probe: is h proportional to the boost generator B?")
    print("=" * 72)
    print("mpmath precision: %d digits; N: %s" % (mp.mp.dps, ", ".join(map(str, NS))))

    A_list, rho_list, R_list, spread_list = [], [], [], []
    lam_mean, lam_spread, lam_by_f = [], [], {}
    profile_last = None

    print("\n[1] 近邻权重比例（局部形状）")
    print("     N |  A_N (最小二乘) | rho_N (中心键) | R_N (残差) | 半区局部比例相对散布")
    for N in NS:
        h = build_h(N)
        w = [abs(h[i, i + 1]) for i in range(N - 1)]
        beta = beta_bond(N)
        uw, ub = w[TRIM:len(w) - TRIM], beta[TRIM:len(beta) - TRIM]
        A = mp.fsum(ub[i] * uw[i] for i in range(len(uw))) / mp.fsum(b * b for b in ub)
        R = mp.sqrt(mp.fsum((uw[i] - A * ub[i]) ** 2 for i in range(len(uw)))) / mp.sqrt(
            mp.fsum(u * u for u in uw)
        )
        c = N // 2
        rho = w[c] / beta[c]
        lo, hi = N // 4, 3 * N // 4
        scales = [w[b] / beta[b] for b in range(lo, hi)]
        mean_s = mp.fsum(scales) / len(scales)
        spread = mp.sqrt(mp.fsum((s - mean_s) ** 2 for s in scales) / len(scales)) / mean_s
        A_list.append(A); rho_list.append(rho); R_list.append(R); spread_list.append(spread)
        print("    %2d | %s | %s | %s | %s"
              % (N, fmt(A, 10), fmt(rho, 10), fmt(R, 10), fmt(spread, 10)))
        if N == NS[-1]:
            profile_last = [(b, w[b] / beta[b]) for b in range(N - 1)]

    print("\n[2] 二次型尺度 lambda_f = q_h(f)/q_B(f)（连续极限的正确比较层）")
    names = list(test_functions(NS[0]).keys())
    print("     N | " + " | ".join("%-14s" % k for k in names) + " | 跨 f 相对散布")
    for N in NS:
        h = build_h(N)
        B = boost_matrix(N)
        funcs = test_functions(N)
        lam_N = []
        row = []
        for name in names:
            f = funcs[name]
            lam = quad(h, f, N) / quad(B, f, N)
            lam_N.append(lam)
            lam_by_f.setdefault(name, []).append(lam)
            row.append("%s" % fmt(lam, 10))
        mean_lam = mp.fsum(lam_N) / len(lam_N)
        spr = mp.sqrt(mp.fsum((x - mean_lam) ** 2 for x in lam_N) / len(lam_N)) / abs(mean_lam)
        lam_mean.append(mean_lam); lam_spread.append(spr)
        print("    %2d | " % N + " | ".join("%-14s" % v for v in row) + " | %s" % fmt(spr, 8))

    print("\n[3] 1/N 外推（y = a + b/N）")
    inv = [mp.mpf(1) / N for N in NS]
    aA, _ = fit_linear(inv, A_list)
    aRho, _ = fit_linear(inv, rho_list)
    aR, _ = fit_linear(inv, R_list)
    aSp, _ = fit_linear(inv, spread_list)
    aLam, _ = fit_linear(inv, lam_mean)
    aLS, _ = fit_linear(inv, lam_spread)
    print("    A_N        -> %s" % fmt(aA, 10))
    print("    rho_N      -> %s" % fmt(aRho, 10))
    print("    R_N        -> %s" % fmt(aR, 10))
    print("    spread     -> %s" % fmt(aSp, 10))
    print("    lambda_N   -> %s" % fmt(aLam, 10))
    print("    lam-spread -> %s" % fmt(aLS, 10))
    for name in names:
        a, _ = fit_linear(inv, lam_by_f[name])
        print("      lambda(%s) -> %s" % (name, fmt(a, 10)))

    print("\n[4] 候选常数对照（对 lambda_inf）")
    cands = [
        ("2 pi", 2 * mp.pi),
        ("pi (=2 pi / v_F, v_F=2)", mp.pi),
        ("pi^2/3 (R15 的经验值)", mp.pi ** 2 / 3),
        ("4 pi / 3", 4 * mp.pi / 3),
    ]
    for name, val in cands:
        print("    %-26s = %s | lambda_inf 相对差 %s" % (name, fmt(val, 10), fmt(abs(abs(aLam) - val) / val, 6)))

    print("\n[5] 最大 N 的局部比例 w_b / beta_b 剖面（含镜像对）")
    N = NS[-1]
    n_b = N - 1
    for b in (n_b // 10, n_b // 4, n_b // 2, 3 * n_b // 4):
        mirror = n_b - 1 - b
        print("    b=%3d (b/l=%s): w/beta=%s   |   mirror b=%3d: %s"
              % (b, fmt(mp.mpf(b) / (N - 1), 4), fmt(profile_last[b], 10),
                 mirror, fmt(profile_last[mirror], 10)))

    print("\n[6] 判定（报告发现，不预设期待值）")
    lam_abs = abs(aLam)
    ests = [("A_N 最小二乘", aA), ("rho_N 中心键", aRho), ("|lambda_N| 二次型", lam_abs)]
    for name, val in ests:
        print("    %-16s -> %s" % (name, fmt(val, 8)))
    lo = min(v for _, v in ests)
    hi = max(v for _, v in ests)
    check("三个估计器的外推常数都落在 [3.2, 3.6]",
          lo > mp.mpf("3.2") and hi < mp.mpf("3.6"),
          "范围 [%s, %s]" % (fmt(lo, 8), fmt(hi, 8)))
    check("估计器之间的相对散布小于 4%",
          (hi - lo) / lo < mp.mpf("0.04"),
          "(hi-lo)/lo=%s" % fmt((hi - lo) / lo, 6))
    check("2 pi 被排除：2 pi / max(估计器) > 1.7",
          (2 * mp.pi) / hi > mp.mpf("1.7"),
          "2pi/%s=%s" % (fmt(hi, 8), fmt((2 * mp.pi) / hi, 6)))
    check("R15 的 pi^2/3 与最小二乘估计器一致到 1% 以内",
          abs(aA - mp.pi ** 2 / 3) / (mp.pi ** 2 / 3) < mp.mpf("0.01"),
          "A_inf=%s, pi^2/3=%s" % (fmt(aA, 10), fmt(mp.pi ** 2 / 3, 10)))
    check("pi^2/3 与二次型估计器相差超过 2%（常数是估计器依赖的）",
          abs(lam_abs - mp.pi ** 2 / 3) / (mp.pi ** 2 / 3) > mp.mpf("0.02"),
          "|lambda_inf|=%s" % fmt(lam_abs, 10))
    check("形状残差 R_N 总体下降并趋于 ~1% 平台",
          R_list[-1] < R_list[0] and R_list[-1] < mp.mpf("0.02"),
          "R: %s" % " -> ".join(fmt(r, 6) for r in R_list))
    check("二次型跨光滑测试函数的散布小于 2%",
          max(abs(s) for s in lam_spread) < mp.mpf("0.02"),
          "lam-spread: %s" % " -> ".join(fmt(s, 6) for s in lam_spread))

    print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
    print("  断言 %d 项，不符项：%d" % (ASSERTIONS, FAILURES))
    if FAILURES:
        print("  未通过")
        sys.exit(1)
    print("  全部通过")
    sys.exit(0)


if __name__ == "__main__":
    main()
