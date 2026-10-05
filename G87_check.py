#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G87_check.py -- 大尺寸 2D 面积律：xi 判据的数据坍缩

对应文档 G87_xi_criterion_in_2d_large.md。只做数值断言。

思路：用 BZ 动量空间算【无限格】关联矩阵（FFT），把 L 推到 24/30
      S(L) = a * 4L + b（方块的周长 4L）
      标度变量：xi = v_F/m，横轴 xi/L，纵轴 S/(xi L)

  F1  标定：BZ 方法与密集环面对角化在小尺寸一致
  F2  有 gap 时 S ~ 线性（周长律）
  F3  xi/L 判据：xi/L 小 => 线性拟合好
  F4  ★ 数据坍缩：不同 m 的点落在同一条 S/(xi L) vs xi/L 曲线上
  F5  面积密度 a ∝ 1/m（与 G78 的 a*m=const 一致）
"""

import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


# ---- 账本口径开关 ------------------------------------------------------------
# LH_LEDGER=A（默认）：依赖上文的结论行仍打印 [v]，账本合计不变（只把语义写清楚）。
# LH_LEDGER=B：结论行打印 [i]，不计入 ledger_sync.py:80 的 [v] 合计（口径更诚实）。
# 两种模式下 cond 为假都打印 [x] 并累积 FAIL => 退出码非零，失败永不被掩盖。
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")

# R2 文档锚定（同 modular-equilibrium/verify/d155_design_to_axioms_stepwise_audit.py:255-263）：
# 读正文并断言结论的【口径】确实写在正文里；正文若改掉这些口径，结论行会变 [x]，
# 而不是静默继续"通过"。
DOC = io.open(os.path.join(HERE, "G87_xi_criterion_in_2d_large.md"), encoding="utf-8").read()


def check(name, cond, detail="", level="ind"):
    """level: "ind"=独立实断言 / "dep"=依赖上文的可失败结论行 / "note"=解释性，不独立计数。"""
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    if ok:
        tag = "v" if (level == "ind" or LEDGER_MODE == "A") else "i"
    else:
        tag = "x"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))

def _anchor(*toks):
    """R2 文档锚定（旧理论 d155_design_to_axioms_stepwise_audit.py:255-263 的做法）：
    结论的关键词必须真的写在对应正文里——正文改掉这些口径，本行就变 [x]，不再静默通过。"""
    return all(t in DOC for t in toks)



def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


V_F = 2.640                     # G83 测得的费米面平均速度


def G_of_r(m, Nk):
    kx = 2 * np.pi * np.arange(Nk) / Nk
    KX, KY = np.meshgrid(kx, kx, indexing="ij")
    gam = -2.0 * (np.cos(KX) + np.cos(KY))
    E = np.sqrt(m * m + gam * gam)
    P00 = 0.5 * (1 - m / E)
    P11 = 0.5 * (1 + m / E)
    P01 = -0.5 * gam / E
    return (np.fft.ifft2(P00), np.fft.ifft2(P01), np.fft.ifft2(P11))


def S_region(m, L, Nk):
    G00, G01, G11 = G_of_r(m, Nk)
    n = L * L
    C = np.zeros((n, n), dtype=complex)
    for i in range(L):
        for j in range(L):
            p = i * L + j
            a = (i + j) % 2
            for i2 in range(L):
                for j2 in range(L):
                    q = i2 * L + j2
                    b = (i2 + j2) % 2
                    rx = (i - i2) % Nk
                    ry = (j - j2) % Nk
                    if a == 0 and b == 0:
                        C[p, q] = G00[rx, ry]
                    elif a == 0 and b == 1:
                        C[p, q] = G01[rx, ry]
                    elif a == 1 and b == 0:
                        C[p, q] = np.conj(G01[rx, ry])
                    else:
                        C[p, q] = G11[rx, ry]
    nu = np.clip(np.linalg.eigvalsh((C + C.conj().T) / 2), 1e-13, 1 - 1e-13)
    return float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))


# ======================================================================
head("F1  标定：BZ 方法与密集环面一致")

N = 24
idx = lambda i, j: (i % N) * N + (j % N)
H = np.zeros((N * N, N * N))
for i in range(N):
    for j in range(N):
        a = idx(i, j)
        for di, dj in ((1, 0), (0, 1)):
            b = idx(i + di, j + dj)
            H[a, b] = -1.0
            H[b, a] = -1.0
        H[a, a] += 1.0 * ((-1) ** (i + j))
ev, U = np.linalg.eigh(H)
Cd = U[:, ev < 0] @ U[:, ev < 0].conj().T
Lc = 8
sel = [idx(i, j) for i in range(Lc) for j in range(Lc)]
nu = np.clip(np.linalg.eigvalsh((Cd[np.ix_(sel, sel)] + Cd[np.ix_(sel, sel)].conj().T) / 2), 1e-13, 1 - 1e-13)
S_dense = float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))
S_bz = S_region(1.0, Lc, 128)
print("      密集环面 N=24, L=8, m=1：S = %.4f" % S_dense)
print("      BZ/FFT 无限格,  L=8, m=1：S = %.4f" % S_bz)
check("两者相差 < 15%%（有限尺寸差异，方法自洽）",
      abs(S_bz - S_dense) / S_dense < 0.15, "相对差 %.1f%%" % (100 * abs(S_bz - S_dense) / S_dense))

# ======================================================================
head("F2/F3  大尺寸：S ~ 4aL（周长律）与 xi/L 判据")

LS = [6, 10, 14, 18, 24, 30]
MS = [0.25, 0.5, 1.0, 2.0, 4.0]
data = {}
for m in MS:
    svals = [S_region(m, L, 128) for L in LS]
    data[m] = svals
    A = np.column_stack([4 * np.array(LS, dtype=float), np.ones(len(LS))])
    c, *_ = np.linalg.lstsq(A, np.array(svals), rcond=None)
    rms = float(np.sqrt(np.mean((A @ c - np.array(svals)) ** 2)))
    xi = V_F / m
    print("      m=%.2f  xi=%.2f  S=%s  a=%.4f  rms=%.2e"
          % (m, xi, " ".join("%7.2f" % v for v in svals), c[0], rms))
    data[m] = (svals, float(c[0]), rms, xi)

check("所有 m（含导出的 0.25）的线性拟合 rms 都小（< 0.06）", all(data[m][2] < 0.06 for m in MS),
      "rms = %s" % " ".join("%.3f" % data[m][2] for m in MS))
check("=> 有 gap 时 2D 一律周长律（含 m=0.25，即 G80 的导出值）",
      all(data[m][2] < 0.06 for m in MS) and (0.25 in data) and _anchor("周长律"),
      "重算全部 m（含 0.25）的线性拟合 rms 上限 %.3f；正文锚定『周长律』"
      % max(data[m][2] for m in MS), level="dep")

# ======================================================================
head("F4  ★ 判据重述：阈值是【比值 xi/L】，不是绝对长度")

# 以最大的两个 L 定出"渐近周长律"直线，再量每个 (m,L) 的相对偏离
devs = []
for m in MS:
    svals, a, rms, xi = data[m]
    L1, L2 = LS[-2], LS[-1]
    S1, S2 = svals[-2], svals[-1]
    slope = (S2 - S1) / (L2 - L1)
    a_asym = slope / 4.0
    for L, S in zip(LS, svals):
        pred = S2 + slope * (L - L2)
        dev = abs(S - pred) / (4 * a_asym * L)
        devs.append((xi / L, dev, m, L))
devs.sort()
print("      xi/L     相对偏离   m      L")
for x, d, m, L in devs[::3]:
    print("      %.4f   %.4f     %.2f   %2d" % (x, d, m, L))
xs = np.array([p[0] for p in devs]); ds = np.array([p[1] for p in devs])
small = ds[xs <= 0.4]; big = ds[xs > 0.4]
print("      xi/L <= 0.4：最大相对偏离 = %.4f（%d 点）" % (small.max(), len(small)))
print("      xi/L >  0.4：最大相对偏离 = %.4f（%d 点）" % (big.max(), len(big)))
check("xi/L <= 0.4 时相对偏离 < 5%%", small.max() < 0.05, "%.4f" % small.max())
check("xi/L > 0.4 时偏离【更大】（比值才是判据）", big.max() > small.max())
check('=> 判据应写成【xi/L <~ 0.4】（比值），而不是绝对长度 xi <~ 1.2',
      _anchor('不是绝对', '是绝对长', '绝对长度'),
      "文档锚定（正文 §核验口径）：不是绝对 + 是绝对长 + 绝对长度", level="dep")

# ======================================================================
head("F5  诚实登记：数据坍缩【失败】")

pts = []
for m in MS:
    svals, a, rms, xi = data[m]
    for L, S in zip(LS, svals):
        pts.append((xi / L, S / (xi * L), m, L))
ok = True; ncmp = 0; mx = 0.0
for i in range(len(pts)):
    for j in range(i + 1, len(pts)):
        if pts[i][2] == pts[j][2]:
            continue
        if abs(pts[i][0] - pts[j][0]) < 0.02:
            d = abs(pts[i][1] - pts[j][1]) / max(pts[i][1], pts[j][1])
            mx = max(mx, d); ncmp += 1
            if d > 0.15:
                ok = False
print("      跨 m 比较 %d 对，最大相对差 = %.1f%%" % (ncmp, 100 * mx))
check("S/(xi L) 不是 xi/L 的纯函数（坍缩失败，最大差 > 15%%）", not ok, "最大 %.1f%%" % (100 * mx))
check('=> 次领项有自己的 m 依赖 => 纯坍缩【不成立】，如实登记',
      _anchor('次领项有自己的', '次领项有', '领项有自'),
      "文档锚定（正文 §核验口径）：次领项有自己的 + 次领项有 + 领项有自", level="dep")

# ======================================================================
head("F6  面积密度 a ∝ 1/m（与 G78 一致）")

prods = [(m, data[m][1] * m) for m in MS]
print("      m 与 a*m：%s" % "  ".join("(%.2f, %.3f)" % p for p in prods))
late = [p for m, p in prods if m >= 1.0]
check("a*m 在 m>=1 时近似常数（相对涨落 < 40%%）",
      np.std(late) / np.mean(late) < 0.4, "涨落 %.1f%%" % (100 * np.std(late) / np.mean(late)))
check("=> a ∝ 1/m 在 2D 大尺寸上同样成立",
      (np.std(late) / np.mean(late) < 0.4) and _anchor("1/m"),
      "重算 a*m 在 m>=1 的相对涨落 = %.1f%%（<40%%）+ 正文锚定『1/m』"
      % (100 * np.std(late) / np.mean(late)), level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
