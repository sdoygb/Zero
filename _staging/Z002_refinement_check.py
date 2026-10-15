#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z002_refinement_check.py —— Z002 §4.3 的复现与核验（细化收敛）

【核的是什么】
§4.3 的"有效度规"是截断闭游走计数（路线二）：
    k = n（被逼出的尺度）
    w_k(A,k)_ij = Σ_{m=2..k} m · A_ij · (A^{m-1})_ij        （注意是 A^{m-1}，同 §4.1 的 (4.2)）
本脚本核两件事：
  (1) 复现旧库 G46_check.py 的表（图 = path(80)，边取 0-based (1,2) 与 (30,31)）；
  (2) 判"收敛与否"：谱分解给出极限（Perron 形式），并给出偏离的衰减长度 τ ∝ N²，
      以及极限的非局域性（只改远端即改变极限）。

【断言（本脚本核验的）】
  A1  复现 G46 表：k=8,16,32,64 → 0.8248587571 / 0.5492835950 / 0.3193997518 / 0.1727912564
  A2  谱路径与精确矩阵幂路径逐位一致（k ≤ 128，差 < 1e-12）
  A3  比值单调下降且 → Perron 比 v_i v_j /(v_k v_l) = 0.0101924429（k=16384 时差 < 1e-12）
  A4  衰减长度 τ ≈ 2N²/3π²（渐近式；N ≥ 80 时相对误差 < 3%）
  A5  极限是全局量（非局域）：N: 80→82 改变极限比 > 1%

用法：python3 _staging/Z002_refinement_check.py     （退出码 0 = 全过）
"""
import math
import sys

import numpy as np

FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    ok = bool(cond)
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("   " + detail) if detail else ""))
    if not ok:
        FAIL += 1


def head(t):
    print("\n" + "=" * 72)
    print(t)
    print("=" * 72)


# ────────────────────── 图与权重（与 G46_check.py 逐字一致） ──────────────────────
def path(N):
    A = np.zeros((N, N))
    for i in range(N - 1):
        A[i, i + 1] = A[i + 1, i] = 1.0
    return A


def w_k(A, k):
    """长度 ≤ k 的闭合游走对每条边的总穿越数（G46_check.py 的定义，逐字照抄）。

    注意 P 的更新次序使 P = A^{m-1}，故权重是 Σ_m m·A_ij·(A^{m-1})_ij（= §4.1 的 (4.2)）。
    """
    W = np.zeros_like(A, dtype=float)
    P = np.eye(A.shape[0])
    for m in range(2, k + 1):
        P = P @ A
        W += m * A * P
    return W


def ratio_matrix(N, k, e1=(1, 2), e2=(30, 31)):
    """精确路径（矩阵幂）。k 大时 2^k 溢出，故只用于小 k。"""
    W = w_k(path(N), k)
    return float(W[e1] / W[e2])


# ────────────────────── 谱路径（可推到 k = 1e4） ──────────────────────
def spectrum(N):
    """P_N 的谱（解析）：λ_j = 2cos(jπ/(N+1))，v^(j)_i = sqrt(2/(N+1)) sin(jπ(i+1)/(N+1))。"""
    lam = np.array([2.0 * math.cos(j * math.pi / (N + 1)) for j in range(1, N + 1)])
    V = np.array([[math.sqrt(2.0 / (N + 1)) * math.sin(j * math.pi * (i + 1) / (N + 1))
                   for i in range(N)] for j in range(1, N + 1)])
    return lam, V


def gt_direct(L, lam):
    """g̃_L(λ) = Σ_{m=2..L} m λ^{m-1}（直接求和；L ≤ 200 时不溢出）。"""
    return sum(m * lam ** (m - 1) for m in range(2, L + 1))


def gt_ratio_big(L, lam, lam1):
    """g̃_L(λ)/g̃_L(λ1) 的大 L 近似。

    精确分子 = λ^{L}(Lλ-(L+1)) + 1（因 g̃ = g/λ，g 的分子 = λ^{L+1}(Lλ-(L+1)) + λ）；
    L ≥ 300 且 |λ|~2 时末项可忽略。g̃ = g/λ 引入的额外因子 λ1/λ 一并写进去。
    """
    if abs(lam) < 1e-13:
        return 0.0
    c = ((L * lam - (L + 1)) / (L * lam1 - (L + 1))) * ((1 - lam1) ** 2) / ((1 - lam) ** 2)
    return (lam / lam1) ** (L + 1) * (lam1 / lam) * c


def ratio_spectral(N, L, e1=(1, 2), e2=(30, 31)):
    """谱分解：W_ij = Σ_j (P_j)_ij g̃_L(λ_j)，比值以 g̃_L(λ1) 归一。"""
    lam, V = spectrum(N)
    P1 = np.outer(V[0], V[0])

    def num(pair):
        if L <= 200:
            return sum(V[j][pair[0]] * V[j][pair[1]] * gt_direct(L, lam[j]) for j in range(N))
        tot = P1[pair]
        for j in range(1, N):
            r = gt_ratio_big(L, lam[j], lam[0])
            if abs(r) < 1e-18:
                continue
            tot += V[j][pair[0]] * V[j][pair[1]] * r
        return tot

    return float(num(e1) / num(e2))


# ══════════════════════════════════ 核验 ══════════════════════════════════
head("A1  复现 G46 表（图 = path(80)，边 = 0-based (1,2) 与 (30,31)）")
G46 = {8: 0.8248587571, 16: 0.5492835950, 32: 0.3193997518, 64: 0.1727912564}
print("      k      本脚本              G46_check.py")
for k, want in G46.items():
    got = ratio_matrix(80, k)
    print("      %-7d %.10f      %.10f" % (k, got, want))
    check("A1 k=%d 复现" % k, abs(got - want) < 1e-9, "差 %.2e" % abs(got - want))

head("A2  两条路径逐位一致（谱 + 闭式 vs 精确矩阵幂）")
for k in (8, 16, 32, 64, 128):
    a, b = ratio_matrix(80, k), ratio_spectral(80, k)
    check("A2 k=%-4d" % k, abs(a - b) < 1e-12, "矩阵=%.10f 谱=%.10f 差=%.1e" % (a, b, abs(a - b)))

head("A3  收敛：极限 = Perron 比 v_i v_j /(v_k v_l)")
lam, V = spectrum(80)
per = float(np.outer(V[0], V[0])[(1, 2)] / np.outer(V[0], V[0])[(30, 31)])
print("      λ1 = %.10f   λ2/λ1 = %.10f   1/(-log(λ2/λ1)) = %.1f 步"
      % (lam[0], lam[1] / lam[0], 1.0 / -math.log(lam[1] / lam[0])))
print("      Perron 极限比 = %.10f" % per)
seq = []
for L in (64, 256, 1024, 2048, 4096, 8192, 16384):
    r = ratio_spectral(80, L)
    seq.append(r)
    print("      k=%-6d 比值=%.10f  与极限之差=%.3e" % (L, r, abs(r - per)))
check("A3 单调下降", all(seq[i] > seq[i + 1] for i in range(len(seq) - 1)))
check("A3 k=4096 差 < 1e-5", abs(ratio_spectral(80, 4096) - per) < 1e-5)
check("A3 k=16384 差 < 1e-12", abs(ratio_spectral(80, 16384) - per) < 1e-12)

head("A4  衰减长度 τ ≈ 2N²/3π²（细化使收敛更慢；该式为渐近式）")
for N in (20, 40, 80, 160, 320, 640):
    l1 = 2 * math.cos(math.pi / (N + 1))
    l2 = 2 * math.cos(2 * math.pi / (N + 1))
    tau = 1.0 / -math.log(l2 / l1)
    pred = 2 * N * N / (3 * math.pi ** 2)
    rel = abs(tau - pred) / tau
    print("      N=%-4d τ=%.1f    2N²/3π²=%.1f    相对差=%.2f%%" % (N, tau, pred, 100 * rel))
    tol = 0.03 if N >= 80 else 0.10
    check("A4 N=%-3d τ≈2N²/3π²（渐近，容差 %.0f%%）" % (N, 100 * tol), rel < tol)

head("A5  极限是全局量（非局域）：只改远端就改变极限比")
r = {}
for N in (80, 82, 100):
    Vv = spectrum(N)[1]
    p = np.outer(Vv[0], Vv[0])
    r[N] = float(p[(1, 2)] / p[(30, 31)])
    print("      N=%-4d 极限比 = %.10f" % (N, r[N]))
check("A5 N:80→82 改变极限比 > 1%", abs(r[80] - r[82]) / r[80] > 0.01,
      "%.5f → %.5f" % (r[80], r[82]))

head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
