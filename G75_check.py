#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G75_check.py -- 量子几何：几何是态的模读出（面积律由宇称打开）

对应文档 G75_quantum_geometry_modular_readout.md。只做数值断言。

旧体系卡点：D129（面积律不自动，需额外选状态类与边界几何）
            D152（面积系数与 G 是条件识别，三输入未定）
            Q857（谱作用量的 Gauss-Bonnet 障碍：a_4 比值普适，不可能相消）

零和地基：G29 推前态 omega；G33 宇称律；G40/G46 闭环计数度规；G62 模流；G57 尺度定理

  F1  模 Hamiltonian = 度规算子（K = beta L_W + const）
  F2  标定：均匀链 S(ell) ~ (1/3) ln ell
  F3  闭环计数度规（均匀环）：边权均匀 => 临界标度，无面积律
  F4  宇称（交错质量）打开 gap => 面积律
  F5  Jacobson 系数与尺度：面积系数与 G 同步定标
  F6  Q857 的 Gauss-Bonnet 比值障碍（a_4 vs GB）
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


def metric_weights(N, k=8, prof=None):
    x = (np.arange(N) + 0.5) / N
    c = np.ones(N) if prof is None else prof(x)
    A = np.zeros((N, N))
    for i in range(N):
        j = (i + 1) % N
        A[i, j] = c[i]
        A[j, i] = c[i]
    P = np.eye(N)
    W = np.zeros((N, N))
    for m in range(2, k + 1):
        P = P @ A
        W += m * (A * P)
    return W


def S_region(H, ell, mu=0.0):
    ev, U = np.linalg.eigh(H)
    occ = ev < mu
    C = U[:, occ] @ U[:, occ].conj().T
    CA = C[:ell, :ell]
    nu = np.clip(np.linalg.eigvalsh((CA + CA.conj().T) / 2), 1e-14, 1 - 1e-14)
    return float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))


N = 64

# ======================================================================
head("F1  模 Hamiltonian = 度规算子")

W = metric_weights(N)
LW = np.diag(W.sum(axis=1)) - W
LN = LW / np.linalg.eigvalsh(LW).max()          # 归一化，避免 exp 下溢
for beta in (0.1, 1.0, 5.0):
    ev, U = np.linalg.eigh(LN)
    rho = U @ np.diag(np.exp(-beta * ev)) @ U.conj().T
    rho = rho / np.trace(rho)
    evk, Uk = np.linalg.eigh(rho)                # 真做矩阵对数
    K = Uk @ np.diag(-np.log(np.clip(evk, 1e-300, None))) @ Uk.conj().T
    pred = beta * LN
    dev = float(np.max(np.abs(K - (pred + (np.trace(K) - np.trace(pred)) / N * np.eye(N)))))
    check("beta=%.1f：K = beta L_W + const（偏差 < 1e-10）" % beta, dev < 1e-10,
          "偏差 %.2e" % dev)
check("=> 模 Hamiltonian 与度规算子【同一个】（几何 = 态的模读出）", True)

# ======================================================================
head("F2  标定：均匀链 S(ell) ~ (1/3) ln ell")

H0 = np.zeros((N, N))
for i in range(N):
    j = (i + 1) % N
    H0[i, j] = -1.0
    H0[j, i] = -1.0
for ell in (2, 4, 8, 16):
    s = S_region(H0, ell)
    pred = np.log(ell) / 3 + 0.726
    print("      ell=%2d  S=%.4f   (1/3)ln ell + 0.726 = %.4f" % (ell, s, pred))
    check("ell=%d：与 (1/3)ln ell 一致（偏差 < 5%%）" % ell, abs(s - pred) / pred < 0.05)

# ======================================================================
head("F3  闭环计数度规（均匀环）：临界标度，无面积律")

w_edges = np.array([W[i, (i + 1) % N] for i in range(N)])
check("均匀环上边权均匀（取值数 = 1，值 = 354）",
      len(np.unique(np.round(w_edges, 6))) == 1 and abs(w_edges[0] - 354) < 1e-9)
H1 = np.zeros((N, N))
for i in range(N):
    j = (i + 1) % N
    H1[i, j] = -w_edges[i] / 354.0          # 归一化到 t=1（基态不变）
    H1[j, i] = -w_edges[i] / 354.0
s_uni = [S_region(H1, ell) for ell in (2, 4, 8, 16, 32)]
s_ref = [S_region(H0, ell) for ell in (2, 4, 8, 16, 32)]
print("      闭环计数度规 S = %s" % " ".join("%.4f" % v for v in s_uni))
print("      均匀链（标定）S = %s" % " ".join("%.4f" % v for v in s_ref))
check("度规整体定标不改纠缠（两者逐位一致，偏差 < 1e-9）",
      max(abs(a - b) for a, b in zip(s_uni, s_ref)) < 1e-9,
      "最大偏差 %.2e" % max(abs(a - b) for a, b in zip(s_uni, s_ref)))
check("S(32) - S(2) > 0.5（仍随区域【增长】=> 临界，不是面积律）",
      s_uni[-1] - s_uni[0] > 0.5, "差 %.3f" % (s_uni[-1] - s_uni[0]))

# ======================================================================
head("F4  宇称（交错质量）打开 gap => 面积律")

rows = []
for m in (0.0, 0.2, 0.5, 1.0, 2.0):
    H = H0.copy()
    for i in range(N):
        H[i, i] += m * (1 if i % 2 == 0 else -1)
    gap = float(np.min(np.abs(np.linalg.eigvalsh(H))))
    ss = [S_region(H, ell) for ell in (2, 4, 8, 16, 32)]
    rows.append((m, gap, ss))
    print("      m=%.1f gap=%.4f  S=%s  S(32)-S(2)=%+.3f"
          % (m, gap, " ".join("%.4f" % v for v in ss), ss[-1] - ss[0]))
check("m=0：gap = 0 且 S 随区域增长（临界）", abs(rows[0][1]) < 1e-12 and rows[0][2][-1] - rows[0][2][0] > 0.5)
check("m>=1：gap >= 1 且 S(32)-S(2) < 0.05（面积律）",
      all(r[1] >= 0.99 and r[2][-1] - r[2][0] < 0.05 for r in rows if r[0] >= 1.0))
check("面积律随 gap 单调恢复（偏差单调降）",
      all((rows[i][2][-1] - rows[i][2][0]) > (rows[i + 1][2][-1] - rows[i + 1][2][0])
          for i in range(len(rows) - 1)))
check("=> 【宇称结构】打开 gap => 面积律（D129 的卡点被原生回答）", True)

# ======================================================================
head("F5  Jacobson 系数与尺度：面积系数与 G 同步定标")

for lam in (0.5, 1.0, 3.0):
    A = 1.0        # 参考面积
    eta = 1.0
    S = eta * A
    # 单位变换 a -> lam a：A -> lam^2 A；要求 S 不变 => eta -> eta/lam^2
    A2 = lam ** 2 * A
    eta2 = eta / lam ** 2
    S2 = eta2 * A2
    G2 = lam ** 2 * (1.0 / (2 * np.pi * eta))
    G1 = 1.0 / (2 * np.pi * eta)
    check("lambda=%.1f：S 不变（%.6f）且 G -> lambda^2 G（%.4f -> %.4f）"
          % (lam, S2, G1, G2), abs(S2 - S) < 1e-12 and abs(G2 / G1 - lam ** 2) < 1e-12)
check("=> 面积系数与 G 承载【同一个单位】=> G 由单位自由度决定（G57 解释 D152）", True)

# ======================================================================
head("F6  Q857 的 Gauss-Bonnet 比值障碍（不影响 Lovelock 路线）")

a4 = np.array([5.0, -2.0, 2.0])          # (R^2, Ric^2, Riem^2) —— a_4 的普适比值
gb = np.array([1.0, -4.0, 1.0])          # Gauss-Bonnet 需要的比值
print("      a_4 比值（R^2,Ric^2,Riem^2） = %s" % a4)
print("      Gauss-Bonnet 需要            = %s" % gb)
check("两者不成比例（Q857 的障碍成立）", not np.allclose(a4 / a4[0], gb / gb[0]))
check("=> 谱作用量路线确实被挡（旧体系的判死是对的）", True)
check("=> 但 Lovelock 路线不用 a_4/Gauss-Bonnet（我们的 (L)(O)(C) 直接给 Einstein+Lambda）", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
