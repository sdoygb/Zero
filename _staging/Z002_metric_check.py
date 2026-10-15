#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z002_metric_check.py —— Z002 §4.1／§4.2 的复现与核验（闭环计数度规）

【核的是什么】
路线一/二的权重：w_ij ∝ φ_i φ_j A_ij（φ = Perron 向量），
其来源是长度 k 的闭合游走对边 (i,j) 的总穿越数
    N^(k)_ij = k · A_ij · (A^{k-1})_ij .
本脚本复现 §4.1／§4.2 的四个数字组：
  A1  §4.1 F2：归一化穿越数 → φφA 的收敛偏差（k=6,12,24,48,96）
  A2  §4.1 F3/F4：正则图上边权精确均匀；不规则图上非均匀
  A3  §4.1 非局域性：远端扰动对中部边权比的影响（N=60 链）
  A4  §4.1 非局域性：影响不随距离衰减
  A5  §4.2 (O)/(C)：加权 Laplacian 二阶性（max|L·1| 小）与恰有 1 个零本征值
  A6  §4.2 k 扫描：截断超过直径后权重开始随 k 变（非局域的起始）

约定与旧库 G40/G41_check.py 一致（图集、perron、cw_weight、wlap 逐字照抄）。

用法：python3 _staging/Z002_metric_check.py      （退出码 0 = 全过）
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


# ────────────────────── 工具（照 G40/G41_check.py） ──────────────────────
def perron(A):
    w, v = np.linalg.eigh(A)
    p = np.abs(v[:, -1])
    return p / np.linalg.norm(p), float(w[-1])


def cw_weight(A):
    phi = perron(A)[0]
    return phi[:, None] * phi[None, :] * A


def path(N, pert=None):
    A = np.zeros((N, N))
    for i in range(N - 1):
        e = 1.0 if (pert is None or i != pert[0]) else pert[1]
        A[i, i + 1] = A[i + 1, i] = e
    return A


def wlap(W):
    return np.diag(W.sum(axis=1)) - W


def graphs():
    G = {}
    G["K4"] = np.ones((4, 4)) - np.eye(4)
    m = 6
    C = np.zeros((m, m))
    for i in range(m):
        C[i, (i + 1) % m] = C[(i + 1) % m, i] = 1
    G["C6"] = C
    A = np.zeros((5, 5))
    for i, j in [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4)]:
        A[i, j] = A[j, i] = 1
    G["irreg"] = A
    S = np.zeros((5, 5))
    for i in range(1, 5):
        S[0, i] = S[i, 0] = 1
    G["star"] = S
    return G


G = graphs()

head("A1  §4.1 F2：归一化穿越数 → φ_i φ_j A_ij 的收敛偏差（k=6,12,24,48,96）")
print("      N^(k)_ij = k·A_ij·(A^{k-1})_ij，按全矩阵和归一化；偏差 = 边上 |N̂ - φ̂φ̂A| 的最大值")
dev_tab = {}
for name, A in G.items():
    phi, lam = perron(A)
    pred = phi[:, None] * phi[None, :] * A
    pred = pred / pred.sum()
    devs = []
    for k in (6, 12, 24, 48, 96):
        Nk = k * A * np.linalg.matrix_power(A, k - 1)
        Nk = Nk / Nk.sum()
        devs.append(float(np.abs(Nk - pred)[A > 0].max()))
    dev_tab[name] = devs
    print("      %-6s λ=%.3f  偏差 %s" % (name, lam, " ".join("%.4f" % d for d in devs)))
check("A1 K4 全为 0（机器精度）", all(d < 1e-12 for d in dev_tab["K4"]))
check("A1 C6 全为 0（机器精度）", all(d < 1e-12 for d in dev_tab["C6"]))
check("A1 irreg 序列 ≈ 0.0201/0.0049/0.0002/0/0",
      all(abs(a - b) < 5e-5 for a, b in zip(dev_tab["irreg"], [0.0201, 0.0049, 0.0002, 0.0, 0.0])),
      "实测 %s" % ["%.4f" % d for d in dev_tab["irreg"]])
check("A1 irreg 的 λ ≈ 2.214", abs(perron(G["irreg"])[1] - 2.214) < 5e-3,
      "λ=%.3f" % perron(G["irreg"])[1])

head("A2  §4.1 F3/F4：正则图上精确均匀；不规则图上非均匀")
for name in ("K4", "C6", "star"):
    A = G[name]
    W = cw_weight(A)
    W = W / W.sum()
    vals = sorted(set(np.round(W[A > 0], 10).tolist()))
    print("      %-6s 边权取值 %s" % (name, [round(v, 6) for v in vals]))
    check("A2 %s 只有 1 个取值（精确均匀）" % name, len(vals) == 1)
check("A2 K4/C6 权 = 0.083333", abs(cw_weight(G["K4"]).sum() > 0 and
      (cw_weight(G["K4"]) / cw_weight(G["K4"]).sum())[G["K4"] > 0][0] - 0.0833333) < 1e-6)
check("A2 星图权 = 0.125",
      abs((cw_weight(G["star"]) / cw_weight(G["star"]).sum())[G["star"] > 0][0] - 0.125) < 1e-6)
A = G["irreg"]
Wi = cw_weight(A) / cw_weight(A).sum()
vals = sorted(set(np.round(Wi[A > 0], 8).tolist()))
print("      irreg  边权取值 %s" % [round(v, 4) for v in vals])
check("A2 irreg 非均匀且 ≈ 0.0239/0.0934/0.1116/0.1355",
      len(vals) == 4 and all(abs(a - b) < 5e-4 for a, b in
                             zip(sorted(vals), [0.0239, 0.0934, 0.1116, 0.1355])),
      "实测 %s" % [round(v, 4) for v in sorted(vals)])

head("A3  §4.1 非局域性：扰动远端边 (58,59) 的权 1→1+k，看中部边权比 r = w[29,30]/w[30,31]")
NN = 60
i1, i2 = NN // 2 - 1, NN // 2
print("      k       闭环计数 r            相对变化        均匀权相对变化")
rows = []
for k in (0.0, 0.02, 0.10, 0.50, 2.0):
    A0, A1 = path(NN), path(NN, pert=(NN - 2, 1.0 + k))
    W0, W1 = cw_weight(A0), cw_weight(A1)
    r0 = W0[i1, i1 + 1] / W0[i2, i2 + 1]
    r1 = W1[i1, i1 + 1] / W1[i2, i2 + 1]
    u0 = A0[i1, i1 + 1] / A0[i2, i2 + 1]
    u1 = A1[i1, i1 + 1] / A1[i2, i2 + 1]
    dcw, dun = abs(r1 - r0) / abs(r0), abs(u1 - u0) / abs(u0)
    rows.append((k, r1, dcw, dun))
    print("      %-7.2f %-22.10f %-15.3e %.3e" % (k, r1, dcw, dun))
check("A3 均匀权一列精确为 0（严格局域）", all(r[3] < 1e-15 for r in rows))
check("A3 闭环计数一列 ≈ 0 / 2.299e-4 / 1.428e-3 / 2.025e-1 / 8.753e-1",
      all(abs(a - b) < 2e-3 * max(1.0, b) for a, b in
          zip([r[2] for r in rows], [0.0, 2.299e-4, 1.428e-3, 2.025e-1, 8.753e-1])),
      "实测 %s" % ["%.3e" % r[2] for r in rows])
check("A3 影响随扰动幅度单调增长", all(rows[i][2] < rows[i + 1][2] for i in range(1, 4)))

head("A4  §4.1 非局域性：影响不随距离衰减")
A0, A1 = path(NN), path(NN, pert=(NN - 2, 1.5))
W0, W1 = cw_weight(A0), cw_weight(A1)
dec = []
for i in (NN - 3, NN - 8, NN - 18, NN - 28, NN // 2 - 1):
    d = (NN - 2) - i
    rel = abs(W1[i, i + 1] - W0[i, i + 1]) / W0[i, i + 1]
    dec.append((d, rel))
    print("      被测边 %-4d 距扰动 d=%-3d 相对变化 %.3e" % (i, d, rel))
far = [r for d, r in dec if d >= 20]
check("A4 远场（d≥20）仍为 O(1)（不衰减）", all(r > 0.1 for r in far),
      "远场 %s" % ["%.3f" % r for r in far])

head("A5  §4.2 (O)/(C)：截断权重（N=80，k=n）的加权 Laplacian")
A80 = path(80)


def w_trunc(k, N=80):
    A = path(N)
    W = np.zeros_like(A)
    P = np.eye(N)
    for m in range(2, k + 1):
        P = P @ A
        W += m * A * P
    return W


Wk = w_trunc(64)
L = wlap(Wk / np.abs(Wk).max())
rO = float(np.abs(L @ np.ones(80)).max())
ev = np.linalg.eigvalsh(L)
nz = int(np.sum(np.abs(ev) < 1e-9 * max(1.0, np.abs(ev).max())))
print("      max|L·1| = %.3e   （二阶性 (O)）" % rO)
print("      零本征值个数 = %d   （守恒源 (C)）" % nz)
check("A5 (O) max|L·1| ≤ 1e-15", rO <= 1e-15, "%.2e" % rO)
check("A5 (C) 恰有 1 个零本征值", nz == 1, "个数 = %d" % nz)

head("A6  §4.2 局域性：远端扰动在 k ≤ 直径时对中部权的影响为精确 0")
A80p = path(80)
A80q = path(80, pert=(58, 1.5))
print("      扰动远端边 (58,59) 1→1.5，看中部边 (1,2) 的权相对变化")
zero_ok = True
for k in (8, 16, 32, 64, 79):
    W0 = w_trunc(k)
    Wq = np.zeros_like(A80q)
    P = np.eye(80)
    for m in range(2, k + 1):
        P = P @ A80q
        Wq += m * A80q * P
    rel = abs(Wq[1, 2] - W0[1, 2]) / W0[1, 2]
    print("      k=%-4d 相对变化 = %.3e" % (k, rel))
    zero_ok = zero_ok and rel == 0.0
check("A6 k ≤ 79（＝直径）时影响为精确 0（(L') 在此范围内成立）", zero_ok)

head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
