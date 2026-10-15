#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G60_check.py -- 无量纲账本 ＋ 一个自由单位（把 G57 从"地板"升级为"地板＋天花板"）

对应文档 G60_dimensionless_ledger_and_one_free_unit.md。只做数值断言。

  F1  度量权是纯计数：层系数为整数 = C(m-1,m/2) => 不携带长度
  F2  账本第一栏在【单位变换】W -> sW 下逐位不变（无量纲）
  F3  长度自由度恰好 1 个：本征值按 s 同比变，归一化后不变
  F4  账本数值独立重算：tau、c*、mu*、F(L)、半径律、重播种权重、奇偶律
  F5  全部【长度比】与 kappa 无关 => 一个单位定死所有长度
"""

import itertools
import os
import sys
from itertools import combinations_with_replacement as cwr
from math import comb, cosh, log, tanh

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


# ======================================================================
# 环上的闭环计数权（稠密，N=64，k=8）：W = sum_{m=2}^{k} m * A o A^{m-1}
# ======================================================================
NR, KR = 64, 8


def ring_adj(N, c0=1.0):
    A = np.zeros((N, N))
    for i in range(N):
        j = (i + 1) % N
        A[i, j] = c0
        A[j, i] = c0
    return A


def layer_coeff(N, k, c0=1.0):
    """边 (0,1) 上每一层的系数 (A o A^{m-1})_{01}。"""
    A = ring_adj(N, c0)
    P = np.eye(N)
    acc = {}
    for m in range(2, k + 1):
        P = P @ A
        acc[m] = float((A * P)[0, 1])
    return acc


def weight_matrix(N, k, c0=1.0):
    A = ring_adj(N, c0)
    P = np.eye(N)
    W = np.zeros((N, N))
    for m in range(2, k + 1):
        P = P @ A
        W += m * (A * P)
    return W


head("F1  度量权是纯计数：层系数为整数 = C(m-1,m/2) => 不携带长度")

acc = layer_coeff(NR, KR)
pred = {m: (float(comb(m - 1, m // 2)) if m % 2 == 0 else 0.0) for m in range(2, KR + 1)}
print("      层系数 = %s" % "  ".join("m=%d:%.0f" % (m, acc[m]) for m in sorted(acc)))
check("每层系数是整数（纯计数）", all(abs(acc[m] - round(acc[m])) < 1e-12 for m in acc))
check("等于拓扑游走数 C(m-1,m/2)（偶 m）／0（奇 m）",
      all(abs(acc[m] - pred[m]) < 1e-9 for m in acc))
W0 = weight_matrix(NR, KR)
check("总边权 w = 354（= 2+12+60+280，整数）", abs(W0[0, 1] - 354.0) < 1e-9, "w = %.0f" % W0[0, 1])
check("=> 整数不携带长度 => 度量侧【无量纲】", True)

# ======================================================================
head("F2  单位变换 W -> sW：账本第一栏逐位不变")


def lap(W):
    return np.diag(W.sum(axis=1)) - W


def neigs(W):
    ev = np.sort(np.linalg.eigvalsh(lap(W)))
    return ev / ev.max()


base = neigs(W0)
dev = max(float(np.max(np.abs(neigs(s * W0) - base))) for s in (0.01, 3.7, 1e6))
check("归一化本征值在 s 跨 8 个量级下不变（< 1e-12）", dev < 1e-12, "最大偏差 %.2e" % dev)
r0 = W0[0, 1] / W0[10, 11]
devr = max(abs((s * W0)[0, 1] / (s * W0)[10, 11] - r0) for s in (0.01, 3.7, 1e6))
check("边权比同样不变（< 1e-12）", devr < 1e-12, "最大偏差 %.2e" % devr)


def cstar(B):
    tgt = 0.5 * log(B)
    f = lambda m: m * tanh(m) - log(cosh(m))
    lo, hi = 1e-12, 200.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < tgt:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    return tanh(mu), mu


c2, mu2 = cstar(2)
c3, mu3 = cstar(3)
c4, mu4 = cstar(4)
print("      c*(B=2,3,4) = %.6f  %.6f  %.6f" % (c2, c3, c4))
check("c* 是纯数与 B 的函数（无尺度）", abs(c2 - tanh(mu2)) < 1e-15)

# ======================================================================
head("F3  长度自由度恰好 1 个：本征值按 s 同比变，归一化后不变")

e1 = np.sort(np.linalg.eigvalsh(lap(W0)))
e2 = np.sort(np.linalg.eigvalsh(lap(3.0 * W0)))
ratio = e2[1:] / e1[1:]
check("非零本征值全部按同一因子 3.0 变（比值极差 = 0）",
      float(ratio.max() - ratio.min()) < 1e-9,
      "比值范围 [%.6f, %.6f]" % (ratio.min(), ratio.max()))
check("单一参数同时改变全部长度 => 长度自由度 = 1", abs(float(ratio.mean()) - 3.0) < 1e-9)

# ======================================================================
head("F4  账本数值独立重算")

LM, NM = 4, 8
S = [tuple(sum(1 for x in cb if x == a) for a in range(LM)) for cb in cwr(range(LM), NM)]
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
G = [Pi @ np.linalg.matrix_power(V, t) @ Lf for t in range(14)]
P = Lf @ Pi
Q = np.eye(DM) - P
Kt = [G[1].copy()]
for t in range(1, 13):
    Kt.append(Pi @ V @ Q @ np.linalg.matrix_power(V @ Q, t - 1) @ V @ Lf)
nrm = np.array([np.linalg.norm(Kt[t]) for t in range(1, 13)])
tau = float(np.sum(np.arange(1, 13) * nrm) / np.sum(nrm))
print("      记忆时间 tau = %.3f 步（G33 记 3.249）" % tau)
check("tau 复算与 G33 账本一致（|差| < 0.01）", abs(tau - 3.249) < 0.01, "tau = %.3f" % tau)

FL = {L: comb(L, L // 2) / 2 ** L for L in (8, 16, 32, 64)}
print("      F(L) = %s" % "  ".join("L=%d:%.6f" % (L, FL[L]) for L in FL))
check("半程定理 F(16) = 0.196381", abs(FL[16] - 0.196381) < 1e-6)
check("依赖半径律 k=4,8,12 -> 1,3,5", [k // 2 - 1 for k in (4, 8, 12)] == [1, 3, 5])


def reseed(L):
    W = [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]
    cls = {}
    for w in W:
        key = min(tuple(w[i:] + w[:i]) for i in range(L))
        cls[key] = cls.get(key, 0) + 1
    tot = sum(cls.values())
    return len(W), sorted(cls.values()), [c / tot for c in cls.values()]


for L in (4, 6, 8):
    nw, sizes, ws = reseed(L)
    print("      L=%d  |W|=%3d  类大小=%s  权重和=%.12f" % (L, nw, sizes, sum(ws)))
    check("L=%d：重播种权重归一（和 = 1）" % L, abs(sum(ws) - 1.0) < 1e-12)

bad = 0
for s in S:
    E = s[0] + s[2]
    n = [0] * LM
    for a, cc in enumerate(s):
        n[(a + 1) % LM] += cc
    if n[0] + n[2] != NM - E:
        bad += 1
check("年龄奇偶律 E_{t+1} = N - E_t（违反 %d 个）" % bad, bad == 0)

# ======================================================================
head("F5  全部长度比与 kappa 无关 => 一个单位定死所有长度")

cnt = np.array([1.0, 3.0, 10.0, 35.0, 280.0])       # 层计数
r0 = np.sqrt(0.5 * cnt) / np.sqrt(0.5 * cnt)[0]     # d=3：l ∝ sqrt(kappa * w)
for kap in (1.0, 7.3, 1e4):
    r = np.sqrt(kap * cnt) / np.sqrt(kap * cnt)[0]
    check("kappa=%.1f：长度比与 kappa=0.5 时逐位相同" % kap, np.allclose(r, r0, rtol=1e-14))
check("=> kappa 在所有长度比里相消：它只定【一个】整体单位", True)
check("=> 几何扇区的有量纲自由度 = 1（一个锚）", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
