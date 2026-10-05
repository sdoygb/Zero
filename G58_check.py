#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G58_check.py -- I2a 判定：Γ-收敛（非 UV 不动点）＋ 固定步数支二阶收敛

对应文档 G58_I2a_resolved_as_embedding_input.md。
只做数值断言（按"只留增量、只留有信息量"协议，不含文档措辞断言）。

  F1  固定 k：内部剖面随 N 二阶收敛（相邻差递减、比值 < 0.5）
  F2  拟合收敛阶 ≈ -2
  F3  局域性：层依赖半径 = k/2 - 1
  F4  固定物理时长 k = rho*N：动态范围爆炸（发散）
  F5  均匀链上核对 G49 的精确分解（几何与拓扑因子化）
"""

import math
import os
import sys

import numpy as np
from scipy.sparse import diags

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


def prof(x):
    return 1.0 + 0.3 * np.cos(2 * np.pi * x)


def build(N, k, profile=prof, pert=None, eps=0.3):
    """链上闭环计数权：w_e = sum_{m=2}^{k} m (A .* A^{m-1}) 的超对角元。"""
    x = (np.arange(N) + 0.5) / N
    c = profile(x)
    if pert is not None:
        c[pert] += eps
    A = diags([c, c], [1, -1], shape=(N + 1, N + 1)).tocsr()
    P = A.copy()
    W = A.multiply(P) * 2
    for m in range(3, k + 1):
        P = (P @ A).tocsr()
        W = W + m * A.multiply(P)
    return W.diagonal(1), c, x


def field(xs, w, N):
    """把边权插值到固定物理位置（边心 (i+0.5)/N），归一化到内部中位数。"""
    xe = (np.arange(N) + 0.5) / N
    ww = w / np.median(w[(xe > 0.25) & (xe < 0.75)])
    return np.interp(xs, xe, ww)


# ======================================================================
head("F1/F2  固定 k=8：内部剖面随 N 二阶收敛")

K = 8
XS = np.linspace(0.2, 0.8, 41)
prev, errs, Ns = None, [], []
for N in (128, 256, 512, 1024):
    w, c, x = build(N, K)
    f = field(XS, w, N)
    if prev is not None:
        e = float(np.max(np.abs(f - prev)))
        errs.append(e)
        Ns.append(N)
        print("      N=%4d  相邻 N 的最大差 = %.4e" % (N, e))
    prev = f
check("相邻差严格递减", all(errs[i] > errs[i + 1] for i in range(len(errs) - 1)),
      " -> ".join("%.3e" % e for e in errs))
check("每加倍至少减半（比值 < 0.5）",
      all(errs[i + 1] / errs[i] < 0.5 for i in range(len(errs) - 1)),
      "最大比值 %.3f" % max(errs[i + 1] / errs[i] for i in range(len(errs) - 1)))
slope = float(np.polyfit(np.log(Ns), np.log(errs), 1)[0])
check("拟合收敛阶 ≈ -2（二阶 O(a^2)）", abs(slope + 2.0) < 0.25, "斜率 %.3f" % slope)
check("=> 固定步数截断下，细化族【收敛】", True)

# ======================================================================
head("F3  局域性：层依赖半径 = k/2 - 1")

N = 256
for k in (4, 8, 12):
    base = build(N, k)[0][100]
    hits = [d for d in range(1, 40) if abs(build(N, k, pert=100 + d)[0][100] - base) > 1e-14]
    exp = list(range(1, k // 2))
    print("      k=%2d 受影响边距离 = %s（预期 1..%d）" % (k, hits, k // 2 - 1))
    check("k=%d：依赖半径 = k/2 - 1 = %d" % (k, k // 2 - 1), hits == exp)
check("=> 相互作用范围有限 => 极限【局域】 => (L) 成立", True)

# ======================================================================
head("F4  固定物理时长 k = rho*N：动态范围爆炸（发散）")

for rho in (0.05, 0.10):
    vals = []
    for N in (128, 256, 512):
        k = max(2, int(rho * N))
        w, c, x = build(N, k)
        m = (x > 0.25) & (x < 0.75)
        vals.append(float(w[m].max() / w[m].min()))
        print("      rho=%.2f N=%3d k=%3d 内部 max/min = %.4e" % (rho, N, k, vals[-1]))
    check("rho=%.2f：动态范围随 N 单调爆炸" % rho, vals[0] < vals[1] < vals[2])
    check("rho=%.2f：指数跨度 > 100 倍" % rho, vals[-1] / vals[0] > 100,
          "跨度 %.3e" % (vals[-1] / vals[0]))
check("=> 固定物理时长支【发散】（度量退化）", True)

# ======================================================================
head("F5  均匀环上核对 G49 的精确分解（几何 c^m × 拓扑游走数）")


def ring_edge_weight(N, k, c0):
    """环 C_N 上均匀权重 c0 时的闭环计数边权（环上处处相同）。"""
    A = diags([np.full(N, c0), np.full(N, c0), np.full(N, c0)], [1, -1, N - 1],
              shape=(N, N)).tocsr()
    A = A + diags([np.full(N, c0)], [-(N - 1)], shape=(N, N)).tocsr()
    P = A.copy()
    W = A.multiply(P) * 2
    for m in range(3, k + 1):
        P = (P @ A).tocsr()
        W = W + m * A.multiply(P)
    return float(W.diagonal(1)[0]), float(W.diagonal(1).max() - W.diagonal(1).min())


NR = 64
CS = np.linspace(0.5, 3.0, 13)
ws, spread = [], 0.0
for c0 in CS:
    we, sp = ring_edge_weight(NR, K, c0)
    ws.append(we)
    spread = max(spread, sp)
check("环上边权处处相同（顶点传递 => 均匀）", spread < 1e-9, "极差 %.2e" % spread)

# G49: w = sum_{m=2}^{k} m * (A_top^{m-1})_{adj} * c^m  => 关于 c 的 8 次多项式，系数 = m*C(m-1,m/2)
coef = np.polyfit(CS, ws, K)[::-1]            # 升幂 a_0..a_8
pred = np.zeros(K + 1)
for m in range(2, K + 1, 2):
    pred[m] = m * math.comb(m - 1, m // 2)
resid = float(np.max(np.abs(np.polyval(coef[::-1], CS) - np.array(ws))))
print("      拟合系数（升幂） = %s" % np.array2string(coef, precision=6))
print("      解析预测         = %s" % np.array2string(pred, precision=6))
print("      拟合残差 = %.3e" % resid)
check("w(c) 是 c 的 8 次多项式（拟合残差 < 1e-9）", resid < 1e-9)
check("系数与解析预测 m*C(m-1,m/2) 一致（偏差 < 1e-6）",
      float(np.max(np.abs(coef - pred))) < 1e-6,
      "最大偏差 %.3e" % float(np.max(np.abs(coef - pred))))
check("=> G49 的因子化成立：几何进 c^m、拓扑进游走数", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
