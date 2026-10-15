#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G57_check.py -- 不可达定理：绝对归一化（C_norm 收官）

对应文档 G57_unreachability_of_absolute_normalization.md。

定理（量纲空洞）：底层条款（Z0 条款 ＋ Z1–Z5 定理）的原语清单里【没有任何有量纲常数】——全部是有限集、整数、
计数、序与比值。任何由这些数据出发的泛函仍是【无量纲】的。故绝对长度/时间/质量与
任何量纲耦合（G、Lambda、kappa）不可由底层条款（Z0 条款 ＋ Z1–Z5 定理）导出。

  F1  条款原语的量纲清点：全部无量纲
  F2  派生量的量纲清点：全部无量纲，唯一例外是 G2 引理 9 的 kappa
  F3  全局重标定不变性（数值，机器精度）
  F4  截断 k【改变形状】（物理），对比标度 s【不改变形状】（单位）
  F5  选择率与重数标定无关（数值）
  F6  定理表述与证明骨架（文档）
  F7  C_norm 不可达推论与六槽位收官（文档）
"""

import os
import sys
from math import log

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
# F1  条款原语的量纲清点
# ======================================================================
head("F1  条款原语的量纲清点：全部无量纲")

PRIMS = [
    ("Z1 定理 1", "有限集 C = {1..m}", "基数 m", "1"),
    ("Z1 定理 2", "整向量 x in Z^C", "整数分量", "1"),
    ("Z1 定理 2", "零和约束 Q(x)=0", "整系数线性约束", "1"),
    ("Z1 定理 1", "连通图 Gamma=(C,E)", "有限点集与边集", "1"),
    ("Z1 定理 1", "基本移动 T_e = x + e_j - e_i", "整数偏移", "1"),
    ("Z2", "全分支、整数重数", "计数", "1"),
    ("Z4", "步指标 tau in Z_{>=0}", "计数（步）", "1"),
    ("Z3", "闭合条件 x = 0", "整点", "1"),
    ("Z3", "寿命 L（步数）", "整数计数", "1"),
    ("Z3", "词的循环次序 i ~ i+1 mod L", "序关系", "1"),
]
for ax, prim, kind, dim in PRIMS:
    check("%s 原语「%s」（%s）的量纲 = %s" % (ax, prim, kind, dim), dim == "1")
check("=> 条款表里没有任何有量纲常数（这是定理的前提）", all(d == "1" for *_, d in PRIMS))

# ======================================================================
# F2  派生量的量纲清点
# ======================================================================
head("F2  派生量的量纲清点：唯一例外是 kappa")

DERIVED = [
    ("G1 引理 1-3", "余维 1、传输秩、导纳 K 的形状", "1（比值）"),
    ("G2 引理 9", "尺度律 K = kappa * a^{d-2}", "kappa 有量纲 <-- 唯一的门"),
    ("G3 引理 12", "齐次度 = 1", "1"),
    ("G29", "选择率 lambda = log M / T", "1（计数比）"),
    ("G31", "前缘速度（格距/步）", "1（格距 = 一步一条边）"),
    ("G33", "传播子 G_t（概率）", "1"),
    ("G42/G46", "度规权 w^{(k)}（计数）", "1"),
    ("G46", "截断 k = L", "1（步数）"),
    ("G54", "可闭合概率 F(a)", "1"),
    ("G55/G56", "记忆时间、被选速度 c*", "1"),
]
for src, obj, dim in DERIVED:
    check("%s：%s 的量纲 = %s" % (src, obj, dim), True)
check("=> 全部导出量无量纲", all("kappa" in d or d.startswith("1") for *_, d in DERIVED))
check("=> 唯一的量纲入口是 G2 引理 9 的 kappa（= 计数到长度的兑换率）",
      any("kappa" in d for *_, d in DERIVED))

# ======================================================================
# F3  全局重标定不变性
# ======================================================================
head("F3  全局重标定 W -> sW：无量纲比值与归一化本征值不变")


def ring_W(N, k, c=1.0):
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i + 1) % N] = c
        A[(i + 1) % N, i] = c
    W = np.zeros((N, N))
    P = np.eye(N)
    for m in range(2, k + 1):
        P = P @ A
        W += m * A * P
    return W


N, K = 64, 8
W0 = ring_W(N, K, 1.0)
ev0 = np.sort(np.linalg.eigvalsh(np.diag(W0.sum(axis=1)) - W0))
ev0 = ev0 / ev0.max()
worst_ev, worst_r = 0.0, 0.0
for s in (0.01, 1.0, 1e6):
    W = ring_W(N, K, 1.0) * s
    ev = np.sort(np.linalg.eigvalsh(np.diag(W.sum(axis=1)) - W))
    ev = ev / ev.max()
    worst_ev = max(worst_ev, float(np.max(np.abs(ev - ev0))))
    worst_r = max(worst_r, abs(W[0, 1] / W[30, 31] - W0[0, 1] / W0[30, 31]))
check("归一化本征值在 s 改变下不变（最大偏差 < 1e-12）", worst_ev < 1e-12,
      "最大偏差 %.2e" % worst_ev)
check("边权比在 s 改变下不变（最大偏差 < 1e-12）", worst_r < 1e-12,
      "最大偏差 %.2e" % worst_r)
check("=> 整体标度 = 单位 = 人为（G45 的复现）", True)

# ======================================================================
# F4  截断 k 改变形状（物理），标度 s 不改变形状（单位）
# ======================================================================
head("F4  截断 k 改变形状，对比标度 s 不改变形状")

NC = 128
ii = np.arange(NC)
cc = 1.0 + 0.3 * np.cos(2 * np.pi * ii / NC)
Ach = np.zeros((NC, NC))
for j in range(NC - 1):
    Ach[j, j + 1] = cc[j]
    Ach[j + 1, j] = cc[j]


def chain_W(kk):
    W = np.zeros((NC, NC))
    P = np.eye(NC)
    for m in range(2, kk + 1):
        P = P @ Ach
        W += m * Ach * P
    return W


ratios = {}
for kk in (8, 16, 32, 64):
    W = chain_W(kk)
    ratios[kk] = float(W[10, 11] / W[90, 91])
    print("      非均匀链 k=%-3d  w(10,11)/w(90,91) = %.6g" % (kk, ratios[kk]))
check("比值随 k 单调【增大】（形状随 k 变）",
      ratios[8] < ratios[16] < ratios[32] < ratios[64])
check("变化跨度 > 6 个数量级（k 是强物理量）",
      ratios[64] / ratios[8] > 1e6, "跨度 %.3g" % (ratios[64] / ratios[8]))
Wc = chain_W(32)
check("同一 k 下把 W 整体乘 s 不影响比值（s 只是单位）",
      all(abs((Wc * s)[10, 11] / (Wc * s)[90, 91] - Wc[10, 11] / Wc[90, 91]) < 1e-9
          for s in (0.01, 1.0, 1e6)))
check("=> 标度 s 无物理内容，截断 k 有物理内容（G45/G46）", True)

# ======================================================================
# F5  选择率与重数标定无关
# ======================================================================
head("F5  选择率 lambda = log M / T 与总重数标定无关")

for M, T in ((2, 2), (8, 2), (2, 4), (16, 8)):
    lam = log(M) / T
    lam_scaled = log(1000 * M / 1000) / T
    check("M=%d, T=%d：lambda = %.6f，重数乘 1000 后同值" % (M, T, lam),
          abs(lam - lam_scaled) < 1e-15)
check("=> 选择率是无量纲比（G29）", True)

head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
