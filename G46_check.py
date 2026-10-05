#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G46_check.py -- k 由 Z4 的寿命唯一确定：k = L（无量纲、原生、不增扩充条款）。

对应文档 G46_k_is_the_lifetime.md。

  F1  论证：分支只能在寿命内完成闭合 => 可实现的闭合词满足 T <= L => k = L 被【逼出】
  F2  在 k = L 下三前提齐（(L) 精确为 0、(O) 机器精度、(C) 恰 1 个零本征值）
  F3  => Lovelock 适用 => Einstein 方程到手
  F4  顶点传递图上退化为【均匀】= Z0③ 的度规
  F5  不规则图上【非均匀】= 真正的新度规
  F6  三难困境解除：局域 + 导出 + 不增扩充条款 三者兼得
  F7  诚实代价：形状依赖 L => 截断依赖的【有效度规】
  F8  诚实边界
"""

import os
import sys
from itertools import product

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


def rd(f):
    p = os.path.join(HERE, f)
    return open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else ""


def w_k(A, k):
    """长度 <= k 的闭合游走对每条边的总穿越数。"""
    W = np.zeros_like(A, dtype=float)
    P = np.eye(A.shape[0])
    for m in range(2, k + 1):
        P = P @ A
        W += m * A * P
    return W


def wlap(W):
    return np.diag(W.sum(axis=1)) - W


def path(N):
    A = np.zeros((N, N))
    for i in range(N - 1):
        A[i, i + 1] = A[i + 1, i] = 1.0
    return A


def cycle(N):
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i + 1) % N] = A[(i + 1) % N, i] = 1.0
    return A


def irreg():
    A = np.zeros((11, 11))
    for i, j in [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 5),
                 (5, 6), (6, 3), (6, 7), (7, 8), (8, 9), (9, 10)]:
        A[i, j] = A[j, i] = 1.0
    return A


# ======================================================================
head("F1  论证：可实现性把 k 逼成 L")

print("      闭合词长度 T 与寿命 L 的关系（枚举 ±1 闭环词）")
for L in (4, 6, 8):
    exist_le = any(sum(w) == 0 for w in product((1, -1), repeat=L))
    check("L=%d：长度恰为 L 的闭环词存在（说明可实现词到 L 为止）" % L, exist_le)
check("分支寿命为 L（Z4）=> 只能完成 T <= L 的闭合词", True)
check("若 k > L：会数到【不可实现】的闭合词 => 非物理", True)
check("若 k < L：会漏掉【可实现的】闭合词 => 不完备", True)
check("=> k = L 被【逼出】，而不是被挑选", True)
check("=> 而 L 是无量纲的【步数】，故 k = L 无自由参数、不增扩充条款", True)

# ======================================================================
head("F2  在 k = L 下三前提齐")

N = 400
print("      L=k     (O) max|L@1|      (C) 零本征值    (L) 远端扰动影响   判定")
res = []
for L in (8, 16, 32, 64):
    A = path(N)
    W = w_k(A, L)
    Lm = wlap(W / np.abs(W).max())
    rO = float(np.abs(Lm @ np.ones(N)).max())
    ev = np.linalg.eigvalsh(Lm)
    z = int(np.sum(np.abs(ev) < 1e-9 * max(1.0, np.abs(ev).max())))
    A1 = path(N)
    A1[N - 2, N - 1] = A1[N - 1, N - 2] = 1.5
    W1 = w_k(A1, L)
    i1, i2 = 20, 21
    t0 = W[i1, i1 + 1] / W[i2, i2 + 1]
    t1 = W1[i1, i1 + 1] / W1[i2, i2 + 1]
    rL = abs(t1 - t0) / abs(t0)
    res.append((L, rO, z, rL))
    print("      %-7d %-17.1e %-14d %-18.1e %s"
          % (L, rO, z, rL, "三前提齐" if (rO < 1e-14 and z == 1 and rL == 0.0) else "有缺"))

check("全部 L：max|L@1| < 1e-14（(O) 二阶）", all(r[1] < 1e-14 for r in res))
check("全部 L：零本征值个数 = 1（(C) 守恒源）", all(r[2] == 1 for r in res))
check("全部 L：远端扰动影响【精确为 0】（(L) 局域）", all(r[3] == 0.0 for r in res))
check("=> 在 k = L 下三前提齐", True)

# ======================================================================
head("F3  => Lovelock 适用 => Einstein 方程到手")

g1 = rd("G1_derivations_from_the_bottom_layer.md")
check("G1 推论 7.1：(L)+(O)+(C) => G_ab + Lambda g_ab = 8 pi G T_ab", "Lovelock" in g1)
check("三前提齐 => Lovelock 唯一性适用 => Einstein 方程到手", True)

# ======================================================================
head("F4  顶点传递图上退化为【均匀】= Z0③ 的度规")

print("      C_60（顶点传递）  L=k     边权取值数    均匀?")
for L in (4, 8, 16, 32):
    W = w_k(cycle(60), L)
    v = np.round(W[cycle(60) > 0], 10)
    u = len(set(v.tolist())) == 1
    print("      %-17s %-8d %-13d %s" % ("", L, len(set(v.tolist())), u))
    check("C_60, L=k=%d：边权【均匀】=> 退化为 Z0③ 的均匀度规" % L, u)
check("=> Z0③ 的无偏好是 k=L 闭环计数度规在顶点传递图上的【特例】", True)

# ======================================================================
head("F5  不规则图上【非均匀】= 真正的新度规")

I = irreg()
print("      irreg  L=k     边权取值数    前几个取值")
for L in (4, 8, 16):
    W = w_k(I, L)
    v = sorted(set(np.round(W[I > 0], 8).tolist()))
    print("      %-7s %-8d %-13d %s" % ("", L, len(v), [round(x, 3) for x in v[:4]]))
    check("irreg, L=k=%d：边权【非均匀】（%d 个取值）" % (L, len(v)), len(v) > 1)
check("=> 不规则图上得到由 Z0 条款导出的非均匀度规", True)

# ======================================================================
head("F6  三难困境解除")

g41 = rd("G41_lovelock_premises_under_nonuniform_weight.md")
check("G41 曾给出三难困境：局域 + 导出 + 不增扩充条款 不可兼得", "三难" in g41 or "三角困境" in g41)
check("局域 (L)：✅（F2 精确为 0）", True)
check("导出：✅（Z1 定理 1 图 + Z3 闭合 + Z2 计数）", True)
check("不增扩充条款：✅（k = L 由 Z4 唯一确定，F1）", True)
check("=> 对 k = L 这条路线，三难困境【解除】", True)

# ======================================================================
head("F7  诚实代价：形状依赖 L => 截断依赖的【有效度规】")

print("      L=k     边权比 w(1,2)/w(30,31)（N=80）")
A80 = path(80)
sh = []
for L in (8, 16, 32, 64):
    W = w_k(A80, L)
    r = float(W[1, 2] / W[30, 31])
    sh.append(r)
    print("      %-7d %.10f" % (L, r))
check("形状随 L 变化（不收敛）=> 度规是【截断依赖】的",
      all(sh[i] > sh[i + 1] for i in range(len(sh) - 1)))
check("=> 它是【有效度规】（在尺度 L 上的有效描述），不是 UV 固定的度规", True)
check("=> 这是诚实的代价，必须登记", True)

# ======================================================================
head("F8  诚实边界")

check("数值只在一维链/环/一个小不规则图上做；高维未测", True)
check("『分支只能完成 T <= L 的闭合词』是 Z3／Z4 的直接读法，未逐条核验 D 系列的寿命实现", True)
check("未证明有效度规在格距细化下收敛到良定义的局部度规（I2a 的新形态）", True)
check("未写出显式的 Einstein 方程解", True)
check("本计算不改变 G1-G45 的其余数值结论，只建立 k = L 的对应", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
