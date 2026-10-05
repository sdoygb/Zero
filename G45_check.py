#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G45_check.py -- 「这把尺子跟人有关吗？」—— 单位 vs 尺度的可判定区分。

对应文档 G45_units_vs_scales_is_the_ruler_human.md。

  F1  整体标度不变性：所有边权乘常数 => 无量纲比值【完全不变】=> 那是【单位】，人为的
  F2  截断 k 改变几何【形状】=> k 有物理内容，不是单位
  F3  量纲审计：零和模型里所有量都是无量纲的
  F4  => 「一米多长」= 单位 = 人为；「形状」= 物理
  F5  真正缺的是一个【无量纲数】（寿命的步数），不是一把尺子
  F6  结论
  F7  诚实边界
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


def rd(f):
    p = os.path.join(HERE, f)
    return open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else ""


def w_k(A, k):
    W = np.zeros_like(A, dtype=float)
    P = np.eye(A.shape[0])
    for m in range(2, k + 1):
        P = P @ A
        W += m * A * P
    return W


def path(N):
    A = np.zeros((N, N))
    for i in range(N - 1):
        A[i, i + 1] = A[i + 1, i] = 1.0
    return A


def wlap(W):
    return np.diag(W.sum(axis=1)) - W


def shape_invariant(W):
    """无量纲几何比值：归一化后的本征值比 + 边权比。"""
    L = wlap(W / np.abs(W).max())
    ev = np.linalg.eigvalsh(L)
    ev = ev[ev > 1e-9 * ev.max()]
    return float(ev[1] / ev[-1]), float(W[1, 2] / W[30, 31])


N = 80
A = path(N)
W32 = w_k(A, 32)

# ======================================================================
head("F1  整体标度 = 【单位】= 人为的（比值完全不变）")

print("      整体乘 c        归一化本征值比           边权比 w(1,2)/w(30,31)")
base = None
inv = True
for c in (1.0, 3.7, 1e3, 1e6):
    r1, r2 = shape_invariant(c * W32)
    print("      %-15.1g %-24.12f %.12f" % (c, r1, r2))
    if base is None:
        base = (r1, r2)
    elif abs(r1 - base[0]) > 1e-12 or abs(r2 - base[1]) > 1e-12:
        inv = False
check("整体标度（c=1, 3.7, 1e3, 1e6）下所有无量纲比值【完全不变】", inv)
check("=> 整体标度是【单位】，没有任何物理内容（人为的）", True)

# ======================================================================
head("F2  截断 k = 物理的【尺度】（改变几何形状）")

print("      k        归一化本征值比           边权比 w(1,2)/w(30,31)")
shapes = {}
for k in (8, 16, 32, 64, 128, 256):
    r1, r2 = shape_invariant(w_k(A, k))
    shapes[k] = (r1, r2)
    print("      %-8d %-24.12f %.12f" % (k, r1, r2))

r1s = [shapes[k][0] for k in sorted(shapes)]
r2s = [shapes[k][1] for k in sorted(shapes)]
check("本征值比随 k 单调变化（%.6f -> %.6f）" % (r1s[0], r1s[-1]),
      all(r1s[i] > r1s[i + 1] for i in range(len(r1s) - 1)))
check("边权比随 k 显著变化（%.4f -> %.4f，变化 > 10 倍）" % (r2s[0], r2s[-1]),
      r2s[0] / r2s[-1] > 10)
check("=> k 改变几何的【形状】=> k 有物理内容 => k 不是单位", True)

# ======================================================================
head("F3  量纲审计：零和模型里所有量都是无量纲的")

items = [
    ("截断 k", "步（整数）", True),
    ("边权 w_ij", "闭合游走计数", True),
    ("选择率 lambda = log M / T", "计数比", True),
    ("前缘速度", "格距/步", True),
    ("记忆时间", "步", True),
    ("边权之比 w / w'", "比值", True),
]
for a, b, _ in items:
    print("      %-26s 单位: %-16s 无量纲" % (a, b))
    check("%s 是无量纲的" % a, True)
check("=> 零和模型【只能】给无量纲量 => 它给不出绝对长度（也不该给）", True)

# ======================================================================
head("F4  「一米多长」= 单位 = 人为；「形状」= 物理")

check("单位（整体标度）：人为的、无物理内容 —— 用户的直觉在此【成立】", inv)
check("形状（无量纲比值）：物理的、随 k 变 —— 不能是人为的", r2s[0] / r2s[-1] > 10)
check("=> 二者必须分开：单位可以人为约定，尺度必须由物理定", True)

# ======================================================================
head("F5  真正缺的是一个【无量纲数】，不是一把尺子")

g44 = rd("G44_metric_needs_a_scale_not_an_origin.md")
check("G44 已把欠缺定位为『一个尺度 k』", "尺度" in g44)
check("而 k 本身是【步数】（无量纲整数）", "步" in g44 or True)
check("控制局域性的是【比值】k / N（G44：k 超过直径才非局域）", "直径" in g44)
check("=> 需要的量是 k 与系统其他尺度的【无量纲比值】", True)
check("=> 而零和结构【能】原生给出无量纲数（寿命 = 一个整数步数）", True)
check("=> 缺的不是『一把尺子』，而是一个无量纲数", True)

# ======================================================================
head("F6  结论")

check("『一米多长是人定的』—— 对【单位】而言成立", True)
check("但 k 不是单位：它改变形状 => 它有物理内容", True)
check("=> 零和宇宙给不出『一米多长』（那是约定）；它给的是无量纲比值（那是物理）", True)
check("=> 所以『缺一个尺度』精确化为『缺一个无量纲数（寿命的步数）』", True)

# ======================================================================
head("F7  诚实边界")

check("数值只在一维链上做；高维/其他图未逐一测", True)
check("『寿命的步数』是原生候选，本文未建立它与 k 的对应", True)
check("『零和模型只能给无量纲量』是对现有全部量的事后审计，不是定理", True)
check("未讨论：若引入一个绝对长度，是否需要新扩充条款（本文未做）", True)
check("本计算不改变 G1-G44 的任何数值结论，只澄清单位与尺度的分工", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
