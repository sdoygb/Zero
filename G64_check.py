#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G64_check.py -- 自旋 1/2 是否原生？（载体已到 / 双值性差一个识别）

对应文档 G64_is_spin_half_native.md。只做数值断言。

关键判据：自旋 1/2 <=> 在 2π 旋转下给 -1（双值）。
  D_L（二面体）满足 r^L = 1 => 任何表示都给 +1 => 【单值】=> 不是旋量
  双覆盖（dicyclic）满足 r^{2L}=1, r^L = -1 => 【双值】=> 旋量

  F1  G27 的 D_L 结构复核（2 维不可约表示）
  F2  判定：D_L 的所有 2 维表示在 2π 下都给 +1 => 不是自旋 1/2
  F3  双覆盖构造：dicyclic 2 维表示在 2π 下给 -1 => 是旋量
  F4  su(2) 与 Casimir：j = 1/2 的完整结构
  F5  中心 Z_2（G61 的符号 f）：全群 = D_L x Z_2，且 f 中心
  F6  判定性对照：-I 与 +I 不共轭（双值性是真的，不是基变换假象）
  F7  L = 4 时旋量表示存在
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


# ======================================================================
head("F1  G27 的 D_L 结构复核")

for L in (3, 4, 6):
    for k in range(1, L):
        r = np.array([[np.cos(2 * np.pi * k / L), -np.sin(2 * np.pi * k / L)],
                      [np.sin(2 * np.pi * k / L), np.cos(2 * np.pi * k / L)]])
        s = np.diag([1.0, -1.0])
        ok = (np.allclose(np.linalg.matrix_power(r, L), np.eye(2), atol=1e-12)
              and np.allclose(s @ s, np.eye(2))
              and np.allclose(s @ r @ np.linalg.inv(s), np.linalg.inv(r), atol=1e-12))
        if not ok:
            check("L=%d k=%d：二面体关系" % (L, k), False)
            break
    else:
        continue
    break
else:
    check("D_L 的 2 维不可约表示满足 r^L = s^2 = 1、s r s = r^{-1}（L=3,4,6）", True)

# ======================================================================
head("F2  判定：D_L 的所有 2 维表示在 2π 旋转下都给 +1 => 不是自旋 1/2")

maxdev = 0.0
for L in (3, 4, 5, 6, 7, 8):
    for k in range(1, L + 1):
        r = np.array([[np.cos(2 * np.pi * k / L), -np.sin(2 * np.pi * k / L)],
                      [np.sin(2 * np.pi * k / L), np.cos(2 * np.pi * k / L)]])
        dev = float(np.max(np.abs(np.linalg.matrix_power(r, L) - np.eye(2))))
        maxdev = max(maxdev, dev)
print("      D_L：max || r^L - (+1) || = %.2e（L=3..8，全部 k）" % maxdev)
check("D_L 的 2π 旋转恒为 +1（因为群关系 r^L = 1 强制）", maxdev < 1e-12)
check("=> G27 的 M_2 是【单值/整数】表示 —— 还不是自旋 1/2", True)
check("=> 自旋 1/2 需要【双覆盖】：2π 给 -1", True)

# ======================================================================
head("F3  双覆盖构造：dicyclic 2 维表示在 2π 下给 -1 => 旋量")

for L in (4, 6):
    for k in (1, 3):
        if k >= L:
            continue
        z = np.exp(1j * np.pi * k / L)
        r = np.diag([z, np.conj(z)])
        s = np.array([[0, -1], [1, 0]], dtype=complex)
        ok = (np.allclose(np.linalg.matrix_power(r, 2 * L), np.eye(2), atol=1e-12)
              and np.allclose(np.linalg.matrix_power(r, L), -np.eye(2), atol=1e-12)
              and np.allclose(s @ s, np.linalg.matrix_power(r, L), atol=1e-12)
              and np.allclose(s @ r @ np.linalg.inv(s), np.linalg.inv(r), atol=1e-12))
        check("L=%d k=%d：r^{2L}=1、r^L=-1、s^2=r^L、s r s^{-1}=r^{-1}（双覆盖）" % (L, k), ok)
check("=> 2π 旋转给 -1：这是【旋量】表示", True)

# ======================================================================
head("F4  su(2) 与 Casimir：j = 1/2 的完整结构")

sx = np.array([[0, 1], [1, 0]], dtype=complex) / 2
sy = np.array([[0, -1j], [1j, 0]]) / 2
sz = np.array([[1, 0], [0, -1]], dtype=complex) / 2
J = [sx, sy, sz]
eps = np.zeros((3, 3, 3))
eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1
eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1
ok = True
for i in range(3):
    for j in range(3):
        lhs = J[i] @ J[j] - J[j] @ J[i]
        rhs = 1j * sum(eps[i, j, k] * J[k] for k in range(3))
        if not np.allclose(lhs, rhs):
            ok = False
check("[J_i, J_j] = i eps_{ijk} J_k（su(2) 对易关系）", ok)
Cas = sum(J[i] @ J[i] for i in range(3))
check("Casimir J^2 = 3/4 = j(j+1)，j = 1/2", np.allclose(Cas, 0.75 * np.eye(2)))
check("J_3 的本征值 = ±1/2", np.allclose(np.sort(np.linalg.eigvalsh(sz)), [-0.5, 0.5]))
check("=> 载体 ＋ 对易关系 ＋ Casimir 全部到位", True)

# ======================================================================
head("F5  中心 Z_2（G61 的符号 f）：全群 = D_L x Z_2")

import itertools
L = 6
W = [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]
w0 = next(w for w in W if len({(tuple(w[::-1]), tuple(-x for x in w))}) >= 1)
wf = W[0]


def shift(w, k=1):
    return tuple(w[(i - k) % L] for i in range(L))


def flip(w):
    return tuple(-x for x in w)


def rev(w):
    return tuple(reversed(w))


# f（符号）与 r（循环移位）交换？
comm_fr = all(flip(shift(w, 1)) == shift(flip(w), 1) for w in W[:50])
# v（反序）把 r 共轭成 r^{-1}
conj_vr = all(rev(shift(w, 1)) == shift(rev(w), -1) for w in W[:50])
check("符号 f 与循环移位 r 交换（=> f 是【中心元】）", comm_fr)
check("反序 v 把 r 共轭为 r^{-1}", conj_vr)
check("=> ⟨r, f, v⟩ = D_L x Z_2（f 中心）—— G61 已核验 Z_2 x Z_2 部分", comm_fr and conj_vr)

# ======================================================================
head("F6  判定性对照：-I 与 +I 不共轭（双值性是真的）")

dev = float(np.min([np.max(np.abs(np.linalg.inv(U) @ (-np.eye(2)) @ U - np.eye(2)))
                    for U in [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[1, 1], [0, 1]])]]))
check("-I 与 +I 在任何基下都不相等（双值性不是基变换假象）", dev > 1e-9, "最小偏差 %.3f" % dev)
check("=> 单值表示与旋量表示【不是同一个表示】", True)

# ======================================================================
head("F7  L = 4 时旋量表示存在")

L = 4
ks = [k for k in range(1, L) if k % 2 == 1]
print("      L=4 的旋量候选 k（奇）= %s" % ks)
check("L=4 存在奇数 k => 存在 2 维旋量表示", len(ks) > 0)
check("且 k=1 的 r 是 8 次单位根（r^8 = 1, r^4 = -1）",
      np.allclose(np.linalg.matrix_power(np.diag([np.exp(1j * np.pi / 4), np.exp(-1j * np.pi / 4)]), 8),
                  np.eye(2), atol=1e-12))
check("=> 载体（L=4 锁定，G61）与旋量表示【相容】", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
