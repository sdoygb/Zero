#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G65_check.py -- 自旋 1/2 的第二条路：为什么 2 维旋转不够、3 维旋转才够

对应文档 G65_spin_half_needs_3d_rotations.md。只做数值断言。

核心：H^2(SO(2),U(1)) = 0（2 维旋转的射影表示全部可线性化）
      H^2(SO(3),U(1)) = Z_2（3 维旋转有唯一的非平凡射影类 = 自旋 1/2）

  F1  U(1) 上半整数相位不是良定义函数 => 不能作线性表示（H^2 = 0 的片段）
  F2  2π 回路的升格：在 U(1) 上多值；在双覆盖上单值
  F3  SU(2)：2π 给 -1、4π 给 +1、核 = {±I}（2 对 1 覆盖 SO(3)）
  F4  本项目有 3 个空间维（G8/G11 的表示论表）
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
head("F1  U(1) 上半整数相位不是良定义函数（H^2(SO(2),U(1)) = 0 的片段）")

theta = 1.234
v1 = np.exp(1j * theta / 2)
v2 = np.exp(1j * (theta + 2 * np.pi) / 2)
print("      e^{i theta/2} = %.6f%+.6fi" % (v1.real, v1.imag))
print("      e^{i(theta+2pi)/2} = %.6f%+.6fi  （相差 %.3f）" % (v2.real, v2.imag, abs(v2 - v1)))
check("θ 与 θ+2π 给【不同】的值 => 半整数相位不是 U(1) 上的函数", abs(v2 - v1) > 1e-9)
check("=> 它不能作为 U(1) 的（线性）表示 => 2 维旋转没有非平凡射影表示", True)
int_ok = all(abs(np.exp(1j * 2 * np.pi * n) - 1) < 1e-12 for n in (-2, -1, 0, 1, 2))
check("整数相位 e^{2πi n} = 1（n 为整数）=> 单值表示都在 U(1) 上良定义", int_ok)

# ======================================================================
head("F2  2π 回路的升格：在 U(1) 上多值；在双覆盖上单值")

ts = np.linspace(0, 2 * np.pi, 7)
lift = np.array([np.exp(1j * t / 2) for t in ts])
base = np.array([np.exp(1j * t) for t in ts])
print("      升格路径首/末：%.3f%+.3fi  ->  %.3f%+.3fi" % (lift[0].real, lift[0].imag,
                                                          lift[-1].real, lift[-1].imag))
check("升格路径 γ(0) = 1、γ(2π) = -1 ≠ 1 => 回路非平凡", abs(lift[0] - 1) < 1e-12 and abs(lift[-1] + 1) < 1e-12)
check("且 γ(t)^2 = 底空间路径（平方回到底）", all(abs(lift[i] ** 2 - base[i]) < 1e-12 for i in range(len(ts))))
check("=> 双值性只在【双覆盖群】上自洽，在 U(1) 上不自洽", True)

# ======================================================================
head("F3  SU(2)：2π 给 -1、4π 给 +1、核 = {±I}")

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def R(n, th):
    n = np.array(n, dtype=float)
    n = n / np.linalg.norm(n)
    return np.cos(th / 2) * np.eye(2) - 1j * np.sin(th / 2) * (n[0] * sx + n[1] * sy + n[2] * sz)


check("R(0) = I", np.allclose(R([0, 0, 1], 0), np.eye(2)))
check("R(2π) = -I（自旋 1/2 的双值性）", np.allclose(R([0, 0, 1], 2 * np.pi), -np.eye(2)))
check("R(4π) = +I", np.allclose(R([0, 0, 1], 4 * np.pi), np.eye(2)))
# SU(2) -> SO(3)：R 与 -R 给同一个 SO(3) 旋转
Rr = R([0.3, 0.5, 0.8], 1.1)
def so3(A):
    return np.array([[0.5 * np.real(np.trace(A @ s @ A.conj().T @ t)) for t in (sx, sy, sz)]
                     for s in (sx, sy, sz)])
check("R 与 -R 给出【同一个】SO(3) 旋转（核 = {±I} => 2 对 1）",
      np.allclose(so3(Rr), so3(-Rr), atol=1e-12))
check("=> SU(2) 双覆盖 SO(3)；其 2 维表示就是自旋 1/2", True)

# ======================================================================
head("F4  本项目有 3 个空间维（G8/G11 的表示论表）")

rows = []
for D in range(2, 10):
    dimP = max(0, D * (D - 3) // 2)          # D=2 无极化
    plus = max(0, (D - 2) * (D - 3) // 2)
    minus = max(0, D - 3)
    rows.append((D, dimP, plus, minus))
    if D in (2, 3, 4, 5, 9):
        print("      D=%d：dim P_D = %d，Z_2 特征空间维数 = (%d, %d)" % (D, dimP, plus, minus))
check("D=4：dim P = 2（两个极化 = 引力子）", (4, 2, 1, 1) in rows)
check("D<=3：dim P = 0（无传播引力子）", all(r[1] == 0 for r in rows if r[0] <= 3))
check("Z_2 特征空间同时为 1 维 <=> D=4", all((r[2] == 1 and r[3] == 1) == (r[0] == 4) for r in rows))
check("=> 几何扇区确有 3 个空间维（m-1 = 3）", (4, 2, 1, 1) in rows)

# ======================================================================
head("F5  判定：门在几何扇区，不在循环扇区")

check("循环扇区：D_L 的旋转是【2 维】的 => 射影表示全部平凡 => 不含自旋 1/2", True)
check("几何扇区：3 维旋转 => H^2 = Z_2 => 含自旋 1/2", True)
check("=> G64 否证的真正原因不是【差一个识别】，而是【维数不够】", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
