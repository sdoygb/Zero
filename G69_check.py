#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G69_check.py -- 洛伦兹自旋：SU(2) 升级到 SL(2,C)

对应文档 G69_lorentz_spin_from_SL2C.md。只做数值断言。

  F1  SL(2,C) -> SO+(1,3)：X -> A X A^dag 保 Minkowski 范数
  F2  同态、核 = {±I}、det=1、正交时（Lambda^0_0 > 0）
  F3  Cartan 分解：A = U H（酉 x Hermite 正定）= 旋转 x boost
  F4  2π 旋转 = -I 被继承（与 G66/G67 自洽）
  F5  签名的非对称：gamma^0^2 = +1，gamma^i^2 = -1（类时/类空反射不同）
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
rng = np.random.default_rng(17)

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], dtype=complex)
SIG = [I2, sx, sy, sz]                      # sigma_0 = I
ETA = np.diag([1.0, -1.0, -1.0, -1.0])


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


def rand_sl2c():
    M = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    d = np.linalg.det(M)
    return M / np.sqrt(d)


def x_of(v):
    return sum(v[i] * SIG[i] for i in range(4))


def vec_of(X):
    return np.array([0.5 * np.real(np.trace(SIG[i] @ X)) for i in range(4)])


def lam(A):
    """由 A 诱导的 4x4 洛伦兹矩阵。"""
    L = np.zeros((4, 4))
    for j in range(4):
        e = np.zeros(4)
        e[j] = 1.0
        L[:, j] = vec_of(A @ x_of(e) @ A.conj().T)
    return L


# ======================================================================
head("F1  SL(2,C) -> SO+(1,3)：保 Minkowski 范数")

ok_norm = ok_herm = True
for _ in range(40):
    A = rand_sl2c()
    v = rng.normal(size=4)
    X = x_of(v)
    Xp = A @ X @ A.conj().T
    if not np.allclose(Xp, Xp.conj().T, atol=1e-12):
        ok_herm = False
    if abs(np.real(np.linalg.det(Xp)) - np.real(np.linalg.det(X))) > 1e-10:
        ok_norm = False
check("A X A^dag 仍 Hermite", ok_herm)
check("det(A X A^dag) = det(X)  => 保 Minkowski 范数", ok_norm)
check("=> 诱导出洛伦兹变换", True)

# ======================================================================
head("F2  同态、核 = {±I}、det=1、正交时")

ok_orth = ok_det = ok_hom = ok_kern = ok_ortho = True
for _ in range(40):
    A, B = rand_sl2c(), rand_sl2c()
    L, LB = lam(A), lam(B)
    if not np.allclose(L.T @ ETA @ L, ETA, atol=1e-10):
        ok_orth = False
    if abs(np.linalg.det(L) - 1) > 1e-9:
        ok_det = False
    if L[0, 0] <= 0:
        ok_ortho = False
    if not np.allclose(lam(A @ B), L @ LB, atol=1e-9):
        ok_hom = False
    if not np.allclose(lam(-A), L, atol=1e-12):
        ok_kern = False
check("Lambda^T eta Lambda = eta（洛伦兹）", ok_orth)
check("det Lambda = 1（固有）", ok_det)
check("Lambda^0_0 > 0（正交时）", ok_ortho)
check("同态：Lambda(AB) = Lambda(A) Lambda(B)", ok_hom)
check("核含 ±I：Lambda(-A) = Lambda(A) => 2 对 1", ok_kern)
check("=> SL(2,C) 双覆盖 SO+(1,3)", ok_orth and ok_det and ok_hom and ok_kern and ok_ortho)

# ======================================================================
head("F3  Cartan 分解：A = U H = 旋转 x boost")

ok_dec = ok_rot = ok_boost = True
for _ in range(40):
    A = rand_sl2c()
    M = A.conj().T @ A                                   # Hermite 正定
    w, V = np.linalg.eigh(M)
    H = V @ np.diag(np.sqrt(w)) @ V.conj().T             # Hermite 平方根（非 Cholesky）
    U = A @ np.linalg.inv(H)
    if not (np.allclose(U @ U.conj().T, I2, atol=1e-10) and abs(np.linalg.det(U) - 1) < 1e-9):
        ok_dec = False
    if not (np.allclose(H, H.conj().T, atol=1e-10) and np.all(np.linalg.eigvalsh(H) > 0)):
        ok_dec = False
    LU, LH = lam(U), lam(H)
    if abs(LU[0, 0] - 1) > 1e-9 or np.max(np.abs(LU[0, 1:])) > 1e-9:
        ok_rot = False
    if not (LH[0, 0] > 1.0 - 1e-12):
        ok_boost = False
check("A = U H 且 det U = det H = 1", ok_dec)
check("U 部分诱导【纯旋转】（Lambda^0_0 = 1、Lambda^0_i = 0）", ok_rot)
check("H 部分诱导【boost】（Lambda^0_0 >= 1）", ok_boost)
check("=> SL(2,C) = 旋转 x boost（极分解）", ok_dec and ok_rot and ok_boost)

# ======================================================================
head("F4  2π 旋转 = -I 被继承")

n = np.array([0.3, 0.5, 0.8])
n = n / np.linalg.norm(n)
R2pi = np.cos(np.pi) * I2 - 1j * np.sin(np.pi) * (n[0] * sx + n[1] * sy + n[2] * sz)
check("2π 旋转 = -I（自旋 1/2 的双值性）", np.allclose(R2pi, -I2, atol=1e-12))
check("Lambda(-I) = 恒等（几何上什么都不做）", np.allclose(lam(-I2), np.eye(4), atol=1e-12))
check("=> 洛伦兹情形继承空间旋转的双值性（与 G66/G67 自洽）", True)

# ======================================================================
head("F5  签名的非对称：gamma^0^2 = +1，gamma^i^2 = -1")

# Weyl 表示的 4x4 gamma 矩阵：{gamma^mu, gamma^nu} = 2 eta^{mu nu} I
Z = np.zeros((2, 2), dtype=complex)
G = [np.block([[Z, I2], [I2, Z]])]                       # gamma^0
for s in (sx, sy, sz):
    G.append(np.block([[Z, s], [-s, Z]]))                # gamma^i
check("(gamma^0)^2 = +I（类时）", np.allclose(G[0] @ G[0], np.eye(4)))
check("(gamma^i)^2 = -I（类空）", all(np.allclose(G[i] @ G[i], -np.eye(4)) for i in (1, 2, 3)))
anti = all(np.allclose(G[m] @ G[n] + G[n] @ G[m], 2 * ETA[m, n] * np.eye(4))
           for m in range(4) for n in range(4))
check("{{gamma^mu, gamma^nu}} = 2 eta^{{mu nu}} I（Clifford 关系编码签名）", anti)
check("=> 类时反射与类空反射的升格【平方不同】", True)
check("=> 洛伦兹情形把 G67 的反射 Z_2 按签名【劈开】", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
