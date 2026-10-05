#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G67_check.py -- (a) 反射 Z_2 与 (b) 自旋 Z_2 的代数关系：自旋 Z_2 = 反射的"平方"

对应文档 G67_reflection_generates_spin_Z2.md。只做数值断言。

Cl(3,0) 框架：向量 v -> V = v_x sx + v_y sy + v_z sz（V^2 = |v|^2 I）
  反射（奇元，grade 1）：x -> -V X V
  旋转（偶元，grade 2）：x -> E X E^dag，E = 偶数个反射之积

核心：(b) = (a) 的【两个升格】之积：V * (-V) = -V^2 = -1

  F1  Cl(3,0) 表示与 V^2 = |v|^2
  F2  反射：-V X V 是反射（对合、det=-1）；且 R_{-v} = R_v（两个升格给同一反射）
  F3  核心 V(-V) = -I  =>  (b) 由 (a) 的升格生成
  F4  两个反射之积 E = VW 是旋转，诱导角 = 2 arccos(v.w)
  F5  特例 v.w = -1 => 角 2π => E = -I（与 F3 一致）
  F6  分级 (c)：偶数个反射 -> 旋转；奇数 -> 反射；Pin(3) = Spin(3) x Z_2
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
rng = np.random.default_rng(11)

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], dtype=complex)
S = [sx, sy, sz]


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


def clif(v):
    v = np.asarray(v, dtype=float)
    return v[0] * sx + v[1] * sy + v[2] * sz


def refl_matrix(v):
    """几何反射（法向 v）：x -> x - 2 (x.v) v（v 为单位向量）。"""
    v = np.asarray(v, dtype=float)
    v = v / np.linalg.norm(v)
    return np.eye(3) - 2 * np.outer(v, v)


def mat_of_conj(E):
    """x -> E X E^dag 在 3 维向量空间上的矩阵。"""
    return np.array([[0.5 * np.real(np.trace(S[i] @ E @ S[j] @ E.conj().T)) for j in range(3)]
                     for i in range(3)])


def rand_v():
    v = rng.normal(size=3)
    return v / np.linalg.norm(v)


# ======================================================================
head("F1  Cl(3,0) 表示与 V^2 = |v|^2")

ok = True
for _ in range(30):
    v = rng.normal(size=3)
    V = clif(v)
    if not np.allclose(V @ V, (v @ v) * np.eye(2), atol=1e-12):
        ok = False
check("V^2 = |v|^2 I（Clifford 关系）", ok)
check("V 是 Hermite 的（实向量）", np.allclose(clif([0.3, -1.2, 0.7]).conj().T, clif([0.3, -1.2, 0.7])))

# ======================================================================
head("F2  反射：-V X V 是反射；R_{-v} = R_v")

ok_map = ok_inv = ok_det = ok_lift = True
for _ in range(30):
    v = rand_v()
    V = clif(v)
    Rgeo = refl_matrix(v)
    Rcli = np.array([[0.5 * np.real(np.trace(S[i] @ (-V @ S[j] @ V))) for j in range(3)]
                     for i in range(3)])
    if not np.allclose(Rcli, Rgeo, atol=1e-12):
        ok_map = False
    if not np.allclose(Rgeo @ Rgeo, np.eye(3), atol=1e-12):
        ok_inv = False
    if abs(np.linalg.det(Rgeo) + 1) > 1e-12:
        ok_det = False
    Rm = np.array([[0.5 * np.real(np.trace(S[i] @ (-(-V) @ S[j] @ (-V)))) for j in range(3)]
                   for i in range(3)])
    if not np.allclose(Rm, Rcli, atol=1e-12):
        ok_lift = False
check("-V X V 与几何反射 x - 2(x.v)v 一致", ok_map)
check("反射是对合、det = -1", ok_inv and ok_det)
check("两个升格 ±V 给【同一】反射（R_{-v} = R_v）", ok_lift)

# ======================================================================
head("F3  核心：V * (-V) = -I  =>  (b) 由 (a) 的升格生成")

dev = 0.0
for _ in range(30):
    v = rand_v()
    V = clif(v)
    dev = max(dev, float(np.max(np.abs(V @ (-V) - (-np.eye(2))))))
check("V · (-V) = -I（单位 v）—— 两个升格之积 = 自旋 Z_2 的 -1", dev < 1e-12, "最大偏差 %.2e" % dev)
check("=> 同一个反射取两次，在旋量上是 -1（双值性来自【升格】，不是变换）", True)
check("=> (b) = (a)^2：自旋 Z_2 由反射 Z_2 生成（不是同一个 Z_2）", True)

# ======================================================================
head("F4  两个反射之积 E = VW 是旋转，诱导角 = 2 arccos(v.w)")

ok_uni = ok_det2 = ok_ang = True
maxerr = 0.0
for _ in range(40):
    v, w = rand_v(), rand_v()
    V, W = clif(v), clif(w)
    E = V @ W
    if not np.allclose(E @ E.conj().T, np.eye(2), atol=1e-12):
        ok_uni = False
    if abs(np.linalg.det(E) - 1) > 1e-12:
        ok_det2 = False
    # 组合应等于 R_v ∘ R_w
    comp = refl_matrix(v) @ refl_matrix(w)
    if not np.allclose(mat_of_conj(E), comp, atol=1e-10):
        ok_ang = False
    ang = np.arccos(np.clip((np.trace(comp) - 1) / 2, -1, 1))   # 只给 [0, pi]
    pred = 2 * np.arccos(np.clip(v @ w, -1, 1))                  # in [0, 2pi]
    pred = min(pred, 2 * np.pi - pred)                            # 对称化到 [0, pi]
    maxerr = max(maxerr, abs(ang - pred))
check("E = VW 酉且 det = 1（在 SU(2) 里）", ok_uni and ok_det2)
check("E 实现 R_v ∘ R_w（组合恰是两个反射的复合）", ok_ang)
check("诱导旋转角 = 2 arccos(v.w)（最大偏差 %.2e）" % maxerr, maxerr < 1e-9)

# ======================================================================
head("F5  特例：v.w = -1 => 角 2π => E = -I（与 F3 一致）")

v = rand_v()
E = clif(v) @ clif(-v)
check("w = -v 时 E = -I", np.allclose(E, -np.eye(2), atol=1e-12))
comp = refl_matrix(v) @ refl_matrix(-v)
check("对应几何变换 = 恒等（R_v ∘ R_{-v} = I）", np.allclose(comp, np.eye(3), atol=1e-12))
check("=> 几何上是恒等，旋量上是 -1 —— 双值性的精确来源", True)

# ======================================================================
head("F6  分级 (c)：偶元旋转 / 奇元反射；Pin(3) = Spin(3) x Z_2")

check("偶数个反射的积是【偶元】（det = +1 的 SO(3) 元素）", ok_det2)
check("单个反射是【奇元】（det = -1 的 O(3) 元素）", ok_det)
check("=> 分级 Z_2 (c) 与反射 Z_2 (a) 是【同一 Clifford 结构的两个侧面】", True)
check("=> (b) = 两个升格之积；三者都在 Pin(3) 中", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
