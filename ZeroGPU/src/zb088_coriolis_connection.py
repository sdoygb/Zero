#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Zb088 · 旋转读出：科氏项是连接项，不是源。

五组判据（对应正文 §2–§3、§6）：
  S1 一阶恒等式   (∂_t - ω∂_φ) v = ∂_t u，其中 v(y,t)=u(R(ωt)y,t)，∂_φ=(Jy)·∇。
  S2 判据方向     去掉 ω∂_φ 项（"把旋转当没发生"）：残差立刻是 ω|k|r 量级。
  S3 二阶（波动） (∂_t - ω∂_φ)² v = c²Δv：残差→0；去掉 ω 两项的残差 ~ 4ωσ|k|r。
  S4 离散一阶     C8 循环图、P 为循环移位：(Δ_τ+I-P⁻¹)v = P^{-(τ+1)}Δ_τ u（定理 6.1）。
  S5 离散二阶     同一图、u 满足 Δ_τ²u=-c²Lu(τ+1)：D_τ²v = -c²LP⁻¹v(τ+1)（定理 6.2）。

导数全部用中心差分（解析场采样，不做 PDE 时间推进）。
运行：python3 ZeroGPU/src/zb088_coriolis_connection.py
输出：ZeroGPU/results/zb088_coriolis_connection.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "zb088_coriolis_connection.json"

OMEGA, C, T0, DT, H = 1.0, 1.0, 0.3, 1e-5, 2.4 / 256
K = np.array([1.0, 0.5])
SIGMA = C * float(np.linalg.norm(K))


def R(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])



def kappa_dot(t):
    return -OMEGA * (J @ kappa(t))


def kappa(t):
    return R(OMEGA * t).T @ K


def Jmat():
    return np.array([[0.0, -1.0], [1.0, 0.0]])


J = Jmat()


def v_field(P, t):
    return np.cos(P @ kappa(t) - SIGMA * t)


def phase(P, t):
    return P @ kappa(t) - SIGMA * t


def analytic(P, t):
    """v 及其解析导数：v=cosθ, θ=κ(t)·y-σt。"""
    th = phase(P, t)
    c, s = np.cos(th), np.sin(th)
    k, kd, kdd = kappa(t), kappa_dot(t), -OMEGA ** 2 * kappa(t)
    A = P @ np.array([k[1], -k[0]])          # y₁κ₂ - y₂κ₁
    B = P @ kd - SIGMA                        # dθ/dt
    Bp = P @ kdd                              # dB/dt
    dv = -B * s
    dphi = A * s * -1.0                       # (Jy)·∇v = A·(-sinθ)
    dphi = -A * s
    dtt = -Bp * s - B ** 2 * c
    dphi_dt = (P @ np.array([kd[1], -kd[0]])) * (-s) + (-A) * (B * c)
    d2phi = (P @ k) * s - A ** 2 * c
    return {"v": c, "dt": dv, "dphi": dphi, "dtt": dtt, "dphi_dt": dphi_dt, "d2phi": d2phi}


def fd_check(P, t, h):
    """用中心差分近似 ∂_t v 与 ∂_φ v，检验二阶收敛。"""
    dv = (v_field(P, t + h) - v_field(P, t - h)) / (2 * h)
    gx = (v_field(P + np.array([h, 0.0]), t) - v_field(P - np.array([h, 0.0]), t)) / (2 * h)
    gy = (v_field(P + np.array([0.0, h]), t) - v_field(P - np.array([0.0, h]), t)) / (2 * h)
    dphi = -P[:, 1] * gx + P[:, 0] * gy
    return dv, dphi


def grid():
    ax = np.linspace(-1.2, 1.2, 256)
    X, Y = np.meshgrid(ax, ax, indexing="ij")
    P = np.stack([X.ravel(), Y.ravel()], axis=1)
    r = np.linalg.norm(P, axis=1)
    return P[(r > 0.45) & (r < 1.05)]


def _cycle_laplacian(n):
    A = np.zeros((n, n))
    for i in range(n):
        A[i, (i + 1) % n] = A[(i + 1) % n, i] = 1.0
    return np.diag(A.sum(1)) - A


def s4_discrete_first_order(n=8, T=10, seed=41):
    P = np.roll(np.eye(n), 1, axis=0)          # (Pu)_i = u_{i+1}
    L = _cycle_laplacian(n)
    rng = np.random.default_rng(seed)
    U = rng.normal(size=(T + 1, n))

    def Pinv(k):
        return np.linalg.matrix_power(P.T, k)

    def v(k):
        return Pinv(k) @ U[k]

    res = ctrl = 0.0
    for k in range(T):
        lhs = v(k + 1) - Pinv(1) @ v(k)
        rhs = Pinv(k + 1) @ (U[k + 1] - U[k])
        res = max(res, float(np.max(np.abs(lhs - rhs))))
        ctrl = max(ctrl, float(np.max(np.abs(v(k + 1) - v(k) - rhs))))
    return {"n": n, "T": T,
            "max_abs_residual": float(res),
            "control_without_connection": float(ctrl),
            "commutator_PLPt_minus_L": float(np.max(np.abs(P @ L @ P.T - L))),
            "pass": bool(res < 1e-12 and ctrl > 1e-2)}


def s5_discrete_second_order(n=8, T=10, c=0.5, seed=42):
    P = np.roll(np.eye(n), 1, axis=0)
    L = _cycle_laplacian(n)
    rng = np.random.default_rng(seed)
    U = np.zeros((T + 1, n))
    U[0], U[1] = rng.normal(size=n), rng.normal(size=n)
    for k in range(T - 1):                      # 离散波动方程：Δ_τ²u = -c²Lu(τ+1)
        U[k + 2] = 2.0 * U[k + 1] - U[k] - c * c * L @ U[k + 1]

    def Pinv(k):
        return np.linalg.matrix_power(P.T, k)

    def v(k):
        return Pinv(k) @ U[k]

    def Dv(k):
        return v(k + 1) - Pinv(1) @ v(k)

    res = ctrl = expn = alg = 0.0
    for k in range(T - 1):
        lhs = Dv(k + 1) - Pinv(1) @ Dv(k)
        rhs = -c * c * L @ (Pinv(1) @ v(k + 1))
        res = max(res, float(np.max(np.abs(lhs - rhs))))
        # 去掉连接项后的对照：只剩 Δ_τ²v
        ctrl = max(ctrl, float(np.max(np.abs((v(k + 2) - 2 * v(k + 1) + v(k)) - rhs))))
        # 展开 (Δ+B)² = Δ² + 2BΔ + B²
        dv = v(k + 1) - v(k)
        expanded = ((v(k + 2) - 2 * v(k + 1) + v(k))
                    + 2 * (dv - Pinv(1) @ dv)
                    + (v(k) - 2 * Pinv(1) @ v(k) + Pinv(2) @ v(k)))
        expn = max(expn, float(np.max(np.abs(lhs - expanded))))
        # 代数恒等式：D²v = P^{-(k+2)}Δ²u
        alg = max(alg, float(np.max(np.abs(
            lhs - Pinv(k + 2) @ (U[k + 2] - 2 * U[k + 1] + U[k])))))
    return {"n": n, "T": T, "c": c,
            "max_abs_residual": float(res),
            "max_abs_expansion": float(expn),
            "max_abs_algebraic": float(alg),
            "control_without_connection": float(ctrl),
            "pass": bool(max(res, expn, alg) < 1e-12 and ctrl > 1e-2)}


def main():
    P = grid()
    t = T0
    a = analytic(P, t)
    rhs = SIGMA * np.sin(phase(P, t))

    # S1：一阶恒等式（解析导数）
    res1 = float(np.max(np.abs(a["dt"] - OMEGA * a["dphi"] - rhs)))
    ctrl1 = float(np.max(np.abs(a["dt"] - rhs)))

    # S1b：有限差分的二阶收敛（步长减半，误差应 ÷4）
    dv_h, dphi_h = fd_check(P, t, 1e-3)
    err_h = float(np.max(np.abs(dv_h - a["dt"])))
    dv_h2, _ = fd_check(P, t, 5e-4)
    err_h2 = float(np.max(np.abs(dv_h2 - a["dt"])))

    # S3：二阶（解析）
    lhs2 = a["dtt"] - 2 * OMEGA * a["dphi_dt"] + OMEGA ** 2 * a["d2phi"]
    lap = -(K @ K) * a["v"]
    res2 = float(np.max(np.abs(lhs2 - C ** 2 * lap)))
    ctrl2 = float(np.max(np.abs(a["dtt"] - C ** 2 * lap)))

    res = {
        "grid_points": int(P.shape[0]),
        "params": {"omega": OMEGA, "c": C, "k": K.tolist(), "sigma": SIGMA},
        "S1_first_order": {"max_abs_residual": res1, "pass": bool(res1 < 1e-10),
                           "control_without_connection": ctrl1},
        "S1b_fd_convergence": {"err_h": err_h, "err_h_over_2": err_h2,
                               "ratio": float(err_h / max(err_h2, 1e-300))},
        "S2_direction": {"ratio_control_over_exact": float(ctrl1 / max(res1, 1e-300)),
                         "fails_as_expected": bool(ctrl1 > 1e-2)},
        "S3_second_order": {"max_abs_residual": res2, "pass": bool(res2 < 1e-10),
                            "control_without_two_terms": ctrl2},
        "S4_discrete_first_order": s4_discrete_first_order(),
        "S5_discrete_second_order": s5_discrete_second_order(),
    }
    res["PASS"] = all([res["S1_first_order"]["pass"],
                       res["S2_direction"]["fails_as_expected"],
                       res["S3_second_order"]["pass"],
                       res["S4_discrete_first_order"]["pass"],
                       res["S5_discrete_second_order"]["pass"]])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(res, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
