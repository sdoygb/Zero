#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R34_mixing_generator_probe.py -- P2 探针：原生代数里能否实现 boost（混合生成元）？

问法（R33 §6 的第一击）：
  (Q1) 原生 M_2(C) 里是否存在 J,K 使 so(1,3) 关系 [J,J]=iJ、[J,K]=iK、[K,K]=-iJ 成立？
  (Q2) 若存在，K 是否 Hermitian（即能否生成**酉**流）？
  (Q3) 原生的模流 sigma_t = Ad(e^{itK_omega}) 能否等于 boost 流？

输出：R34_mixing_generator_results.json（供 R34 文档引用）
"""

from __future__ import annotations

import io
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [sx, sy, sz]
EPS = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1,
       (1, 0, 2): -1, (2, 1, 0): -1, (0, 2, 1): -1}


def comm(A, B):
    return A @ B - B @ A


def eps(i, j, k):
    return EPS.get((i, j, k), 0)


out = {}

# ---------------------------------------------------------------- Q1
# J_i = sigma_i/2（Hermitian），K_i = i sigma_i/2（anti-Hermitian）
J = [s / 2 for s in PAULI]
K = [1j * s / 2 for s in PAULI]

q1_jj = all(np.allclose(comm(J[i], J[j]), 1j * sum(eps(i, j, k) * J[k] for k in range(3)))
            for i in range(3) for j in range(3))
q1_jk = all(np.allclose(comm(J[i], K[j]), 1j * sum(eps(i, j, k) * K[k] for k in range(3)))
            for i in range(3) for j in range(3))
q1_kk = all(np.allclose(comm(K[i], K[j]), -1j * sum(eps(i, j, k) * J[k] for k in range(3)))
            for i in range(3) for j in range(3))
out["Q1_so13_realized"] = bool(q1_jj and q1_jk and q1_kk)
out["Q1_detail"] = {"J_J": bool(q1_jj), "J_K": bool(q1_jk), "K_K": bool(q1_kk)}

# ---------------------------------------------------------------- Q2
herm_J = all(np.allclose(J[i].conj().T, J[i]) for i in range(3))
herm_K = all(np.allclose(K[i].conj().T, K[i]) for i in range(3))
out["Q2_J_hermitian"] = bool(herm_J)
out["Q2_K_hermitian"] = bool(herm_K)

# Q2b: 是否存在**实数** c 使 K_i = c*sigma_i 满足 [K,K] = -i J？（Hermitian 解不存在）
real_c_works = []
for c in np.linspace(-5, 5, 2001):
    if abs(c) < 1e-12:
        continue
    Kc = [c * s for s in PAULI]
    if all(np.allclose(comm(Kc[i], Kc[j]), -1j * sum(eps(i, j, k) * J[k] for k in range(3)))
           for i in range(3) for j in range(3)):
        real_c_works.append(float(c))
out["Q2b_real_hermitian_boost_exists"] = bool(real_c_works)
out["Q2b_required_c_squared"] = -0.25  # 2c^2 = -1/2

# Q2c: e^{i theta K} 是否酉？取 K_1 = i sigma_1/2，则 i theta K = -theta sigma_1/2（实指数）
th = 1.0
U = np.array([[np.cosh(th / 2), -np.sinh(th / 2)], [-np.sinh(th / 2), np.cosh(th / 2)]])
out["Q2c_boost_exp_unitary"] = bool(np.allclose(U.conj().T @ U, np.eye(2)))
out["Q2c_boost_exp_norm"] = float(np.linalg.norm(U, 2))
out["Q2c_rotation_exp_unitary"] = bool(np.allclose(
    (lambda R: R.conj().T @ R)(np.array([[np.cos(0.5), -1j * np.sin(0.5)],
                                         [-1j * np.sin(0.5), np.cos(0.5)]])) , np.eye(2)))

# ---------------------------------------------------------------- Q3
# 原生模 Hamiltonian：K_omega = -log omega（G72 §1 的推前；用整数计数构造）
counts = np.array([1.0, 2.0, 10.0, 50.0])
omega = counts / counts.sum()
K_omega = np.diag(-np.log(omega))
out["Q3_K_omega_spectrum"] = [float(x) for x in np.diag(K_omega)]
out["Q3_K_omega_span"] = float(np.diag(K_omega).max() - np.diag(K_omega).min())
out["Q3_K_omega_hermitian"] = bool(np.allclose(K_omega.conj().T, K_omega))

# 模流的闭包紧性：有限谱 ⇒ 准周期 ⇒ 闭包是环面（紧）
evs = np.diag(K_omega)
periods = [2 * np.pi / abs(a - b) for i, a in enumerate(evs) for b in evs[i + 1:] if abs(a - b) > 1e-12]
out["Q3_modular_flow_compact"] = True
out["Q3_pair_periods_count"] = len(periods)

# 对比：boost 流的"闭包"非紧——单参数实指数族 ||e^{-theta sigma/2}|| -> inf
norms = [float(np.linalg.norm(np.array([[np.cosh(t / 2), -np.sinh(t / 2)],
                                        [-np.sinh(t / 2), np.cosh(t / 2)]]), 2)) for t in (0, 1, 5, 10)]
out["Q3_boost_norm_growth"] = norms
out["Q3_boost_flow_noncompact"] = bool(norms[-1] > 10 * norms[0])

# ---------------------------------------------------------------- Q4
# Killing 形式：so(4)（紧）负定 vs so(1,3)（非紧）不定
def so_basis(signature):
    """signature=(1,1,1,1) -> so(4)；signature=(1,-1,-1,-1) -> so(1,3)"""
    eta = np.diag(signature).astype(float)
    basis = []
    for a in range(4):
        for b in range(a + 1, 4):
            M = np.zeros((4, 4))
            M[a, b] = 1.0
            M[b, a] = -1.0
            # 保度量条件：M^T eta + eta M = 0 的无穷小形式
            gen = M @ eta
            basis.append(gen)
    return basis


def killing_signature(basis):
    n = len(basis)
    # 结构常数：c[k][i][j] with [e_i,e_j] = sum_k c^k_ij e_k
    B = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            ad_i = np.zeros((n, n))
            for k in range(n):
                C = basis[i] @ basis[k] - basis[k] @ basis[i]
                # 解 C = sum_l ad_i[l][k] basis[l]
                for l in range(n):
                    denom = np.vdot(basis[l], basis[l]).real
                    if denom > 1e-12 and np.allclose(C, (np.vdot(basis[l], C) / denom) * basis[l], atol=1e-9):
                        ad_i[l, k] = (np.vdot(basis[l], C) / denom).real
            ad_j = np.zeros((n, n))
            for k in range(n):
                C = basis[j] @ basis[k] - basis[k] @ basis[j]
                for l in range(n):
                    denom = np.vdot(basis[l], basis[l]).real
                    if denom > 1e-12 and np.allclose(C, (np.vdot(basis[l], C) / denom) * basis[l], atol=1e-9):
                        ad_j[l, k] = (np.vdot(basis[l], C) / denom).real
            B[i, j] = np.trace(ad_i @ ad_j)
    ev = np.linalg.eigvalsh((B + B.T) / 2)
    return ev


ev_so4 = killing_signature(so_basis((1, 1, 1, 1)))
ev_so13 = killing_signature(so_basis((1, -1, -1, -1)))
out["Q4_killing_eigs_so4"] = [float(x) for x in ev_so4]
out["Q4_killing_eigs_so13"] = [float(x) for x in ev_so13]
out["Q4_so4_negative_definite"] = bool(np.all(ev_so4 < 1e-9))
out["Q4_so13_indefinite"] = bool(np.any(ev_so13 > 1e-9) and np.any(ev_so13 < -1e-9))

# ---------------------------------------------------------------- 结论
out["verdict"] = {
    "mixing_generator_exists_algebraically": out["Q1_so13_realized"],
    "boost_can_be_hermitian": bool(out["Q2_K_hermitian"] or out["Q2b_real_hermitian_boost_exists"]),
    "boost_flow_unitary": bool(out["Q2c_boost_exp_unitary"]),
    "modular_flow_is_compact_at_finite_T": True,
    "so13_noncompact": out["Q4_so13_indefinite"],
    "finite_T_boost_as_modular_flow_possible": False,
    "target_moves_to": "细化极限须给出 type III（几何模流需 type III_1，BW）",
}

with io.open(os.path.join(HERE, "R34_mixing_generator_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)

print(json.dumps(out, ensure_ascii=False, indent=1))
