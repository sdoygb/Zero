#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立复核 Zb041–Zb045（外部审校报告 + 文档自检）。
只做数值/代数核验，不修改正文。
"""
import json
import math
import numpy as np

out = {}

# ---------------------------------------------------------------- Zb041
# (A1) 实循环矩阵 vs 厄米循环矩阵：本征值形态
def circulant(first_row):
    n = len(first_row)
    return np.array([[first_row[(j - i) % n] for j in range(n)] for i in range(n)], dtype=complex)

a0, a1 = 1.0, 0.6 * np.exp(1j * 0.7)
C_herm = circulant([a0, a1, np.conj(a1)])
C_real = circulant([1.0, 0.5, -0.3])          # 一般实循环矩阵
eig_herm_num = np.sort(np.linalg.eigvalsh(C_herm).real)
eig_herm_formula = np.sort([a0 + 2 * abs(a1) * np.cos(np.angle(a1) + 2 * np.pi * (n - 1) / 3) for n in [1, 2, 3]])
out["Zb041_hermitian_circulant"] = {
    "eig_num": eig_herm_num.tolist(),
    "eig_formula": eig_herm_formula.tolist(),
    "max_abs_err": float(np.max(np.abs(eig_herm_num - eig_herm_formula))),
    "real_circulant_eig_imag_max": float(np.max(np.abs(np.imag(np.linalg.eigvals(C_real))))),
    "real_circulant_eig": [[float(z.real), float(z.imag)] for z in np.linalg.eigvals(C_real)],
}

# (A2) 命题 2.5 表：b=sqrt2
b = np.sqrt(2)
tab = []
for th_deg in [0, 12.73, 60, 90, 180]:
    th = np.radians(th_deg)
    brackets = [1 + b * np.cos(th + 2 * np.pi * (n - 1) / 3) for n in [1, 2, 3]]
    sqrt_m = [abs(x) for x in brackets]          # sqrt(m_n) = |A * bracket|
    m = [x ** 2 for x in sqrt_m]
    Q_num = sum(m) / (sum(sqrt_m) ** 2)
    Q_f = (1 + b ** 2 / 2) / 3
    tab.append({"theta_deg": th_deg, "brackets": brackets, "Q_num": Q_num, "Q_formula": Q_f})
out["Zb041_prop25_table"] = tab

# (A3) 轻子 DFT 读数
m_e, m_mu, m_tau = 0.51099895, 105.6583755, 1776.86
s = np.array([np.sqrt(m_e), np.sqrt(m_mu), np.sqrt(m_tau)])
s_hat_0 = np.mean(s)
k = np.array([0, 1, 2])
s_hat_1 = np.sum(s * np.exp(-2j * np.pi * k / 3)) / 3.0
arg1 = np.angle(s_hat_1)
if arg1 < 0:
    arg1 += 2 * np.pi
b_val = 2 * abs(s_hat_1) / abs(s_hat_0)
Q_val = (1 + b_val ** 2 / 2) / 3
out["Zb041_DFT"] = {
    "s_hat_0": float(s_hat_0),
    "arg_s_hat_1": float(arg1),
    "arg_minus_2pi3": float(arg1 - 2 * np.pi / 3),
    "two_ninths": 2 / 9,
    "abs_diff": float(arg1 - 2 * np.pi / 3 - 2 / 9),
    "rel_diff": float((arg1 - 2 * np.pi / 3 - 2 / 9) / (2 / 9)),
    "b": float(b_val),
    "Q": float(Q_val),
    "Q_direct": float(np.sum(s ** 2) / np.sum(s) ** 2),
}

# (A4) 汤川恒等式
v = 246.22e3  # MeV
y = np.sqrt(2) * np.array([m_e, m_mu, m_tau]) / v
K = float(np.sum(np.sqrt(y)) ** 2 / np.sum(y))
c_l = 313.856 / v
u_s = float(np.sum(y) / (np.sqrt(2) * c_l * 6))
out["Zb041_yukawa"] = {
    "y": y.tolist(),
    "sum_y": float(np.sum(y)),
    "K": K,
    "K_rel": float((K - 1.5) / 1.5),
    "c_l": float(c_l),
    "u_s": u_s,
    "1_minus_u_s": float(1 - u_s),
}

# ---------------------------------------------------------------- Zb042
# (B1) 均匀球：球心势 / 表面势 / 组装积分
R = 1.0
Q = 1.0
# 球心势：V(0) = (1/(4 pi eps0)) * (3/5? no) -> 3Q/(8 pi eps0 R)
# 用 eps0=1 归一：V(0)=3Q/(8 pi R), V(R)=Q/(4 pi R)
V_center = 3 * Q / (8 * np.pi * R)
V_surface = Q / (4 * np.pi * R)
# 组装积分：dW = V_surface(r) dq, q=Q r^3/R^3
r = np.linspace(1e-9, R, 200001)
q = Q * (r / R) ** 3
dqdr = 3 * Q * r ** 2 / R ** 3
E_num = getattr(np, 'trapezoid', np.trapz)((q / (4 * np.pi * r)) * dqdr, r)
E_exact = 3 * Q ** 2 / (20 * np.pi * R)  # 3/5 * Q^2/(4 pi eps0 R)
out["Zb042_potential"] = {
    "V_center_over_Q": float(V_center),
    "V_surface_over_Q": float(V_surface),
    "ratio_center_surface": float(V_center / V_surface),
    "E_num": float(E_num),
    "E_exact": float(E_exact),
    "rel_err_integral": float(abs(E_num - E_exact) / E_exact),
}

# (B2) 每槽能量与质量
v_gev = 246.22
A2_exact = v_gev * 1e3 / 784.5
A2_rounded = 313.856
theta_kappa = np.pi / 12 - 2 / 9
c_N = np.sqrt(5) * theta_kappa ** 2
for tag, A2 in [("exact", A2_exact), ("rounded", A2_rounded)]:
    q_slot = A2 * c_N
    mp = 3 * A2 - 3 * q_slot
    out[f"Zb042_slot_{tag}"] = {
        "A2": float(A2),
        "q_slot": float(q_slot),
        "2q": float(2 * q_slot),
        "8q": float(8 * q_slot),
        "2q_rel_vs_2.2246": float((2 * q_slot - 2.2246) / 2.2246),
        "8q_rel_vs_8.7946": float((8 * q_slot - 8.7946) / 8.7946),
        "mp": float(mp),
        "mp_rel_vs_938.2721": float((mp - 938.2721) / 938.2721),
        "theta_kappa": float(theta_kappa),
        "c_N": float(c_N),
    }
    out[f"Zb042_surface_sym_{tag}"] = {
        "aS_A": float(48.75 * 0.4 * q_slot),
        "aS_sets": float(38.7 * 0.4 * q_slot),
        "aA": float(20 * q_slot),
        "aS_over_aA": float(48.75 * 0.4 / 20),
    }

# (B3) 质子磁矩与质量比
cos2_15 = np.cos(np.pi / 12) ** 2
g_p = 2 * 3 * cos2_15 * (1 - np.sqrt(2) * theta_kappa ** 2)
M_e = 1 + np.sqrt(2) * np.cos(2 * np.pi / 3 + 2 / 9)
ratio = M_e ** 2 / (3 * (1 - np.sqrt(5) * theta_kappa ** 2))
m_e_exp, m_p_exp = 0.51099895, 938.27208816
exp_ratio = m_e_exp / m_p_exp
out["Zb042_gp_ratio"] = {
    "cos2_15": float(cos2_15),
    "g_p": float(g_p),
    "g_p_rel_vs_5.58569": float((g_p - 5.58569) / 5.58569),
    "M_e": float(M_e),
    "ratio": float(ratio),
    "exp_ratio": float(exp_ratio),
    "ratio_rel_vs_exp": float((ratio - exp_ratio) / exp_ratio),
    "ratio_rel_vs_doc_5.446166": float((ratio - 5.446166e-4) / 5.446166e-4),
}

# ---------------------------------------------------------------- Zb044
# (C1) 斜投影反例
E1 = np.array([[1.0, 0.0], [1.0, 0.0]])
E2 = np.array([[1.0, 0.0], [2.0, 0.0]])
out["Zb044_oblique_counterexample"] = {
    "E1_idempotent": bool(np.allclose(E1 @ E1, E1)),
    "E2_idempotent": bool(np.allclose(E2 @ E2, E2)),
    "E2E1_eq_E2": bool(np.allclose(E2 @ E1, E2)),
    "E1E2_eq_E2": bool(np.allclose(E1 @ E2, E2)),
    "E1E2": (E1 @ E2).tolist(),
    "E2": E2.tolist(),
}

# (C2) 正交投影下命题 3.3 的随机核验；并核验斜投影下会失败
rng = np.random.default_rng(20261031)
def random_orth_proj(dim):
    X = rng.normal(size=(dim, dim))
    Qm, _ = np.linalg.qr(X)
    r = int(rng.integers(1, dim))
    U = Qm[:, :r]
    return U @ U.T

orth_fail = 0
orth_trials = 300
for _ in range(orth_trials):
    d = 4
    P1 = random_orth_proj(d)
    P2 = random_orth_proj(d)
    # 只测试一个方向：若 P2 P1 = P2，则 ran P2 ⊆ ran P1
    if np.allclose(P2 @ P1, P2):
        # P1 P2 = P2 当且仅当 ran P2 ⊆ ran P1
        if not np.allclose(P1 @ P2, P2):
            orth_fail += 1
out["Zb044_orthogonal_prop33"] = {"trials": orth_trials, "failures": orth_fail}

# (C3) 斜投影下定理 3.1（乘积为投影 iff 交换）是否也需要自伴
# 随机搜索 2x2 幂等矩阵：E1E2 幂等但 E1E2 != E2E1？若存在，定理 3.1 的陈述必须限定正交投影。
def random_idempotent():
    # 构造 E = u v^T / (v^T u)，u,v 随机且 v^T u != 0
    u = rng.normal(size=2)
    v = rng.normal(size=2)
    denom = v @ u
    if abs(denom) < 1e-8:
        return None
    return np.outer(u, v) / denom

oblique_theorem31 = None
for _ in range(200000):
    A = random_idempotent()
    B = random_idempotent()
    if A is None or B is None:
        continue
    AB = A @ B
    if np.allclose(AB @ AB, AB) and not np.allclose(A @ B, B @ A):
        oblique_theorem31 = {"A": A.tolist(), "B": B.tolist(),
                             "AB": AB.tolist(),
                             "commutator": (A @ B - B @ A).tolist()}
        break
out["Zb044_oblique_theorem31_counterexample"] = oblique_theorem31

# ---------------------------------------------------------------- Zb043
cap = []
for T in range(2, 9):
    cT = int(math.comb(T, T // 2))
    cap.append({"T": T, "c_T": cT, "4(T+1)": 4 * (T + 1), "2(T+1)": 2 * (T + 1),
                "pass1": cT <= 4 * (T + 1), "pass2": cT <= 2 * (T + 1)})
out["Zb043_capacity"] = cap

print(json.dumps(out, ensure_ascii=False, indent=2))
with open("_staging/Zb041-045_复核_20261031.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
