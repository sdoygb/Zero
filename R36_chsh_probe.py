#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R36_chsh_probe.py -- CHSH 探针：Zero 现有结构能否违反 Bell 不等式？

判据（Horodecki 1995，两比特）：S_max = 2*sqrt(u1^2 + u2^2)，
  u1,u2 = 关联矩阵 T_ij = Tr(rho sigma_i ⊗ sigma_j) 的两个最大奇异值。
  S_max > 2  <=>  纠缠（对两比特，CHSH 违反 = 纠缠判据）。

输出：R36_chsh_results.json
"""

from __future__ import annotations

import io
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
PAULI = [sx, sy, sz]


def correlation_matrix(rho):
    T = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            T[i, j] = np.trace(rho @ np.kron(PAULI[i], PAULI[j])).real
    return T


def s_max(rho):
    u = np.linalg.svd(correlation_matrix(rho), compute_uv=False)
    return 2.0 * float(np.sqrt(u[0] ** 2 + u[1] ** 2)), u


def bell_state():
    psi = np.zeros(4, dtype=complex)
    psi[0], psi[3] = 1 / np.sqrt(2), 1 / np.sqrt(2)   # (|00> + |11>)/sqrt2
    return np.outer(psi, psi.conj())


def werner(p):
    return p * bell_state() + (1 - p) * np.eye(4, dtype=complex) / 4


def product(a, b):
    return np.kron(a, b)


out = {}

# ---------------------------------------------------------------- S1 机制校验
out["S1_validation"] = {}
s_bell, _ = s_max(bell_state())
out["S1_validation"]["bell_state_S"] = s_bell
out["S1_validation"]["bell_expected"] = 2 * np.sqrt(2)
out["S1_validation"]["bell_ok"] = bool(abs(s_bell - 2 * np.sqrt(2)) < 1e-9)

prod = product(np.array([[1, 0], [0, 0]], dtype=complex), np.array([[1, 0], [0, 0]], dtype=complex))
s_prod, _ = s_max(prod)
out["S1_validation"]["product_state_S"] = s_prod

class_corr = 0.5 * (product(np.array([[1, 0], [0, 0]], dtype=complex), np.array([[1, 0], [0, 0]], dtype=complex))
                    + product(np.array([[0, 0], [0, 1]], dtype=complex), np.array([[0, 0], [0, 1]], dtype=complex)))
s_cc, _ = s_max(class_corr)
out["S1_validation"]["classically_correlated_S"] = s_cc

# Werner 阈值：p > 1/sqrt(2)
ps = np.linspace(0.0, 1.0, 2001)
thr = None
for p in ps:
    s, _ = s_max(werner(p))
    if s > 2 + 1e-9:
        thr = float(p)
        break
out["S1_validation"]["werner_threshold"] = thr
out["S1_validation"]["werner_threshold_expected"] = 1 / np.sqrt(2)

# ---------------------------------------------------------------- S2 G68 的多路径设定
a68 = np.array([0.6 + 0.3j, 0.5 - 0.2j, 0.2 - 0.1j, -0.35 + 0.15j])
rho_path = np.outer(a68, a68.conj())          # 单系统四条路径的纯态
rho_path = rho_path / np.trace(rho_path)
merge = [(0, 1), (2, 3)]
cross = float(sum(abs(a68[i] + a68[j]) ** 2 - (abs(a68[i]) ** 2 + abs(a68[j]) ** 2) for i, j in merge))
out["S2_G68_setting"] = {
    "n_paths": int(len(a68)),
    "norm_sq": float(np.sum(np.abs(a68) ** 2)),
    "merged_total_sq": float(sum(abs(a68[i] + a68[j]) ** 2 for i, j in merge)),
    "cross_total": cross,
    "bipartite": False,
    "why": "单系统多路径：两条『party』不是动力学独立且类空分离的子系统，故不构成 Bell 检验的场地",
}
# 强行把它当两比特（路径 {0,1} vs {2,3}）时的最大 CHSH（用 2x2 截断）
psi_pairs = np.zeros(4, dtype=complex)
psi_pairs[0] = a68[0]           # |0>_A ⊗ |0>_B
psi_pairs[1] = a68[1]           # |0>_A ⊗ |1>_B
psi_pairs[2] = a68[2]           # |1>_A ⊗ |0>_B
psi_pairs[3] = a68[3]           # |1>_A ⊗ |1>_B
psi_pairs = psi_pairs / np.linalg.norm(psi_pairs)
rho_pairs = np.outer(psi_pairs, psi_pairs.conj())
s_pair, u_pair = s_max(rho_pairs)
out["S2_G68_setting"]["forced_two_party_S"] = s_pair
out["S2_G68_setting"]["forced_two_party_singular_values"] = [float(x) for x in u_pair]

# ---------------------------------------------------------------- S3 G82 的两个独立 Z2
# 两个 Z2 的乘积态（独立约束 => 直积）：年龄奇偶 × 词的取向
z_age = np.array([1, 0], dtype=complex)      # 年龄奇偶的 +1 本征态
z_ori = np.array([1, 0], dtype=complex)      # 取向的 +1 本征态
prod_z2 = product(np.outer(z_age, z_age.conj()), np.outer(z_ori, z_ori.conj()))
s_z2, u_z2 = s_max(prod_z2)
# 若两 Z2 有一般关联（最大混合 + 关联 r）：rho = I/4 + r*(sigma_z⊗sigma_z)/4 等
def correlated_zz(rzz):
    rho = np.eye(4, dtype=complex) / 4
    rho = rho + rzz * np.kron(sz, sz) / 4
    return rho
s_zz, u_zz = s_max(correlated_zz(1.0))
out["S3_G82_two_Z2"] = {
    "product_state_S": s_z2,
    "product_singular_values": [float(x) for x in u_z2],
    "max_correlated_zz_S": s_zz,
    "max_correlated_zz_singular_values": [float(x) for x in u_zz],
    "note": "两个 Z2 由『独立动机』选出（G82: B=2x2），故联合态跨该二分割是直积；单方向关联 (Z⊗Z) 只给 rank-1 的 T，S 至多 2",
}

# ---------------------------------------------------------------- S4 违反所需的定量门槛
# 需要 T 的 rank>=2 且 u1^2+u2^2>1
def two_direction(r1, r2):
    """构造 T = diag(r1, r2, 0) 的态：rho = I/4 + (r1 sx⊗sx + r2 sy⊗sy)/4"""
    rho = np.eye(4, dtype=complex) / 4
    rho = rho + r1 * np.kron(sx, sx) / 4 + r2 * np.kron(sy, sy) / 4
    return rho

grid = []
for r in (0.5, 0.8, 0.9, 1.0):
    s, u = s_max(two_direction(r, r))
    grid.append({"r1=r2": r, "S": s, "u": [float(x) for x in u]})
out["S4_threshold"] = {
    "criterion": "S>2 iff u1^2+u2^2>1（T 需秩>=2 且两个方向都有足够关联）",
    "examples": grid,
    "werner_like_target": 1 / np.sqrt(2),
}

out["verdict"] = {
    "G68_setting_is_bipartite": False,
    "G82_two_Z2_state_entangled": False,
    "product_state_S": s_z2,
    "any_native_bipartition_violates": False,
    "conclusion": ("现有可检验的二分割上 CHSH<=2（Bell 局域）；"
                   "违反需要跨二分割的纠缠，而仓库的候选（空间二分割，G76/G78 的面积律）"
                   "尚无张量分解（D31 列为未恢复），故 CHSH 无法在此形式体系内被检验为‘违反’"),
    "quantitative_target": "若要违反，须存在二分割使 u1^2+u2^2>1（两比特等价于纠缠）；参考 Werner 门槛 p>1/sqrt(2)=0.7071",
}

with io.open(os.path.join(HERE, "R36_chsh_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)

print(json.dumps(out, ensure_ascii=False, indent=1))
