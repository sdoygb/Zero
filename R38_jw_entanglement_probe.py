#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R38_jw_entanglement_probe.py -- 路线 B′：纠缠能否来自"共同起因 ＋ 未记录的位点自由度"？

设定（Z15 的 Jordan-Wigner 构造）：
  链上位点 j 带自旋 1/2；c_j = (prod_{k<j} sigma_k^z) sigma_j^-，c_j† = (prod_{k<j} sigma_k^z) sigma_j^+
  单粒子扇区：{|0>, c_j†|0>}（共同起因 = "存在一个费米子" = 共享的闭合记录）
判据（R36）：S_max = 2*sqrt(u1^2+u2^2)，S>2 <=> 两站点间纠缠

输出：R38_jw_entanglement_results.json
"""

from __future__ import annotations

import io
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

I2 = np.eye(2, dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)      # |0>=空=自旋下 -> sz|0> = +|0>
sp = np.array([[0, 0], [1, 0]], dtype=complex)       # sigma^+ = |1><0|
sm = sp.conj().T
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
PAULI = [sx, sy, sz]
N = 4


def op(single, site):
    """在第 site 位作用 single，其余为恒等（kron 顺序：site 0 在最左）"""
    out = np.array([[1.0 + 0j]])
    for k in range(N):
        out = np.kron(out, single if k == site else I2)
    return out


def c_dag(j):
    """c_j† = (prod_{k<j} sz_k) sp_j"""
    out = np.array([[1.0 + 0j]])
    for k in range(N):
        if k < j:
            out = np.kron(out, sz)
        elif k == j:
            out = np.kron(out, sp)
        else:
            out = np.kron(out, I2)
    return out


vac = np.zeros(2 ** N, dtype=complex)
vac[0] = 1.0                                        # 所有位点 |0>（空）

basis = [("vacuum", vac)] + [("1p_site%d" % j, c_dag(j) @ vac) for j in range(N)]
superpos = (c_dag(0) + c_dag(1)) @ vac / np.sqrt(2)  # 跨位点 0,1 的相干叠加
mix = 0.5 * np.outer(c_dag(0) @ vac, (c_dag(0) @ vac).conj()) \
    + 0.5 * np.outer(c_dag(1) @ vac, (c_dag(1) @ vac).conj())


def reduce_two_sites(psi):
    """把位点 2..N-1 缩并掉，得到位点 0,1 的 2 比特密度矩阵"""
    T = psi.reshape([2] * N)
    rho = np.tensordot(T, T.conj(), axes=([list(range(2, N)), list(range(2, N))]))
    return rho.reshape(4, 4)


def s_max(rho):
    Tm = np.array([[np.trace(rho @ np.kron(PAULI[i], PAULI[j])).real for j in range(3)] for i in range(3)])
    u = np.linalg.svd(Tm, compute_uv=False)
    return 2.0 * float(np.sqrt(u[0] ** 2 + u[1] ** 2)), u


out = {}
rows = []
for name, psi in basis + [("superposition_sites0_1", superpos)]:
    psi = psi / np.linalg.norm(psi)
    rho = reduce_two_sites(psi)
    S, u = s_max(rho)
    rows.append({"state": name, "S_max": S, "singular_values": [float(x) for x in u],
                 "entangled": bool(S > 2 + 1e-9)})

rho_mix = reduce_two_sites(np.zeros(2 ** N, dtype=complex))  # 占位，下面直接算
# 混合态：直接由两站点构造（缩并对角混合）
rho_mix2 = 0.5 * reduce_two_sites(basis[1][1]) + 0.5 * reduce_two_sites(basis[2][1])
S_mix, u_mix = s_max(rho_mix2)
rows.append({"state": "classical_mixture_(which_site_recorded)", "S_max": S_mix,
             "singular_values": [float(x) for x in u_mix], "entangled": bool(S_mix > 2 + 1e-9)})

out["states"] = rows
out["bell_expected"] = 2 * np.sqrt(2)
sup_row = next(r for r in rows if r["state"] == "superposition_sites0_1")
out["verdict"] = {
    "single_particle_superposition_across_two_sites_is_entangled":
        bool(sup_row["entangled"]),
    "S_superposition": sup_row["S_max"],
    "classical_mixture_S": S_mix,
    "which_site_record_decides": bool(S_mix <= 2 + 1e-9),
    "shared_origin": "同一个费米子（费米子数/宇称为共享的闭合记录）",
    "unrecorded_dof": "位点指标（Jordan-Wigner 字符串把这个自由度变成非局域）",
    "conclusion": ("共同起因（同一费米子）＋ 未记录的位点自由度 ⇒ 跨位点纠缠（S=2√2）；"
                   "一旦位点被记录（经典混合）⇒ S=2"),
}

with io.open(os.path.join(HERE, "R38_jw_entanglement_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)

print(json.dumps(out, ensure_ascii=False, indent=1))
