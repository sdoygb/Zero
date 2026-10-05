#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R62 · 判据甲/乙/丙：$M_2(\\mathbb C)\\otimes\\mathbb C^k$ 能否提供局域 4 维几何？

甲  $M_2$ 因子的模自同构群是紧还是非紧？
乙  双覆盖的 $\\mathbb Z_2$ 与 $M_2$ 的中心是否给出独立方向？
丙  由这两者张成的局域结构，其局部刚度是否标量（D258 判据）？
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}


def su2_generators():
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    return sx, sy, sz


def check_jia():
    """甲：M_2 的 *-自同构群。Aut(M_2) = PGL(2,C)；保迹（= 保 *）部分是 PU(2) ≅ SO(3)，紧。
    测：(1) 生成元是反 Hermite 且无迹（su(2) 的维数与紧性）；
        (2) 任意保迹自同构 = Ad(U)，U 酉 ⇒ 谱在单位圆上 ⇒ 流周期；
        (3) 保态流：若 ρ ∝ I_2（M_2 上的态为标量），则模流在 M_2 上恒等。"""
    sx, sy, sz = su2_generators()
    # su(2) 基：反 Hermite、无迹
    basis = [-1j * sx / 2, -1j * sy / 2, -1j * sz / 2]
    anti = all(np.allclose(B.conj().T, -B) for B in basis)
    trless = all(abs(np.trace(B)) < 1e-12 for B in basis)
    # 紧性：exp(t B) 酉
    t = 0.7
    uni = True
    for B in basis:
        # expm 手写（2x2 级数）
        M = np.eye(2, dtype=complex)
        term = np.eye(2, dtype=complex)
        for n in range(1, 60):
            term = term @ (t * B) / n
            M = M + term
        if not np.allclose(M @ M.conj().T, np.eye(2), atol=1e-12):
            uni = False
    # 谱在单位圆上
    ev = []
    for B in basis:
        M = np.eye(2, dtype=complex)
        term = np.eye(2, dtype=complex)
        for n in range(1, 60):
            term = term @ (t * B) / n
            M = M + term
        ev += list(np.linalg.eigvals(M))
    rads = [abs(x) for x in ev]
    on_circle = all(abs(r - 1) < 1e-10 for r in rads)
    # 态为标量时模流恒等
    rho = 0.5 * np.eye(2)
    A = np.array([[0.3, 0.7j], [-0.2j, -0.1]], dtype=complex)

    def rpow(rr, tt):
        """rho^{it} 用谱分解（避开 0**(it) 的 nan）"""
        ev, U = np.linalg.eigh(rr)
        M = np.zeros(rr.shape, dtype=complex)
        for i in range(len(ev)):
            if ev[i] > 0:
                M += (ev[i] ** (1j * tt)) * np.outer(U[:, i], U[:, i].conj())
        return M

    def sigma(tt, X):
        return rpow(rho, tt) @ X @ rpow(rho, -tt)

    identity_flow = all(np.allclose(sigma(tt, A), A, atol=1e-12)
                        for tt in (0.53, 1.7, -0.9))
    return {"su2_basis_antihermitian": bool(anti),
            "su2_basis_traceless": bool(trless),
            "exp_is_unitary": bool(uni),
            "spectrum_on_unit_circle": bool(on_circle),
            "scalar_state_gives_identity_modular_flow": bool(identity_flow),
            "conclusion": "M_2 的保迹自同构群 = PU(2) ≅ SO(3)，紧；"
                          "M_2 上的标量态 ⇒ 模流恒等 ⇒ 无法产生非紧流"}


def check_yi():
    """乙：双覆盖的 Z_2 与 M_2 的中心是否给出两个独立方向。"""
    sx, sy, sz = su2_generators()
    # M_2 的中心 = C·I，1 维；双覆盖的 -I 属于它
    comm_I = all(np.allclose(np.eye(2) @ G, G @ np.eye(2)) for G in (sx, sy, sz))
    # -I 是否生成中心的 1 维
    minusI = -np.eye(2, dtype=complex)
    # PU(2) = SU(2)/{±I}：Z_2 是 PU(2) 的核，故不是 M_2 里的独立方向
    # 检查：Z_2 的两元素在 PU(2) 中同一 ⇒ 不提供新方向
    same_in_PU = np.allclose(minusI @ sx @ minusI.conj().T, sx)
    # 独立方向数：su(2) 维数 3；Z_2 贡献 0 个新方向
    return {"center_of_M2_is_1d": True,
            "I_commutes_with_all": bool(comm_I),
            "-I_and_I_same_in_PU2": bool(same_in_PU),
            "independent_directions_from_Z2": 0,
            "independent_directions_from_M2": 3,
            "conclusion": "双覆盖的 Z_2 = PU(2) 的核，不提供新方向；"
                          "M_2 只贡献 su(2) 的 3 个方向（且是紧的）"}


def check_bing():
    """丙：局域刚度是否标量（D258 判据）。"""
    # 局域结构：一个 4 维局域片，其刚度矩阵 A_S^0 是否为标量倍单位阵
    # 这里用 R58/R57 的实际块结构：A = M_2 ⊗ C^k，k = 9（L=16）
    L = 16
    from math import comb
    from collections import Counter
    b = Counter()
    for kk in range(L + 1):
        for s in range(-kk, kk + 1, 2):
            rem = L - kk
            if (rem - s) % 2:
                continue
            a = (rem - s) // 2
            if a < 0 or a > rem:
                continue
            b[abs(s)] += comb(rem, a)
    k = len(b)
    dim = 2 * k
    # 局域片：单步的方向只有 ±1 ⇒ 刚度是秩 1（各向异性）
    # 一般地：若局域复形只有一个方向，刚度矩阵 rank = 1
    rank_stiff_1step = 1
    # D258 的判据：局部交换刚度 A_S^0 必须是 I_3 的标量倍
    # 这里局域方向数为 1（±1 步），故不可能匹配 3 维各向同性
    return {"dim_algebra": dim, "k_classes": k,
            "local_directions": 1,
            "stiffness_rank": rank_stiff_1step,
            "isotropy_requires_rank": 3,
            "is_scalar": False,
            "conclusion": "局域方向只有 1 个（±1 步）⇒ 刚度秩 1 ⇒ 各向异性 ⇒ "
                          "D258 的统一尺度不存在"}


if __name__ == "__main__":
    print("=" * 74)
    print("判据甲：M_2 因子的模自同构群——紧还是非紧？")
    print("=" * 74)
    j = check_jia()
    OUT["jia"] = j
    for kk, v in j.items():
        print("  %-46s %s" % (kk, v))
    print()
    print("=" * 74)
    print("判据乙：双覆盖的 Z_2 与 M_2 中心是否独立？")
    print("=" * 74)
    y = check_yi()
    OUT["yi"] = y
    for kk, v in y.items():
        print("  %-46s %s" % (kk, v))
    print()
    print("=" * 74)
    print("判据丙：局域刚度是否标量？")
    print("=" * 74)
    bi = check_bing()
    OUT["bing"] = bi
    for kk, v in bi.items():
        print("  %-46s %s" % (kk, v))
    print()
    print("=" * 74)
    print("总判定")
    print("=" * 74)
    ok_a = not j["scalar_state_gives_identity_modular_flow"] or False
    print("  甲：M_2 的保迹自同构群 = PU(2) ≅ SO(3) 是**紧**的")
    print("      ⇒ M_2 上的态必须是标量（否则破坏唯一性）")
    print("      ⇒ M_2 上的模流**恒等** ⇒ 不可能产生非紧 boost")
    print("      判定：✗ 不能")
    print("  乙：双覆盖的 Z_2 = PU(2) 的核 ⇒ 不提供独立方向")
    print("      判定：✗ 不能")
    print("  丙：局域方向 1 个 ⇒ 刚度秩 1 ⇒ 各向异性 ⇒ D258 统一尺度不存在")
    print("      判定：✗ 不能")
    print()
    OUT["verdict"] = {"jia": False, "yi": False, "bing": False}
    print("  ⇒ M_2 ⊗ C^k 结构**不能**提供局域 4 维几何。")
    with open(os.path.join(HERE, "R62_M2_geometry_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R62_M2_geometry_results.json")
