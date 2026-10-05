#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R78 · 物质侧：从补偿移动构造微观 T_{ab}

构造（乙）：能量型
    补偿移动 T_{ex}: x -> x + e_j - e_i 由通道对 (i,j) 指标化。
    在 K_{D+1} 的 1-骨架上，每条通道对给一个**方向**。
    设 K_{ij} 为该对上的转移率（由测度/续接数给）。
    定义微观应力张量：
        T^{ab} = Σ_{(i,j)} K_{ij} · u^a_{(ij)} u^b_{(ij)}
    其中 u_{(ij)} 是通道对 (i,j) 在通道空间中的方向向量。

检验：
  A  守恒性：由逐边配平（Z1 定理 1），Σ_j T^{ab} 的散度是否为零
  B  形状：在 Stab(0)=S_D 下，T^{ab} 的非对角部分是否为零（完美流体预言）
  C  迹与偏迹：trace 是否对应单守恒荷；偏迹是否为零
  D  方程状态：w = p/ρ 是否由结构定出
"""
from __future__ import annotations

import json
import math
import os
from fractions import Fraction
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}
D4 = 4
NV = D4 + 1
VERTS = list(range(NV))
EDGES = list(combinations(VERTS, 2))


def direction_vector(i, j, nv):
    """通道对 (i,j) 的方向：在 R^{nv} 中取 e_j - e_i"""
    v = np.zeros(nv)
    v[i] = -1.0
    v[j] = 1.0
    return v


def build_T(weights=None, nv=NV):
    """T^{ab} = Σ_{(i,j)} K_{ij} u u^T"""
    T = np.zeros((nv, nv))
    for idx, (i, j) in enumerate(EDGES):
        K = 1.0 if weights is None else weights[idx]
        u = direction_vector(i, j, nv)
        T += K * np.outer(u, u)
    return T


if __name__ == "__main__":
    print("=" * 74)
    print("A · 均匀权重下的 T^{ab}（K_5，D=4）")
    print("=" * 74)
    T = build_T()
    np.set_printoptions(precision=4, suppress=True, linewidth=120)
    print("  T^{ab} =")
    print(T)
    print()
    ev = np.linalg.eigvalsh(T)
    print("  本征值 = %s" % np.round(ev, 6))
    print("  rank(T) = %d" % int(np.sum(np.abs(ev) > 1e-10)))
    print()
    print("  在 Stab(0)=S_4 下（置换顶点 1..4，固定 0）检验 S_4 不变性：")
    from itertools import permutations
    base = T.copy()
    maxdev = 0.0
    for p in permutations([1, 2, 3, 4]):
        perm = [0] + list(p)
        P = np.zeros((NV, NV))
        for a in range(NV):
            P[perm[a], a] = 1.0
        Tp = P @ base @ P.T
        maxdev = max(maxdev, float(np.abs(Tp - base).max()))
    print("    max|P T P^T − T| = %.3e" % maxdev)
    OUT["A"] = {"eigenvalues": [float(x) for x in ev],
                "stabilizer_invariance": maxdev}

    print()
    print("=" * 74)
    print("B · 形状：非对角（各向异性）部分")
    print("=" * 74)
    off = T - np.diag(np.diag(T))
    print("  非对角部分 max|T_off| = %.3e" % np.abs(off).max())
    print("  对角部分 = %s" % np.round(np.diag(T), 6))
    # 在“通道空间”里，各向异性指: T 在 S_4 不可约表示下的分解
    # S_4 作用在 4 维（顶点 1..4）上：分解为 trivial + standard(3维)
    # 拆掉根 0 那一维后，看 4x4 块
    Tsub = T[1:, 1:]
    tr = np.trace(Tsub) / 4.0
    iso = tr * np.eye(4)
    aniso = Tsub - iso
    print()
    print("  限制到 4 个非根通道（4x4 块）：")
    print("    迹/4 = %.6f" % tr)
    print("    各向异性部分 max|aniso| = %.3e" % np.abs(aniso).max())
    print("    ⇒ 各向异性 = %s" % ("零（完美流体型）" if np.abs(aniso).max() < 1e-10
                                 else "非零"))
    OUT["B"] = {"aniso_max": float(np.abs(aniso).max()),
                "trace_over_4": float(tr)}

    print()
    print("=" * 74)
    print("C · 守恒性：散度是否为零")
    print("=" * 74)
    # 逐边配平：对每个通道 i，Σ_j (x_i - x_j) = 0 当且仅当 x 为常量
    # T 的散度 ∂_a T^{ab}：在通道空间即 Σ_a T^{ab}
    div = T.sum(axis=0)
    print("  Σ_a T^{ab} = %s" % np.round(div, 10))
    print("  max|div| = %.3e" % np.abs(div).max())
    print()
    print("  物理读法：均匀转移率下，每个通道收到的净动量流相消")
    print("  ⇒ 守恒（对应于 Z1 定理 1 的逐边配平）")
    print("  ⇒ 散度为零 = %s" % ("是" if np.abs(div).max() < 1e-10 else "否"))
    OUT["C"] = {"div_max": float(np.abs(div).max())}

    print()
    print("=" * 74)
    print("D · 方程状态 w = p/ρ")
    print("=" * 74)
    # 单守恒荷给出的是"迹"部分；各向同性给出压力
    # 取通道 0（根）为时间方向的候选
    rho = T[0, 0]
    p_spatial = np.mean(np.diag(T)[1:])
    print("  T^{00}（根通道）      = %.6f" % rho)
    print("  空间对角平均          = %.6f" % p_spatial)
    print("  ⇒ w = p/ρ = %.6f" % (p_spatial / rho if rho else float('nan')))
    print()
    # 换一个自然读法：把 D 个非根通道当空间，根当时间
    # T^{00} = Σ K_{0j} * 1  （根与每个通道配对）
    # 空间对角 T^{jj} = K_{0j} + Σ_{k≠0,j} K_{jk}
    t00 = sum(1 for (a, b) in EDGES if a == 0 or b == 0)
    print("  更精确：T^{00} = 与根相连的通道对数 = %d" % t00)
    tjj = []
    for j in range(1, NV):
        cnt = sum(1 for (a, b) in EDGES if (a == j) ^ (b == j))
        tjj.append(cnt)
    print("  空间对角 T^{jj} = 与通道 j 相连的对数 = %s" % tjj)
    print("  ⇒ 空间各向同性（全相等）= %s" % (len(set(tjj)) == 1))
    w = (sum(tjj) / len(tjj)) / t00
    print("  ⇒ w = p/ρ = %.6f" % w)
    OUT["D"] = {"T00": int(t00), "Tjj": tjj, "w": float(w)}

    print()
    print("=" * 74)
    print("总判定")
    print("=" * 74)
    print("  A  S_4 不变性（Stab(0)）：%s" % ("成立" if maxdev < 1e-10 else "不成立"))
    print("  B  各向异性为零（完美流体预言）：%s"
          % ("成立" if np.abs(aniso).max() < 1e-10 else "不成立"))
    print("  C  散度为零（守恒）：%s"
          % ("成立" if np.abs(div).max() < 1e-10 else "不成立"))
    print("  D  w = p/ρ = %.4f" % w)
    OUT["verdict"] = {"stabilizer": maxdev < 1e-10,
                      "isotropic": bool(np.abs(aniso).max() < 1e-10),
                      "conserved": bool(np.abs(div).max() < 1e-10),
                      "w": float(w)}

    with open(os.path.join(HERE, "R78_stress_tensor_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R78_stress_tensor_results.json")
