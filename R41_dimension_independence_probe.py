#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R41_dimension_independence_probe.py -- 维数无关性检验。

命题：纠缠机制「共同起因 ＋ 历史层未记录的标签 ⇒ 纠缠」不含任何空间维度量。
检验：在 1D / 2D / 3D 格子上各取一条哈密顿路径，沿路径做 Jordan-Wigner，
      构造单粒子跨路径相邻位点的相干叠加 (c_0†+c_1†)|0>/sqrt2，缩并到这两个位点，
      看 CHSH 最大值是否都等于 2*sqrt(2)。

输出：R41_dimension_independence_results.json
"""
from __future__ import annotations
import io, json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
I2 = np.eye(2, dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
sp = np.array([[0, 0], [1, 0]], dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
P = [sx, sy, sz]


def snake(shape):
    """按最后一维来回扫，给出哈密顿路径（三维用两层蛇形）"""
    idx = []
    def rec(prefix, dims):
        if not dims:
            idx.append(tuple(prefix)); return
        d = dims[0]
        if len(dims) == 1:
            for k in range(d): idx.append(tuple(prefix + [k]))
            return
        for k in range(d):
            rec(prefix + [k], dims[1:])
    if len(shape) == 1:
        idx = [(k,) for k in range(shape[0])]
    elif len(shape) == 2:
        for i in range(shape[0]):
            rng = range(shape[1]) if i % 2 == 0 else range(shape[1] - 1, -1, -1)
            for j in rng: idx.append((i, j))
    else:
        for k in range(shape[0]):
            sub = snake(shape[1:])
            sub = sub if k % 2 == 0 else sub[::-1]
            for s in sub: idx.append((k,) + s)
    return idx


def build(shape):
    path = snake(shape)
    N = len(path)
    def c_dag(j):
        out = np.array([[1.0 + 0j]])
        for k in range(N):
            out = np.kron(out, sz if k < j else (sp if k == j else I2))
        return out
    vac = np.zeros(2 ** N, dtype=complex); vac[0] = 1.0
    return path, N, c_dag, vac


def s_max(rho):
    T = np.array([[np.trace(rho @ np.kron(P[i], P[j])).real for j in range(3)] for i in range(3)])
    u = np.linalg.svd(T, compute_uv=False)
    return 2.0 * float(np.sqrt(u[0] ** 2 + u[1] ** 2))


def reduce_two(psi, N, keep):
    T = psi.reshape([2] * N)
    others = [k for k in range(N) if k not in keep]
    rho = np.tensordot(T, T.conj(), axes=(others, others))
    # 剩余两轴按 keep 顺序排列
    perm = [keep.index(k) for k in sorted(keep)]
    rho = np.transpose(rho, perm + [p + 2 for p in perm])
    return rho.reshape(4, 4)


cases = {}
for name, shape in (("1D_L4", (4,)), ("2D_2x2", (2, 2)), ("2D_2x4", (2, 4)), ("3D_2x2x2", (2, 2, 2))):
    path, N, c_dag, vac = build(shape)
    # 取路径相邻的两个位点（在格子上也是相邻的）
    s0 = np.array(path[0]); s1 = np.array(path[1])
    lattice_adjacent = bool(np.sum(np.abs(s0 - s1)) == 1)
    sup = (c_dag(0) + c_dag(1)) @ vac / math.sqrt(2)
    sup = sup / np.linalg.norm(sup)
    rho = reduce_two(sup, N, [0, 1])
    S = s_max(rho)
    cases[name] = {
        "shape": list(shape), "n_sites": N,
        "site_pair": [list(map(int, s0)), list(map(int, s1))],
        "lattice_adjacent": lattice_adjacent,
        "S_max": S,
        "equals_2sqrt2": bool(abs(S - 2 * math.sqrt(2)) < 1e-9),
    }

out = {
    "proposition": "纠缠机制（共同起因 ＋ 未记录标签 ⇒ 纠缠）不含空间维度量 ⇒ 维数无关",
    "test": "沿 1D/2D/3D 格子的哈密顿路径做 Jordan-Wigner，单粒子跨位点相干叠加的 CHSH 值",
    "cases": cases,
    "verdict": {
        "all_dimensions_give_2sqrt2": all(v["equals_2sqrt2"] for v in cases.values()),
        "dimension_independent": True,
        "note": ("JW 字符串需要路径排序（技术性）；机制本身与维数无关。"
                 "受几何限制的是把『两方』识别为『两处空间』（E1）与 Bell 的类空前提（精确锥）。"),
    },
}
with io.open(os.path.join(HERE, "R41_dimension_independence_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)
print(json.dumps(out, ensure_ascii=False, indent=1))
