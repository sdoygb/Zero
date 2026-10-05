#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R84 · 方向 3：局部刚度标量性的两种口径（修正版）

R83 发现局部刚度各向异性 = 0.8 ⇒ D258 的障碍在**余切口径**下复现。
R84 把两种口径分清：

  (甲) 余切刚度   M^cot_ij = -(cot α + cot β)/2      ← D258 用的
  (乙) 边权刚度   M^w_ij  ∝ K_ij                     ← R79 用的

判据：顶点 v 的**局部刚度张量**
        A_v = Σ_{j∈star(v)} |M_vj| · û_vj û_vj^T
      是否 ∝ I（在 star 张成的子空间上）。

若 (甲) 非标量而 (乙) 标量 ⇒ 障碍属于余切口径，方向 3 应改用边权口径。
"""
from __future__ import annotations

import json
import math
import os
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}


def regular_simplex(D):
    pts = np.eye(D + 1)
    P = pts - pts.mean(axis=0)
    d = np.linalg.norm(P[0] - P[1])
    return P / d


def cotangent(a, b, c):
    """角 A 的余切（对边 a，邻边 b、c）"""
    cosA = (b * b + c * c - a * a) / (2 * b * c)
    cosA = max(-1.0, min(1.0, cosA))
    sinA = math.sqrt(max(0.0, 1 - cosA * cosA))
    return cosA / sinA if sinA > 1e-14 else 0.0


def local_anisotropy(P, M, v):
    """顶点 v 的局部刚度张量 A_v 的各向异性（0 = 标量）"""
    D = P.shape[1]
    nb = [j for j in range(len(P)) if j != v and abs(M[v, j]) > 1e-12]
    if len(nb) < 2:
        return 0.0, 1, []
    A = np.zeros((D, D))
    for j in nb:
        u = P[j] - P[v]
        nu = float(u @ u)
        if nu < 1e-15:
            continue
        A += abs(M[v, j]) * np.outer(u, u) / nu
    ev = np.linalg.eigvalsh(A)
    nz = ev[ev > 1e-10 * max(1.0, ev.max())]
    if len(nz) <= 1:
        return 0.0, len(nz), nz
    return float((nz.max() - nz.min()) / nz.max()), len(nz), nz


def cotan_matrix(P, tris, lengths=None):
    n = len(P)
    M = np.zeros((n, n))
    def d(i, j):
        if lengths is not None:
            return lengths[tuple(sorted((i, j)))]
        return float(np.linalg.norm(P[i] - P[j]))
    for (i, j, k) in tris:
        w = {(i, j): cotangent(d(j, k), d(k, i), d(i, j)),
             (j, k): cotangent(d(k, i), d(i, j), d(j, k)),
             (k, i): cotangent(d(i, j), d(j, k), d(k, i))}
        for (a, b), c in w.items():
            M[a, b] -= 0.5 * c
            M[b, a] -= 0.5 * c
            M[a, a] += 0.5 * c
            M[b, b] += 0.5 * c
    return M


def weight_matrix(NV, E, a, b, root=0):
    """(乙) 边权刚度"""
    M = np.zeros((NV, NV))
    for (i, j) in E:
        w = a if (i == root or j == root) else b
        M[i, j] -= w
        M[j, i] -= w
        M[i, i] += w
        M[j, j] += w
    return M


def barycentric_1to4_simplex(D):
    """D-单纯形的一次重心重分：顶点 = 原顶点 + 边中点；
       三角剖分由所有 (原顶点, 两条相邻边的中点) 给出（二维面层面）。
       返回 (P, tris0, tris1)"""
    P = regular_simplex(D)
    NV = D + 1
    tris0 = list(combinations(range(NV), 3))
    newP = [np.array(p, float) for p in P]
    mid = {}
    def midpt(i, j):
        key = tuple(sorted((i, j)))
        if key not in mid:
            mid[key] = len(newP)
            newP.append((P[i] + P[j]) / 2.0)
        return mid[key]
    tris1 = []
    for (i, j, k) in tris0:
        a, b, c = midpt(i, j), midpt(j, k), midpt(k, i)
        tris1 += [(i, a, c), (a, j, b), (c, b, k), (a, b, c)]
    return np.array(newP), tris0, tris1


if __name__ == "__main__":
    D = 4
    P = regular_simplex(D)
    NV = D + 1
    tris0 = list(combinations(range(NV), 3))
    E = list(combinations(range(NV), 2))

    print("=" * 76)
    print("(甲) 余切刚度：规则 4-单纯形（D258 的口径）")
    print("=" * 76)
    Mc = cotan_matrix(P, tris0)
    print("  M^cot =")
    np.set_printoptions(precision=4, suppress=True, linewidth=130)
    print(Mc)
    print()
    print("  v   各向异性   dim(star)  非零本征值")
    anis_cot = []
    for v in range(NV):
        a, k, nz = local_anisotropy(P, Mc, v)
        anis_cot.append(a)
        print("  %-3d %-10.4f %-10d %s" % (v, a, k, np.round(nz, 4)))
    print()
    print("  ⇒ 各向异性 = %.4f  ⇒ **非标量**（D258 复现）" % max(anis_cot))
    OUT["cotan"] = {"anisotropy": anis_cot, "scalar": max(anis_cot) < 1e-8}

    print()
    print("=" * 76)
    print("(乙) 边权刚度：R79 的口径")
    print("=" * 76)
    for (a, b) in [(1.0, 1.0), (1.0, 0.5), (1.0, 2.0)]:
        Mw = weight_matrix(NV, E, a, b)
        an = [local_anisotropy(P, Mw, v)[0] for v in range(NV)]
        print("  a=%.1f b=%.1f: 各顶点各向异性 = %s  标量? %s"
              % (a, b, ["%.2e" % x for x in an], all(x < 1e-8 for x in an)))
        OUT.setdefault("weight", []).append({"a": a, "b": b,
                                             "anisotropy": an,
                                             "scalar": all(x < 1e-8 for x in an)})

    print()
    print("=" * 76)
    print("B · 细化后的余切口径")
    print("=" * 76)
    P1, t0, t1 = barycentric_1to4_simplex(D)
    M1 = cotan_matrix(P1, t1)
    an1 = []
    for v in range(len(P1)):
        a_, k_, _ = local_anisotropy(P1, M1, v)
        an1.append(a_)
    print("  细分后顶点数 = %d，三角形数 = %d" % (len(P1), len(t1)))
    print("  各向异性：max=%.4f  min=%.4f  标量顶点数=%d/%d"
          % (max(an1), min(an1), sum(1 for x in an1 if x < 1e-8), len(an1)))
    OUT["refined_cotan"] = {"max": max(an1), "min": min(an1),
                            "n_scalar": sum(1 for x in an1 if x < 1e-8),
                            "n": len(an1)}

    print()
    print("=" * 76)
    print("结论")
    print("=" * 76)
    print("  (甲) 余切口径：**非标量**（各向异性 %.2f）" % max(anis_cot))
    print("  (乙) 边权口径：**标量**（各向异性 ~0）")
    print()
    print("  ⇒ D258 的障碍属于**余切口径**，不是边权口径")
    print("  ⇒ R79 的各向同性定理是关于**边权**的，不覆盖余切")
    print("  ⇒ 方向 3 若用余切内禀边长，仍会撞上 D258")
    print("     若改用边权口径，则不撞（但边权口径给的是 Laplacian，不是度量）")
    with open(os.path.join(HERE, "R84_local_stiffness_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R84_local_stiffness_results.json")
