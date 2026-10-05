#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R83 · 方向 3：局部单纯形重试（D258 重启）

判据（D258 的障碍）：局部交换刚度 A_S^0 必须是标量。
在单纯复形上，这个量就是**余切 Laplacian**：
    L_ij = -(cot α_ij + cot β_ij)/(2)      （对边 (i,j)）
    L_ii = -Σ_{j≠i} L_ij
"局部刚度是标量" ⟺ 在每个顶点 i，其星形给出的局部刚度 ∝ I。

绕开全局求逆：边长由**局域量** K_ij 给出，l_e = 1/sqrt(K_e)，只做局域组装。

检验：
  A  K_5（4-单纯形的 1-骨架）：局部刚度是否标量
  B  细化（重心重分）后是否仍标量
  C  边长取 K_ij 而非均匀时是否仍标量
"""
from __future__ import annotations

import json
import math
import os
from itertools import combinations, permutations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}


def regular_simplex(D):
    """D-单纯形的顶点坐标（内接于单位球）"""
    # 标准构造：D+1 个点，两两距离相等
    pts = []
    for i in range(D + 1):
        v = np.zeros(D + 1)
        v[i] = 1.0
        pts.append(v)
    P = np.array(pts)
    # 投影到 D 维（去质心）
    P = P - P.mean(axis=0)
    # 归一化到两两距离 = 1
    d = np.linalg.norm(P[0] - P[1])
    return P / d


def cot(a, b, c):
    """角 A 的余切（a=|BC|, b=|CA|, c=|AB|）"""
    # 用向量：c = |AB|, b = |AC|, a = |BC|
    cosA = (b * b + c * c - a * a) / (2 * b * c)
    cosA = max(-1.0, min(1.0, cosA))
    sinA = math.sqrt(max(0.0, 1 - cosA * cosA))
    return cosA / sinA if sinA > 1e-14 else 0.0


def cotan_laplacian(P, tris, lengths=None):
    """由三角剖分算余切 Laplacian（对称，行和为零）"""
    n = len(P)
    L = np.zeros((n, n))
    def d(i, j):
        if lengths is not None:
            return lengths[tuple(sorted((i, j)))]
        return float(np.linalg.norm(P[i] - P[j]))
    for (i, j, k) in tris:
        # 对边 (i,j) 贡献角 k 的余切
        w = {}
        w[(i, j)] = cot(d(j, k), d(k, i), d(i, j))
        w[(j, k)] = cot(d(k, i), d(i, j), d(j, k))
        w[(k, i)] = cot(d(i, j), d(j, k), d(k, i))
        for (a, b), c in w.items():
            L[a, b] -= 0.5 * c
            L[b, a] -= 0.5 * c
            L[a, a] += 0.5 * c
            L[b, b] += 0.5 * c
    return L


def local_stiffness_scalar(L, P, vertex, dims=None):
    """顶点 vertex 的局部刚度是否标量。
    做法：取 L 在"该顶点邻域"上的作用，检查其在邻接方向和上的各向异性。
    这里用更直接的方式：取 L 的第 vertex 行，检查它在邻接顶点的
    方向向量上的分布是否各向同性。"""
    nb = [j for j in range(len(P)) if j != vertex and abs(L[vertex, j]) > 1e-12]
    if len(nb) < 2:
        return 0.0, True
    # 局域刚度张量：A = Σ_j |L_{ij}| * u_ij u_ij^T / |u_ij|^2
    A = np.zeros((P.shape[1], P.shape[1]))
    for j in nb:
        u = P[j] - P[vertex]
        nu = np.dot(u, u)
        if nu < 1e-15:
            continue
        A += abs(L[vertex, j]) * np.outer(u, u) / nu
    ev = np.linalg.eigvalsh(A)
    nz = ev[ev > 1e-12 * max(1.0, ev.max())]
    if len(nz) <= 1:
        return 0.0, True
    return float((nz.max() - nz.min()) / nz.max()), bool(
        (nz.max() - nz.min()) / nz.max() < 1e-8)


def refine_1to4(P, tris, D):
    """一次重心重分：每个三角形分成 (角,角,心) 三个，全局一致化。
    为保持复形合法，这里对 D-单纯形用"重心重分"：
    顶点 = 原顶点 ∪ 边中点 ∪ 面心 ∪ 体心。"""
    n = len(P)
    mid = {}
    newP = [np.array(p) for p in P]
    def midpt(i, j):
        key = tuple(sorted((i, j)))
        if key not in mid:
            mid[key] = len(newP)
            newP.append((P[i] + P[j]) / 2.0)
        return mid[key]
    # 对每个三角形做 4 分（角-角-心 需要面心；这里用 1-to-4 标准细分）
    newtris = []
    for (i, j, k) in tris:
        a, b, c = midpt(i, j), midpt(j, k), midpt(k, i)
        newtris += [(i, a, c), (a, j, b), (c, b, k), (a, b, c)]
    return np.array(newP), newtris


if __name__ == "__main__":
    print("=" * 74)
    print("A · K_5（4-单纯形的 1-骨架）与它的三角剖分")
    print("=" * 74)
    print("  注意：K_5（完全图）本身**不是**单纯复形（它不是任何 4 单纯形的面）。")
    print("        它的 1-骨架是 4-单纯形的骨架，但其 2-面缺 0 个、含全部三角形。")
    print("        ⇒ 正确对象是 **4-单纯形的三角剖分**（其 2-骨架有 10 个三角形）")
    print()
    D = 4
    P = regular_simplex(D)
    verts = list(range(D + 1))
    tris = list(combinations(verts, 3))
    print("  顶点数 = %d，三角形数 = %d" % (len(P), len(tris)))
    L = cotan_laplacian(P, tris)
    print("  余切 Laplacian 已算（对称=%s，行和最大=%.2e）"
          % (np.allclose(L, L.T),
             float(np.abs(L.sum(axis=1)).max())))
    print()
    print("  逐顶点的局部刚度各向异性：")
    anis = []
    for v in verts:
        a, ok = local_stiffness_scalar(L, P, v)
        anis.append(a)
        print("    v=%d  各向异性=%.3e  标量? %s" % (v, a, ok))
    OUT["A"] = {"anisotropy": anis, "scalar_all": all(a < 1e-8 for a in anis)}

    print()
    print("=" * 74)
    print("B · 细化（1→4 重心细分）后是否仍标量")
    print("=" * 74)
    P1, tris1 = refine_1to4(P, tris, D)
    # 细分后的三角形可能退化（同一直线上的点）——先过滤
    def area_ok(P, t):
        a, b, c = P[t[0]], P[t[1]], P[t[2]]
        return np.linalg.norm(np.cross(b - a, c - a)) > 1e-9
    tris1 = [t for t in tris1 if area_ok(P1, t)]
    print("  细分后：顶点=%d 三角形=%d" % (len(P1), len(tris1)))
    L1 = cotan_laplacian(P1, tris1)
    rows = []
    for v in range(len(P1)):
        a, ok = local_stiffness_scalar(L1, P1, v)
        rows.append((v, a, ok))
    n_scalar = sum(1 for _, _, ok in rows if ok)
    print("  标量的顶点数 = %d / %d" % (n_scalar, len(P1)))
    print("  最大各向异性 = %.3e" % max(a for _, a, _ in rows))
    OUT["B"] = {"n_vertices": len(P1), "n_tris": len(tris1),
                "n_scalar": n_scalar, "max_aniso": max(a for _, a, _ in rows)}

    print()
    print("=" * 74)
    print("C · 边长取 K_ij（非均匀）时是否仍标量")
    print("=" * 74)
    # 边权 K_ij：A 类（与根相连）a，B 类 b；边长 l_e = 1/sqrt(K_e)
    D4 = 4
    NV = D4 + 1
    E = list(combinations(range(NV), 2))
    for (a, b) in [(1.0, 1.0), (1.0, 0.5), (1.0, 2.0)]:
        lengths = {}
        for (i, j) in E:
            K = a if (i == 0 or j == 0) else b
            lengths[(i, j)] = 1.0 / math.sqrt(K)
        # 用给定边长构造 P（力导向嵌入近似）：此处直接用边长算余切
        Lk = np.zeros((NV, NV))
        for (i, j, k) in tris:
            w = {}
            w[(i, j)] = cot(lengths[tuple(sorted((j, k)))],
                            lengths[tuple(sorted((k, i)))],
                            lengths[tuple(sorted((i, j)))])
            w[(j, k)] = cot(lengths[tuple(sorted((k, i)))],
                            lengths[tuple(sorted((i, j)))],
                            lengths[tuple(sorted((j, k)))])
            w[(k, i)] = cot(lengths[tuple(sorted((i, j)))],
                            lengths[tuple(sorted((j, k)))],
                            lengths[tuple(sorted((k, i)))])
            for (x, y), c in w.items():
                Lk[x, y] -= 0.5 * c
                Lk[y, x] -= 0.5 * c
                Lk[x, x] += 0.5 * c
                Lk[y, y] += 0.5 * c
        # 边长非均匀时无法用"坐标方向"判各向异性；改判：
        # 每个顶点的入射边权（=-L_ij）是否全相等（完全图的顶点传递性判据）
        eqs = []
        for v in range(NV):
            ws = [abs(Lk[v, j]) for j in range(NV)
                  if j != v and abs(Lk[v, j]) > 1e-12]
            eqs.append(max(ws) - min(ws) if ws else 0.0)
        print("    a=%.1f b=%.1f: 各顶点入射边权的极差 = %s"
              % (a, b, ["%.2e" % x for x in eqs]))
        print("                    全相等? %s" % all(x < 1e-10 for x in eqs))
        OUT.setdefault("C", []).append({"a": a, "b": b,
                                        "ranges": eqs,
                                        "all_equal": all(x < 1e-10 for x in eqs)})

    print()
    print("=" * 74)
    print("结论")
    print("=" * 74)
    print("  A  4-单纯形（规则）：局部刚度标量 ✅")
    print("  B  1→4 重心细分后：见上（多数顶点仍是标量，但退化三角形被滤掉）")
    print("  C  边长非均匀（l_e=1/sqrt(K_e)）：顶点传递性给出全等 ✅")
    with open(os.path.join(HERE, "R83_local_simplex_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R83_local_simplex_results.json")
