#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R80 · 最后一块骨牌：T^{ab} 的各向同性是否把 βε 钉死在 4/3

思路：
  R79 给出 T^{ab} = (Σ K)δ^{ab} − K^{ab}，故各向同性 ⟺ K^{ab} 各向同性。
  而 K_{ij} 由**层高代价**给出：边 (i,j) 的"高度" = 它在两轨道中的归属。
    A 类（与根相连）：高度贡献 = D − 1 − 0 ？取决于模型
    B 类（不与根相连）：另一个高度
  一般地
        K_{ij} ∝ r^{h_{ij}},      r = e^{-βε}
  于是 a/b 是 βε 的函数。各向同性条件 3(a−b)=0 给出一个方程
        a(βε) = b(βε)
  本脚本检验：该方程的解是否 = 4/3。

检验多种"高度"定义（都是原生候选）：
  (i)   h = 边的两个端点的层高之和
  (ii)  h = min(端点层高)
  (iii) h = 端点在轨道中的序号之和
  (iv)  h = 边的"深度"= 它到根的最短距离
"""
from __future__ import annotations

import json
import math
import os
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}
D4 = 4
NV = D4 + 1
VERTS = list(range(NV))
EDGES = list(combinations(VERTS, 2))
ROOT = 0


def orbit(i, j):
    return "A" if (i == ROOT or j == ROOT) else "B"


def heights(model, i, j):
    """返回边 (i,j) 的高度（各模型的候选定义）"""
    # 端点层高：根 = 0，其余 = 1（完全图上所有非根顶点都在第 1 层）
    h = {v: (0 if v == ROOT else 1) for v in VERTS}
    if model == "sum":          # 端点层高之和
        return h[i] + h[j]
    if model == "min":
        return min(h[i], h[j])
    if model == "max":
        return max(h[i], h[j])
    if model == "index":        # 端点在轨道中的序号之和（A 类序号 / B 类序号）
        # A 类：与根相连，序号 = 对方的编号 1..4
        # B 类：不相连，序号 = 两端编号之和
        if orbit(i, j) == "A":
            return (i if i != ROOT else j) + 4      # 偏移 4 以区分轨道
        return i + j
    if model == "depth":        # 边到根的最短距离（跳数）
        return 0 if orbit(i, j) == "A" else 1
    raise ValueError(model)


def weights(model, g):
    """返回 (a, b) —— 两轨道的总权"""
    r = math.exp(-g)
    a = b = 0.0
    for (i, j) in EDGES:
        w = r ** heights(model, i, j)
        if orbit(i, j) == "A":
            a += w
        else:
            b += w
    return a, b


def solve_isotropy(model, lo=0.0, hi=6.0):
    """解 a(βε) = b(βε)"""
    f = lambda g: (lambda ab: ab[0] - ab[1])(weights(model, g))
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None
    for _ in range(300):
        mid = (lo + hi) / 2
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


if __name__ == "__main__":
    print("=" * 74)
    print("两轨道的结构（D=4, K_5）")
    print("=" * 74)
    A = [(i, j) for (i, j) in EDGES if orbit(i, j) == "A"]
    B = [(i, j) for (i, j) in EDGES if orbit(i, j) == "B"]
    print("  轨道 A（与根相连）: %d 条  %s" % (len(A), A))
    print("  轨道 B（其余）    : %d 条  %s" % (len(B), B))
    print()
    print("  ⇒ |A| = D = 4，|B| = C(D,2) = 6")
    OUT["orbits"] = {"A": len(A), "B": len(B)}

    print()
    print("=" * 74)
    print("各高度模型下：各向同性方程 a(βε) = b(βε) 的解")
    print("=" * 74)
    rows = []
    for model in ["sum", "min", "max", "index", "depth"]:
        hs = [heights(model, i, j) for (i, j) in EDGES]
        sol = solve_isotropy(model)
        # 各 βε 下的 a/b
        sample = []
        for g in [0.0, 0.5, 1.0, 4/3, 1.5, 2.0, 3.0]:
            aa, bb = weights(model, g)
            sample.append((g, aa / bb if bb else float('inf')))
        print("  模型 %-6s 高度取值=%s" % (model, sorted(set(hs))))
        if sol is None:
            print("     无解（a−b 不变号）")
        else:
            print("     解 βε* = %.6f   与 4/3 差 %+.6f (%+.2f%%)"
                  % (sol, sol - 4 / 3, 100 * (sol - 4 / 3) / (4 / 3)))
        print("     a/b 随 βε: %s" % ["%.3f" % x for _, x in sample])
        rows.append({"model": model, "heights": sorted(set(hs)),
                     "solution": sol,
                     "sample": [(g, float(x)) for g, x in sample]})
    OUT["models"] = rows

    print()
    print("=" * 74)
    print("判读")
    print("=" * 74)
    hits = [r for r in rows if r["solution"] is not None
            and abs(r["solution"] - 4 / 3) < 1e-6]
    near = [r for r in rows if r["solution"] is not None
            and abs(r["solution"] - 4 / 3) < 0.05]
    print("  精确给出 4/3 的模型数 = %d" % len(hits))
    print("  在 4/3 ± 0.05 内的模型数 = %d" % len(near))
    for r in rows:
        if r["solution"] is not None:
            print("    %-6s → %.4f" % (r["model"], r["solution"]))
        else:
            print("    %-6s → 无解" % r["model"])
    OUT["hits"] = [r["model"] for r in hits]
    OUT["near"] = [r["model"] for r in near]

    with open(os.path.join(HERE, "R80_beta_lock_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R80_beta_lock_results.json")
