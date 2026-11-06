#!/usr/bin/env python3
"""Zb078 -- L(K4) = octahedron: line graph spectrum, perfect matchings, antipodal axes.

Checks:
  (1) K4 line graph: A_L = B^T B - 2I; L(K4) is isomorphic to the octahedron K_{2,2,2}.
  (2) Spectra: Q(K4) = D + A, line graph, octahedron.
  (3) Perfect matchings of K4 = 3 = antipodal pairs of the octahedron.
  (4) Complement of the octahedron = 3K2 (Zb075 C01).

Output: ZeroGPU/results/zb078_k4_line_graph.json
"""
import itertools
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "ZeroGPU", "results", "zb078_k4_line_graph.json")


def k4():
    n = 4
    A = np.ones((n, n)) - np.eye(n)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    B = np.zeros((n, len(edges)))
    for k, (i, j) in enumerate(edges):
        B[i, k] = B[j, k] = 1.0
    return A, edges, B


def octahedron():
    n = 6
    grp = [0, 0, 1, 1, 2, 2]
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j and grp[i] != grp[j]:
                A[i, j] = 1.0
    return A


def is_iso(A, B):
    n = A.shape[0]
    for p in itertools.permutations(range(n)):
        P = np.zeros((n, n))
        for i, pi in enumerate(p):
            P[i, pi] = 1.0
        if np.allclose(P @ A @ P.T, B):
            return True, p
    return False, None


def perfect_matchings(n=4):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    out = []
    for e1, e2 in itertools.combinations(edges, 2):
        if len({e1[0], e1[1], e2[0], e2[1]}) == 4:
            out.append((e1, e2))
    return out


def main():
    A, edges, B = k4()
    A_L = B.T @ B - 2 * np.eye(len(edges))
    A_line = np.zeros((len(edges), len(edges)))
    for k, (a, b) in enumerate(edges):
        for l, (c, d) in enumerate(edges):
            if k != l and (a == c or a == d or b == c or b == d):
                A_line[k, l] = 1.0
    A_oct = octahedron()
    iso, perm = is_iso(A_L, A_oct)
    Q = B @ B.T
    pms = perfect_matchings(4)
    # non-adjacent pairs in L(K4): edges of K4 that are disjoint
    disjoint = [(k, l) for k in range(len(edges)) for l in range(k + 1, len(edges))
                if len(set(edges[k]) | set(edges[l])) == 4]
    # complement of octahedron
    comp = np.ones((6, 6)) - np.eye(6) - A_oct
    data = {
        "k4": {"n": 4, "edges": len(edges), "spec_A": np.round(np.linalg.eigvalsh(A), 9).tolist()},
        "Q_K4": np.round(np.linalg.eigvalsh(Q), 9).tolist(),
        "line_graph": {
            "A_L_eq_BtB_minus_2I": bool(np.allclose(A_L, A_line)),
            "spec": np.round(np.linalg.eigvalsh(A_L), 9).tolist(),
            "isomorphic_to_octahedron": bool(iso),
            "permutation": [int(x) for x in perm] if perm else None,
        },
        "octahedron": {"spec": np.round(np.linalg.eigvalsh(A_oct), 9).tolist()},
        "perfect_matchings": {"count": len(pms), "pairs": [[list(e1), list(e2)] for e1, e2 in pms]},
        "disjoint_edge_pairs": {"count": len(disjoint), "pairs": [list(x) for x in disjoint]},
        "complement_octahedron": {
            "spec": np.round(np.linalg.eigvalsh(comp), 9).tolist(),
            "is_3K2": bool(np.allclose(comp @ comp, np.eye(6))),
        },
    }
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
