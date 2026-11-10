#!/usr/bin/env python3
"""Zb081 -- structure conjugation, concrete models for C13/C14/C24.

C13 matroid dual: linear matroid of the 2-core Jacobian + K4 graphic matroid.
C14 poset dual: 2-layer / 4-layer chains, ideals <-> filters by complement.
C24 layer loop: record->seed writeback K and its Moore-Penrose inverse K+;
    K K+ K = K (Zb045 general closure), quotient dim = rank K.

Output: ZeroGPU/results/zb081_conj_models.json
"""
import json
import os
from itertools import combinations

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "ZeroGPU", "results", "zb081_conj_models.json")


def c13_matroid():
    # 2-core: variables (theta_M, theta_C, theta_I, kappa0, eta0);
    # equations: angle sum (1,1,1,0,0), product kappa0*eta0 = const.
    # linearize at a point with kappa0=eta0=1: row2 = (0,0,0,1,1)
    B = np.array([[1.0, 1, 1, 0, 0], [0, 0, 0, 1, 1]])
    r = int(np.linalg.matrix_rank(B))
    # K4 graphic matroid: 6 edges, rank V-1 = 3
    n, m = 4, 6
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    incidence = np.zeros((n, m))
    for k, (i, j) in enumerate(edges):
        incidence[i, k] = -1.0
        incidence[j, k] = 1.0
    rK4 = int(np.linalg.matrix_rank(incidence))
    return {
        "two_core": {"n_vars": 5, "n_eqs": 2, "rank": r,
                     "dual_rank": 5 - r, "dof": 5 - r,
                     "B": B.astype(int).tolist()},
        "K4_graphic": {"n_edges": m, "rank": rK4, "dual_rank": m - rK4},
    }


def ideals_of_chain(k):
    # chain 0<1<...<k-1; ideals are initial segments
    return [frozenset(range(i)) for i in range(k + 1)]


def filters_of_chain(k):
    # filters are final segments
    return [frozenset(range(i, k)) for i in range(k + 1)]


def c14_poset():
    out = {}
    for k in [2, 4]:
        I = ideals_of_chain(k)
        F = filters_of_chain(k)
        comp = {frozenset(range(k)) - s for s in I}
        out[f"chain_{k}"] = {"ideals": len(I), "filters": len(F),
                             "complement_bijection": comp == set(F)}
    return out


def c24_layer():
    # record -> seed writeback K (rank 2, kernel dim 1)
    K = np.array([[1.0, 0, 1], [0, 1, 1], [1, 0, 1]])
    Kp = np.linalg.pinv(K)
    r = int(np.linalg.matrix_rank(K))
    return {
        "K": K.tolist(),
        "rank": r, "kernel_dim": 3 - r,
        "KKpK_eq_K": bool(np.allclose(K @ Kp @ K, K)),
        "quotient_dim": r,
        "pseudoinverse": np.round(Kp, 6).tolist(),
    }


def main():
    data = {"C13": c13_matroid(), "C14": c14_poset(), "C24": c24_layer()}
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
