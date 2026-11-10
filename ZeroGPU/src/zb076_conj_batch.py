#!/usr/bin/env python3
"""Zb076 -- structure conjugation, first theorem batch.

Matrix batch: transpose/Hodge, line graph, bipartite sign, Seidel switching, Schur.
Word batch: reversal and negation.
Constraint/group batch: constraint transpose, stabilizer conjugation.

Output: ZeroGPU/results/zb076_conj_batch.json
"""
import itertools
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "ZeroGPU", "results", "zb076_conj_batch.json")


# ---------- C01 补图 ----------
def c01_complement():
    n = 6
    grp = [0, 0, 1, 1, 2, 2]
    A = np.array([[1.0 if (i != j and grp[i] != grp[j]) else 0.0
                   for j in range(n)] for i in range(n)])
    J, I = np.ones((n, n)), np.eye(n)
    C = J - I - A
    return {
        "spec_A": np.round(np.linalg.eigvalsh(A), 9).tolist(),
        "spec_C": np.round(np.linalg.eigvalsh(C), 9).tolist(),
        "involution": bool(np.allclose(J - I - C, A)),
    }


# ---------- C02 转置 / Hodge ----------
def c02_transpose():
    rng = np.random.default_rng(7)
    B = rng.integers(-1, 2, size=(5, 9)).astype(float)
    sv = np.linalg.svd(B, compute_uv=False)
    svT = np.linalg.svd(B.T, compute_uv=False)
    return {
        "rank": int(np.linalg.matrix_rank(B)),
        "rank_T": int(np.linalg.matrix_rank(B.T)),
        "sv_match": bool(np.allclose(np.sort(sv), np.sort(svT))),
        "ker_dim_B": 9 - int(np.linalg.matrix_rank(B)),
        "ker_dim_BT": 5 - int(np.linalg.matrix_rank(B)),
    }


# ---------- C03 线图 ----------
def c03_line_graph():
    # 随机简单图 n=6
    rng = np.random.default_rng(3)
    A = np.zeros((6, 6))
    edges = [(i, j) for i in range(6) for j in range(i + 1, 6) if rng.random() < 0.5]
    for i, j in edges:
        A[i, j] = A[j, i] = 1.0
    E = len(edges)
    B = np.zeros((6, E))
    for k, (i, j) in enumerate(edges):
        B[i, k] = B[j, k] = 1.0
    A_L = B.T @ B - 2 * np.eye(E)
    Q = B @ B.T
    spec_AL = np.linalg.eigvalsh(A_L)
    spec_Q = np.linalg.eigvalsh(Q)
    r = int(np.linalg.matrix_rank(B))
    pred = np.sort(np.concatenate([spec_Q[spec_Q > 1e-9] - 2,
                                   -2 * np.ones(E - r)]))
    return {
        "n": 6, "edges": E, "rank_B": r,
        "spec_line": np.round(spec_AL, 9).tolist(),
        "spec_pred": np.round(pred, 9).tolist(),
        "match": bool(np.allclose(np.sort(spec_AL), pred)),
    }


# ---------- C04 二部图 Laplacian = 符号 Laplacian ----------
def c04_bipartite():
    rng = np.random.default_rng(11)
    nx, ny = 4, 4
    A = np.zeros((nx + ny, nx + ny))
    for i in range(nx):
        for j in range(ny):
            if rng.random() < 0.5:
                A[i, nx + j] = A[nx + j, i] = 1.0
    d = A.sum(1)
    D = np.diag(d)
    L = D - A
    Q = D + A
    sig = np.array([1.0] * nx + [-1.0] * ny)
    S = np.diag(sig)
    return {
        "DAD_eq_minusA": bool(np.allclose(S @ A @ S, -A)),
        "spec_L": np.round(np.linalg.eigvalsh(L), 9).tolist(),
        "spec_Q": np.round(np.linalg.eigvalsh(Q), 9).tolist(),
        "spec_match": bool(np.allclose(np.linalg.eigvalsh(L), np.linalg.eigvalsh(Q))),
    }


# ---------- C05 Seidel switching ----------
def c05_seidel():
    rng = np.random.default_rng(5)
    n = 7
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                A[i, j] = A[j, i] = 1.0
    S = np.ones((n, n)) - np.eye(n) - 2 * A
    sig = rng.choice([-1.0, 1.0], size=n)
    D = np.diag(sig)
    A_sw = np.where(np.outer(sig, sig) > 0, A, 1.0 - A)
    np.fill_diagonal(A_sw, 0.0)
    S_sw = np.ones((n, n)) - np.eye(n) - 2 * A_sw
    return {
        "seidel_equivariance": bool(np.allclose(S_sw, D @ S @ D)),
        "spec_match": bool(np.allclose(np.sort(np.linalg.eigvalsh(S_sw)),
                                       np.sort(np.linalg.eigvalsh(S)))),
    }


# ---------- C06/C07 词反转与取反 ----------
def height(w):
    h, out = 0, [0]
    for x in w:
        h += x
        out.append(h)
    return out


def area(w):
    return sum(height(w)[:-1])


def runs(w):
    return 1 + sum(1 for a, b in zip(w, w[1:]) if a != b) if w else 0


def sector(w):
    # 环上 +1 段的个数
    if all(x == w[0] for x in w):
        return 1 if w[0] == 1 else 0
    return sum(1 for i in range(len(w)) if w[i] == 1 and w[i - 1] == -1)


def c0607_words(L=6):
    words = [w for w in itertools.product([-1, 1], repeat=L) if sum(w) == 0]
    rev_ok = neg_ok = True
    for w in words:
        wr = tuple(reversed(w))
        wn = tuple(-x for x in w)
        if area(wr) != -area(w) or runs(wr) != runs(w) or sector(wr) != sector(w):
            rev_ok = False
        if area(wn) != -area(w) or runs(wn) != runs(w) or sector(wn) != sector(w):
            neg_ok = False
        if sum(wn) != 0 or sum(1 for x in wn if x == 1) != L - sum(1 for x in w if x == 1):
            neg_ok = False
    return {
        "L": L, "balanced_words": len(words),
        "reversal_area_odd": rev_ok,
        "negation_area_odd": neg_ok,
        "example": {"w": list(words[0]), "area": area(words[0]),
                    "area_rev": area(tuple(reversed(words[0]))),
                    "area_neg": area(tuple(-x for x in words[0])),
                    "runs": runs(words[0]), "sector": sector(words[0])},
    }


# ---------- C08 约束转置 ----------
def c08_constraint():
    rng = np.random.default_rng(13)
    B = rng.integers(-1, 2, size=(7, 11)).astype(float)
    r = int(np.linalg.matrix_rank(B))
    return {
        "rank": r, "rank_T": int(np.linalg.matrix_rank(B.T)),
        "dof_var": 11 - r, "dof_eq": 7 - r,
        "rank_match": bool(r == np.linalg.matrix_rank(B.T)),
    }


# ---------- C09 Schur 补 ----------
def c09_schur():
    rng = np.random.default_rng(17)
    A = rng.normal(size=(3, 3))
    B = rng.normal(size=(3, 2))
    C = rng.normal(size=(2, 3))
    Dm = rng.normal(size=(2, 2))
    A = A @ A.T + np.eye(3)
    M = np.block([[A, B], [C, Dm]])
    S = Dm - C @ np.linalg.solve(A, B)
    return {
        "det_M": float(np.linalg.det(M)),
        "det_A_det_S": float(np.linalg.det(A) * np.linalg.det(S)),
        "match": bool(np.allclose(np.linalg.det(M), np.linalg.det(A) * np.linalg.det(S))),
    }


# ---------- C10 子群共扼 ----------
def c10_subgroup():
    from itertools import permutations
    def compose(p, q):  # p o q
        return tuple(p[q[i]] for i in range(3))
    e = (0, 1, 2)
    g = (1, 0, 2)          # (0 1): 0 -> 1
    h = (0, 2, 1)          # (1 2): fixes 0
    H = {e, h}             # stabilizer of 0
    ginv = tuple(np.argsort(g))
    gH = {compose(compose(g, x), ginv) for x in H}
    stab0 = {p for p in permutations(range(3)) if p[0] == 0}
    stab1 = {p for p in permutations(range(3)) if p[1] == 1}
    return {
        "stab0": sorted(stab0), "stab1": sorted(stab1),
        "gH": sorted(gH),
        "conj_match": bool(gH == stab1),
    }


# ---------- C11 划分转置 ----------
def c11_partition():
    def parts(n, mx=None):
        if mx is None:
            mx = n
        if n == 0:
            yield ()
        else:
            for k in range(min(n, mx), 0, -1):
                for rest in parts(n - k, k):
                    yield (k,) + rest
    def transpose(lam):
        if not lam:
            return ()
        return tuple(sum(1 for x in lam if x >= j) for j in range(1, lam[0] + 1))
    ok = all(transpose(transpose(l)) == l and sum(transpose(l)) == sum(l) for l in parts(6))
    return {"partitions_of_6": len(list(parts(6))), "involution_size": ok}


# ---------- C12 码对偶 ----------
def c12_code():
    # [7,4,3] Hamming
    H = np.array([[1,0,1,0,1,0,1],
                  [0,1,1,0,0,1,1],
                  [0,0,0,1,1,1,1]], dtype=int)
    # generator: null space of H over GF(2), via RREF
    M = H.copy() % 2
    m, n = M.shape
    pivots, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, m) if M[i, c]), None)
        if p is None:
            continue
        M[[r, p]] = M[[p, r]]
        for i in range(m):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        pivots.append(c)
        r += 1
        if r == m:
            break
    free = [c for c in range(n) if c not in pivots]
    G = []
    for f in free:
        v = np.zeros(n, dtype=int)
        v[f] = 1
        for i, c in enumerate(pivots):
            v[c] = M[i, f]
        G.append(v)
    G = np.array(G, dtype=int)
    prod = (G @ H.T) % 2
    return {"dim_C": int(G.shape[0]), "dim_dual": int(H.shape[0]), "sum": int(n),
            "GHT_zero": bool(np.all(prod == 0)),
            "rank_G": int(np.linalg.matrix_rank(G))}


def main():
    data = {
        "C01_complement": c01_complement(),
        "C02_transpose": c02_transpose(),
        "C03_line_graph": c03_line_graph(),
        "C04_bipartite": c04_bipartite(),
        "C05_seidel": c05_seidel(),
        "C06C07_words": c0607_words(),
        "C08_constraint": c08_constraint(),
        "C09_schur": c09_schur(),
        "C10_subgroup": c10_subgroup(),
        "C11_partition": c11_partition(),
        "C12_code": c12_code(),
    }
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
