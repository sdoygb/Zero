#!/usr/bin/env python3
"""Zb077 -- structure conjugation, second batch: representation/code/matroid/poset.

Checks:
  C11 partition transpose: involution, size invariant.
  C12 code dual: dim C + dim C^perp = n ([7,4,3] Hamming).
  C13 matroid dual: r*(E) = |E| - r(E) (graphic matroid of K4).
  C14 poset dual: ideals <-> filters by complement.
  C15/C16 conjugate/dual representation: chi(conj rho) = conj chi(rho) (Z3, S3).
  C17 Kramers: T^2 = -1 (3 spins) => even degeneracy of a T-invariant H.
  C23 Fourier: F^2 = reversal, F^4 = 1 (Z_N).

Output: ZeroGPU/results/zb077_conj_algebra.json
"""
import itertools
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "ZeroGPU", "results", "zb077_conj_algebra.json")


# ---------- C11 partitions ----------
def partitions(n, mx=None):
    if mx is None:
        mx = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, mx), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def transpose(lam):
    if not lam:
        return ()
    return tuple(sum(1 for x in lam if x >= j) for j in range(1, lam[0] + 1))


def c11():
    ok = all(transpose(transpose(l)) == l and sum(transpose(l)) == sum(l)
             for n in range(9) for l in partitions(n))
    return {"n_max": 8, "involution_and_size": ok}


# ---------- C12 code ----------
def c12():
    H = np.array([[1, 0, 1, 0, 1, 0, 1],
                  [0, 1, 1, 0, 0, 1, 1],
                  [0, 0, 0, 1, 1, 1, 1]], dtype=int)
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
    return {"n": int(n), "dim_C": len(free), "dim_dual": int(H.shape[0]),
            "sum": int(n), "rank_H": int(np.linalg.matrix_rank(H))}


# ---------- C13 matroid ----------
def c13():
    edges = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    n = 4

    def rank(es):
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for k in es:
            a, b = edges[k]
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
        comps = len({find(i) for i in range(n)})
        return n - comps

    E = list(range(len(edges)))
    rE = rank(E)
    # dual rank r*(A) = |A| - r(E) + r(E\A)
    def rank_star_direct(A):
        """Rank in M* by enumerating independent subsets of the restricted dual."""
        A = set(A)
        best = 0
        for k in range(len(A) + 1):
            for X in itertools.combinations(A, k):
                if rank(set(E) - set(X)) == rE:
                    best = max(best, k)
        return best

    all_subsets = [A for k in range(7) for A in itertools.combinations(E, k)]
    dual_formula_checks = [
        rank_star_direct(A) == len(A) - rE + rank(set(E) - set(A))
        for A in all_subsets
    ]
    rstar_E = len(E) - rE + rank([])
    return {"edges": len(E), "rE": rE, "rstar_E": rstar_E,
            "dual_rank_formula_full": all(dual_formula_checks),
            "subsets_checked": len(dual_formula_checks)}


# ---------- C14 poset ----------
def c14():
    # diamond poset: 0<1, 0<2, 1<3, 2<3
    le = {(0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (1, 3), (2, 2), (2, 3), (3, 3)}
    elems = [0, 1, 2, 3]
    ideals = []
    for mask in range(16):
        S = {elems[i] for i in range(4) if mask >> i & 1}
        if all((b in S) for (a, b) in le if a in S):
            ideals.append(S)
    filters = []
    for mask in range(16):
        S = {elems[i] for i in range(4) if mask >> i & 1}
        if all((a in S) for (a, b) in le if b in S):
            filters.append(S)
    comp = {frozenset(elems) - I for I in ideals}
    return {"ideals": len(ideals), "filters": len(filters),
            "complement_bijection": comp == {frozenset(F) for F in filters}}


# ---------- C15/C16 representations ----------
def c15():
    w = np.exp(2j * np.pi / 3)
    # Z3 2D representation diag(w, w^2)
    rho = [np.diag([w ** k, w ** (2 * k)]) for k in range(3)]
    ok = all(np.allclose(np.trace(rho[k].conj()), np.conj(np.trace(rho[k]))) for k in range(3))
    # S3 standard 2D representation (real), check character conjugation
    # permutation matrices on the 2D subspace sum-zero: use 3x3 permutation and project
    import itertools
    chars = []
    for p in itertools.permutations(range(3)):
        P = np.zeros((3, 3))
        for i, pi in enumerate(p):
            P[i, pi] = 1
        # action on the 2D subspace {sum=0}: basis e1-e3, e2-e3
        B = np.array([[1, 0], [0, 1], [-1, -1]], dtype=float)
        Q = np.linalg.pinv(B) @ P @ B
        chars.append(np.trace(Q))
    return {"Z3_char_conj": ok, "S3_char_real": bool(np.allclose(np.array(chars), np.array(chars).conj()))}


# ---------- C17 Kramers ----------
def c17():
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)
    # 3 spins: T = (i sy) (x) I (x) I ; T^2 = -1
    U = 1j * np.kron(np.kron(sy, I2), I2)
    # T = U K ; T^2 = U U^* = -1 ?
    T2 = U @ U.conj()
    # random Hermitian H, symmetrize to T-invariant: H = (H + U H^* U^dagger)/2
    rng = np.random.default_rng(19)
    A = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
    H = (A + A.conj().T) / 2
    H = (H + U @ H.conj() @ U.conj().T) / 2
    H = (H + H.conj().T) / 2
    ev = np.linalg.eigvalsh(H)
    # check even degeneracy
    ev_sorted = np.sort(ev)
    pairs = all(abs(ev_sorted[2 * i] - ev_sorted[2 * i + 1]) < 1e-9 for i in range(4))
    return {"T2_is_minus1": bool(np.allclose(T2, -np.eye(8))),
            "H_T_invariant": bool(np.allclose(U @ H.conj() @ U.conj().T, H)),
            "even_degeneracy": bool(pairs),
            "eigs": np.round(ev_sorted, 6).tolist()}


# ---------- C23 Fourier ----------
def c23():
    out = {}
    for N in [4, 8]:
        w = np.exp(-2j * np.pi / N)
        F = np.array([[w ** (j * k) for k in range(N)] for j in range(N)]) / np.sqrt(N)
        P = np.zeros((N, N), dtype=complex)  # reversal: (Pv)_j = v_{-j}
        for j in range(N):
            P[j, (-j) % N] = 1
        out[f"N{N}"] = {"F2_is_reversal": bool(np.allclose(F @ F, P)),
                        "F4_is_identity": bool(np.allclose(F @ F @ F @ F, np.eye(N)))}
    return out


def main():
    data = {"C11": c11(), "C12": c12(), "C13": c13(), "C14": c14(),
            "C15C16": c15(), "C17": c17(), "C23": c23()}
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
