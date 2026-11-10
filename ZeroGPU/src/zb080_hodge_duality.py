#!/usr/bin/env python3
"""Zb080 -- 4d Hodge duality on 2-forms: star, self-dual decomposition, EM duality.

Checks:
  (1) Euclidean star on the 6 basis 2-forms (12,13,14,23,24,34):
      *(12)=34, *(13)=-24, *(14)=23, *(23)=14, *(24)=-13, *(34)=12; star^2 = 1.
  (2) eigenvalues: +/-1, each multiplicity 3; projectors P+/- = (1 +/- star)/2.
  (3) F ^ F = |F+|^2 - |F-|^2 for random F (with F^F = 2(f12 f34 - f13 f24 + f14 f23) vol).
  (4) Lorentzian model: a real 6x6 complex structure J with J^2 = -1; eigenvalues +/- i.
  (5) 4-cube direction-pair subspace: 01<->23, 02<->13, 03<->12; star^2 = 1.
      The full 4-cube has 24 square 2-faces, not 6.

Output: ZeroGPU/results/zb080_hodge_duality.json
"""
import json
import itertools
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "ZeroGPU", "results", "zb080_hodge_duality.json")

BASIS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
IDX = {b: i for i, b in enumerate(BASIS)}


def euclidean_star():
    S = np.zeros((6, 6), dtype=int)
    pairs = {((1, 2), (3, 4), 1), ((1, 3), (2, 4), -1), ((1, 4), (2, 3), 1),
             ((2, 3), (1, 4), 1), ((2, 4), (1, 3), -1), ((3, 4), (1, 2), 1)}
    for a, b, s in pairs:
        S[IDX[a], IDX[b]] = s
    return S


def main():
    S = euclidean_star()
    P_plus = (np.eye(6) + S) / 2
    P_minus = (np.eye(6) - S) / 2
    ev = np.linalg.eigvalsh(S)
    # F ^ F = 2(f12 f34 - f13 f24 + f14 f23)
    W = np.zeros((6, 6))
    for a, b, s in [((1, 2), (3, 4), 1), ((1, 3), (2, 4), -1), ((1, 4), (2, 3), 1)]:
        W[IDX[a], IDX[b]] = s
        W[IDX[b], IDX[a]] = s
    rng = np.random.default_rng(23)
    f = rng.normal(size=6)
    wedge = 2 * (f[0] * f[5] + f[1] * f[4] * (-1) + f[2] * f[3])
    norm2 = float(f @ f)
    star_f = S @ f
    # |F+|^2 = (1/2)(|F|^2 + F.SF), |F-|^2 = (1/2)(|F|^2 - F.SF)
    fp2 = 0.5 * (norm2 + float(f @ star_f))
    fm2 = 0.5 * (norm2 - float(f @ star_f))
    # Lorentzian complex structure
    J2 = np.zeros((6, 6))
    for i in range(3):
        J2[2 * i, 2 * i + 1] = -1
        J2[2 * i + 1, 2 * i] = 1
    # 4-cube direction-pair classes: (01<->23, 02<->13, 03<->12).
    cube_pairs = {((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))}
    all_square_faces = {
        (mu, nu, fixed)
        for mu, nu in itertools.combinations(range(4), 2)
        for fixed in itertools.product((0, 1), repeat=2)
    }
    data = {
        "star": S.tolist(),
        "star_sq_is_identity": bool(np.array_equal(S @ S, np.eye(6, dtype=int))),
        "eigenvalues": np.round(ev, 9).tolist(),
        "projectors_orthogonal": bool(np.allclose(P_plus @ P_minus, 0)),
        "projectors_complete": bool(np.allclose(P_plus + P_minus, np.eye(6))),
        "projectors_idempotent": bool(np.allclose(P_plus @ P_plus, P_plus)
                                      and np.allclose(P_minus @ P_minus, P_minus)),
        "wedge_example": {"F": f.tolist(), "F_wedge_F": float(wedge),
                          "Fplus2_minus_Fminus2": float(fp2 - fm2),
                          "identity": bool(abs(wedge - (fp2 - fm2)) < 1e-9)},
        "complex_structure": {"J_sq_is_minus_identity": bool(np.allclose(J2 @ J2, -np.eye(6))),
                              "eigenvalues_re_im": [[round(float(z.real), 9), round(float(z.imag), 9)] for z in np.linalg.eigvals(J2)]},
        "cube_direction_pair_classes": [[list(a), list(b)] for a, b in cube_pairs],
        "cube_direction_pair_count": len(cube_pairs),
        "cube_square_face_count": len(all_square_faces),
    }
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=2)
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
