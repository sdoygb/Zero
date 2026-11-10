#!/usr/bin/env python3
"""Zb082 -- three generations from the 3 perfect matchings / 3 octahedron axes.

(1) Explicit map: matchings M1,M2,M3 -> octahedron axes -> generations n=1,2,3;
    the K4 3-cycle (1 2 3) cyclically permutes the matchings.
(2) Z3-circulant mass matrix: sqrt(m_n) = A[1 + b cos(theta + 2 pi (n-1)/3)];
    fit the charged leptons; Q = (1 + b^2/2)/3, b = sqrt(2) => Q = 2/3.
(3) Mixing: the Z3 limit gives the Fourier matrix (sin^2 th13 = 1/3,
    sin^2 th12 = sin^2 th23 = 1/2); the complement path
    M(t) = (1-t) M_Z3 + t diag(m) gives a computable family.

Output: ZeroGPU/results/zb082_three_generations.json
"""
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "ZeroGPU", "results", "zb082_three_generations.json")

MATCHINGS = {
    "M1": [(0, 1), (2, 3)],
    "M2": [(0, 2), (1, 3)],
    "M3": [(0, 3), (1, 2)],
}
GENERATIONS = {"M1": 1, "M2": 2, "M3": 3}


def cycle_vertex(v, cyc=(1, 2, 3)):
    # 3-cycle 1->2->3->1, 0 fixed
    if v == 0:
        return 0
    return cyc[(cyc.index(v) + 1) % 3]


def cycle_matching(edges):
    return tuple(tuple(sorted((cycle_vertex(a), cycle_vertex(b)))) for a, b in edges)


def z3_check():
    out = {}
    for name, edges in MATCHINGS.items():
        img = cycle_matching(edges)
        for other, other_edges in MATCHINGS.items():
            if set(img) == set(tuple(sorted(e)) for e in other_edges):
                out[name] = other
    return out


def fit_leptons():
    m = np.array([0.51099895, 105.6583755, 1776.86])  # MeV
    y = np.sqrt(m)
    A = float(y.mean())
    w = np.exp(-2j * np.pi / 3)
    z = (2 / 3) * sum(y[n] * w ** (-n) for n in range(3))
    C = abs(z)
    b = C / A
    theta = float(np.angle(z))
    Q = float(m.sum() / y.sum() ** 2)
    pred = (A * (1 + b * np.cos(theta + 2 * np.pi * np.arange(3) / 3))) ** 2
    return {"m": m.tolist(), "A": A, "C": C, "b": b, "theta_deg": float(np.degrees(theta)),
            "Q": Q, "pred": pred.tolist(), "sqrt2": float(np.sqrt(2)), "two_thirds": 2 / 3}


def mixing_table():
    m = np.array([0.51099895, 105.6583755, 1776.86])
    y = np.sqrt(m)
    w = np.exp(-2j * np.pi / 3)
    F = np.array([[w ** (j * k) for k in range(3)] for j in range(3)]) / np.sqrt(3)
    S = F @ np.diag(y) @ F.conj().T
    M0 = S @ S
    D = np.diag(m)
    rows = []
    for t in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0]:
        M = (1 - t) * M0 + t * D
        ev, U = np.linalg.eigh(M)
        s13 = abs(U[0, 2]) ** 2
        s12 = abs(U[0, 1]) ** 2 / (1 - s13) if 1 - s13 > 1e-12 else 0.0
        s23 = abs(U[1, 2]) ** 2 / (1 - s13) if 1 - s13 > 1e-12 else 0.0
        rows.append({"t": t, "s13_sq": float(s13), "s12_sq": float(s12),
                     "s23_sq": float(s23), "m_eigs": ev.tolist()})
    return rows



OBS = np.array([0.307, 0.5, 0.022])  # approximate observed s12^2, s23^2, s13^2


def _angles(M):
    ev, V = np.linalg.eigh(M)
    s13 = abs(V[0, 2]) ** 2
    s12 = abs(V[0, 1]) ** 2 / (1 - s13) if 1 - s13 > 1e-14 else 0.0
    s23 = abs(V[1, 2]) ** 2 / (1 - s13) if 1 - s13 > 1e-14 else 0.0
    return np.array([s12, s23, s13]), ev


def tbm_analysis():
    """TBM basis + Z3-breaking perturbations B12, B13, B23."""
    from scipy.optimize import least_squares
    m = np.array([0.51099895, 105.6583755, 1776.86])
    U = np.array([[np.sqrt(2 / 3), np.sqrt(1 / 3), 0],
                  [-np.sqrt(1 / 6), np.sqrt(1 / 3), np.sqrt(1 / 2)],
                  [np.sqrt(1 / 6), -np.sqrt(1 / 3), np.sqrt(1 / 2)]])
    M0 = U @ np.diag(m) @ U.T
    B12 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0.0]])
    B13 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0.0]])
    B23 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0.0]])

    def fit(ops, x0):
        def resid(p):
            M = M0 + sum(x * m.max() * B for x, B in zip(p, ops))
            ang, ev = _angles(M)
            # Enforce PSD on the first (smallest) eigenvalue relative to m_e.
            pen = 1e3 * max(0.0, -ev[0]) / m[0]
            return np.concatenate([ang - OBS, [pen]])
        res = least_squares(resid, x0, bounds=([-0.5] * len(x0), [0.5] * len(x0)))
        M = M0 + sum(x * m.max() * B for x, B in zip(res.x, ops))
        ang, ev = _angles(M)
        return {"params": [float(x) for x in res.x], "angles": [float(x) for x in ang],
                "residual": [float(x) for x in ang - OBS], "masses": [float(x) for x in ev],
                "min_eigenvalue": float(ev[0]), "psd_valid": bool(ev[0] >= 0.0)}

    return {"obs": [float(x) for x in OBS], "tbm_zero": [1 / 3, 1 / 2, 0.0],
            "two_direction": fit([B12, B13], [0.03, 0.03]),
            "three_direction": fit([B12, B13, B23], [0.05, 0.1, 0.03])}



def separation_model():
    """Mass-mixing separation: masses fixed, mixing = R12(a) R13(b) U_TBM.
    Fit (a,b) to s12^2 and s13^2; s23^2 is the prediction."""
    from scipy.optimize import least_squares
    m = np.array([0.51099895, 105.6583755, 1776.86])
    U = np.array([[np.sqrt(2 / 3), np.sqrt(1 / 3), 0],
                  [-np.sqrt(1 / 6), np.sqrt(1 / 3), np.sqrt(1 / 2)],
                  [np.sqrt(1 / 6), -np.sqrt(1 / 3), np.sqrt(1 / 2)]], dtype=complex)

    def R12(t):
        c, s = np.cos(t), np.sin(t)
        return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]], dtype=complex)

    def R13(t):
        c, s = np.cos(t), np.sin(t)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=complex)

    def angles(M):
        P = np.abs(M) ** 2
        s13 = P[0, 2]
        s12 = P[0, 1] / (1 - s13) if 1 - s13 > 1e-14 else 0.0
        s23 = P[1, 2] / (1 - s13) if 1 - s13 > 1e-14 else 0.0
        return np.array([s12, s23, s13])

    obs = np.array([0.307, 0.5, 0.022])

    def resid(p):
        ang = angles(R12(p[0]) @ R13(p[1]) @ U)
        return np.array([ang[0] - obs[0], ang[2] - obs[2]])

    res = least_squares(resid, [0.05, 0.1], bounds=([-1.5, -1.5], [1.5, 1.5]))
    Uf = R12(res.x[0]) @ R13(res.x[1]) @ U
    ang = angles(Uf)
    M = Uf @ np.diag(m) @ Uf.conj().T
    return {"params_rad": [float(x) for x in res.x], "angles": [float(x) for x in ang],
            "obs": [float(x) for x in obs], "s23_residual": float(ang[1] - obs[1]),
            "mass_eigs": [float(x) for x in np.linalg.eigvalsh(M)]}



def cp_phase_test():
    """Insert diag(1, e^{i delta}, 1) between R12 and R13; check angles and J."""
    a, b = 0.08524777813883505, 0.12539983068284222
    U = np.array([[np.sqrt(2 / 3), np.sqrt(1 / 3), 0],
                  [-np.sqrt(1 / 6), np.sqrt(1 / 3), np.sqrt(1 / 2)],
                  [np.sqrt(1 / 6), -np.sqrt(1 / 3), np.sqrt(1 / 2)]], dtype=complex)

    def R12(t):
        c, s = np.cos(t), np.sin(t)
        return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]], dtype=complex)

    def R13(t):
        c, s = np.cos(t), np.sin(t)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=complex)

    def angles(M):
        P = np.abs(M) ** 2
        s13 = P[0, 2]
        s12 = P[0, 1] / (1 - s13) if 1 - s13 > 1e-14 else 0.0
        s23 = P[1, 2] / (1 - s13) if 1 - s13 > 1e-14 else 0.0
        return np.array([s12, s23, s13])

    rows = []
    for d in [0.0, np.pi / 4, np.pi / 2, np.pi]:
        Uf = R12(a) @ np.diag([1, np.exp(1j * d), 1]) @ R13(b) @ U
        ang = angles(Uf)
        J = np.imag(Uf[0, 0] * Uf[1, 1] * np.conj(Uf[0, 1]) * np.conj(Uf[1, 0]))
        rows.append({"delta": float(d), "angles": [float(x) for x in ang], "J": float(J)})

    from scipy.optimize import least_squares
    obs = np.array([0.307, 0.5, 0.022])

    def fun(p):
        return R12(p[0]) @ np.diag([1, np.exp(1j * p[2]), 1]) @ R13(p[1]) @ U

    best = None
    for a0 in [-0.3, 0.0, 0.085, 0.3]:
        for b0 in [-0.3, 0.0, 0.125, 0.3]:
            for d0 in [0.0, 0.5, 1.0, 2.0, 3.0]:
                res = least_squares(lambda p: angles(fun(p)) - obs, [a0, b0, d0],
                                    bounds=([-1.5, -1.5, -np.pi], [1.5, 1.5, np.pi]))
                err = float(np.sum((angles(fun(res.x)) - obs) ** 2))
                if best is None or err < best[0]:
                    Uf = fun(res.x)
                    J = float(np.imag(Uf[0, 0] * Uf[1, 1] * np.conj(Uf[0, 1]) * np.conj(Uf[1, 0])))
                    best = (err, [float(x) for x in res.x],
                            [float(x) for x in angles(Uf)], J)
    def branch(a, b, d):
        Uf = R12(a) @ np.diag([1, np.exp(1j * d), 1]) @ R13(b) @ U
        J = float(np.imag(Uf[0, 0] * Uf[1, 1] * np.conj(Uf[0, 1]) * np.conj(Uf[1, 0])))
        return {"params": [float(a), float(b), float(d)],
                "angles": [float(x) for x in angles(Uf)], "J": J,
                "angle_cost": float(np.sum((angles(Uf) - obs) ** 2))}

    # Two rephasing-equivalent representatives both have J=0.
    zero_j_branches = [branch(0.083488628, 0.12407255, 0.0),
                       branch(-0.083488628, 0.12407255, np.pi)]
    nonzero_j_near = branch(-0.14833067, -0.06185573, 0.8)
    return {"scan": rows,
            "fit": {"err": best[0], "params": best[1], "angles": best[2], "J": best[3]},
            "zero_J_branches": zero_j_branches,
            "nonzero_J_near_solution": nonzero_j_near,
            "J_not_fixed_by_angles_alone": True}


def main():
    data = {
        "mapping": {"matchings": MATCHINGS, "generations": GENERATIONS,
                    "z3_action": z3_check()},
        "lepton_fit": fit_leptons(),
        "mixing": mixing_table(),
        "z3_limit_angles": {"sin2_th13": 1 / 3, "sin2_th12": 1 / 2, "sin2_th23": 1 / 2},
        "tbm": tbm_analysis(),
        "separation": separation_model(),
        "cp_phase": cp_phase_test(),
    }
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
