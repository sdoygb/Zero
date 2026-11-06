#!/usr/bin/env python3
"""Zb075 -- structure conjugation: complement spectral transport and defect zero modes.

(1) Graph complement C(A) = J - I - A on K4 and on the octahedron K_{2,2,2}.
    Check: C is an involution; for a connected k-regular n-vertex graph,
      spec(C(A)) = {n-1-k} u {-1-lam : lam in spec(A), lam != k}.
(2) 1D Dirac/BdG H = -i v tau_z d_x + m(x) tau_x.
    Sign-changing m: normalizable zero mode, width ~ v/|m|.
    Uniform m: spectral gap 2|m|, no zero mode.
(3) s-wave SNS lattice: min|E| at phase difference pi vs delta0 (finite-size trend
    toward the Andreev zero E(pi)=0 of a short junction).

Outputs
  ZeroGPU/results/zb075_conj_defect.json
  figures/zb075_zero_mode.png
"""
import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# ---------------------------------------------------------------- graphs
def octahedron():
    """K_{2,2,2}: 6 vertices, 3 pairs; edges join different pairs."""
    n = 6
    grp = [0, 0, 1, 1, 2, 2]
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j and grp[i] != grp[j]:
                A[i, j] = 1.0
    return A


def complement(A):
    n = A.shape[0]
    return np.ones((n, n)) - np.eye(n) - A


def is_connected(A):
    n = A.shape[0]
    seen = {0}
    stack = [0]
    while stack:
        i = stack.pop()
        for j in range(n):
            if A[i, j] != 0 and j not in seen:
                seen.add(j)
                stack.append(j)
    return len(seen) == n


def graph_block():
    out = {}
    for name, A in [("K4", np.ones((4, 4)) - np.eye(4)), ("octahedron", octahedron())]:
        n = A.shape[0]
        k = float(A.sum(axis=1)[0])
        lam = np.linalg.eigvalsh(A)
        comp = complement(A)
        lam_c = np.linalg.eigvalsh(comp)
        pred = np.sort(np.concatenate([[n - 1 - k],
                                       [-1.0 - x for x in lam if abs(x - k) > 1e-9]]))
        out[name] = {
            "n": n,
            "k": k,
            "connected": bool(is_connected(A)),
            "spec_A": [round(float(x), 9) for x in lam],
            "spec_complement": [round(float(x), 9) for x in lam_c],
            "spec_formula": [round(float(x), 9) for x in pred],
            "formula_match": bool(np.allclose(pred, lam_c)),
            "involution": bool(np.allclose(complement(comp), A)),
        }
    return out


# ---------------------------------------------------------------- Dirac
def dirac_1d(N, m, v=1.0, a=1.0):
    tz = np.array([[1, 0], [0, -1]], dtype=complex)
    tx = np.array([[0, 1], [1, 0]], dtype=complex)
    H = np.zeros((2 * N, 2 * N), dtype=complex)
    for j in range(N):
        H[2 * j:2 * j + 2, 2 * j:2 * j + 2] = m[j] * tx
    for j in range(N - 1):
        blk = -1j * v / (2 * a) * tz
        H[2 * j:2 * j + 2, 2 * (j + 1):2 * (j + 1) + 2] = blk
        H[2 * (j + 1):2 * (j + 1) + 2, 2 * j:2 * j + 2] = blk.conj().T
    return H


def dirac_block(N=400, v=1.0):
    x = np.arange(N)
    out = {}
    profiles = {}
    for label, m in [
        ("uniform", 0.5 * np.ones(N)),
        ("wall", 0.5 * np.tanh((x - N / 2) / 2.0)),
        ("wall_weak", 0.2 * np.tanh((x - N / 2) / 2.0)),
    ]:
        H = dirac_1d(N, m, v=v)
        w, vec = np.linalg.eigh(H)
        idx = np.argsort(np.abs(w))[:2]
        k = idx[0]
        psi = vec[:, k]
        prob = np.abs(psi[0::2]) ** 2 + np.abs(psi[1::2]) ** 2
        prob = prob / prob.sum()
        pr = 1.0 / float(np.sum(prob ** 2))
        out[label] = {
            "min_abs_E": float(np.abs(w).min()),
            "lowest": [float(w[i]) for i in idx],
            "participation_ratio": pr,
            "peak_site": int(np.argmax(prob)),
        }
        profiles[label] = prob
    return out, profiles, x


# ---------------------------------------------------------------- SNS
def bdg_s(N, delta, phi, t=1.0, mu=0.0):
    h = np.zeros((N, N), dtype=complex)
    for j in range(N - 1):
        h[j, j + 1] = h[j + 1, j] = -t
    h += -mu * np.eye(N)
    D = np.diag(np.asarray(delta, dtype=complex) * np.exp(1j * np.asarray(phi)))
    return np.block([[h, D], [D.conj(), -h.T]])


def sns_block(N=400, Ln=10):
    x = np.arange(N)
    j1 = N // 2 - Ln // 2
    j2 = j1 + Ln
    out = {}
    for delta0 in [0.5, 0.2, 0.1]:
        delta = np.zeros(N)
        delta[:j1] = delta0
        delta[j2:] = delta0
        row = {}
        for phi in [0.0, math.pi / 2, math.pi]:
            ph = np.zeros(N)
            ph[j2:] = phi
            H = bdg_s(N, delta, ph)
            w = np.linalg.eigvalsh(H)
            row[f"phi_{phi:.2f}"] = float(np.abs(w).min())
        out[f"delta0_{delta0}"] = row
    return out


def main():
    graph = graph_block()
    dirac, profiles, x = dirac_block()
    sns = sns_block()
    data = {"graph": graph, "dirac": dirac, "sns": sns}
    res_path = os.path.join(ROOT, "ZeroGPU", "results", "zb075_conj_defect.json")
    os.makedirs(os.path.dirname(res_path), exist_ok=True)
    with open(res_path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(json.dumps(data, ensure_ascii=False, indent=2))

    # figure: zero-mode profiles + SNS trend
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), dpi=160)
        ax = axes[0]
        for label, style in [("uniform", "-"), ("wall", "--"), ("wall_weak", ":")]:
            ax.plot(x, profiles[label], style, lw=1.2, label=label)
        ax.set_xlabel("site")
        ax.set_ylabel(r"$|\psi(x)|^2$")
        ax.set_title("1D Dirac: sign-changing mass -> zero mode")
        ax.legend(fontsize=7)
        ax = axes[1]
        picks = [0.5, 0.2, 0.1]
        vals = [sns[f"delta0_{d}"]["phi_3.14"] for d in picks]
        ax.loglog(picks, vals, "o-")
        ax.set_xlabel(r"$\Delta_0$")
        ax.set_ylabel(r"$\min|E|$ at $\varphi=\pi$")
        ax.set_title("SNS lattice: trend toward $E(\\pi)=0$")
        fig.tight_layout()
        fig_path = os.path.join(ROOT, "figures", "zb075_zero_mode.png")
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        fig.savefig(fig_path)
        print(f"saved {fig_path}")
    except Exception as e:  # noqa: BLE001
        print(f"figure skipped: {e}")


if __name__ == "__main__":
    main()
