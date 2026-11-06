#!/usr/bin/env python3
"""Zb079 -- planar Ising: Kramers-Wannier duality and the self-dual point.

Checks (open LxL square lattice, dual graph built from faces):
  (1) high-T expansion: Z = 2^V (cosh K)^E * sum_{A even} (tanh K)^|A|
  (2) low-T expansion:  Z = 2 e^{K E} * sum_{B cut} e^{-2K|B|}
  (3) exact finite planar KW for connected G:
      Z_G(K) = 2^{1-|V*|+E/2} (sinh 2K)^{E/2} Z_G*(K*),
      with sinh(2K) sinh(2K*) = 1
  (4) L=2 cut/cycle size histograms agree
  (5) self-dual point Kc = (1/2) ln(1+sqrt 2), Tc = 1/Kc (J/kB)

Output: ZeroGPU/results/zb079_ising_duality.json
"""
import itertools
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "ZeroGPU", "results", "zb079_ising_duality.json")


def square_open(L):
    V = [(i, j) for i in range(L) for j in range(L)]
    idx = {v: k for k, v in enumerate(V)}
    E = []
    for i in range(L):
        for j in range(L):
            if i + 1 < L:
                E.append((idx[(i, j)], idx[(i + 1, j)]))
            if j + 1 < L:
                E.append((idx[(i, j)], idx[(i, j + 1)]))
    return V, E


def connected_graph(nv, edges):
    adj = [[] for _ in range(nv)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    seen = {0}
    stack = [0]
    while stack:
        a = stack.pop()
        for b in adj[a]:
            if b not in seen:
                seen.add(b)
                stack.append(b)
    return len(seen) == nv


def dual_edges(L):
    faces = [(a, b) for a in range(L - 1) for b in range(L - 1)] + ["out"]
    fidx = {f: k for k, f in enumerate(faces)}

    def face(a, b):
        return (a, b) if 0 <= a <= L - 2 and 0 <= b <= L - 2 else "out"

    DE = []
    for i in range(L - 1):
        for j in range(L):
            DE.append((fidx[face(i, j - 1)], fidx[face(i, j)]))
    for i in range(L):
        for j in range(L - 1):
            DE.append((fidx[face(i - 1, j)], fidx[face(i, j)]))
    return faces, DE


def Z_graph(nv, edges, K):
    z = 0.0
    for s in itertools.product([-1, 1], repeat=nv):
        e = sum(1 if s[a] == s[b] else -1 for a, b in edges)
        z += math.exp(K * e)
    return z


def even_subsets_edges(nv, edges):
    out = []
    for mask in range(1 << len(edges)):
        deg = [0] * nv
        for k, (a, b) in enumerate(edges):
            if mask >> k & 1:
                deg[a] += 1
                deg[b] += 1
        if all(d % 2 == 0 for d in deg):
            out.append(bin(mask).count("1"))
    return out


def even_subsets(L):
    V, E = square_open(L)
    out = []
    for mask in range(1 << len(E)):
        deg = [0] * len(V)
        for k, (a, b) in enumerate(E):
            if mask >> k & 1:
                deg[a] += 1
                deg[b] += 1
        if all(d % 2 == 0 for d in deg):
            out.append(bin(mask).count("1"))
    return out


def cuts(L):
    V, E = square_open(L)
    seen = set()
    for mask in range(1 << len(V)):
        B = 0
        for k, (a, b) in enumerate(E):
            if ((mask >> a) & 1) != ((mask >> b) & 1):
                B |= 1 << k
        seen.add(B)
    return [bin(b).count("1") for b in seen]


def hist(xs):
    out = {}
    for x in xs:
        out[x] = out.get(x, 0) + 1
    return dict(sorted(out.items()))


def check(L, K=0.3):
    V, E = square_open(L)
    faces, DE = dual_edges(L)
    nv, ne, nvs = len(V), len(E), len(faces)
    ZG = Z_graph(nv, E, K)
    t = math.exp(-2 * K)
    Ks = math.atanh(t)
    ZS = Z_graph(nvs, DE, Ks)
    H = sum(math.tanh(K) ** k for k in even_subsets(L))
    Ls = sum(math.exp(-2 * K * k) for k in cuts(L))
    high_T = 2 ** nv * math.cosh(K) ** ne * H
    low_T = 2 * math.exp(K * ne) * Ls
    pref = 2 ** (1 - nvs + ne / 2) * (math.sinh(2 * K)) ** (ne / 2)
    return {
        "L": L, "V": nv, "E": ne, "Vstar": nvs,
        "connected": connected_graph(nv, E),
        "Z": ZG, "high_T": high_T, "low_T": low_T,
        "high_T_match": abs(ZG - high_T) < 1e-9 * max(1.0, abs(ZG)),
        "low_T_match": abs(ZG - low_T) < 1e-9 * max(1.0, abs(ZG)),
        "KW_ratio": ZG / (pref * ZS),
        "sinh_product": math.sinh(2 * K) * math.sinh(2 * Ks),
        "cut_hist": hist(cuts(L)),
        "even_hist": hist(even_subsets(L)),
        "dual_cycle_hist": hist(even_subsets_edges(nvs, DE)),
    }


def main():
    data = {
        "L2": check(2),
        "L3": check(3),
        "Kc": 0.5 * math.log(1 + math.sqrt(2)),
        "Tc_over_J_over_kB": 1.0 / (0.5 * math.log(1 + math.sqrt(2))),
    }
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
