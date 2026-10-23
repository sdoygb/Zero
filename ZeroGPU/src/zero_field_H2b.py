#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_field_H2b.py — Hartree 修好（4πa² + 1s 锚）后重跑 H2。

Hartree 核：V_H = c * [4π a² (G*N) - (4π a² G0 - ONSITE a²) N]
锚：精确 1s 密度给 E_H=0.3125（c 由它定）。
核势仍用路线 A：离源 -1/r、在源 -ONSITE/a。
输出：results/zero_field_H2b.json
"""
from __future__ import annotations
import json, pathlib, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.fft import fftn, ifftn

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "zero_field_H2b.json"
ONSITE = 2.3800756243710857


def lap(n):
    k = 2 * np.pi * np.fft.fftfreq(n)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    return 2 * (3 - np.cos(KX) - np.cos(KY) - np.cos(KZ))


def lattice(n):
    idx = {}; coords = []
    for x in range(n):
        for y in range(n):
            for z in range(n):
                idx[(x, y, z)] = len(coords); coords.append((x, y, z))
    N = len(coords); rows = []; cols = []; vals = []
    for (x, y, z), i in idx.items():
        nb = 0
        for dx, dy, dz in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)):
            j = idx.get((x+dx, y+dy, z+dz))
            if j is not None:
                nb += 1; rows.append(i); cols.append(j); vals.append(-1.0)
        rows.append(i); cols.append(i); vals.append(float(nb))
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N)), idx, coords


def Vn(xyz, ps, a):
    V = np.zeros(len(xyz))
    for p in ps:
        r = np.linalg.norm(xyz - p, axis=1)
        on = r < 1e-9; r = np.where(on, a, r)
        V += np.where(on, -ONSITE / a, -1.0 / r)
    return V


def main():
    a = 0.25; L = 14.0; Ne = 2
    n = int(round(L / a)) + 1
    Lm, idx, coords = lattice(n)
    xyz = np.array(coords, float) * a
    lam = lap(n); inv = np.zeros_like(lam); m = lam > 1e-12; inv[m] = 1 / lam[m]
    G0 = np.real(ifftn(inv))[0, 0, 0]; vol = a ** 3
    c0 = xyz[idx[(n//2, n//2, n//2)]]

    def VH_raw(Nsite1):
        n3 = Nsite1.reshape(n, n, n)
        v = 4 * np.pi * a ** 2 * np.real(ifftn(inv * fftn(n3))) - (4 * np.pi * a ** 2 * G0 - ONSITE * a ** 2) * n3
        return v.reshape(-1)

    # 1s 锚：E_H[1s]=0.3125
    r = np.linalg.norm(xyz - c0, axis=1)
    rho1s = (1 / np.pi) * np.exp(-2 * r)
    N1s = rho1s * vol
    V1 = VH_raw(N1s)
    EH_raw = 0.5 * np.sum(N1s * V1)
    c = 0.3125 / EH_raw
    print(f"1s 锚: E_H_raw={EH_raw:.4f}, c={c:.4f}")

    def run(R):
        nuc = [c0 - np.array([0, 0, R/2]), c0 + np.array([0, 0, R/2])]
        Vext = Vn(xyz, nuc, a)
        r0 = np.linalg.norm(xyz - c0, axis=1)
        rho = np.exp(-r0**2 / 2); rho *= Ne / (rho.sum() * vol)
        iters = 12
        for it in range(iters):
            Ns = rho * vol
            VH = c * VH_raw(Ns)
            Veff = Vext + VH * (1 - 1.0/Ne)
            H = 0.5 / a**2 * Lm + sp.diags(Veff)
            eps, vec = spl.eigsh(H, k=1, which="SA")
            phi = vec[:, 0]; phi /= np.sqrt(vol * np.sum(phi**2))
            rho = 0.5 * rho + 0.5 * Ne * phi**2
        Ns = rho * vol; VH = c * VH_raw(Ns)
        phi3 = phi.reshape(n, n, n)
        T = 0.5 * np.sum(sum((np.roll(phi3, -1, ax) - phi3)**2 for ax in (0,1,2))) * vol / a**2
        Ttot = Ne * T; Een = np.sum(Ns * Vext); EH = 0.5 * np.sum(Ns * VH)
        Eelec = Ttot + Een + EH * (1 - 1.0/Ne)
        return {"R": R, "E_tot": Eelec + 1/R, "binding": Eelec + 1/R + 1,
                "T": Ttot, "Een": Een, "EH": EH}

    rows = []
    for R in (1.2, 1.4, 1.6):
        t = time.time(); rr = run(R); rr["time_s"] = round(time.time() - t, 1)
        rows.append(rr)
        print(f"  R={R}: E_tot={rr['E_tot']:+.4f}  bind={rr['binding']:+.4f}  T={rr['T']:.3f} Een={rr['Een']:+.3f} EH={rr['EH']:.3f} ({rr['time_s']}s)")
    OUT.write_text(json.dumps({"a": a, "c": c, "EH_raw_1s": EH_raw, "rows": rows}, ensure_ascii=False, indent=2))
    print("->", OUT)


if __name__ == "__main__":
    main()
