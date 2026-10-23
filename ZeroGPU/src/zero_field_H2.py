#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_field_H2.py — 净荷 FA 泛函：H2（2 电子）单轨道 SCF。

核 = 路线 A 的立方胞平均：离源 -1/r、在源 -2.3800756/a。
V_eff = V_n + V_H(1-1/N)（N=2 → V_H/2）；E_elec = T + V_en + E_H(1-1/N)；
总能量 = E_elec + 1/R。对照：H2 R=1.4 精确总能量约 -1.174 Ha（结合 0.174）。

运行：/usr/local/bin/python3 ZeroGPU/src/zero_field_H2.py
输出：results/zero_field_H2.json
"""
from __future__ import annotations
import json, pathlib, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.fft import fftn, ifftn

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "zero_field_H2.json"
ONSITE = 2.3800756243710857
OK, BAD, ROWS = [], [], []


def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def lap(n):
    k = 2 * np.pi * np.fft.fftfreq(n)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    return 2 * (3 - np.cos(KX) - np.cos(KY) - np.cos(KZ))


def lattice(n):
    idx = {}
    coords = []
    for x in range(n):
        for y in range(n):
            for z in range(n):
                idx[(x, y, z)] = len(coords)
                coords.append((x, y, z))
    N = len(coords)
    rows, cols, vals = [], [], []
    for (x, y, z), i in idx.items():
        nb = 0
        for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            j = idx.get((x + dx, y + dy, z + dz))
            if j is not None:
                nb += 1
                rows.append(i); cols.append(j); vals.append(-1.0)
        rows.append(i); cols.append(i); vals.append(float(nb))
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N)), idx, coords


def Vn(xyz, ps, a):
    V = np.zeros(len(xyz))
    for p in ps:
        r = np.linalg.norm(xyz - p, axis=1)
        on = r < 1e-9
        r = np.where(on, a, r)
        V += np.where(on, -ONSITE / a, -1.0 / r)
    return V


def run(a=0.25, R=1.4, L=14.0, Ne=2, iters=8):
    t0 = time.time()
    n = int(round(L / a)) + 1
    Lm, idx, coords = lattice(n)
    xyz = np.array(coords, float) * a
    c = xyz[idx[(n // 2, n // 2, n // 2)]]
    nuc = [c - np.array([0, 0, R / 2]), c + np.array([0, 0, R / 2])]
    Vext = Vn(xyz, nuc, a)
    lam = lap(n); inv = np.zeros_like(lam); m = lam > 1e-12; inv[m] = 1 / lam[m]
    G0 = np.real(ifftn(inv))[0, 0, 0]
    vol = a ** 3
    # initial density: gaussian on the bond centre
    r0 = np.linalg.norm(xyz - c, axis=1)
    rho = np.exp(-r0 ** 2 / (2 * 1.0 ** 2)); rho *= Ne / (rho.sum() * vol)
    def hartree(rho1):
        n3 = (rho1 * vol).reshape(n, n, n)
        v3 = 4 * np.pi * np.real(ifftn(inv * fftn(n3))) / a - (4 * np.pi * G0 - ONSITE) / a * n3
        return v3.reshape(-1), n3
    for it in range(iters):
        VH, _ = hartree(rho)
        Veff = Vext + VH * (1 - 1.0 / Ne)
        H = 0.5 / a ** 2 * Lm + sp.diags(Veff)
        eps, vec = spl.eigsh(H, k=1, which="SA")
        phi = vec[:, 0]
        phi = phi / np.sqrt(vol * np.sum(phi ** 2))
        rho_new = Ne * phi ** 2
        rho = 0.5 * rho + 0.5 * rho_new
    VH, n3 = hartree(rho)
    phi3 = phi.reshape(n, n, n)
    T = 0.5 * np.sum(sum((np.roll(phi3, -1, ax) - phi3) ** 2 for ax in (0, 1, 2))) * vol / a ** 2
    Ttot = Ne * T
    Een = np.sum(rho * vol * Vext)
    EH = 0.5 * np.sum(n3 * VH.reshape(n, n, n))
    Eelec = Ttot + Een + EH * (1 - 1.0 / Ne)
    Etot = Eelec + 1.0 / R
    return {"a": a, "n": n, "R": R, "E_elec": Eelec, "E_tot": Etot,
            "binding": Etot + 1.0, "T": Ttot, "Een": Een, "EH": EH,
            "eps": float(eps[0]), "time_s": round(time.time() - t0, 1)}


def main():
    print("== 净荷 FA：H2 单轨道 SCF ==")
    for R in (1.2, 1.4, 1.6):
        r = run(R=R)
        print(f"  R={R} a={r['a']} n={r['n']}  E_elec={r['E_elec']:+.4f}  E_tot={r['E_tot']:+.4f}  "
              f"bind={r['binding']:+.4f}  T={r['T']:.3f} Een={r['Een']:+.3f} EH={r['EH']:.3f} ({r['time_s']}s)")
    out = {"rows": [], "checks": ROWS, "pass": len(OK), "fail": len(BAD)}
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
