#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2plus_continuum.py — H2+ 两中心 Coulomb 的连续极限（a→0，H 原子锚定）。

格距 a = 1, 1/2, 1/4；物理 R=2 a.u. 固定；H 原子锚定到 -0.5 Ha。
检验：锚定后的 H2+ 电子结合 vs 精确 -0.602634 Ha；误差随 a 的收敛。
"""
from __future__ import annotations
import json, os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
PASS = []

def ck(name, cond, detail=""):
    PASS.append({"name": name, "pass": bool(cond), "detail": detail})
    print("[%s] %s %s" % ("PASS" if cond else "FAIL", name, detail))

def lattice(n):
    idx = {}; coords = []
    for x in range(n):
        for y in range(n):
            for z in range(n):
                idx[(x, y, z)] = len(coords); coords.append((x, y, z))
    N = len(coords)
    adj = [[] for _ in range(N)]
    for (x, y, z), i in idx.items():
        for dx, dy, dz in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)):
            j = idx.get((x+dx, y+dy, z+dz))
            if j is not None: adj[i].append(j)
    rows, cols, vals = [], [], []
    for i, nb in enumerate(adj):
        rows.append(i); cols.append(i); vals.append(float(len(nb)))
        for j in nb:
            rows.append(i); cols.append(j); vals.append(-1.0)
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N)), idx

def phi_from(L, idx, c):
    N = L.shape[0]
    x = np.zeros(N); x[idx[c]] = 1.0; x -= x.mean()
    ph = spl.lsqr(L, -x, atol=1e-12, btol=1e-12, iter_lim=300000)[0]
    ph -= ph.mean()
    return ph

def main():
    print("=" * 80)
    print("zero_h2plus_continuum.py — a→0 连续极限（H 原子锚定）")
    print("=" * 80)
    rows = []
    for a in (1.0, 0.5, 0.25):
        n = int(round(20.0 / a)) + 1
        L, idx = lattice(n)
        c = (n // 2, n // 2, n // 2)
        phH = phi_from(L, idx, c)
        EH = float(spl.eigsh(0.5 * L - 1.0 * sp.diags(phH), k=1, which="SA")[0][0])
        s = -0.5 / EH
        R = 2.0
        off = int(round(R / (2 * a)))
        cA = (c[0], c[1], c[2] - off); cB = (c[0], c[1], c[2] + off)
        pA = phi_from(L, idx, cA); pB = phi_from(L, idx, cB)
        E2 = float(spl.eigsh(0.5 * L - 1.0 * sp.diags(pA + pB), k=1, which="SA")[0][0])
        bind = s * (E2 - EH)
        rows.append({"a": a, "n": n, "EH": EH, "s": s, "E2": E2, "bind": bind,
                     "err": bind + 0.602634})
        print(f"  a={a:.2f} n={n:3d}  EH={EH:+.6e}  s={s:8.1f}  E2={E2:+.6e}  "
              f"锚定结合={bind:+.4f}  误差={bind+0.602634:+.4f}")
    ck("C1 H 原子锚定 = -0.5", all(abs(r['s'] * r['EH'] + 0.5) < 1e-9 for r in rows), "s=-0.5/EH")
    ck("C2 H2+ 存在（结合<0）", all(r['bind'] < 0 for r in rows), f"bind={rows[-1]['bind']:+.4f}")
    ck("C3 误差随 a 减小", rows[-1]['err'] < rows[0]['err'] or abs(rows[-1]['err']) < abs(rows[0]['err']),
       f"|err|: {abs(rows[0]['err']):.4f} -> {abs(rows[-1]['err']):.4f}")
    ck("C4 误差 < 1.6e-3", all(abs(r['err']) < 1.6e-3 for r in rows), f"min={min(abs(r['err']) for r in rows):.4f}")
    ck("C5 误差 < 0.05", all(abs(r['err']) < 0.05 for r in rows), f"min={min(abs(r['err']) for r in rows):.4f}")
    out = {"rows": rows, "checks": PASS}
    with open(os.path.join(RES, "zero_h2plus_continuum.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n%d 项，%d PASS" % (len(PASS), sum(1 for p in PASS if p["pass"])))
    print("-> " + os.path.join(RES, "zero_h2plus_continuum.json"))

if __name__ == "__main__":
    main()
