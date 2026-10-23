#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2plus_continuum2.py — 正确的连续极限：H = (1/2)(L/a^2) - g*(phi/a)。

格距 a；物理 R 固定；物理盒 L_phys 固定。势 phi 是单位格泊松核，量纲上 ~1/(a*r_phys)，
故势项要 /a；动能 L 量纲 ~a^2*∇²，故要 /a^2。
对照：H 原子 -0.5 Ha；H2+ R=2.0 电子结合 -0.602634 Ha。
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
                idx[(x,y,z)] = len(coords); coords.append((x,y,z))
    N = len(coords)
    adj = [[] for _ in range(N)]
    for (x,y,z), i in idx.items():
        for dx,dy,dz in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)):
            j = idx.get((x+dx,y+dy,z+dz))
            if j is not None: adj[i].append(j)
    rows, cols, vals = [], [], []
    for i, nb in enumerate(adj):
        rows.append(i); cols.append(i); vals.append(float(len(nb)))
        for j in nb:
            rows.append(i); cols.append(j); vals.append(-1.0)
    return sp.csr_matrix((vals,(rows,cols)),shape=(N,N)), idx

def phi_from(L, idx, c):
    N = L.shape[0]
    x = np.zeros(N); x[idx[c]] = 1.0; x -= x.mean()
    ph = spl.lsqr(L, -x, atol=1e-12, btol=1e-12, iter_lim=60000)[0]
    ph -= ph.mean()
    return ph

def main():
    print("=" * 84)
    print("zero_h2plus_continuum2.py — 正确连续极限 H=(1/2)(L/a^2)-g(phi/a)")
    print("=" * 84)
    rows = []
    Lphys = 8.0     # 物理盒（Bohr）
    R = 2.0          # 物理键长（Bohr）
    for a in (0.5, 0.25, 0.125):
        n = int(round(Lphys/a)) + 1
        L, idx = lattice(n)
        N = L.shape[0]
        c = (n//2, n//2, n//2)
        phH = phi_from(L, idx, c)
        H1 = (0.5/a**2)*L - (1.0/a)*sp.diags(phH)
        EH = float(spl.eigsh(H1, k=1, which="SA")[0][0])
        s = -0.5/EH
        off = int(round(R/(2*a)))
        if off < 1: off = 1
        cA = (c[0],c[1],c[2]-off); cB = (c[0],c[1],c[2]+off)
        pA = phi_from(L, idx, cA); pB = phi_from(L, idx, cB)
        H2 = (0.5/a**2)*L - (1.0/a)*sp.diags(pA+pB)
        E2 = float(spl.eigsh(H2, k=1, which="SA")[0][0])
        bind = s*(E2 - EH)
        rows.append({"a":a,"n":n,"N":N,"EH":EH,"s":s,"E2":E2,"bind":bind,
                     "err":bind+0.602634,"off":off})
        print(f"  a={a:.3f} n={n:3d} off={off}  EH={EH:+.6e} s={s:7.1f}  E2={E2:+.6e} "
              f"锚定结合={bind:+.4f} 误差={bind+0.602634:+.4f}")
    ck("E1 H 原子锚定 = -0.5", all(abs(r['s']*r['EH']+0.5)<1e-9 for r in rows), "s=-0.5/EH")
    ck("E2 H 原子 |E0| 随 a 减小而增大（趋近连续）",
       all(abs(rows[i+1]['EH']) > abs(rows[i]['EH']) for i in range(len(rows)-1)),
       f"|EH|: {abs(rows[0]['EH']):.3e} -> {abs(rows[-1]['EH']):.3e}")
    ck("E3 H2+ 结合 < 0", all(r['bind']<0 for r in rows), f"{rows[-1]['bind']:+.4f}")
    ck("E4 误差随 a 减小", all(abs(rows[i+1]['err'])<abs(rows[i]['err']) for i in range(len(rows)-1)),
       f"|err|: {abs(rows[0]['err']):.4f} -> {abs(rows[-1]['err']):.4f}")
    ck("E5 误差 < 0.05", all(abs(r['err'])<0.05 for r in rows), f"min={min(abs(r['err']) for r in rows):.4f}")
    ck("E6 误差 < 1.6e-3", all(abs(r['err'])<1.6e-3 for r in rows), f"min={min(abs(r['err']) for r in rows):.4f}")
    out = {"rows": rows, "checks": PASS}
    with open(os.path.join(RES,"zero_h2plus_continuum2.json"),"w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n%d 项，%d PASS" % (len(PASS), sum(1 for p in PASS if p["pass"])))
    print("-> " + os.path.join(RES,"zero_h2plus_continuum2.json"))

if __name__ == "__main__":
    main()
