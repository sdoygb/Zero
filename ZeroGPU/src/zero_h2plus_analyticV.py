#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2plus_analyticV.py — 路线 A：Zero 动能（图拉普拉斯，1/a^2）+ 解析库仑形状 1/r。

Zero 给：d=3（Z021）、幅度 A=1/(σ ℓ_b)（Z147，立方=1.0000）、动能=图拉普拉斯（Z033）。
外部给：1/r 的函数形式（解析）。
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
    return sp.csr_matrix((vals,(rows,cols)),shape=(N,N)), idx, coords

def main():
    print("=" * 84)
    print("zero_h2plus_analyticV.py — 路线 A：图拉普拉斯动能 + 解析 1/r")
    print("=" * 84)
    Lphys = 8.0; R = 2.0; A = 1.0   # 原子单位库仑 -1/r（e^2/(4πε0)=1）；Zero 的 4π 归一化由 g=4π 吸收
    rows = []
    for a in (0.5, 0.25, 0.125):
        n = int(round(Lphys/a)) + 1
        L, idx, coords = lattice(n)
        N = L.shape[0]
        c = (n//2, n//2, n//2)
        # 物理坐标
        xyz = np.array(coords, float)*a
        ctr = xyz[idx[c]]
        def V_of(protons):
            V = np.zeros(N)
            for p in protons:
                r = np.linalg.norm(xyz - p, axis=1)
                # 立方格胞平均 1/r 的奇点正则化：在源点用 2.3801/a（标准立方自能常数）
                onsite = r < 1e-9
                r = np.where(onsite, a, r)
                V -= A/r
                V[onsite] = -A*2.3801/a
            return V
        VH = V_of([ctr])
        H1 = (0.5/a**2)*L + sp.diags(VH)
        EH = float(spl.eigsh(H1, k=1, which="SA")[0][0])
        cA = ctr - np.array([0,0,R/2]); cB = ctr + np.array([0,0,R/2])
        VAB = V_of([cA, cB])
        H2 = (0.5/a**2)*L + sp.diags(VAB)
        E2 = float(spl.eigsh(H2, k=1, which="SA")[0][0])
        bind = E2 - EH
        rows.append({"a":a,"n":n,"N":N,"EH":EH,"E2":E2,"bind":bind,"err":bind+0.602634})
        print(f"  a={a:.2f} n={n:3d} N={N:6d}  EH={EH:+.8f}  E2={E2:+.8f}  "
              f"电子结合={bind:+.6f}  误差={bind+0.602634:+.6f}")
    ck("G1 H 原子 E0 -> -0.5", all(abs(r['EH']+0.5)<0.05 for r in rows), f"EH(最细)={rows[-1]['EH']:+.6f}")
    ck("G2 H 原子误差随 a 减小", all(abs(rows[i+1]['EH']+0.5)<abs(rows[i]['EH']+0.5) for i in range(len(rows)-1)),
       f"|err|: {abs(rows[0]['EH']+0.5):.5f} -> {abs(rows[-1]['EH']+0.5):.5f}")
    ck("G3 H2+ 结合 < 0", all(r['bind']<0 for r in rows), f"{rows[-1]['bind']:+.6f}")
    ck("G4 H2+ 误差随 a 减小", all(abs(rows[i+1]['err'])<abs(rows[i]['err']) for i in range(len(rows)-1)),
       f"|err|: {abs(rows[0]['err']):.5f} -> {abs(rows[-1]['err']):.5f}")
    ck("G5 H2+ 误差 < 0.05", all(abs(r['err'])<0.05 for r in rows), f"min={min(abs(r['err']) for r in rows):.5f}")
    ck("G6 H2+ 误差 < 1.6e-3", all(abs(r['err'])<1.6e-3 for r in rows), f"min={min(abs(r['err']) for r in rows):.5f}")
    out = {"rows": rows, "checks": PASS}
    with open(os.path.join(RES,"zero_h2plus_analyticV.json"),"w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n%d 项，%d PASS" % (len(PASS), sum(1 for p in PASS if p["pass"])))
    print("-> " + os.path.join(RES,"zero_h2plus_analyticV.json"))

if __name__ == "__main__":
    main()
