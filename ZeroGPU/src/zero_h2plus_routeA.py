#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2plus_routeA.py — 路线 A 定量验证：Zero 动能 + Zero 幅度 + 解析 1/r 形状。

Zero 给：d=3（Z021）、幅度 A（Z147：A=1/(σℓ_b)，立方舞台归一化到物理库仑）、
        动能 = 舞台图拉普拉斯（Z033，量纲 1/a^2）。
外部给：1/r 的函数形式（解析；Z147 说形状由 d=3 给，此处只借它的显式函数）。
正则化：立方胞平均 1/r 在源点 = 2.3800756/a（数值核验）。
对照：H 原子 -0.5 Ha；H2+ R=2.0 电子结合 -0.602634 Ha；化学精度 1.6e-3 Ha。
输出：results/zero_h2plus_routeA.json
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

ONSITE = 2.3800756243710857   # 立方胞平均 1/r（a=1）

def run(a, Lphys, R=2.0, A=1.0):
    n = int(round(Lphys/a)) + 1
    L, idx, coords = lattice(n)
    N = L.shape[0]
    xyz = np.array(coords, float)*a
    c = (n//2, n//2, n//2)
    ctr = xyz[idx[c]]
    def V_of(ps):
        V = np.zeros(N)
        for p in ps:
            r = np.linalg.norm(xyz - p, axis=1)
            on = r < 1e-9
            r = np.where(on, a, r)
            V -= A/r
            V[on] = -A*ONSITE/a
        return V
    H1 = (0.5/a**2)*L + sp.diags(V_of([ctr]))
    EH = float(spl.eigsh(H1, k=1, which="SA")[0][0])
    H2 = (0.5/a**2)*L + sp.diags(V_of([ctr - np.array([0,0,R/2]), ctr + np.array([0,0,R/2])]))
    E2 = float(spl.eigsh(H2, k=1, which="SA")[0][0])
    return {"a":a,"Lphys":Lphys,"n":n,"N":N,"EH":EH,"E2":E2,"bind":E2-EH,
            "errH":EH+0.5,"errB":E2-EH+0.602634}

def main():
    print("=" * 84)
    print("zero_h2plus_routeA.py — Zero 动能 + 解析 1/r 形状（定量）")
    print("=" * 84)
    rows = []
    for a, Lp in [(0.25, 14.0), (0.125, 14.0)]:
        r = run(a, Lp)
        rows.append(r)
        print(f"  a={a:.3f} L={Lp:.0f} n={r['n']:3d} N={r['N']:7d}  "
              f"EH={r['EH']:+.8f} (err {r['errH']:+.2e})  "
              f"bind={r['bind']:+.8f} (err {r['errB']:+.2e})")
    ck("R1 H 原子误差 < 1.6e-3", all(abs(r['errH'])<1.6e-3 for r in rows),
       f"min={min(abs(r['errH']) for r in rows):.2e}")
    ck("R2 H2+ 结合误差 < 1.6e-3", all(abs(r['errB'])<1.6e-3 for r in rows),
       f"min={min(abs(r['errB']) for r in rows):.2e}")
    ck("R3 误差随 a 减小", all(abs(rows[i+1]['errB'])<abs(rows[i]['errB']) for i in range(len(rows)-1)),
       f"|errB|: {abs(rows[0]['errB']):.2e} -> {abs(rows[-1]['errB']):.2e}")
    ck("R4 H2+ 成键（bind<0）", all(r['bind']<0 for r in rows), f"{rows[-1]['bind']:+.6f}")
    out = {"rows": rows, "checks": PASS}
    with open(os.path.join(RES,"zero_h2plus_routeA.json"),"w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n%d 项，%d PASS" % (len(PASS), sum(1 for p in PASS if p["pass"])))
    print("-> " + os.path.join(RES,"zero_h2plus_routeA.json"))

if __name__ == "__main__":
    main()
