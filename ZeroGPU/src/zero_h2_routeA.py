#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2_routeA.py — H2 两电子（路线 A：Zero 动能 + Zero 幅度 + 解析 1/r）。

一体：H1 = (1/2a^2)L - A/r_A - A/r_B（已验证到化学精度）
两体：V_abcd = Σ_ij φ_a,i φ_b,i K_ij φ_c,j φ_d,j，K_ij = A/|r_i-r_j|（解析库仑，源点立方胞平均）
反对称：Slater（singlet 空间对称）
对照：H2 R=2.0 a.u. 总能量 -1.0678 Ha（两 H 原子 -1.0；结合 -0.0678）
输出：results/zero_h2_routeA.json
"""
from __future__ import annotations
import json, os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
PASS = []
ONSITE = 2.3800756243710857

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
    print("zero_h2_routeA.py — H2 两电子（路线 A，直接核）")
    print("=" * 84)
    a = 0.5; Lphys = 8.0; R = 2.0; A = 1.0
    n = int(round(Lphys/a)) + 1
    L, idx, coords = lattice(n); N = L.shape[0]
    xyz = np.array(coords, float)*a
    c = (n//2, n//2, n//2); ctr = xyz[idx[c]]
    pA = ctr - np.array([0,0,R/2]); pB = ctr + np.array([0,0,R/2])
    def Vnuc(ps):
        V = np.zeros(N)
        for p in ps:
            r = np.linalg.norm(xyz-p, axis=1)
            on = r < 1e-9; r = np.where(on, a, r)
            V -= A/r; V[on] = -A*ONSITE/a
        return V
    H1 = (0.5/a**2)*L + sp.diags(Vnuc([pA, pB]))
    M = 6
    ev, evec = spl.eigsh(H1, k=M, which="SA")
    o = np.argsort(ev); ev = ev[o]; evec = evec[:, o]
    print("  H2+ 单电子最低 %d: %s" % (M, np.array2string(ev, precision=5)))
    # 核 K（直接，N^2）
    D = np.linalg.norm(xyz[:,None,:] - xyz[None,:,:], axis=2)
    np.fill_diagonal(D, 1e9)
    K = A/D
    np.fill_diagonal(K, A*ONSITE/a)
    # 两体：V_abcd = (φ_aφ_b)^T K (φ_cφ_d)
    def two_body(aa, bb, cc, dd):
        return float((evec[:,aa]*evec[:,bb]) @ K @ (evec[:,cc]*evec[:,dd]))
    pairs = [(x,y) for x in range(M) for y in range(x,M)]
    print("  CI 维度 = %d" % len(pairs))
    H = np.zeros((len(pairs), len(pairs)))
    for I,(x,y) in enumerate(pairs):
        for J,(cc,dd) in enumerate(pairs):
            v1 = 0.0
            if y == dd: v1 += ev[x]
            if x == cc: v1 += ev[y]
            v2 = 0.0
            for (x1,y1) in ((x,y),(y,x)):
                for (c1,d1) in ((cc,dd),(dd,cc)):
                    v2 += 0.25*two_body(x1,y1,c1,d1)
            H[I,J] = v1+v2
    w, vec = np.linalg.eigh(H)
    E0 = float(w[0])
    H1H = (0.5/a**2)*L + sp.diags(Vnuc([ctr]))
    eH = float(spl.eigsh(H1H, k=1, which="SA")[0][0])
    E2atoms = 2*eH
    bind = E0 - E2atoms
    print("  H2 两电子基态 = %+.6f" % E0)
    print("  两个独立 H 原子 = %+.6f" % E2atoms)
    print("  结合能 = %+.6f（精确 -0.0678）" % bind)
    print("  σ_g^2 单组态能量 = %+.6f" % H[0,0])
    ck("S1 H2 两电子基态存在", E0 < 0, f"E0={E0:+.6f}")
    ck("S2 结合能 < 0（成键）", bind < 0, f"bind={bind:+.6f}")
    ck("S3 结合能误差 < 0.05", abs(bind+0.0678) < 0.05, f"err={bind+0.0678:+.6f}")
    ck("S4 结合能误差 < 1.6e-3", abs(bind+0.0678) < 1.6e-3, f"err={bind+0.0678:+.6f}")
    out = {"a":a,"Lphys":Lphys,"R":R,"n":n,"N":N,"eps":ev.tolist(),"E0":E0,
           "E2atoms":float(E2atoms),"bind":float(bind),"sigma_g2":float(H[0,0]),
           "dim":len(pairs),"checks":PASS}
    with open(os.path.join(RES,"zero_h2_routeA.json"),"w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n%d 项，%d PASS" % (len(PASS), sum(1 for p in PASS if p["pass"])))
    print("-> " + os.path.join(RES,"zero_h2_routeA.json"))

if __name__ == "__main__":
    main()
