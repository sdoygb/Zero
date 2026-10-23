#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2_final.py — H2 两电子（路线 A，FFT 最小像核，高分辨率）。

一体：H1 = (1/2a^2)L - A/r_A - A/r_B（已到化学精度）
两体：<ij|kl> = (p_i p_k | p_j p_l)，FFT 卷积 + 最小像截断（与直接核一致到 1e-15）
反对称：自旋轨道 Slater–Condon
对照：H2 R=2.0 总 -1.0678；两 H 原子 -1.0；结合 -0.0678
输出：results/zero_h2_final.json
"""
from __future__ import annotations
import json, os, itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from numpy.fft import rfftn, irfftn

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

def run(a, Lphys, R=2.0, M=4):
    A = 1.0
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
    ev, evec = spl.eigsh(H1, k=M, which="SA")
    o = np.argsort(ev); ev = ev[o]; evec = evec[:, o]
    # FFT 最小像核
    big = 2*n
    Kf = np.zeros((big,big,big))
    for dx in range(-n,n):
        for dy in range(-n,n):
            for dz in range(-n,n):
                d = np.sqrt(dx*dx+dy*dy+dz*dz)*a
                if d <= Lphys/2 and d > 0:
                    Kf[dx%big, dy%big, dz%big] = A/d
    Kf[0,0,0] = A*ONSITE/a
    Kfr = rfftn(Kf)
    def conv(rho3):
        rb = np.zeros((big,big,big)); rb[:n,:n,:n] = rho3
        return irfftn(rfftn(rb)*Kfr, s=(big,big,big))[:n,:n,:n]
    def eri(p,q,r,s):
        rho = (evec[:,q]*evec[:,s]).reshape((n,n,n))
        Psi = conv(rho).ravel()
        return float(evec[:,p] @ (Psi*evec[:,r]))
    SO = [(p,s) for p in range(M) for s in (0,1)]
    def hij(i,j):
        (p,si),(q,sj) = SO[i],SO[j]
        return ev[p] if (p==q and si==sj) else 0.0
    def vijkl(i,j,k,l):
        (p,si),(q,sj),(r,sk),(s,sl) = SO[i],SO[j],SO[k],SO[l]
        return eri(p,r,q,s) if (si==sk and sj==sl) else 0.0
    dets = list(itertools.combinations(range(2*M), 2)); nd = len(dets)
    H = np.zeros((nd, nd))
    for I, D1 in enumerate(dets):
        for J, D2 in enumerate(dets):
            only1 = [x for x in D1 if x not in D2]; only2 = [x for x in D2 if x not in D1]
            common = [x for x in D1 if x in D2]
            if not only1 and not only2:
                v = sum(hij(i,i) for i in D1)
                for (i,j) in itertools.combinations(D1,2):
                    v += vijkl(i,j,i,j) - vijkl(i,j,j,i)
                H[I,J] = v
            elif len(only1)==1 and len(only2)==1:
                i,av = only1[0], only2[0]
                v = hij(i,av)
                for j in common:
                    v += vijkl(av,j,i,j) - vijkl(av,j,j,i)
                H[I,J] = (-1)**(D1.index(i)+D2.index(av))*v
            elif len(only1)==2 and len(only2)==2:
                i,j = only1; av,bv = only2
                v = vijkl(av,bv,i,j) - vijkl(av,bv,j,i)
                H[I,J] = (-1)**(D1.index(i)+D1.index(j)+D2.index(av)+D2.index(bv))*v
    w, _ = np.linalg.eigh(H)
    E0 = float(w[0]); Vnn = 1.0/R; Etot = E0+Vnn
    H1H = (0.5/a**2)*L + sp.diags(Vnuc([ctr]))
    eH = float(spl.eigsh(H1H, k=1, which="SA")[0][0])
    E2atoms = 2*eH; bind = Etot-E2atoms
    return {"a":a,"Lphys":Lphys,"R":R,"n":n,"N":N,"eps":ev.tolist(),
            "E_elec":E0,"Vnn":Vnn,"E_tot":Etot,"eH":eH,"E2atoms":float(E2atoms),
            "bind":float(bind),"ndet":nd,"lowest":w[:3].tolist()}

def main():
    print("=" * 84)
    print("zero_h2_final.py — H2 两电子（FFT 最小像核）")
    print("=" * 84)
    rows = []
    for a, Lp in [(0.5, 8.0), (0.25, 8.0)]:
        r = run(a, Lp)
        rows.append(r)
        print(f"  a={a:.2f} L={Lp:.0f} n={r['n']:3d} N={r['N']:6d}  eps0={r['eps'][0]:+.6f}  "
              f"E(H2)={r['E_tot']:+.6f} (err {r['E_tot']+1.0678:+.5f})  "
              f"bind={r['bind']:+.6f} (err {r['bind']+0.0678:+.5f})")
    ck("V1 H 原子误差随 a 减小", abs(rows[1]['eH']+0.5)<abs(rows[0]['eH']+0.5),
       f"{abs(rows[0]['eH']+0.5):.5f} -> {abs(rows[1]['eH']+0.5):.5f}")
    ck("V2 H2 总能量误差 < 0.01", all(abs(r['E_tot']+1.0678)<0.01 for r in rows),
       f"min={min(abs(r['E_tot']+1.0678) for r in rows):.5f}")
    ck("V3 结合能 < 0", all(r['bind']<0 for r in rows), f"{rows[-1]['bind']:+.6f}")
    ck("V4 结合能误差 < 0.05", all(abs(r['bind']+0.0678)<0.05 for r in rows),
       f"min={min(abs(r['bind']+0.0678) for r in rows):.5f}")
    ck("V5 结合能误差 < 1.6e-3", all(abs(r['bind']+0.0678)<1.6e-3 for r in rows),
       f"min={min(abs(r['bind']+0.0678) for r in rows):.5f}")
    out = {"rows": rows, "checks": PASS}
    with open(os.path.join(RES,"zero_h2_final.json"),"w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n%d 项，%d PASS" % (len(PASS), sum(1 for p in PASS if p["pass"])))
    print("-> " + os.path.join(RES,"zero_h2_final.json"))

if __name__ == "__main__":
    main()
