#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2_ci3.py — H2 两电子（路线 A，分块核，更高分辨率）。

一体：H1 = (1/2a^2)L - A/r_A - A/r_B
两体：<ij|kl> = (p_i p_k | p_j p_l) 化学记法，直接分块算（不存 N^2 核）
反对称：自旋轨道 Slater–Condon（已修正记法）
对照：H2 R=2.0 总能量 -1.0678；两 H 原子 -1.0；结合 -0.0678
输出：results/zero_h2_ci3.json
"""
from __future__ import annotations
import json, os, itertools
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
    print("zero_h2_ci3.py — H2 两电子（分块核，高分辨率）")
    print("=" * 84)
    a = 0.25; Lphys = 8.0; R = 2.0; A = 1.0
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
    M = 4
    ev, evec = spl.eigsh(H1, k=M, which="SA")
    o = np.argsort(ev); ev = ev[o]; evec = evec[:, o]
    print("  H2+ 单电子最低 %d: %s" % (M, np.array2string(ev, precision=6)))
    def matvec(rho):
        """Ψ_i = Σ_j A/|r_i-r_j| ρ_j，分块。"""
        out = np.zeros(N)
        CH = 2000
        for s in range(0, N, CH):
            e = min(s+CH, N)
            D = np.linalg.norm(xyz[s:e,None,:] - xyz[None,:,:], axis=2)
            on = D < 1e-9
            D = np.where(on, a, D)
            blk = A/D
            blk[on] = A*ONSITE/a
            out[s:e] = blk @ rho
        return out
    def eri(p,q,r,s):
        rho = evec[:,q]*evec[:,s]        # 化学记法 (p q | r s) = ∫ φ_p(1)φ_q(1) K φ_r(2)φ_s(2)
        Psi = matvec(rho)
        return float(evec[:,p] @ (Psi * evec[:,r]))
    SO = [(p,s) for p in range(M) for s in (0,1)]
    def hij(i,j):
        (p,si),(q,sj) = SO[i],SO[j]
        return ev[p] if (p==q and si==sj) else 0.0
    def vijkl(i,j,k,l):
        (p,si),(q,sj),(r,sk),(s,sl) = SO[i],SO[j],SO[k],SO[l]
        if si==sk and sj==sl:
            return eri(p,r,q,s)
        return 0.0
    dets = list(itertools.combinations(range(2*M), 2)); nd = len(dets)
    print("  自旋轨道 %d，行列式 %d" % (2*M, nd))
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
    w, vec = np.linalg.eigh(H)
    E0 = float(w[0]); Vnn = 1.0/R; Etot = E0+Vnn
    H1H = (0.5/a**2)*L + sp.diags(Vnuc([ctr]))
    eH = float(spl.eigsh(H1H, k=1, which="SA")[0][0])
    E2atoms = 2*eH; bind = Etot - E2atoms
    print("  H2 电子 %+.6f + Vnn %.3f = 总 %+.6f（精确 -1.0678）" % (E0, Vnn, Etot))
    print("  两个独立 H 原子 = %+.6f（精确 -1.0）" % E2atoms)
    print("  结合能 = %+.6f（精确 -0.0678）" % bind)
    ck("U1 H 原子误差 < 0.01", abs(eH+0.5)<0.01, f"eH={eH:+.6f}")
    ck("U2 H2 总能量误差 < 0.01", abs(Etot+1.0678)<0.01, f"err={Etot+1.0678:+.6f}")
    ck("U3 结合能 < 0", bind<0, f"{bind:+.6f}")
    ck("U4 结合能误差 < 0.05", abs(bind+0.0678)<0.05, f"err={bind+0.0678:+.6f}")
    ck("U5 结合能误差 < 1.6e-3", abs(bind+0.0678)<1.6e-3, f"err={bind+0.0678:+.6f}")
    out = {"a":a,"Lphys":Lphys,"R":R,"n":n,"N":N,"eps":ev.tolist(),"E_elec":E0,
           "Vnn":Vnn,"E_tot":Etot,"E2atoms":float(E2atoms),"bind":float(bind),
           "ndet":nd,"checks":PASS}
    with open(os.path.join(RES,"zero_h2_ci3.json"),"w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n%d 项，%d PASS" % (len(PASS), sum(1 for p in PASS if p["pass"])))
    print("-> " + os.path.join(RES,"zero_h2_ci3.json"))

if __name__ == "__main__":
    main()
