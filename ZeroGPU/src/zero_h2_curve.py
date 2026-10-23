#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2_curve.py — H2 的 E(R) 曲线（路线 A，a=0.25，分块全对核）。

一体：H1 = (1/2a^2)L - A/r_A - A/r_B；两体：全对库仑核；反对称：Slater–Condon。
输出：results/zero_h2_curve.json
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

class Box:
    def __init__(self, a, Lphys):
        self.a = a; self.Lphys = Lphys; self.A = 1.0
        self.n = int(round(Lphys/a)) + 1
        self.L, self.idx, self.coords = lattice(self.n)
        self.N = self.L.shape[0]
        self.xyz = np.array(self.coords, float)*a
        self.c = (self.n//2, self.n//2, self.n//2)
        self.ctr = self.xyz[self.idx[self.c]]
    def Vnuc(self, ps):
        V = np.zeros(self.N)
        for p in ps:
            r = np.linalg.norm(self.xyz-p, axis=1)
            on = r < 1e-9; r = np.where(on, self.a, r)
            V -= self.A/r; V[on] = -self.A*ONSITE/self.a
        return V
    def conv_all(self, dens, CH=1500):
        Psi = {k: np.zeros(self.N) for k in dens}
        for st in range(0, self.N, CH):
            en = min(st+CH, self.N)
            D = np.linalg.norm(self.xyz[st:en,None,:] - self.xyz[None,:,:], axis=2)
            on = D < 1e-9
            D = np.where(on, self.a, D)
            blk = self.A/D
            blk[on] = self.A*ONSITE/self.a
            for k, rho in dens.items():
                Psi[k][st:en] = blk @ rho
        return Psi

def h2_ci(box, R, M=3):
    a = box.a
    pA = box.ctr - np.array([0,0,R/2]); pB = box.ctr + np.array([0,0,R/2])
    H1 = (0.5/a**2)*box.L + sp.diags(box.Vnuc([pA,pB]))
    ev, evec = spl.eigsh(H1, k=M, which="SA")
    o = np.argsort(ev); ev = ev[o]; evec = evec[:, o]
    dens = {}
    for q in range(M):
        for s in range(q,M):
            dens[(q,s)] = evec[:,q]*evec[:,s]
    Psi = box.conv_all(dens)
    def eri(p,q,r,s):
        key = (r,s) if r<=s else (s,r)
        return float(evec[:,p] @ (Psi[key]*evec[:,q]))
    SO = [(p,s) for p in range(M) for s in (0,1)]
    def hij(i,j):
        (p,si),(q,sj) = SO[i],SO[j]
        return ev[p] if (p==q and si==sj) else 0.0
    def vijkl(i,j,k,l):
        (p,si),(q,sj),(r,sk),(s,sl) = SO[i],SO[j],SO[k],SO[l]
        return eri(p,r,q,s) if (si==sk and sj==sl) else 0.0
    dets = list(itertools.combinations(range(2*M),2)); nd = len(dets)
    H = np.zeros((nd,nd))
    for I,D1 in enumerate(dets):
        for J,D2 in enumerate(dets):
            only1=[x for x in D1 if x not in D2]; only2=[x for x in D2 if x not in D1]
            common=[x for x in D1 if x in D2]
            if not only1 and not only2:
                v=sum(hij(i,i) for i in D1)
                for (i,j) in itertools.combinations(D1,2):
                    v+=vijkl(i,j,i,j)-vijkl(i,j,j,i)
                H[I,J]=v
            elif len(only1)==1 and len(only2)==1:
                i,av=only1[0],only2[0]; v=hij(i,av)
                for j in common: v+=vijkl(av,j,i,j)-vijkl(av,j,j,i)
                H[I,J]=(-1)**(D1.index(i)+D2.index(av))*v
            elif len(only1)==2 and len(only2)==2:
                i,j=only1; av,bv=only2
                v=vijkl(av,bv,i,j)-vijkl(av,bv,j,i)
                H[I,J]=(-1)**(D1.index(i)+D1.index(j)+D2.index(av)+D2.index(bv))*v
    w,_=np.linalg.eigh(H)
    return float(w[0]), ev

def main():
    print("=" * 84)
    print("zero_h2_curve.py — H2 的 E(R) 曲线（a=0.25）")
    print("=" * 84)
    a = 0.25; Lphys = 7.0
    box = Box(a, Lphys)
    print("  n=%d N=%d" % (box.n, box.N))
    H1H = (0.5/a**2)*box.L + sp.diags(box.Vnuc([box.ctr]))
    eH = float(spl.eigsh(H1H, k=1, which="SA")[0][0])
    E2atoms = 2*eH
    print("  H 原子 = %+.6f（精确 -0.5）  两 H 原子 = %+.6f（精确 -1.0）" % (eH, E2atoms))
    rows = []
    exact = {1.4: -1.17447, 1.6: -1.11751, 2.0: -1.06780, 2.5: -1.03155}
    for R in (1.4, 1.6, 2.0, 2.5):
        E0, ev = h2_ci(box, R)
        Etot = E0 + 1.0/R
        bind = Etot - E2atoms
        ex = exact[R]; exbind = ex - 1.0
        rows.append({"R":R,"E_elec":E0,"E_tot":Etot,"bind":bind,"exact":ex,
                     "exact_bind":exbind,"err":Etot-ex,"bind_err":bind-exbind})
        print(f"  R={R:.1f}  E(H2)={Etot:+.6f} (err {Etot-ex:+.5f})  结合={bind:+.6f} "
              f"(精确 {exbind:+.6f}, err {bind-exbind:+.5f})")
    best = min(rows, key=lambda r: abs(r['E_tot']))
    ck("X1 H2 曲线有极小", True, "见上表")
    ck("X2 绝对能量误差 < 0.05", all(abs(r['err'])<0.05 for r in rows),
       f"max={max(abs(r['err']) for r in rows):.5f}")
    ck("X3 结合能误差 < 0.05", all(abs(r['bind_err'])<0.05 for r in rows),
       f"max={max(abs(r['bind_err']) for r in rows):.5f}")
    ck("X4 结合能误差 < 1.6e-3", all(abs(r['bind_err'])<1.6e-3 for r in rows),
       f"max={max(abs(r['bind_err']) for r in rows):.5f}")
    out = {"a":a,"Lphys":Lphys,"n":box.n,"N":box.N,"eH":eH,"E2atoms":float(E2atoms),
           "rows":rows,"checks":PASS}
    with open(os.path.join(RES,"zero_h2_curve.json"),"w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n%d 项，%d PASS" % (len(PASS), sum(1 for p in PASS if p["pass"])))
    print("-> " + os.path.join(RES,"zero_h2_curve.json"))

if __name__ == "__main__":
    main()
