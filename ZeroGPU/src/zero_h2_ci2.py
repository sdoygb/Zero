#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2_ci2.py — H2 两电子（路线 A，自旋轨道行列式 CI，Slater–Condon）。

一体：H1 = (1/2a^2)L - A/r_A - A/r_B（已验证到化学精度）
两体：V_abcd = Σ_ij φ_a,i φ_b,i K_ij φ_c,j φ_d,j，K_ij = A/|r_i-r_j|（解析库仑）
反对称：自旋轨道 Slater 行列式；标准 Slater–Condon 规则
对照：H2 R=2.0 a.u. 总能量 -1.0678 Ha（两 H 原子 -1.0；结合 -0.0678）
输出：results/zero_h2_ci2.json
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
    print("zero_h2_ci2.py — H2 自旋轨道 CI（Slater–Condon）")
    print("=" * 84)
    a = 0.5; Lphys = 10.0; R = 2.0; A = 1.0
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
    Vn = Vnuc([pA, pB])
    H1 = (0.5/a**2)*L + sp.diags(Vn)
    M = 4
    ev, evec = spl.eigsh(H1, k=M, which="SA")
    o = np.argsort(ev); ev = ev[o]; evec = evec[:, o]
    print("  H2+ 单电子最低 %d: %s" % (M, np.array2string(ev, precision=5)))
    D = np.linalg.norm(xyz[:,None,:] - xyz[None,:,:], axis=2)
    np.fill_diagonal(D, 1e9)
    K = A/D
    np.fill_diagonal(K, A*ONSITE/a)
    def eri(p,q,r,s):
        return float((evec[:,p]*evec[:,q]) @ K @ (evec[:,r]*evec[:,s]))
    # 自旋轨道：2M 个，(p,sigma)
    SO = [(p,s) for p in range(M) for s in (0,1)]
    NS = len(SO)
    def hij(i,j):
        (p,si),(q,sj) = SO[i],SO[j]
        return ev[p] if (p==q and si==sj) else 0.0
    def vijkl(i,j,k,l):
        """物理学家记法 <ij|kl> = ∫ φ_i(1)φ_j(2) (1/r) φ_k(1)φ_l(2)；
        用化学记法 (p q | r s) = ∫ φ_p(1)φ_q(1) (1/r) φ_r(2)φ_s(2) 表示：
        <ij|kl> = (p_i p_k | p_j p_l)，自旋需匹配 σ_i=σ_k, σ_j=σ_l。"""
        (p,si),(q,sj),(r,sk),(s,sl) = SO[i],SO[j],SO[k],SO[l]
        if si==sk and sj==sl:
            return eri(p,r,q,s)
        return 0.0
    dets = list(itertools.combinations(range(NS), 2))
    nd = len(dets)
    print("  自旋轨道 %d 个，行列式 %d 个" % (NS, nd))
    H = np.zeros((nd, nd))
    for I, D1 in enumerate(dets):
        for J, D2 in enumerate(dets):
            common = [x for x in D1 if x in D2]
            only1 = [x for x in D1 if x not in D2]
            only2 = [x for x in D2 if x not in D1]
            if len(only1)==0 and len(only2)==0:
                val = sum(hij(i,i) for i in D1)
                for (i,j) in itertools.combinations(D1,2):
                    val += vijkl(i,j,i,j) - vijkl(i,j,j,i)
                H[I,J] = val
            elif len(only1)==1 and len(only2)==1:
                i, av = only1[0], only2[0]
                val = hij(i,av)
                for j in common:
                    val += vijkl(av,j,i,j) - vijkl(av,j,j,i)
                pos_i = D1.index(i); pos_a = D2.index(av)
                phase = (-1)**(pos_i + pos_a)
                H[I,J] = phase*val
            elif len(only1)==2 and len(only2)==2:
                i,j = only1; av,bv = only2
                val = vijkl(av,bv,i,j) - vijkl(av,bv,j,i)
                pi = D1.index(i)+D1.index(j); pa = D2.index(av)+D2.index(bv)
                H[I,J] = (-1)**(pi+pa)*val
    w, vec = np.linalg.eigh(H)
    E0 = float(w[0])
    # 总能量 = 电子 + 核排斥 1/R
    Vnn = 1.0/R
    Etot = E0 + Vnn
    H1H = (0.5/a**2)*L + sp.diags(Vnuc([ctr]))
    eH = float(spl.eigsh(H1H, k=1, which="SA")[0][0])
    #print("    [debug] min Vnuc([ctr]) = %+.4f  min Vnuc([pA,pB]) = %+.4f" % (Vnuc([ctr]).min(), Vnuc([pA,pB]).min()))
    _V=Vnuc([ctr])
    #print("    [debug] 单核矩阵元数 =", H1H.nnz, " 两核 =", H1.nnz)
    E2atoms = 2*eH
    bind = Etot - E2atoms
    print("  H2 电子基态 = %+.6f   +Vnn(%.3f) = 总 %+.6f" % (E0, Vnn, Etot))
    print("  两个独立 H 原子 = %+.6f" % E2atoms)
    print("  结合能 = %+.6f（精确 -0.0678）" % bind)
    print("  最低 3 个本征值: %s" % np.array2string(w[:3], precision=5))
    ck("T1 H2 总能量 < 0", Etot < 0, f"{Etot:+.6f}")
    ck("T2 结合能 < 0（成键）", bind < 0, f"{bind:+.6f}")
    ck("T3 结合能误差 < 0.05", abs(bind+0.0678) < 0.05, f"err={bind+0.0678:+.6f}")
    ck("T4 结合能误差 < 1.6e-3", abs(bind+0.0678) < 1.6e-3, f"err={bind+0.0678:+.6f}")
    out = {"a":a,"Lphys":Lphys,"R":R,"n":n,"N":N,"eps":ev.tolist(),"E_elec":E0,
           "Vnn":Vnn,"E_tot":Etot,"E2atoms":float(E2atoms),"bind":float(bind),
           "ndet":nd,"lowest":w[:4].tolist(),"checks":PASS}
    with open(os.path.join(RES,"zero_h2_ci2.json"),"w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n%d 项，%d PASS" % (len(PASS), sum(1 for p in PASS if p["pass"])))
    print("-> " + os.path.join(RES,"zero_h2_ci2.json"))

if __name__ == "__main__":
    main()
