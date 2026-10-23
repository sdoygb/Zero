#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2_basis.py — 按对称性+局域性选 CI 基组，重跑 H2 的 E(R)。

诊断（§4.15）：最低 M 个本征态混入有限盒驻波（盒态），M=5 反而更差。
修法：对每个候选本征态算
  · z 镜像宇称（质子中点对称/反对称）—— σg/σu 的格点版
  · 局域性：|ψ|² 在核附近半径 R_loc 内的占比（盒态占比低）
只保留"局域 + 宇称明确"的轨道。
输出：results/zero_h2_basis.json
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
        self.a=a; self.Lphys=Lphys; self.A=1.0
        self.n=int(round(Lphys/a))+1
        self.L,self.idx,self.coords=lattice(self.n)
        self.N=self.L.shape[0]
        self.xyz=np.array(self.coords,float)*a
        self.c=(self.n//2,self.n//2,self.n//2)
        self.ctr=self.xyz[self.idx[self.c]]
        # z 镜像映射：z -> 2*z_mid - z（格点上）
        self.zmid = self.ctr[2]
        self.mirror = np.zeros(self.N, dtype=int)
        for i,(x,y,z) in enumerate(self.coords):
            z2 = int(round(2*self.ctr[2]/a - z))
            self.mirror[i] = self.idx[(x,y,z2)]
    def Vnuc(self, ps):
        V=np.zeros(self.N)
        for p in ps:
            r=np.linalg.norm(self.xyz-p,axis=1)
            on=r<1e-9; r=np.where(on,self.a,r)
            V-=self.A/r; V[on]=-self.A*ONSITE/self.a
        return V
    def conv_all(self, dens, CH=1200):
        Psi={k:np.zeros(self.N) for k in dens}
        for st in range(0,self.N,CH):
            en=min(st+CH,self.N)
            D=np.linalg.norm(self.xyz[st:en,None,:]-self.xyz[None,:,:],axis=2)
            on=D<1e-9; D=np.where(on,self.a,D)
            blk=self.A/D; blk[on]=self.A*ONSITE/self.a
            for k,rho in dens.items(): Psi[k][st:en]=blk@rho
        return Psi

def analyze(box, evec, ps, Rloc=3.0):
    """宇称 + 局域性。"""
    out=[]
    for j in range(evec.shape[1]):
        v = evec[:,j]
        # 宇称：P v = v[mirror]
        p = float(v @ v[box.mirror])
        # 局域性：核附近 |v|^2 占比
        w = v**2
        near = np.zeros(box.N, dtype=bool)
        for q in ps:
            near |= (np.linalg.norm(box.xyz-q, axis=1) < Rloc)
        loc = float(w[near].sum())
        out.append({"parity":p, "loc":loc})
    return out

def main():
    print("=" * 84)
    print("zero_h2_basis.py — 按对称性+局域性选基组")
    print("=" * 84)
    a=0.5; Lphys=14.0; R=1.4
    box=Box(a,Lphys)
    pA=box.ctr-np.array([0,0,R/2]); pB=box.ctr+np.array([0,0,R/2])
    H1=(0.5/a**2)*box.L+sp.diags(box.Vnuc([pA,pB]))
    K=10
    ev,evec=spl.eigsh(H1,k=K,which="SA")
    o=np.argsort(ev); ev=ev[o]; evec=evec[:,o]
    an=analyze(box,evec,[pA,pB])
    print("  候选轨道（%d 个）：" % K)
    print("   j    eps        parity    loc(核附近占比)")
    for j in range(K):
        print("   %d  %+.6f   %+.3f     %.3f" % (j, ev[j], an[j]["parity"], an[j]["loc"]))
    # 选择：|parity| > 0.5（宇称明确）且 loc > 0.5（局域）
    sel=[j for j in range(K) if abs(an[j]["parity"])>0.5 and an[j]["loc"]>0.5]
    print("  选中（|parity|>0.5 且 loc>0.5）：", sel)
    ck("Z1 选中数 >= 2", len(sel)>=2, f"{sel}")
    ck("Z2 选中数 <= 6", len(sel)<=6, f"{len(sel)}")
    out={"a":a,"Lphys":Lphys,"R":R,"K":K,"eps":ev.tolist(),
         "parity":[an[j]["parity"] for j in range(K)],
         "loc":[an[j]["loc"] for j in range(K)],"selected":sel,"checks":PASS}
    with open(os.path.join(RES,"zero_h2_basis.json"),"w") as f:
        json.dump(out,f,ensure_ascii=False,indent=2)
    print("\n%d 项，%d PASS"%(len(PASS),sum(1 for p in PASS if p["pass"])))
    print("-> "+os.path.join(RES,"zero_h2_basis.json"))

if __name__=="__main__":
    main()
