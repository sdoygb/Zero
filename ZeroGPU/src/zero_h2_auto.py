#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2_auto.py — H2 的 E(R)，基组自动按宇称+局域性筛（不依赖固定 M）。

对每个 (a,L,R)：取最低 K=8 个本征态，按 |宇称|>0.5 且 局域占比>0.5 筛出物理轨道，
再用它们做 CI。这修掉 §4.15/§4.16 的"最低 M 个混入盒态"问题。
输出：results/zero_h2_auto.json
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
    idx={};coords=[]
    for x in range(n):
        for y in range(n):
            for z in range(n):
                idx[(x,y,z)]=len(coords);coords.append((x,y,z))
    N=len(coords);adj=[[] for _ in range(N)]
    for (x,y,z),i in idx.items():
        for dx,dy,dz in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)):
            j=idx.get((x+dx,y+dy,z+dz))
            if j is not None: adj[i].append(j)
    rows,cols,vals=[],[],[]
    for i,nb in enumerate(adj):
        rows.append(i);cols.append(i);vals.append(float(len(nb)))
        for j in nb: rows.append(i);cols.append(j);vals.append(-1.0)
    return sp.csr_matrix((vals,(rows,cols)),shape=(N,N)),idx,coords

class Box:
    def __init__(self,a,Lphys):
        self.a=a;self.Lphys=Lphys;self.A=1.0
        self.n=int(round(Lphys/a))+1
        self.L,self.idx,self.coords=lattice(self.n)
        self.N=self.L.shape[0]
        self.xyz=np.array(self.coords,float)*a
        self.c=(self.n//2,self.n//2,self.n//2)
        self.ctr=self.xyz[self.idx[self.c]]
        self.mirror=np.zeros(self.N,dtype=int)
        for i,(x,y,z) in enumerate(self.coords):
            z2=(self.n-1)-z
            self.mirror[i]=self.idx[(x,y,z2)]
    def Vnuc(self,ps):
        """核势：离格核用 r_min=a/2 软化（实测把 H 原子从 -0.642 修到 -0.492）。"""
        V=np.zeros(self.N)
        for p in ps:
            r=np.linalg.norm(self.xyz-p,axis=1)
            r=np.maximum(r, self.a/2.0)
            V-=self.A/r
        return V
    def conv_all(self,dens,CH=1200):
        Psi={k:np.zeros(self.N) for k in dens}
        for st in range(0,self.N,CH):
            en=min(st+CH,self.N)
            D=np.linalg.norm(self.xyz[st:en,None,:]-self.xyz[None,:,:],axis=2)
            on=D<1e-9;D=np.where(on,self.a,D)
            blk=self.A/D;blk[on]=self.A*ONSITE/self.a
            for k,rho in dens.items(): Psi[k][st:en]=blk@rho
        return Psi

def select_orbitals(box,ps,Rloc=3.0,K=8,par_tol=0.5,loc_tol=0.5):
    H1=(0.5/box.a**2)*box.L+sp.diags(box.Vnuc(ps))
    ev,evec=spl.eigsh(H1,k=K,which="SA")
    o=np.argsort(ev);ev=ev[o];evec=evec[:,o]
    sel=[]
    info=[]
    near=np.zeros(box.N,dtype=bool)
    for q in ps:
        near |= (np.linalg.norm(box.xyz-q,axis=1)<Rloc)
    for j in range(K):
        v=evec[:,j]
        par=float(v@v[box.mirror])
        loc=float((v**2)[near].sum())
        info.append({"j":j,"eps":float(ev[j]),"parity":par,"loc":loc})
        if abs(par)>par_tol and loc>loc_tol:
            sel.append(j)
    return ev,evec,sel,info

def h2_ci(box,R):
    a=box.a
    pA=box.ctr-np.array([0,0,R/2]);pB=box.ctr+np.array([0,0,R/2])
    ev,evec,sel,info=select_orbitals(box,[pA,pB])
    M=len(sel)
    if M<2: return None,sel,info
    E=evec[:,sel]
    # 重新正交化（eigh 已正交）
    dens={}
    for qi in range(M):
        for si in range(qi,M): dens[(qi,si)]=E[:,qi]*E[:,si]
    Psi=box.conv_all(dens)
    def eri(p,q,r,s):
        key=(r,s) if r<=s else (s,r)
        return float(E[:,p]@(Psi[key]*E[:,q]))
    SO=[(p,s) for p in range(M) for s in (0,1)]
    def hij(i,j):
        (p,si),(q,sj)=SO[i],SO[j]
        return ev[sel[p]] if (p==q and si==sj) else 0.0
    def vijkl(i,j,k,l):
        (p,si),(q,sj),(r,sk),(s,sl)=SO[i],SO[j],SO[k],SO[l]
        return eri(p,r,q,s) if (si==sk and sj==sl) else 0.0
    dets=list(itertools.combinations(range(2*M),2));nd=len(dets)
    H=np.zeros((nd,nd))
    for I,D1 in enumerate(dets):
        for J,D2 in enumerate(dets):
            only1=[x for x in D1 if x not in D2];only2=[x for x in D2 if x not in D1]
            common=[x for x in D1 if x in D2]
            if not only1 and not only2:
                v=sum(hij(i,i) for i in D1)
                for (i,j) in itertools.combinations(D1,2): v+=vijkl(i,j,i,j)-vijkl(i,j,j,i)
                H[I,J]=v
            elif len(only1)==1 and len(only2)==1:
                i,av=only1[0],only2[0];v=hij(i,av)
                for j in common: v+=vijkl(av,j,i,j)-vijkl(av,j,j,i)
                H[I,J]=(-1)**(D1.index(i)+D2.index(av))*v
            elif len(only1)==2 and len(only2)==2:
                i,j=only1;av,bv=only2
                v=vijkl(av,bv,i,j)-vijkl(av,bv,j,i)
                H[I,J]=(-1)**(D1.index(i)+D1.index(j)+D2.index(av)+D2.index(bv))*v
    w,_=np.linalg.eigh(H)
    return float(w[0]),sel,info

def main():
    print("="*84)
    print("zero_h2_auto.py — H2 E(R)，基组自动筛选")
    print("="*84)
    exact={1.4:-1.174470,1.6:-1.117510,2.0:-1.067800,2.5:-1.031550}
    rows=[]
    for (a,Lp) in [(0.5,14.0),(0.4,14.0)]:
        box=Box(a,Lp)
        H1H=(0.5/a**2)*box.L+sp.diags(box.Vnuc([box.ctr]))
        eH=float(spl.eigsh(H1H,k=1,which="SA")[0][0]);E2a=2*eH
        print("  a=%.2f L=%.0f N=%d  H 原子=%+.6f"%(a,Lp,box.N,eH))
        for R in (1.4,2.0):
            E0,sel,info=h2_ci(box,R)
            if E0 is None:
                print("    R=%.1f  选中<2，跳过"%R); continue
            Etot=E0+1.0/R;bind=Etot-E2a;ex=exact[R];exb=ex+1.0
            rows.append({"a":a,"L":Lp,"N":box.N,"R":R,"nsel":len(sel),"sel":sel,
                         "E":Etot,"err":Etot-ex,"bind":bind,"bind_err":bind-exb})
            print("    R=%.1f 选中%d个%s  E=%+.6f (err %+.5f)  结合=%+.6f (err %+.5f)"%(
                R,len(sel),sel,Etot,Etot-ex,bind,bind-exb))
    ck("A1 自动筛选给出 >=2 个轨道", all(r["nsel"]>=2 for r in rows), str([r["nsel"] for r in rows]))
    ck("A2 绝对能量误差 < 0.05", all(abs(r["err"])<0.05 for r in rows),
       f"max={max(abs(r['err']) for r in rows):.5f}")
    ck("A3 结合能误差 < 0.05", all(abs(r["bind_err"])<0.05 for r in rows),
       f"max={max(abs(r['bind_err']) for r in rows):.5f}")
    ck("A4 结合能误差 < 1.6e-3", all(abs(r["bind_err"])<1.6e-3 for r in rows),
       f"max={max(abs(r['bind_err']) for r in rows):.5f}")
    out={"rows":rows,"checks":PASS}
    with open(os.path.join(RES,"zero_h2_auto.json"),"w") as f:
        json.dump(out,f,ensure_ascii=False,indent=2)
    print("\n%d 项，%d PASS"%(len(PASS),sum(1 for p in PASS if p["pass"])))
    print("-> "+os.path.join(RES,"zero_h2_auto.json"))

if __name__=="__main__":
    main()
