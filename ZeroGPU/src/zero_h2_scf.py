#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2_scf.py — H2 的 Zero 原生自洽（轨道能量极小化，不用 Fock/RDM）。

思路：把占据轨道写成物理候选轨道的线性组合 ψ = Σ c_p φ_p，
直接极小化 E(c) = 2 h_ψψ + (ψψ|ψψ) + V_nn（闭壳层、两电子）。
再在优化轨道 + 其余候选轨道上做 CI（Slater–Condon）。
对照：H2 R=1.4 精确 -1.174470；R=2.0 精确 -1.067800
输出：results/zero_h2_scf.json
"""
from __future__ import annotations
import json, os, itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.optimize import minimize

HERE=os.path.dirname(os.path.abspath(__file__))
RES=os.path.normpath(os.path.join(HERE,"..","results"))
PASS=[]

def ck(name,cond,detail=""):
    PASS.append({"name":name,"pass":bool(cond),"detail":detail})
    print("[%s] %s %s"%("PASS" if cond else "FAIL",name,detail))

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
        g=a/(2*np.sqrt(3.0)); gx=np.array([-g,g])
        GX,GY,GZ=np.meshgrid(gx,gx,gx,indexing='ij')
        self.U=np.stack([GX.ravel(),GY.ravel(),GZ.ravel()],axis=1)
    def Vnuc(self,ps):
        V=np.zeros(self.N)
        for p in ps:
            for i,r in enumerate(self.xyz):
                d=np.linalg.norm(r+self.U-p,axis=1); d=np.maximum(d,1e-6)
                V[i]-=(1.0/d).mean()
        return V
    def conv(self,rho,CH=1000):
        out=np.zeros(self.N)
        for st in range(0,self.N,CH):
            en=min(st+CH,self.N)
            D=np.linalg.norm(self.xyz[st:en,None,:]-self.xyz[None,:,:],axis=2)
            on=D<1e-9;D=np.where(on,self.a,D)
            blk=self.A/D;blk[on]=self.A*2.3800756243710857/self.a
            out[st:en]=blk@rho
        return out

def h2(box,R,M=6,K=10):
    a=box.a
    pA=box.ctr-np.array([0,0,R/2]);pB=box.ctr+np.array([0,0,R/2])
    H1=(0.5/a**2)*box.L+sp.diags(box.Vnuc([pA,pB]))
    ev,evec=spl.eigsh(H1,k=K,which="SA")
    o=np.argsort(ev);ev=ev[o];evec=evec[:,o]
    B=evec[:,:M]                      # 候选轨道（N x M）
    # 两体在候选基里：(pq|rs)
    def eri(p,q,r,s):
        rho=B[:,q]*B[:,s]
        Psi=box.conv(rho)
        return float(B[:,p]@(Psi*B[:,r]))
    def E_closed(c):
        c=c/np.linalg.norm(c)
        psi=B@c
        h=c@(ev[:M]*c)
        rho=psi*psi
        Psi=box.conv(rho)
        J=float(psi@(Psi*psi))
        return 2*h+J
    # 初值：最低轨道
    c0=np.zeros(M); c0[0]=1.0
    res=minimize(E_closed, c0, method="Nelder-Mead",
                 options={"maxiter":4000,"xatol":1e-8,"fatol":1e-10})
    c=res.x/np.linalg.norm(res.x)
    E_hf=E_closed(c)+1.0/R
    # 把优化后的占据轨道与其余候选正交化，得到 CI 基
    psi=B@c
    Q,_=np.linalg.qr(np.column_stack([psi,B]))
    Eci=Q[:,:M]
    # CI（自旋轨道 Slater–Condon）
    def eri2(p,q,r,s):
        rho=Eci[:,q]*Eci[:,s]; Psi=box.conv(rho)
        return float(Eci[:,p]@(Psi*Eci[:,r]))
    h1=[float(Eci[:,p]@( (0.5/a**2)*box.L @ Eci[:,p] ) + Eci[:,p]@(box.Vnuc([pA,pB])*Eci[:,p])) for p in range(M)]
    SO=[(p,s) for p in range(M) for s in (0,1)]
    def hij(i,j):
        (p,si),(q,sj)=SO[i],SO[j]
        return h1[p] if (p==q and si==sj) else 0.0
    def vijkl(i,j,k,l):
        (p,si),(q,sj),(r,sk),(s,sl)=SO[i],SO[j],SO[k],SO[l]
        return eri2(p,r,q,s) if (si==sk and sj==sl) else 0.0
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
    E_ci=float(w[0])+1.0/R
    return E_hf, E_ci, c

def main():
    print("="*84)
    print("zero_h2_scf.py — H2 的 Zero 原生自洽（轨道极小化）")
    print("="*84)
    exact={1.4:-1.174470,2.0:-1.067800}
    a=0.5; Lp=12.0
    box=Box(a,Lp)
    H1H=(0.5/a**2)*box.L+sp.diags(box.Vnuc([box.ctr]))
    eH=float(spl.eigsh(H1H,k=1,which="SA")[0][0]);E2a=2*eH
    print("  a=%.2f L=%.0f N=%d  H 原子=%+.6f"%(a,Lp,box.N,eH))
    rows=[]
    for R in (1.4,2.0):
        Ehf,Eci,c=h2(box,R)
        b_hf=Ehf-E2a; b_ci=Eci-E2a
        ex=exact[R]; exb=ex+1.0
        rows.append({"a":a,"L":Lp,"R":R,"E_hf":Ehf,"E_ci":Eci,
                     "err_hf":Ehf-ex,"err_ci":Eci-ex,
                     "bind_ci":b_ci,"bind_err":b_ci-exb})
        print("    R=%.1f  E_HF=%+.6f (err %+.5f)  E_CI=%+.6f (err %+.5f)  结合=%+.6f (err %+.5f)"%(
            R,Ehf,Ehf-ex,Eci,Eci-ex,b_ci,b_ci-exb))
    ck("D1 CI 绝对能量误差 < 0.05", all(abs(r['err_ci'])<0.05 for r in rows),
       f"max={max(abs(r['err_ci']) for r in rows):.5f}")
    ck("D2 CI 结合能误差 < 0.05", all(abs(r['bind_err'])<0.05 for r in rows),
       f"max={max(abs(r['bind_err']) for r in rows):.5f}")
    ck("D3 CI 结合能误差 < 1.6e-3", all(abs(r['bind_err'])<1.6e-3 for r in rows),
       f"max={max(abs(r['bind_err']) for r in rows):.5f}")
    out={"rows":rows,"checks":PASS}
    with open(os.path.join(RES,"zero_h2_scf.json"),"w") as f:
        json.dump(out,f,ensure_ascii=False,indent=2)
    print("\n%d 项，%d PASS"%(len(PASS),sum(1 for p in PASS if p["pass"])))
    print("-> "+os.path.join(RES,"zero_h2_scf.json"))

if __name__=="__main__":
    main()
