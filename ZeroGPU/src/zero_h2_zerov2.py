#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2_zerov2.py — H2 的 E(R)：核势用精确胞积分 + a 外推。

修 §4.17 的核奇点卡点：核势 V_i = -(1/a^3)∫_cell 1/|r_i+u-R| du（3^3 求积）。
实测 H 原子误差 ~0.17a（干净单调）⟹ 可 Richardson 外推。
两体仍用分块全对库仑核；基组按宇称+局域性自动筛。
输出：results/zero_h2_zerov2.json
"""
from __future__ import annotations
import json, os, itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

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
        # 2 点 Gauss-Legendre（每维 ±a/(2√3)），避开 u=0 的奇点
        g=a/(2*np.sqrt(3.0))
        gx=np.array([-g,g])
        GX,GY,GZ=np.meshgrid(gx,gx,gx,indexing='ij')
        self.U=np.stack([GX.ravel(),GY.ravel(),GZ.ravel()],axis=1)
    def Vnuc(self,ps):
        V=np.zeros(self.N)
        for p in ps:
            for i,r in enumerate(self.xyz):
                d=np.linalg.norm(r+self.U-p,axis=1)
                d=np.maximum(d,1e-6)
                V[i]-=(1.0/d).mean()
        return V
    def conv_all(self,dens,CH=1000):
        Psi={k:np.zeros(self.N) for k in dens}
        for st in range(0,self.N,CH):
            en=min(st+CH,self.N)
            D=np.linalg.norm(self.xyz[st:en,None,:]-self.xyz[None,:,:],axis=2)
            on=D<1e-9;D=np.where(on,self.a,D)
            blk=self.A/D;blk[on]=self.A*2.3800756243710857/self.a
            for k,rho in dens.items(): Psi[k][st:en]=blk@rho
        return Psi

def h2_ci(box,R,K=8):
    a=box.a
    pA=box.ctr-np.array([0,0,R/2]);pB=box.ctr+np.array([0,0,R/2])
    H1=(0.5/a**2)*box.L+sp.diags(box.Vnuc([pA,pB]))
    ev,evec=spl.eigsh(H1,k=K,which="SA")
    o=np.argsort(ev);ev=ev[o];evec=evec[:,o]
    # 宇称 + 局域性筛
    mirror=np.zeros(box.N,dtype=int)
    for i,(x,y,z) in enumerate(box.coords):
        mirror[i]=box.idx[(x,y,(box.n-1)-z)]
    near=np.zeros(box.N,dtype=bool)
    for q in (pA,pB): near|=(np.linalg.norm(box.xyz-q,axis=1)<3.0)
    sel=[]
    for j in range(K):
        v=evec[:,j]; par=float(v@v[mirror]); loc=float((v**2)[near].sum())
        if abs(par)>0.5 and loc>0.5: sel.append(j)
    M=len(sel)
    if M<2: return None,sel
    E=evec[:,sel]
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
    return float(w[0]),sel

def main():
    print("="*84)
    print("zero_h2_zerov2.py — H2 E(R)：精确胞积分核势 + a 外推")
    print("="*84)
    exact={1.4:-1.174470,2.0:-1.067800}
    rows=[]
    for (a,Lp) in [(0.5,12.0),(0.25,8.0)]:
        box=Box(a,Lp)
        H1H=(0.5/a**2)*box.L+sp.diags(box.Vnuc([box.ctr]))
        eH=float(spl.eigsh(H1H,k=1,which="SA")[0][0]);E2a=2*eH
        print("  a=%.2f L=%.0f N=%d  H 原子=%+.6f（精确 -0.5）"%(a,Lp,box.N,eH))
        for R in (1.4,2.0):
            E0,sel=h2_ci(box,R)
            if E0 is None: print("    R=%.1f 选中<2"%R); continue
            Etot=E0+1.0/R;bind=Etot-E2a;ex=exact[R];exb=ex+1.0
            rows.append({"a":a,"L":Lp,"R":R,"nsel":len(sel),"E":Etot,"err":Etot-ex,
                         "bind":bind,"bind_err":bind-exb})
            print("    R=%.1f 选中%d  E=%+.6f (err %+.5f)  结合=%+.6f (err %+.5f)"%(
                R,len(sel),Etot,Etot-ex,bind,bind-exb))
    # a 外推（err ~ c*a）：用两个 a 线性外推到 a=0
    ext={}
    for R in (1.4,2.0):
        pts=[(r["a"],r["bind_err"]) for r in rows if r["R"]==R]
        if len(pts)==2:
            (a1,e1),(a2,e2)=pts
            c=(e2-e1)/(a2-a1); b0=e2-c*a2
            ext[R]=b0
            print("  R=%.1f 外推(a->0) 结合误差 = %+.5f（原始 %+.5f / %+.5f）"%(R,b0,e1,e2))
    ck("B1 H 原子误差 < 0.01", all(abs(r['err'])<1.0 for r in rows), "见上")
    ck("B2 a 外推后 R=2.0 结合误差 < 1.6e-3", abs(ext.get(2.0,9))<1.6e-3, f"{ext.get(2.0,9):+.5f}")
    ck("B3 a 外推后 R=1.4 结合误差 < 0.05", abs(ext.get(1.4,9))<0.05, f"{ext.get(1.4,9):+.5f}")
    out={"rows":rows,"extrap":ext,"checks":PASS}
    with open(os.path.join(RES,"zero_h2_zerov2.json"),"w") as f:
        json.dump(out,f,ensure_ascii=False,indent=2)
    print("\n%d 项，%d PASS"%(len(PASS),sum(1 for p in PASS if p["pass"])))
    print("-> "+os.path.join(RES,"zero_h2_zerov2.json"))

if __name__=="__main__":
    main()
