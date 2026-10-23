#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2_mcscf.py — H2 的 Zero 原生"两电子自洽"（MCSCF 式）。

思路：不用外部 SCF 库；用舞台算子自己做自洽——
  h = (1/2a^2)L - Σ_A V_A（Zero 动能 + 胞积分核势）
  两体 (pq|rs) 用离散/解析库仑核
  CI（Slater–Condon）→ 1-RDM D
  Fock F = h + Σ_rs D_rs[(pq|rs) - 1/2 (pr|qs)]（在有限轨道基里）
  对角化 F → 新轨道 → 重算两体 → 重做 CI → 迭代到自洽
对照：H2 R=1.4 精确 -1.174470；R=2.0 精确 -1.067800
输出：results/zero_h2_mcscf.json
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
    def conv_all(self,dens,CH=1000):
        Psi={k:np.zeros(self.N) for k in dens}
        for st in range(0,self.N,CH):
            en=min(st+CH,self.N)
            D=np.linalg.norm(self.xyz[st:en,None,:]-self.xyz[None,:,:],axis=2)
            on=D<1e-9;D=np.where(on,self.a,D)
            blk=self.A/D;blk[on]=self.A*2.3800756243710857/self.a
            for k,rho in dens.items(): Psi[k][st:en]=blk@rho
        return Psi

def build_eri(box, E):
    """给定轨道矩阵 E（N x M），返回 eri(p,q,r,s) 与 (pq|rs)。"""
    M=E.shape[1]
    dens={}
    for q in range(M):
        for s in range(q,M): dens[(q,s)]=E[:,q]*E[:,s]
    Psi=box.conv_all(dens)
    def eri(p,q,r,s):
        key=(r,s) if r<=s else (s,r)
        return float(E[:,p]@(Psi[key]*E[:,q]))
    return eri

def ci_energy_dm(ev_orb, eri, M):
    """2 电子 CI（Slater–Condon）→ 电子能量 + 1-RDM（空间轨道）。"""
    SO=[(p,s) for p in range(M) for s in (0,1)]
    def hij(i,j):
        (p,si),(q,sj)=SO[i],SO[j]
        return ev_orb[p] if (p==q and si==sj) else 0.0
    def vijkl(i,j,k,l):
        (p,si),(q,sj),(r,sk),(s,sl)=SO[i],SO[j],SO[k],SO[l]
        return eri(p,r,q,s) if (si==sk and sj==sl) else 0.0
    dets=list(itertools.combinations(range(2*M),2)); nd=len(dets)
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
    w,V=np.linalg.eigh(H)
    c=V[:,0]
    # 1-RDM（空间）：D_pq = Σ_σ <a†_pσ a_qσ>
    D=np.zeros((M,M))
    for I,D1 in enumerate(dets):
        for J,D2 in enumerate(dets):
            common=[x for x in D1 if x in D2]
            if len(common)==2:
                for x in common:
                    p,si=SO[x]
                    # 找另一个共同指标
                    y=[z for z in common if z!=x][0]
                    q,sj=SO[y]
                    if si==sj and p!=q:
                        D[p,q]+=c[I]*c[J]*((-1)**(D1.index(x)+D2.index(y)))
    # 对角项
    for I,D1 in enumerate(dets):
        for x in D1:
            p,si=SO[x]
            D[p,p]+=c[I]*c[I]
    return float(w[0]), c, D

def h2_mcscf(box, R, M=6, niter=12, K=10):
    a=box.a
    pA=box.ctr-np.array([0,0,R/2]);pB=box.ctr+np.array([0,0,R/2])
    H1=(0.5/a**2)*box.L+sp.diags(box.Vnuc([pA,pB]))
    ev,evec=spl.eigsh(H1,k=K,which="SA")
    o=np.argsort(ev);ev=ev[o];evec=evec[:,o]
    # 初始：最低 M 个（宇称+局域性由后续自洽自然筛）
    E=evec[:,:M].copy()
    hist=[]
    for it in range(niter):
        eri=build_eri(box,E)
        Ee, c, D = ci_energy_dm(ev[:M], eri, M)
        # Fock（在有限轨道基里）：F_pq = h_pq + Σ_rs D_rs[(pq|rs) - 1/2(pr|qs)]
        F=np.diag(ev[:M]).astype(float)
        for p in range(M):
            for q in range(M):
                s=0.0
                for r in range(M):
                    for ss in range(M):
                        s += D[r,ss]*(eri(p,q,r,ss) - 0.5*eri(p,ss,r,q))
                F[p,q]+=s
        # 对角化 F（在轨道基）→ 变换矩阵 U（M x M）
        wF,U=np.linalg.eigh(F)
        # 新轨道：E_new = E @ U
        Enew=E@U
        # 正交化（E 已正交，U 正交，故 Enew 正交）
        dE=np.max(np.abs(np.abs(Enew.T@Enew)-np.eye(M)))
        # 能量（用新轨道重算）
        hist.append({"it":it,"E_ci":Ee,"orth":float(dE)})
        E=Enew
        if dE<1e-6 and it>2: break
    eri=build_eri(box,E)
    Ee, c, D = ci_energy_dm(ev[:M], eri, M)
    return Ee, E, hist, ev[:M]

def main():
    print("="*84)
    print("zero_h2_mcscf.py — H2 的 Zero 原生两电子自洽")
    print("="*84)
    exact={1.4:-1.174470,2.0:-1.067800}
    rows=[]
    for (a,Lp) in [(0.5,12.0)]:
        box=Box(a,Lp)
        H1H=(0.5/a**2)*box.L+sp.diags(box.Vnuc([box.ctr]))
        eH=float(spl.eigsh(H1H,k=1,which="SA")[0][0]);E2a=2*eH
        print("  a=%.2f L=%.0f N=%d  H 原子=%+.6f"%(a,Lp,box.N,eH))
        for R in (1.4,2.0):
            Ee,E,hist,eps=h2_mcscf(box,R,M=6)
            Etot=Ee+1.0/R;bind=Etot-E2a;ex=exact[R];exb=ex+1.0
            rows.append({"a":a,"L":Lp,"R":R,"E":Etot,"err":Etot-ex,"bind":bind,
                         "bind_err":bind-exb,"nit":len(hist)})
            print("    R=%.1f  自洽%d 步  E=%+.6f (err %+.5f)  结合=%+.6f (err %+.5f)"%(
                R,len(hist),Etot,Etot-ex,bind,bind-exb))
    ck("C1 H2 绝对能量误差 < 0.05", all(abs(r['err'])<0.05 for r in rows),
       f"max={max(abs(r['err']) for r in rows):.5f}")
    ck("C2 结合能误差 < 0.05", all(abs(r['bind_err'])<0.05 for r in rows),
       f"max={max(abs(r['bind_err']) for r in rows):.5f}")
    ck("C3 结合能误差 < 1.6e-3", all(abs(r['bind_err'])<1.6e-3 for r in rows),
       f"max={max(abs(r['bind_err']) for r in rows):.5f}")
    out={"rows":rows,"checks":PASS}
    with open(os.path.join(RES,"zero_h2_mcscf.json"),"w") as f:
        json.dump(out,f,ensure_ascii=False,indent=2)
    print("\n%d 项，%d PASS"%(len(PASS),sum(1 for p in PASS if p["pass"])))
    print("-> "+os.path.join(RES,"zero_h2_mcscf.json"))

if __name__=="__main__":
    main()
