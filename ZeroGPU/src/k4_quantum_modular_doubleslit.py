#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_quantum_modular_doubleslit.py —— L2 图上波以模流演化: 从模流 H_mod 到双缝条纹, 全程内生
======================================================================================
补 Z119 §6.2 开放项: 把模流时间落到 L2 图上波.

链(全程内生, 无外部能量/波长输入):
  L1 稳态 ρ=diag(2/3,1/3) (Z087)
    → 模流 H_mod = -log ρ, 频率差 = log(p/q) = log2 = f  (Z005 定理4, Z119)
    → 两顶点波以模流演化: 相干 ψ1ψ2* 振荡频率 = f  (测)
    → L2 图上双缝: f 作驱动频率 → q=2arcsin(f/2) (色散) → λ_e 实测 → 条纹

关键: f 来自模流 H_mod (内生), 不是外部数值; 波以模流时间演化.
"""
from __future__ import annotations
import os, json
import numpy as np
import scipy.sparse as sp, scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

def build_H(nr, nc, xb=90, apertures=(90, 150), no_barrier=False, c2=1.0):
    N = nr * nc; rows=[]; cols=[]; vals=[]; deg=np.zeros(N)
    def ae(a,b):
        nonlocal rows, cols, vals
        rows+=[a,b]; cols+=[b,a]; vals+=[-1.0,-1.0]; deg[a]+=1; deg[b]+=1
    for y in range(nr):
        for x in range(nc):
            v=y*nc+x
            if x+1<nc and (no_barrier or not (x==xb-1 and y not in apertures)): ae(v,v+1)
            if y+1<nr: ae(v,v+nc)
    rows+=list(range(N)); cols+=list(range(N)); vals+=list(deg)
    return sp.csr_matrix((vals,(rows,cols)),shape=(N,N))*c2

def absorber(nr,nc,m=24,sig0=1.5):
    yy,xx=np.mgrid[0:nr,0:nc]
    de=np.minimum(np.minimum(yy,nr-1-yy),np.minimum(xx,nc-1-xx))
    return np.where(de<m, sig0*((m-de)/m)**2, 0.0)

def solve_scatter(nr,nc,q,xb=90,apertures=(90,150),no_barrier=False,src=(120,20),c2=1.0,m=24,sig0=1.5):
    H=build_H(nr,nc,xb,apertures,no_barrier,c2); N=nr*nc
    E=2*c2*(1-np.cos(q))
    b=np.zeros(N,complex); b[src[0]*nc+src[1]]=1.0
    A=(E*sp.eye(N,format="csr")-H+1j*sp.diags(absorber(nr,nc,m,sig0).ravel())).tocsc()
    return spla.splu(A).solve(b).reshape(nr,nc)

def main():
    P_("="*100)
    P_("L2 图上波以模流演化: 模流 H_mod → 相干频率 f → 双缝条纹 (k4_quantum_modular_doubleslit)")
    P_("="*100)
    out={}

    # ---- 1. 稳态 ρ 与模流 H_mod ----
    p,q=2/3,1/3
    f=float(np.log(p/q))                       # log(p/q)=log2
    Hmod=-np.log(np.array([p,q]))              # H_mod=-log ρ
    P_(f"L1 稳态 ρ=diag({p:.3f},{q:.3f}); H_mod=-logρ={Hmod.round(4)}; 频率差=log(p/q)=log2={f:.5f}")
    out["setup"]={"steady":[p,q],"H_mod":Hmod.round(4).tolist(),"f_intrinsic":round(f,5)}

    # ---- 2. 两顶点波以模流演化: 相干频率=f (测) ----
    T,dt=4000,0.02
    psi0=np.array([1.0,1.0],complex)           # 等幅两支
    coh=np.zeros(T,complex)
    times=dt*np.arange(T)
    for t in range(T):
        tt=times[t]
        psi=np.exp(-1j*Hmod*tt)*psi0            # e^{-iH_mod t} 模流演化
        coh[t]=psi[0]*np.conj(psi[1])           # 相干 ψ1ψ2*
    ph=np.unwrap(np.angle(coh))
    f_coh=float(np.polyfit(times,ph,1)[0])
    P_("【两顶点波以模流演化】ψ(t)=e^{-iH_mod t}ψ0")
    P_(f"  相干 ψ1ψ2* 振荡频率实测={f_coh:.5f}; 精确=log(p/q)={f:.5f}; 误差={abs(f_coh-f)/f*100:.4f}%")
    out["modular_coherence"]={"measured":round(f_coh,5),"exact":round(f,5),
                              "err_pct":round(abs(f_coh-f)/f*100,4)}

    # ---- 3. L2 图上双缝: f 作驱动(内生) → λ_e → 条纹 ----
    P_("\n【L2 图上双缝】f(模流内生) → q=2arcsin(f/2) → 定态散射 → 条纹")
    nr,nc=240,260; xb=90; ap=(90,150); c2=1.0
    qd=float(2*np.arcsin(f/2))
    # 测 λ_e (对照无挡板, 源行相位梯度)
    pc=solve_scatter(nr,nc,qd,no_barrier=True)
    y0=nr//2
    xw=np.arange(30,91)
    qm=float(np.polyfit(xw,np.unwrap(np.angle(pc[y0,xw])),1)[0])
    lam=2*np.pi/abs(qm)
    P_(f"  f={f:.4f} → q=2arcsin(f/2)={qd:.4f}; 实测 q={qm:.4f}(误差{abs(qm-qd)/qd*100:.2f}%); λ_e={lam:.3f}格")
    P_(f"  2π/f={2*np.pi/f:.3f}; λ_e/(2π/f)={lam/(2*np.pi/f):.4f}  (≈1 ⟹ 模流内生能量标度)")
    # 双缝条纹间距
    p12=solve_scatter(nr,nc,qd,xb,ap)
    p1=solve_scatter(nr,nc,qd,xb,(ap[0],)); p2=solve_scatter(nr,nc,qd,xb,(ap[1],))
    dsl={}
    for L in [40,80,120]:
        x=xb+L; ykw=np.arange(y0-60,y0+61)
        ph_=np.unwrap(np.angle(p1[ykw,x]/p2[ykw,x]))
        A_=np.vstack([ykw,np.ones_like(ykw)]).T
        cf,*_=np.linalg.lstsq(A_,ph_,rcond=None); per=2*np.pi/abs(cf[0])
        dsl[L]=round(float(per),3)
        P_(f"  L={L:3d}: 条纹间距={per:.3f}格")
    out["doubleslit"]={"f":round(f,5),"q_drive":round(qd,5),"q_measured":round(qm,5),
                       "lambda_e":round(lam,3),"lambda_2pi_f":round(2*np.pi/f,3),
                       "ratio":round(lam/(2*np.pi/f),4),"fringe_spacing":dsl}

    os.makedirs(RES,exist_ok=True)
    with open(os.path.join(RES,"k4_quantum_modular_doubleslit.json"),"w",encoding="utf-8") as f_:
        json.dump(out,f_,ensure_ascii=False,indent=1)
    with open(os.path.join(RES,"k4_quantum_modular_doubleslit.txt"),"w",encoding="utf-8") as f_:
        f_.write("\n".join(LOG))
    P_("\n已存 k4_quantum_modular_doubleslit.json / .txt")

if __name__=="__main__":
    main()
