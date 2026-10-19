#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_chain.py —— 4→18→6 接起来 + 敏感性 + T=4 继承
================================================================
三步：
  A. 4 势 → 6 边（差）→ x*（质量目标）→ 6 层动力学收敛；
  B. 敏感性：γ,κ,A,T 变化下 6 层是否稳定；
  C. T=4 继承：用 18 引擎的 shadow 6 作驱动（不用外部 sin），看 6 层是否带 T=4。

输出：results/zero_world6_chain.json、zero_world6_chain.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_chain.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
RES=os.path.join(HERE,"..","results")
EDGES=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
A_line=np.zeros((6,6))
for i,(u,v) in enumerate(EDGES):
    for j,(a,b) in enumerate(EDGES):
        if i!=j and len({u,v}&{a,b})>0: A_line[i,j]=1
L_line=np.diag(A_line.sum(1))-A_line
W=np.array([0.111,0.531,0.0,0.358]); P=-8.33
Q=np.array([2.16,4.67,93.4,1270.,4180.,172760.])
def target_from_w(W=W,P=P,C=0.0):
    e=np.array([abs(W[u]-W[v]) for u,v in EDGES])
    m=np.exp(C+P*np.log(e)); m=np.sort(m)
    return m/m.sum()
def run6(xstar, drive=None, gens=20000, dt=0.01, gamma=0.2, kappa=0.05, amp=0.004, T=4, seed=0, every=20):
    rng=np.random.default_rng(seed)
    x=np.log(xstar)+0.05*rng.normal(size=6)
    v=np.array([1.,-1.,1.,-1.,1.,-1.]); v/=np.linalg.norm(v)
    om=2*np.pi/T; traj=[]; totals=[]
    for t in range(gens):
        if drive is None:
            u=amp*np.sin(om*t*dt)*v
        else:
            d=np.asarray(drive)
            if d.ndim==1:
                u=amp*(d[t%len(d)]-d.mean())
            else:
                u=amp*(d[t%len(d)]-d.mean(axis=0))
        dx=-gamma*(x-np.log(xstar))-kappa*(L_line@(x-np.log(xstar)))+u
        x=x+dt*dx
        et=np.exp(x); tot=float(et.sum()); e=et/tot
        if t%every==0:
            traj.append(e.copy()); totals.append(tot)
    return np.array(traj), np.array(totals)
def fft_peak(traj):
    y=traj[:,0]-traj[:,0].mean(); Y=np.abs(np.fft.rfft(y)); k=int(np.argmax(Y[1:])+1)
    return k, len(Y)-1
def main():
    xs=target_from_w()
    tr,tots=run6(xs)
    tail=tr[len(tr)//2:]
    rel=float(np.max(np.abs(tail.mean(0)-xs)/np.maximum(xs,1e-30)))
    # B 敏感性
    grid=[]
    for gamma in (0.05,0.2,0.5):
        for kappa in (0.0,0.05,0.2):
            for amp in (0.001,0.004,0.02):
                for T in (3,4,5,8):
                    t,_=run6(xs,gens=3000,gamma=gamma,kappa=kappa,amp=amp,T=T)
                    tl=t[len(t)//2:]
                    ok=bool(np.isfinite(t).all() and tl.min()>0 and np.max(np.abs(tl.mean(0)-xs)/np.maximum(xs,1e-30))<0.1)
                    grid.append(dict(gamma=gamma,kappa=kappa,amp=amp,T=T,ok=ok,
                                     min=float(tl.min()),rel=float(np.max(np.abs(tl.mean(0)-xs)/np.maximum(xs,1e-30)))))
    # C 18 shadow drive
    from zero_evo import Params
    from zero_evo_18 import ZeroEvo18
    import zero_sim as zs
    base=dict(micro=True,L=4,T=4,retention="top2",eta_n=0.2,eta_h=0.2,eta_psi=0.2,mu=0.5,kappa_p=0.1,
              gamma=0.2,K=2.0,capacity="soft",K_hard=1.5,g_src=1.0,c_h=1.2,Gamma_h=0.2,kappa_h=0.5,
              alpha=0.0,wave_on=True,beta=0.5,g_wave=0.5,kappa_psi=0.5,Gamma_w=0.5)
    e=ZeroEvo18(zs.k4(),Params(**base),seed=np.array([8.,4.,1.5,1.]))
    shad=[]
    for _ in range(8000):
        e.step(); shad.append(e.W18.sum(0).copy())
    shad=np.array(shad); drive=shad/shad.mean(axis=0)
    tr2,tot2=run6(xs,drive=drive,gens=4000,every=1,amp=0.05,dt=0.1,gamma=0.05,kappa=0.02)
    tl2=tr2[len(tr2)//2:]
    rel2=float(np.max(np.abs(tl2.mean(0)-xs)/np.maximum(xs,1e-30)))
    def t4amp(a):
        a=np.asarray(a); ph=[a[t::4].mean(0) for t in range(4)]
        return float(np.max(np.ptp(np.array(ph),axis=0)))
    t4_6=t4amp(tot2); t4_shape=t4amp(tr2); t4_drv=t4amp(drive[:len(tr2)])
    k,n=fft_peak(tr2)
    checks=[
      dict(name="4势→6边→x*：6层收敛（相对误差<0.1）", pass_=rel<0.1, detail=f"{rel:.2e}"),
      dict(name="敏感性：全部 108 点稳定", pass_=all(g["ok"] for g in grid),
           detail=f"{sum(g['ok'] for g in grid)}/{len(grid)}"),
      dict(name="18 shadow 驱动下 6 层稳定/保正", pass_=bool(tl2.min()>0 and np.isfinite(tr2).all()),
           detail=f"min={tl2.min():.2e}, rel={rel2:.2e}"),
      dict(name="6 层继承 T=4（总权重相位 mod 4 振幅 >1e-4）", pass_=t4_6>1e-4,
           detail=f"t4_total={t4_6:.2e}, t4_shape={t4_shape:.2e}, t4_drive={t4_drv:.2e}"),
    ]
    fig,axes=plt.subplots(1,2,figsize=(11,4.0))
    axes[0].plot(np.arange(len(tr))*20,tr[:,:],lw=.8); axes[0].set_yscale("log")
    axes[0].set_title("4→6 target: 6 layers converge"); axes[0].grid(alpha=.3)
    axes[1].plot(np.arange(len(tr2))*20,tr2[:,0],lw=.8,label="layer0 (18 shadow drive)")
    axes[1].plot(np.arange(len(tr2))*20,drive[:len(tr2)],lw=.6,label="drive")
    axes[1].set_title("T=4 inherited from 18 destruction"); axes[1].legend(fontsize=8); axes[1].grid(alpha=.3)
    fig.tight_layout(); plot=os.path.join(RES,"zero_world6_chain.png"); fig.savefig(plot,dpi=130); plt.close(fig)
    out=dict(xstar=xs.tolist(),rel=rel,grid=grid,
             coupled=dict(rel=rel2,min=float(tl2.min()),t4_total=t4_6,t4_shape=t4_shape,t4_drive=t4_drv,fft_peak=k),
             checks=[{**{k:v for k,v in c.items() if k!='pass_'},"pass":bool(c['pass_'])} for c in checks])
    with open(os.path.join(RES,"zero_world6_chain.json"),"w") as f: json.dump(out,f,indent=1)
    print("4→18→6 接起来 + 敏感性 + T=4 继承")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  PASS {sum(c['pass'] for c in out['checks'])}/{len(out['checks'])} → {plot}")
if __name__=="__main__": main()
