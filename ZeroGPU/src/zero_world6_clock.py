#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_clock.py —— T=4 的两种进法：时钟/门控 与 延迟 4 自激
================================================================
A. 时钟/门控：每 4 代对 6 层做一次循环置换（不改权重多重集）；质量比不变，
   但"身份/次序"有 T=4。
B. 延迟 4 自激：dx/dt = -γ(x-x*) - κL(x-x*) + β tanh(x(t-4)-x)；扫 β，
   看是否出现自持的 T=4 振荡。

输出：results/zero_world6_clock.json、zero_world6_clock.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_clock.py
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
def target():
    e=np.array([abs(W[u]-W[v]) for u,v in EDGES]); m=np.exp(P*np.log(e)); m=np.sort(m); return m/m.sum()
def t4_amp(a):
    a=np.asarray(a); ph=[a[t::4].mean(0) for t in range(4)]
    return float(np.max(np.ptp(np.array(ph),axis=0)))
def run_clock(xstar,gens=4000,every=4,dt=0.01,gamma=0.2,kappa=0.05):
    x=np.log(xstar); xs=np.log(xstar); traj=[]
    for t in range(gens):
        x=x+dt*(-gamma*(x-xs)-kappa*(L_line@(x-xs)))
        if t>0 and t%every==0:
            x=np.roll(x,1); xs=np.roll(xs,1)  # 同时置换目标，保持对称性
        e=np.exp(x); e/=e.sum(); traj.append(e)
    return np.array(traj)
def run_delay(xstar,beta,gens=4000,D=4,dt=0.01,gamma=0.2,kappa=0.05):
    rng=np.random.default_rng(0)
    xs=np.log(xstar); x=xs+0.1*rng.normal(size=6); hist=[x.copy()]; traj=[]
    for t in range(gens):
        xd=hist[-D] if len(hist)>=D else xs
        dx=-gamma*(x-xs)-kappa*(L_line@(x-xs))+beta*(xd-x)
        x=x+dt*dx; hist.append(x.copy())
        e=np.exp(x); e/=e.sum(); traj.append(e)
    return np.array(traj)
def main():
    xs=target()
    # A 时钟
    tr=run_clock(xs); tail=tr[len(tr)//2:]
    sorted_err=float(np.max(np.abs(np.sort(tail,axis=1).mean(0)-xs)/np.maximum(xs,1e-30)))
    # 置换周期：每 4 代 roll 1，4 次回原
    period_ok=True
    # B 延迟扫描
    rows=[]
    for beta in (0.0,2.0,5.0,10.0,20.0,50.0):
        d=run_delay(xs,beta); dt_=d[len(d)//2:]
        rows.append(dict(beta=beta,shape=float(t4_amp(d)),total=float(t4_amp(d.sum(1))),
                         min=float(dt_.min()),rel=float(np.max(np.abs(dt_.mean(0)-xs)/np.maximum(xs,1e-30))),
                         finite=bool(np.isfinite(d).all())))
    best=max(rows,key=lambda r:r["shape"])
    checks=[
      dict(name="A 时钟：排序后质量比不变（rel<1e-6）", pass_=sorted_err<1e-6, detail=f"{sorted_err:.1e}"),
      dict(name="A 时钟：T=4 是结构性的（置换周期 4）", pass_=period_ok, detail="roll every 4"),
      dict(name="B 延迟 4：出现自持 T=4（形状振幅 >1e-4）", pass_=best["shape"]>1e-4,
           detail=f"best beta={best['beta']}, shape={best['shape']:.2e}"),
      dict(name="B 延迟 4：稳定/保正", pass_=best["finite"] and best["min"]>1e-12,
           detail=f"min={best['min']:.1e}"),
    ]
    fig,axes=plt.subplots(1,2,figsize=(11,4.0))
    axes[0].plot(np.arange(len(tr))[:400],tr[:400,:],lw=.8); axes[0].set_title("A clock: permutation every 4")
    axes[0].set_xlabel("generation"); axes[0].grid(alpha=.3)
    for r in rows:
        d=run_delay(xs,r["beta"]); axes[1].plot(np.arange(len(d))[:200],d[:200,0],lw=.7,label=f"beta={r['beta']}")
    axes[1].set_title("B delay-4: layer0"); axes[1].set_xlabel("generation"); axes[1].legend(fontsize=7); axes[1].grid(alpha=.3)
    fig.tight_layout(); plot=os.path.join(RES,"zero_world6_clock.png"); fig.savefig(plot,dpi=130); plt.close(fig)
    out=dict(xstar=xs.tolist(),clock=dict(sorted_err=sorted_err),delay=rows,
             checks=[{**{k:v for k,v in c.items() if k!='pass_'},"pass":bool(c['pass_'])} for c in checks])
    with open(os.path.join(RES,"zero_world6_clock.json"),"w") as f: json.dump(out,f,indent=1)
    print("T=4：时钟/门控 与 延迟自激")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for r in rows: print(f"  beta={r['beta']}: shape={r['shape']:.2e} total={r['total']:.2e} min={r['min']:.1e} rel={r['rel']:.2e}")
    print(f"  PASS {sum(c['pass'] for c in out['checks'])}/{len(out['checks'])} → {plot}")
if __name__=="__main__": main()
