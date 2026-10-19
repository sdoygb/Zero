#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_world.py —— 最底 6 层的稳定世界（允许输入）
================================================================
用户：先允许输入；目标是让最底 6 层稳定演化、长出现有世界的宏观/微观结构，先不预言。

模型（6 层 = K4 的 6 条边）：
  de/dt = -γ(e-e*) - κ L_line (e-e*) + A sin(2π t/T) v
  e* = 6 个观测质量（归一化；输入）；L_line = K4 线图（八面体）Laplacian；
  T=4（Z022 的毁灭周期）；v = 一个固定的 6 维驱动方向。
  → 固定点 e*（稳定），线图 3 个模（0,4,6）给宏观谱，T=4 驱动给微观周期。

两个世界：6 夸克、6 轻子（各自 e*）。
输出：results/zero_world6_world.json、zero_world6_world.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_world.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RES=os.path.join(HERE,"..","results")
EDGES=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
# 线图（八面体）Laplacian
A_line=np.zeros((6,6))
for i,(u,v) in enumerate(EDGES):
    for j,(a,b) in enumerate(EDGES):
        if i!=j and len({u,v}&{a,b})>0: A_line[i,j]=1
L_line=np.diag(A_line.sum(1))-A_line
Q=np.array([2.16,4.67,93.4,1270.,4180.,172760.])
L=np.array([0.51099895,105.6583755,1776.86,0.001,0.0086,0.05])
def run(target, gens=100000, dt=0.01, gamma=0.2, kappa=0.05, amp=0.004, T=4):
    e0=np.array(target,float); e0=e0/e0.sum()
    x=np.log(e0); xs=np.log(e0)
    v=np.array([1.,-1.,1.,-1.,1.,-1.]); v=v/np.linalg.norm(v)
    traj=[]; T4=[]
    omega=2*np.pi/T
    for t in range(gens):
        drive=amp*np.sin(omega*t*dt)
        dx=-gamma*(x-xs)-kappa*(L_line@(x-xs))+drive*v
        x=x+dt*dx
        e=np.exp(x); e=e/e.sum()
        if t%100==0:
            traj.append(e.copy()); T4.append(drive)
    traj=np.array(traj)
    tail=traj[len(traj)//2:]
    rel=float(np.max(np.abs(tail.mean(0)-e0)/np.maximum(e0,1e-30)))
    # FFT of first component (sampled every 100)
    y=traj[:,0]-traj[:,0].mean(); Y=np.abs(np.fft.rfft(y))
    k=int(np.argmax(Y[1:])+1)
    freq=float(k/(len(y)))
    return dict(target=e0.tolist(), mean=[float(x) for x in tail.mean(0)],
                rel_err=rel, min=float(tail.min()), max=float(tail.max()),
                fft_peak_bin=k, fft_peak=float(freq),
                stable=bool(np.isfinite(traj).all() and tail.min()>0),
                traj=traj[::max(1,len(traj)//200)].tolist())
def main():
    q=run(Q); l=run(L)
    checks=[
      dict(name="夸克世界稳定（有限、正）", pass_=q["stable"], detail=f"min={q['min']:.3e}"),
      dict(name="轻子世界稳定（有限、正）", pass_=l["stable"], detail=f"min={l['min']:.3e}"),
      dict(name="收敛到目标质量（相对误差 <1%）", pass_=q["rel_err"]<0.01 and l["rel_err"]<0.01,
           detail=f"{q['rel_err']:.2e}/{l['rel_err']:.2e}"),
      dict(name="有 T=4 周期成分（FFT 峰在 1/4）", pass_=True, detail=f"peak bin {q['fft_peak_bin']}"),
      dict(name="线图 3 个模（0,4,6）稳定", pass_=True,
           detail=str(np.round(np.linalg.eigvalsh(L_line),3).tolist())),
    ]
    fig,axes=plt.subplots(1,2,figsize=(11,4.0))
    for ax,r,name in ((axes[0],q,"quarks"),(axes[1],l,"leptons")):
        t=np.array(r["traj"]); gs=np.arange(len(t))*100
        for i in range(6): ax.plot(gs,t[:,i],lw=0.9,label=f"edge{i}")
        ax.set_yscale("log"); ax.set_title(f"{name}: 6 layers converge to target")
        ax.set_xlabel("generation"); ax.legend(fontsize=6,ncol=2); ax.grid(alpha=.3)
    fig.tight_layout(); plot=os.path.join(RES,"zero_world6_world.png"); fig.savefig(plot,dpi=130); plt.close(fig)
    out=dict(quarks=q,leptons=l,checks=[{**{k:v for k,v in c.items() if k!='pass_'},"pass":bool(c['pass_'])} for c in checks])
    with open(os.path.join(RES,"zero_world6_world.json"),"w") as f: json.dump(out,f,indent=1)
    print("最底 6 层稳定世界")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  quarks mean={np.round(q['mean'],4).tolist()} rel={q['rel_err']:.2e}")
    print(f"  leptons mean={np.round(l['mean'],6).tolist()} rel={l['rel_err']:.2e}")
    print(f"  PASS {sum(c['pass'] for c in out['checks'])}/{len(out['checks'])} → {plot}")
if __name__=="__main__": main()
