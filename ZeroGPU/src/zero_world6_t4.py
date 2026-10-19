#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_t4.py —— 把 T=4 传进 6 层：相位/门控 vs 振幅驱动
================================================================
上一轮：shadow 振幅驱动被滤掉（响应 1.9e-5）。本轮试两种"传 T=4"的方式：
  A. **相位/方波**：u = amp * sign(sin(2πt/T)) * v（每 4 代翻转一次）；
  B. **毁灭门控**：每 4 代对 6 层做一次"保留 top2"的毁灭/重播种（照 18 的毁灭规则）。
对照：shadow 振幅驱动（上一轮，失败）。

输出：results/zero_world6_t4.json、zero_world6_t4.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_t4.py
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
def run(xstar, mode, gens=4000, dt=0.01, gamma=0.2, kappa=0.05, amp=0.02, T=4):
    x=np.log(xstar); xs=np.log(xstar); v=np.array([1.,-1.,1.,-1.,1.,-1.]); v/=np.linalg.norm(v)
    traj=[]; tot=[]
    for t in range(gens):
        if mode=="square":
            s=1.0 if (t%T)<T//2 else -1.0; u=amp*s*v
        else: u=np.zeros(6)
        dx=-gamma*(x-xs)-kappa*(L_line@(x-xs))+u
        x=x+dt*dx
        if mode=="destroy" and t>0 and t%T==0:
            e=np.exp(x); e/=e.sum()
            idx=np.argsort(e)[::-1][:2]; mask=np.zeros(6); mask[idx]=1
            e=e*mask; e/=e.sum(); x=np.log(np.maximum(e,1e-300))
        et=np.exp(x); tot.append(float(et.sum())); e=et/et.sum(); traj.append(e)
    return np.array(traj), np.array(tot)
def t4(a):
    a=np.asarray(a); ph=[a[t::4].mean(0) for t in range(4)]
    return float(np.max(np.ptp(np.array(ph),axis=0)))
def main():
    xs=target()
    r={}
    for mode in ("square","destroy"):
        tr,tot=run(xs,mode); tail=tr[len(tr)//2:]
        r[mode]=dict(shape=float(t4(tr)),total=float(t4(tot)),
                     rel=float(np.max(np.abs(tail.mean(0)-xs)/np.maximum(xs,1e-30))),
                     min=float(tail.min()),ok=bool(np.isfinite(tr).all() and tail.min()>0),
                     trace=tr[:800].tolist())
    checks=[
      dict(name="A 相位/方波：6 层出现 T=4（形状振幅 >1e-4）", pass_=r["square"]["shape"]>1e-4,
           detail=f"shape={r['square']['shape']:.2e}, total={r['square']['total']:.2e}"),
      dict(name="B 毁灭门控：6 层出现 T=4（形状振幅 >1e-4）", pass_=r["destroy"]["shape"]>1e-4,
           detail=f"shape={r['destroy']['shape']:.2e}"),
      dict(name="两种方式都稳定/保正", pass_=r["square"]["ok"] and r["destroy"]["ok"],
           detail=f"min={r['square']['min']:.1e}/{r['destroy']['min']:.1e}"),
      dict(name="两种方式都收敛到 x*（相对误差 <0.1）", pass_=r["square"]["rel"]<0.1 and r["destroy"]["rel"]<0.1,
           detail=f"{r['square']['rel']:.2e}/{r['destroy']['rel']:.2e}"),
    ]
    fig,axes=plt.subplots(1,2,figsize=(11,4.0))
    for ax,key in ((axes[0],"square"),(axes[1],"destroy")):
        t=np.array(r[key]["trace"]); gs=np.arange(len(t))
        for i in range(6): ax.plot(gs,t[:,i],lw=0.9,label=f"L{i}")
        ax.set_title(f"{key}: T=4 transmission"); ax.set_xlabel("generation"); ax.legend(fontsize=6,ncol=2); ax.grid(alpha=.3)
    fig.tight_layout(); plot=os.path.join(RES,"zero_world6_t4.png"); fig.savefig(plot,dpi=130); plt.close(fig)
    out=dict(xstar=xs.tolist(),square=r["square"],destroy=r["destroy"],
             checks=[{**{k:v for k,v in c.items() if k!='pass_'},"pass":bool(c['pass_'])} for c in checks])
    with open(os.path.join(RES,"zero_world6_t4.json"),"w") as f: json.dump(out,f,indent=1)
    print("T=4 传进 6 层")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  PASS {sum(c['pass'] for c in out['checks'])}/{len(out['checks'])} → {plot}")
if __name__=="__main__": main()
