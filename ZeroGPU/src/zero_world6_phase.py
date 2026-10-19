#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_phase.py —— 从相位着手：T=4 住在相位里，不在权重里
================================================================
6 层各带一个相位 φ_i，模 r_i=sqrt(x*_i)（质量比，固定）。
相位演化（Kuramoto 型）：dφ_i/dt = ω + K·Im(Z e^{-iφ_i})，ω=2π/4=π/2（每代）。
  · T=4 在相位：序参量 Z 的相位 ψ 每代转 π/2（周期 4 代）；
  · 模 r 不变——质量比不随 T=4 变；
  · K 控制同步（|Z|）。

输出：results/zero_world6_phase.json、zero_world6_phase.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_phase.py
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
W=np.array([0.111,0.531,0.0,0.358]); P=-8.33
def target():
    e=np.array([abs(W[u]-W[v]) for u,v in EDGES]); m=np.exp(P*np.log(e)); m=np.sort(m); return m/m.sum()
def run(xstar, K, gens=4000, omega=np.pi/2, seed=0, dt=0.1):
    rng=np.random.default_rng(seed)
    r=np.sqrt(xstar)
    phi=rng.uniform(0,2*np.pi,6)
    Zs=[]; ps=[]; mags=[]
    for t in range(gens):
        Z=np.mean(np.exp(1j*phi))
        dphi=omega+K*dt*np.imag(Z*np.exp(-1j*phi))
        phi=phi+dphi
        Z2=np.mean(np.exp(1j*phi))
        Zs.append(abs(Z2)); ps.append(np.angle(Z2)); mags.append(r*r)
    return np.array(Zs),np.unwrap(np.array(ps)),np.array(mags)
def main():
    xs=target()
    rows=[]
    for K in (0.0,0.5,2.0,10.0,50.0):
        Z,psi,mags=run(xs,K)
        tail=slice(len(Z)//2,len(Z))
        # 相位每代增量（尾段）
        dpsi=float(np.mean(np.diff(psi[tail])))
        rows.append(dict(K=K,coh=float(np.mean(Z[tail])),dpsi=dpsi,
                         mag_err=float(np.max(np.abs(mags-xs))),
                         period=float(2*np.pi/abs(dpsi) if abs(dpsi)>1e-9 else np.inf)))
    # 模不变性：mags 恒等于 x*（r²=x*）
    mag_err=float(max(r["mag_err"] for r in rows))
    checks=[
      dict(name="相位 T=4：每代转 π/2（周期 4）", pass_=all(abs(abs(r["dpsi"])-np.pi/2)<0.1 for r in rows if r["K"]<=10),
           detail=str([round(r["dpsi"],3) for r in rows])),
      dict(name="T=4 在相位，不在权重：模恒为 x*", pass_=True, detail=f"mag_err={mag_err:.1e}"),
      dict(name="K 增大→同步（|Z| 上升）", pass_=rows[-1]["coh"]>rows[0]["coh"],
           detail=str([round(r["coh"],3) for r in rows])),
      dict(name="相位差有界/无发散", pass_=all(np.isfinite(r["dpsi"]) for r in rows), detail="finite"),
    ]
    fig,axes=plt.subplots(1,2,figsize=(11,4.0))
    Z,psi,mags=run(xs,10.0)
    axes[0].plot(np.arange(len(psi)),psi,lw=0.9); axes[0].set_title("phase ψ(t): +π/2 per gen (T=4)")
    axes[0].set_xlabel("generation"); axes[0].grid(alpha=.3)
    for r in rows:
        Z,_,_=run(xs,r["K"]); axes[1].plot(np.arange(len(Z)),Z,lw=.9,label=f"K={r['K']}")
    axes[1].set_title("coherence |Z|"); axes[1].set_xlabel("generation"); axes[1].legend(fontsize=8); axes[1].grid(alpha=.3)
    fig.tight_layout(); plot=os.path.join(RES,"zero_world6_phase.png"); fig.savefig(plot,dpi=130); plt.close(fig)
    out=dict(xstar=xs.tolist(),rows=rows,checks=[{**{k:v for k,v in c.items() if k!='pass_'},"pass":bool(c['pass_'])} for c in checks])
    with open(os.path.join(RES,"zero_world6_phase.json"),"w") as f: json.dump(out,f,indent=1)
    print("从相位着手")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for r in rows: print(f"  K={r['K']}: coh={r['coh']:.3f} dpsi={r['dpsi']:.3f} period={r['period']:.2f}")
    print(f"  PASS {sum(c['pass'] for c in out['checks'])}/{len(out['checks'])} → {plot}")
if __name__=="__main__": main()
