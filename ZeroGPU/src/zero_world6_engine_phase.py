#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_engine_phase.py —— 用世界模型拟合当前相位（微观、引擎内）
================================================================
模型：6 层质量被 T=4 相对相位调制
    m_i(φ) = m_i^0 · exp( ε cos(2πφ + δ_i) + c )
m_i^0 = 4→6 差映射的预测质量；δ_i 取 3 个匹配的相位 0, 2π/3, 4π/3（每代 2 层）。
拟合 ε、φ、c 到 6 个观测质量（对数 χ²），看 φ 是否被约束。
输出：results/zero_world6_engine_phase.json、zero_world6_engine_phase.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_engine_phase.py
"""
from __future__ import annotations
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"..","results")
EDGES=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
W=np.array([0.111,0.531,0.0,0.358]); P=-8.33
def m0_pred():
    e=np.array([abs(W[u]-W[v]) for u,v in EDGES]); m=np.exp(P*np.log(e)); return np.sort(m)/np.sort(m).sum()
def fit(obs, delta):
    obs=np.sort(np.asarray(obs,float)); obs=obs/obs.sum()
    m0=m0_pred()
    def chi(phi):
        best=None
        for eps0 in (-1.0,-0.5,0.0,0.5,1.0):
            def loss(x):
                eps,c=x
                pred=m0*np.exp(eps*np.cos(2*np.pi*phi+delta)+c)
                pred=pred/pred.sum()
                return float(np.sum((np.log(pred)-np.log(obs))**2))
            r=minimize(loss,[eps0,0.0],method="Nelder-Mead",options={"maxiter":3000,"xatol":1e-10,"fatol":1e-12})
            if best is None or r.fun<best.fun: best=r
        return best.fun,best.x
    phis=np.linspace(0,1,201,endpoint=False); vals=[]
    for phi in phis:
        f,x=chi(phi); vals.append(f)
    vals=np.array(vals); k=int(np.argmin(vals))
    return phis,vals,float(phis[k]),float(vals[k]),float(vals.mean()-vals[k]),m0,obs
def main():
    Q=np.array([2.16,4.67,93.4,1270.,4180.,172760.])
    L=np.array([0.51099895,105.6583755,1776.86,0.001,0.0086,0.05])
    delta=np.array([0,0,2*np.pi/3,2*np.pi/3,4*np.pi/3,4*np.pi/3])
    phq,chq,bq,minq,gapq,m0q,obsq=fit(Q,delta)
    phl,chl,bl,minl,gapl,m0l,obsl=fit(L,delta)
    checks=[
      dict(name="φ 有唯一极小（Δχ²>4）——夸克", pass_=gapq>4, detail=f"φ*={bq:.3f}, Δχ²={gapq:.2f}"),
      dict(name="φ 有唯一极小（Δχ²>4）——轻子", pass_=gapl>4, detail=f"φ*={bl:.3f}, Δχ²={gapl:.2f}"),
      dict(name="两扇区 φ* 一致（|Δφ|<0.1）", pass_=abs(bq-bl)<0.1, detail=f"{bq:.3f} vs {bl:.3f}"),
      dict(name="拟合可接受（min χ² < 6·σ²，σ~0.5 dex）", pass_=(minq<6*0.25 and minl<6*0.25),
           detail=f"χ²={minq:.2f}/{minl:.2f}"),
    ]
    fig,axes=plt.subplots(1,2,figsize=(11,4.0))
    axes[0].plot(phq,chq,lw=1.0); axes[0].set_title("quarks: χ²(φ)"); axes[0].set_xlabel("φ"); axes[0].grid(alpha=.3)
    axes[1].plot(phl,chl,lw=1.0); axes[1].set_title("leptons: χ²(φ)"); axes[1].set_xlabel("φ"); axes[1].grid(alpha=.3)
    fig.tight_layout(); plot=os.path.join(RES,"zero_world6_engine_phase.png"); fig.savefig(plot,dpi=130); plt.close(fig)
    out=dict(delta=delta.tolist(),quarks=dict(phi=bq,chi2=minq,gap=gapq),
             leptons=dict(phi=bl,chi2=minl,gap=gapl),
             checks=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES,"zero_world6_engine_phase.json"),"w") as f: json.dump(out,f,indent=1)
    print("世界模型拟合当前相位")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  quarks φ*={bq:.3f} χ²={minq:.2f} gap={gapq:.2f}; leptons φ*={bl:.3f} χ²={minl:.2f} gap={gapl:.2f}")
if __name__=="__main__": main()
