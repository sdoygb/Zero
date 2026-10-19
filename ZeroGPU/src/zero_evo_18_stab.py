#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_evo_18_stab.py —— 给 gamma>0 加份额负反馈，找稳定化窗口
================================================================
上一轮：gamma>0 长程漂移。本轮在代权重上加减项 -lambda*(share - gen_w)，
扫 lambda ∈ {0,1,5,10,20}，gamma=2，5000 代/点。
输出：results/zero_evo_18_stab.json、zero_evo_18_stab.png
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from zero_evo import Params  # noqa: E402
from zero_evo_18 import ZeroEvo18  # noqa: E402
import zero_sim as zs  # noqa: E402
RES=os.path.join(HERE,"..","results")
BASE=dict(micro=True,L=4,T=4,retention="top2",eta_n=0.2,eta_h=0.2,eta_psi=0.2,mu=0.5,kappa_p=0.1,
          gamma=0.2,K=2.0,capacity="soft",K_hard=1.5,g_src=1.0,c_h=1.2,Gamma_h=0.2,kappa_h=0.5,
          alpha=0.0,wave_on=True,beta=0.5,g_wave=0.5,kappa_psi=0.5,Gamma_w=0.5)
def run(lam,gamma=2.0,gens=5000):
    e=ZeroEvo18(zs.k4(),Params(**BASE),seed=np.array([8.,4.,1.5,1.]))
    e.gamma_gen=gamma; e.lambda_gen=lam
    sh=[]
    for _ in range(gens):
        e.step(); sh.append(e.W18.sum(1).copy())
    sh=np.array(sh); tail=sh[gens//2:]; t=np.arange(len(sh))[gens//2:]
    drift=float(np.polyfit(t,sh[gens//2:,0],1)[0])
    return dict(lam=lam,shares=[float(x) for x in tail.mean(0)],
                relstd=[float(x) for x in tail.std(0)/np.maximum(tail.mean(0),1e-12)],
                drift=drift,finite=bool(np.isfinite(sh).all()),
                trace=sh[::max(1,len(sh)//100)].tolist())
def main():
    rows=[run(l) for l in (0.0,1.0,5.0,10.0,20.0)]
    stable=[r for r in rows if max(r["relstd"])<0.1 and abs(r["drift"])<1e-6]
    checks=[
      dict(name="lambda=0（无反馈）不稳", pass_=max(rows[0]["relstd"])>0.1 or abs(rows[0]["drift"])>1e-6,
           detail=f"relstd={np.round(rows[0]['relstd'],3).tolist()}, drift={rows[0]['drift']:.1e}"),
      dict(name="存在稳定窗口（relstd<0.1 且 |drift|<1e-6）", pass_=bool(stable),
           detail=str([r["lam"] for r in stable])),
      dict(name="lambda=10 稳定", pass_=bool(any(abs(r["lam"]-10)<1e-9 and max(r["relstd"])<0.1 for r in rows)),
           detail=f"relstd={np.round([r for r in rows if r['lam']==10][0]['relstd'],3).tolist()}"),
      dict(name="反馈不破坏有界性", pass_=all(r["finite"] for r in rows), detail="all finite"),
    ]
    fig,ax=plt.subplots(figsize=(8,4.0))
    for r in rows:
        ax.plot(np.arange(len(r["trace"]))*50, np.array(r["trace"])[:,0], lw=0.9, label=f"lam={r['lam']}")
    ax.set_xlabel("generation"); ax.set_ylabel("gen0 share")
    ax.set_title("negative feedback on shares (gamma=2)")
    ax.legend(fontsize=8); ax.grid(alpha=.3)
    fig.tight_layout(); plot=os.path.join(RES,"zero_evo_18_stab.png"); fig.savefig(plot,dpi=130); plt.close(fig)
    out=dict(rows=rows,checks=[{**{k:v for k,v in c.items() if k!='pass_'},"pass":bool(c['pass_'])} for c in checks])
    with open(os.path.join(RES,"zero_evo_18_stab.json"),"w") as f: json.dump(out,f,indent=1)
    print("gamma=2 + 份额负反馈")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for r in rows: print(f"  lam={r['lam']}: shares={np.round(r['shares'],3).tolist()} relstd={np.round(r['relstd'],3).tolist()} drift={r['drift']:.1e}")
    print(f"  PASS {sum(c['pass'] for c in out['checks'])}/{len(out['checks'])} → {plot}")
if __name__=="__main__": main()
