#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_evo_18_long_lam.py —— λ=10, γ=2 的长跑（10^5 代）
输出：results/zero_evo_18_long_lam.json
"""
from __future__ import annotations
import json, os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from zero_evo import Params  # noqa: E402
from zero_evo_18 import ZeroEvo18  # noqa: E402
import zero_sim as zs  # noqa: E402
RES=os.path.join(HERE,"..","results"); GENS=100000
BASE=dict(micro=True,L=4,T=4,retention="top2",eta_n=0.2,eta_h=0.2,eta_psi=0.2,mu=0.5,kappa_p=0.1,
          gamma=0.2,K=2.0,capacity="soft",K_hard=1.5,g_src=1.0,c_h=1.2,Gamma_h=0.2,kappa_h=0.5,
          alpha=0.0,wave_on=True,beta=0.5,g_wave=0.5,kappa_psi=0.5,Gamma_w=0.5)
def main():
    e=ZeroEvo18(zs.k4(),Params(**BASE),seed=np.array([8.,4.,1.5,1.]))
    e.gamma_gen=2.0; e.lambda_gen=10.0
    sh=[]; ns=[]
    for g in range(1,GENS+1):
        e.step()
        if g%200==0:
            sh.append(e.W18.sum(1).copy()); ns.append(float(e.n.std()))
    sh=np.array(sh); tail=sh[len(sh)//2:]
    t=np.arange(len(sh))[len(sh)//2:]
    drift=float(np.polyfit(t,sh[len(sh)//2:,0],1)[0])
    out=dict(gens=GENS,shares=[float(x) for x in tail.mean(0)],
             relstd=[float(x) for x in tail.std(0)/np.maximum(tail.mean(0),1e-12)],
             drift=drift,n_std=float(np.mean(ns[len(ns)//2:])),
             finite=bool(np.isfinite(sh).all()),
             tail_first=[float(x) for x in sh[len(sh)//4]], tail_last=[float(x) for x in sh[-1]])
    with open(os.path.join(RES,"zero_evo_18_long_lam.json"),"w") as f: json.dump(out,f,indent=1)
    print(out)
if __name__=="__main__": main()
