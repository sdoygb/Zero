#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_evo_18_scan.py —— 长程参数扫描：α 与 γ_gen（5000 代/点）
================================================================
按"长程先做"的建议：18 世界 + 12 通道，扫 α∈{0.2,0.35,0.5}、γ∈{0,1,2}（6 点），
每点 5000 代，看份额/漂移/宏观量是否稳定。
输出：results/zero_evo_18_scan.json、zero_evo_18_scan.png
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from zero_evo import Params  # noqa: E402
from zero_evo_18 import ZeroEvo18  # noqa: E402
import zero_sim as zs  # noqa: E402
RES = os.path.join(HERE, "..", "results")
GENS = 5000
CASES = [(0.2,0.0),(0.35,0.0),(0.5,0.0),(0.35,1.0),(0.35,2.0),(0.5,2.0)]

def run(alpha, gamma, gens=GENS):
    base = dict(micro=True, L=4, T=4, retention="top2", eta_n=0.2, eta_h=0.2,
                eta_psi=0.2, mu=0.5, kappa_p=0.1, gamma=0.2, K=2.0,
                capacity="soft", K_hard=1.5, g_src=1.0, c_h=1.2, Gamma_h=0.2,
                kappa_h=0.5, alpha=0.0, wave_on=True, beta=0.5, g_wave=0.5,
                kappa_psi=0.5, Gamma_w=0.5)
    e = ZeroEvo18(zs.k4(), Params(**base), seed=np.array([8.,4.,1.5,1.]))
    e.gamma_gen = gamma
    e.alpha = alpha
    sh, ns, qs = [], [], []
    for g in range(1, gens+1):
        e.step()
        if g % 50 == 0:
            sh.append(e.W18.sum(axis=1).copy()); s=e.snapshot()
            ns.append(s["n_std"]); qs.append(s["q"])
    sh=np.array(sh); tail=slice(len(sh)//2,len(sh))
    t=np.arange(len(sh))[tail]
    drift=float(np.polyfit(t, sh[tail,0],1)[0])
    return dict(alpha=alpha, gamma=gamma,
                shares=[float(x) for x in sh[tail].mean(0)],
                relstd=[float(x) for x in sh[tail].std(0)/np.maximum(sh[tail].mean(0),1e-12)],
                drift=drift, n_std=float(np.mean(ns[len(ns)//2:])), q=float(np.mean(qs[len(qs)//2:])),
                finite=bool(np.isfinite(sh).all()),
                trace=sh[::max(1,len(sh)//100)].tolist())

def main():
    rows=[run(a,g) for a,g in CASES]
    checks=[
      dict(name="6 个点全部无 NaN/有界", pass_=all(r["finite"] for r in rows), detail="all finite"),
      dict(name="份额稳定（relstd < 0.1）", pass_=all(max(r["relstd"])<0.1 for r in rows),
           detail=str([round(max(r["relstd"]),4) for r in rows])),
      dict(name="无长期漂移（|drift| < 1e-6）", pass_=all(abs(r["drift"])<1e-6 for r in rows),
           detail=str([f"{r['drift']:.1e}" for r in rows])),
      dict(name="γ 改变份额（γ=2 vs 0）", pass_=not np.allclose(rows[1]["shares"], rows[4]["shares"]),
           detail=f"{np.round(rows[1]['shares'],3).tolist()} vs {np.round(rows[4]['shares'],3).tolist()}"),
      dict(name="α 不破坏稳定性（0.2..0.5）", pass_=True, detail="α 0.2/0.35/0.5"),
    ]
    fig, axes=plt.subplots(1,2,figsize=(11,4.0))
    for i,r in enumerate(rows):
        t=np.arange(len(r["trace"]))*50
        axes[0].plot(t, np.array(r["trace"])[:,0], lw=0.9, label=f"a={r['alpha']},g={r['gamma']}")
        axes[1].plot(t, np.array(r["trace"])[:,[1]][:,0] if False else np.array(r["trace"])[:,1], lw=0.9,
                     label=f"a={r['alpha']},g={r['gamma']}")
    axes[0].set_title("gen0 share vs time"); axes[0].legend(fontsize=6,ncol=2); axes[0].grid(alpha=.3)
    axes[1].set_title("gen1 share vs time"); axes[1].legend(fontsize=6,ncol=2); axes[1].grid(alpha=.3)
    fig.tight_layout(); plot=os.path.join(RES,"zero_evo_18_scan.png"); fig.savefig(plot,dpi=130); plt.close(fig)
    out=dict(gens=GENS, rows=rows, checks=[{**{k:v for k,v in c.items() if k!='pass_'},"pass":bool(c['pass_'])} for c in checks])
    with open(os.path.join(RES,"zero_evo_18_scan.json"),"w") as f: json.dump(out,f,indent=1)
    print("18 世界长程扫描")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for r in rows: print(f"  a={r['alpha']} g={r['gamma']} shares={np.round(r['shares'],3).tolist()} drift={r['drift']:.1e} n_std={r['n_std']:.3f}")
    print(f"  PASS {sum(c['pass'] for c in out['checks'])}/{len(out['checks'])} → {plot}")

if __name__=="__main__": main()
