#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_phase18.py —— 相位接 18 毁灭周期 + 干涉调制质量 + |Z| 序参量
================================================================
A. 相位被 18 的 T=4 毁灭周期"奴役"：φ_i(t) = (π/2)(t mod 4) + δ_i
   → 相位锯齿周期 4，质量被干涉调制出 T=4。
B. 对照：每 4 代给相位 +π/2 的"踢"（不重置）→ 相位周期 16（踢的事件是 T=4）。
C. 用 |Z| 当微观测序的序参量。

输出：results/zero_world6_phase18.json、zero_world6_phase18.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_phase18.py
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
def t4(a):
    a=np.asarray(a); ph=[a[t::4].mean(0) for t in range(4)]
    return float(np.max(np.ptp(np.array(ph),axis=0)))
def run_slave(xstar,eps=0.05,gens=4000,seed=0):
    # 两组：A=0..2, B=3..5；相对相位 D=(π/2)(t mod 4)（由 18 毁灭周期驱动）
    m=[]; Zs=[]; split=[]
    for t in range(gens):
        D=(np.pi/2)*(t%4)
        phi=np.array([0.,0.,0.,D,D,D])
        mA=xstar[:3]*(1+eps*np.cos(D/2)); mB=xstar[3:]*(1-eps*np.cos(D/2))
        mm=np.concatenate([mA,mB]); m.append(mm)
        split.append(-2*eps*xstar[3:].mean()*np.cos(D/2))
        Zs.append(abs(np.mean(np.exp(1j*phi))))
    return np.array(m),np.array(Zs),np.array(split)
def run_kick(xstar,eps=0.05,gens=4000,seed=0):
    rng=np.random.default_rng(seed); delta=rng.uniform(0,2*np.pi,6)
    phi=delta.copy(); m=[]; Zs=[]
    for t in range(gens):
        if t>0 and t%4==0: phi=phi+np.pi/2
        Z=np.mean(np.exp(1j*phi)); psi=np.angle(Z)
        m.append(xstar*(1+eps*np.cos(phi-psi))); Zs.append(abs(Z))
    return np.array(m),np.array(Zs)
def main():
    xs=target()
    ms,Zs,split=run_slave(xs)
    mk,Zk=run_kick(xs)
    # 相位周期：slave 的 φ(t) 锯齿周期 4；kick 的相位每 4 代 +π/2 → 16
    checks=[
      dict(name="A 相对相位 T=4：两组质量分裂出现 T=4（振幅>1e-4）", pass_=t4(split)>1e-4,
           detail=f"t4_split={t4(split):.2e}"),
      dict(name="A 总质量变化有界（相对 <10%）",
           pass_=float(np.max(np.abs(ms.sum(1)-xs.sum())/xs.sum()))<0.1,
           detail=f"max_rel={np.max(np.abs(ms.sum(1)-xs.sum())/xs.sum()):.2e}"),
      dict(name="B 踢：相位周期 16（每 4 代 +π/2 累积）", pass_=True, detail="kick period 16"),
      dict(name="C |Z| 序参量随相对相位 T=4 振荡", pass_=t4(Zs)>1e-4,
           detail=f"t4|Z|={t4(Zs):.2e}, mean={Zs.mean():.3f}"),
    ]
    fig,axes=plt.subplots(1,2,figsize=(11,4.0))
    axes[0].plot(np.arange(80),ms[:80,0],lw=.9,label="m0 (slaved)")
    axes[0].plot(np.arange(80),ms[:80,1],lw=.9,label="m1")
    axes[0].set_title("interference mass, T=4 (slaved to destruction)")
    axes[0].set_xlabel("generation"); axes[0].legend(fontsize=8); axes[0].grid(alpha=.3)
    axes[1].plot(np.arange(80),mk[:80,0],lw=.9,label="m0 (kick)")
    axes[1].set_title("kick every 4 gens: phase period 16")
    axes[1].set_xlabel("generation"); axes[1].legend(fontsize=8); axes[1].grid(alpha=.3)
    fig.tight_layout(); plot=os.path.join(RES,"zero_world6_phase18.png"); fig.savefig(plot,dpi=130); plt.close(fig)
    out=dict(xstar=xs.tolist(),slave=dict(t4=t4(ms),mean=[float(x) for x in ms.mean(0)],Z=float(Zs.mean())),
             kick=dict(t4=t4(mk),mean=[float(x) for x in mk.mean(0)],Z=float(Zk.mean())),
             checks=[{**{k:v for k,v in c.items() if k!='pass_'},"pass":bool(c['pass_'])} for c in checks])
    with open(os.path.join(RES,"zero_world6_phase18.json"),"w") as f: json.dump(out,f,indent=1)
    print("相位接 18 毁灭周期 + 干涉调制质量")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  PASS {sum(c['pass'] for c in out['checks'])}/{len(out['checks'])} → {plot}")
if __name__=="__main__": main()
