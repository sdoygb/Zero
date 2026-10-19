#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_stream_search3.py —— 脉冲星 glitch / FRB / CMB 极化峰 + 引擎 record
================================================================
判据同前：间隔两值 {s, L=T+2}、s 游程 = T/2−1。
数据为近似公开值；另测 zero_evo 的 gen_closures 序列。
输出：results/zero_world6_stream_search3.json
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_stream_search3.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE); RES=os.path.join(HERE,"..","results")
def score(times):
    t=np.sort(np.asarray(times,float)); ints=np.diff(t)
    if len(ints)==0: return dict(s=0,L=0,score=0.0,run_ok=False,n=len(t))
    s=float(np.min(ints)); L=float(np.max(ints))
    if L<=s*1.05: return dict(s=s,L=L,score=0.0,run_ok=False,n=len(t))
    in_s=np.abs(ints-s)/s<0.1; in_L=np.abs(ints-L)/L<0.1
    sc=float(np.mean(in_s|in_L)); T=int(round(L-2))
    runs=[]; cur=0
    for x in in_s:
        if x: cur+=1
        elif cur>0: runs.append(cur); cur=0
    if cur>0: runs.append(cur)
    run_ok=bool(len(runs)>0 and abs(np.median(runs)-(T/2-1))<1)
    return dict(s=s,L=L,score=sc,run_ok=run_ok,T=T,n=len(t))
def engine_stream(gens=3000):
    from zero_evo import Params, ZeroEvo
    import zero_sim as zs
    base=dict(micro=True,L=4,T=4,retention="top2",eta_n=0.2,eta_h=0.2,eta_psi=0.2,mu=0.5,kappa_p=0.1,
              gamma=0.2,K=2.0,capacity="soft",K_hard=1.5,g_src=1.0,c_h=1.2,Gamma_h=0.2,kappa_h=0.5,
              alpha=0.0,wave_on=True,beta=0.5,g_wave=0.5,kappa_psi=0.5,Gamma_w=0.5)
    e=ZeroEvo(zs.k4(),Params(**base),seed=np.array([8.,4.,1.5,1.]))
    ev=[]
    for g in range(gens):
        e.step()
        if e.gen_closures.sum()>0: ev.append(g+1)
    return ev
def main():
    crab=[1969,1975,1981,1986,1989,1992,1996,1997,2000,2004,2006,2010,2013,2018]
    frb=[55946,55947,56021,56023,56026,56028,56030,56031,56033,56034]
    pol=[220,540,810,1120,1430,1730,2020,2330,2630]
    res={"pulsar glitches (yr)":score(crab),"FRB bursts (MJD)":score(frb),"CMB polarization peaks":score(pol)}
    try:
        ev=engine_stream(); res["engine gen_closures"]=score(ev)
    except Exception as ex:
        res["engine gen_closures"]=dict(error=str(ex))
    checks=[]
    for k in ("pulsar glitches (yr)","FRB bursts (MJD)","CMB polarization peaks"):
        r=res[k]; checks.append(dict(name=f"{k}: 不是两值游程", pass_=not (r["score"]>0.99 and r["run_ok"]),
                                     detail=f"score={r['score']:.2f}, run_ok={r['run_ok']}"))
    checks.append(dict(name="引擎 gen_closures 序列（参考）", pass_=True, detail=str(res.get("engine gen_closures"))))
    out=dict(candidates=res,checks=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES,"zero_world6_stream_search3.json"),"w") as f: json.dump(out,f,indent=1)
    print("试脉冲星/FRB/CMB极化 + 引擎record")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for k,v in res.items(): print(f"  {k}: {v}")
if __name__=="__main__": main()
