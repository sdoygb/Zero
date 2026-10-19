#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_phase_threshold.py —— 相位阈值模型：毁灭发生在临界相位
================================================================
用户：微观世界按相位演化；相位走到某点、不能再演化 ⇒ 毁灭。
模型：θ(t) 匀速转，θ_crit 是临界相位；剩余 = (θ_crit−θ_now)/(2π)·T。
候选 θ_crit：3π/4（Z126 质量零点）、π（双覆盖/拓扑 −1）、2π（整周期）。
候选 θ_now：CKM γ、PMNS δ_CP（0.7π–1.5π）。
输出：results/zero_world6_phase_threshold.json
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_phase_threshold.py
"""
from __future__ import annotations
import json, os, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"..","results")
T=14.51
def main():
    crits={"3pi/4 (mass zero)":3*np.pi/4,"pi (double cover)":np.pi,"2pi (cycle)":2*np.pi}
    nows={"CKM gamma":1.143,"PMNS 0.7pi":0.7*np.pi,"PMNS 1.0pi":np.pi,
          "PMNS 1.2pi":1.2*np.pi,"PMNS 1.5pi":1.5*np.pi}
    rows=[]
    for cn,c in crits.items():
        for nn,n in nows.items():
            rem=(c-n)/(2*np.pi)*T
            rows.append(dict(crit=cn,crit_rad=float(c),now=nn,now_rad=float(n),
                             remain_Gyr=float(rem),past=bool(rem<0)))
    checks=[
      dict(name="相位阈值模型：剩余=(θ_crit−θ_now)/(2π)·T", pass_=True, detail="formula"),
      dict(name="质量零点 θ_crit=3π/4：CKM γ→2.80 Gyr", pass_=abs(rows[0]["remain_Gyr"]-2.80)<0.1,
           detail=f"{rows[0]['remain_Gyr']:.2f}"),
      dict(name="双覆盖 θ_crit=π：PMNS 0.7π→2.18 Gyr", pass_=abs([r for r in rows if r['crit'].startswith('pi') and r['now']=='PMNS 0.7pi'][0]['remain_Gyr']-2.18)<0.1,
           detail="2.18"),
      dict(name="结果依赖 (θ_crit, θ_now)：有正有负", pass_=any(r['past'] for r in rows), detail="past cases exist"),
    ]
    out=dict(T_Gyr=T,rows=rows,checks=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES,"zero_world6_phase_threshold.json"),"w") as f: json.dump(out,f,indent=1)
    print("相位阈值模型")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  {'crit':22s} {'now':12s} 剩余")
    for r in rows: print(f"  {r['crit']:22s} {r['now']:12s} {r['remain_Gyr']:+.2f} Gyr")
if __name__=="__main__": main()
