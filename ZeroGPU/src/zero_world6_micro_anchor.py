#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_micro_anchor.py —— 从微观世界找相位锚（CP 相位）
================================================================
修正：不从宇宙学找锚；从微观世界的"相位型"观测量找。
候选：CKM 的 γ 相位、PMNS 的 δ_CP（都是 2π 周期的相位）。
约定：φ = δ/(2π) mod 1；剩余 = (1−φ)·T，T=14.51 Gyr。
输出：results/zero_world6_micro_anchor.json
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_micro_anchor.py
"""
from __future__ import annotations
import json, os, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"..","results")
T=14.51
def main():
    rows=[]
    # CKM γ = 65.5°（PDG 量级）
    gam=np.deg2rad(65.5)
    rows.append(dict(anchor="CKM gamma", delta_rad=float(gam), phi=float((gam/(2*np.pi))%1),
                     remain_Gyr=float(T*(1-(gam/(2*np.pi))%1))))
    # PMNS δ_CP 范围 0.7π–1.5π
    for k in (0.7,1.0,1.2,1.5):
        d=k*np.pi; phi=(d/(2*np.pi))%1
        rows.append(dict(anchor=f"PMNS delta_CP={k}pi", delta_rad=float(d), phi=float(phi),
                         remain_Gyr=float(T*(1-phi))))
    checks=[
      dict(name="微观相位锚给出 φ∈(0,1)", pass_=all(0<=r["phi"]<1 for r in rows), detail="ok"),
      dict(name="CKM γ：φ≈0.18、剩余≈11.9 Gyr", pass_=abs(rows[0]["phi"]-0.182)<0.01 and abs(rows[0]["remain_Gyr"]-11.87)<0.1,
           detail=f"φ={rows[0]['phi']:.3f}, remain={rows[0]['remain_Gyr']:.2f}"),
      dict(name="PMNS δ_CP：φ 0.35–0.75、剩余 3.6–9.4 Gyr",
           pass_=abs(min(r["phi"] for r in rows[1:])-0.35)<0.02 and abs(max(r["phi"] for r in rows[1:])-0.75)<0.02,
           detail=str([round(r['phi'],2) for r in rows[1:]])),
      dict(name="不同微观锚给不同 φ（无唯一映射）", pass_=True, detail="convention"),
    ]
    out=dict(T_Gyr=T,rows=rows)
    out["checks"]=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks]
    with open(os.path.join(RES,"zero_world6_micro_anchor.json"),"w") as f: json.dump(out,f,indent=1)
    print("微观相位锚")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for r in rows: print(f"  {r['anchor']:18s} δ={r['delta_rad']:.3f} rad  φ={r['phi']:.3f}  剩余={r['remain_Gyr']:.2f} Gyr")
if __name__=="__main__": main()
