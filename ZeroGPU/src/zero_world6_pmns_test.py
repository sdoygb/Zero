#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_pmns_test.py —— 预言 δ_CP<3π/4 与 PMNS 测量对比
================================================================
预言（θ_crit=3π/4）：δ_CP < 135°（0.75π）。
测量（NuFIT 6.0 近似值，NO）：best-fit ≈ 1.2π ≈ 216°；1σ ≈ [0.9π,1.5π]=[162°,270°]。
对照：θ_crit=π（δ<180°）、2π（无约束）。
输出：results/zero_world6_pmns_test.json
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_pmns_test.py
"""
from __future__ import annotations
import json, os
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"..","results")
PRED={"3pi/4 (mass zero)":135.0,"pi (double cover)":180.0,"2pi (cycle)":360.0}
OBS={"NO best-fit":216.0,"NO 1sigma low":162.0,"NO 1sigma high":270.0,
     "IO best-fit":270.0}
def main():
    rows=[]
    for name,b in PRED.items():
        inside = OBS["NO 1sigma low"] <= b <= OBS["NO 1sigma high"]
        gap=max(0.0, OBS["NO 1sigma low"]-b)
        rows.append(dict(pred=name,boundary_deg=b,inside_1sigma=inside,gap_deg=float(gap)))
    checks=[
      dict(name="θ_crit=3π/4：预言 δ<135°", pass_=rows[0]["boundary_deg"]==135.0, detail="135°"),
      dict(name="与 NO 1σ [162°,270°] 不相交（tension）", pass_=not rows[0]["inside_1sigma"],
           detail=f"gap={rows[0]['gap_deg']:.0f}°"),
      dict(name="θ_crit=π：与 1σ 有重叠（边缘）", pass_=rows[1]["inside_1sigma"],
           detail=f"boundary={rows[1]['boundary_deg']}°"),
      dict(name="θ_crit=2π：无约束", pass_=rows[2]["boundary_deg"]==360.0, detail="no constraint"),
    ]
    out=dict(pred=PRED,obs=OBS,rows=rows,checks=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES,"zero_world6_pmns_test.json"),"w") as f: json.dump(out,f,indent=1)
    print("预言 δ_CP<3π/4 vs PMNS")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for r in rows: print(f"  {r['pred']:20s} boundary={r['boundary_deg']:.0f}° inside1σ={r['inside_1sigma']} gap={r['gap_deg']:.0f}°")
if __name__=="__main__": main()
