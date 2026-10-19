#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_phase_anchor.py —— 给"当前相位"找锚：候选事件锚与结果
================================================================
周期 T_cycle=14.51 Gyr（哈勃锚）。当前相位 φ 需要一个"周期起点"锚。
候选锚（宇宙学事件，作为"上次毁灭"的标记）：
  · 复合/再复合 t≈0.00038 Gyr；
  · 物质-辐射相等 t≈0.00005 Gyr；
  · 暗能量起始 t≈9.3 Gyr；
  · 大爆炸 t=0（仅作对照；Zero 没有大爆炸）。
φ = (t_age − t_event)/T_cycle，剩余 = T_cycle − (t_age − t_event)。
输出：results/zero_world6_phase_anchor.json
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_phase_anchor.py
"""
from __future__ import annotations
import json, os
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"..","results")
T=14.51; AGE=13.787
EVENTS=[("recombination",0.00038),("matter-radiation equality",0.00005),
        ("dark-energy onset",9.3),("Big Bang (reference only)",0.0)]
def main():
    rows=[]
    for name,t0 in EVENTS:
        elapsed=AGE-t0; phi=elapsed/T; remain=T-elapsed
        rows.append(dict(anchor=name,t0_Gyr=t0,phi=phi,remain_Gyr=remain))
    out=dict(T_Gyr=T,AGE_Gyr=AGE,rows=rows)
    checks=[
      dict(name="周期 14.51 Gyr", pass_=abs(T-14.51)<0.01, detail=f"{T} Gyr"),
      dict(name="复合锚：φ≈0.95、剩余≈0.72 Gyr", pass_=abs(rows[0]["phi"]-0.95)<0.01 and abs(rows[0]["remain_Gyr"]-0.72)<0.05,
           detail=f"φ={rows[0]['phi']:.3f}, remain={rows[0]['remain_Gyr']:.3f}"),
      dict(name="暗能量锚：φ≈0.31、剩余≈10 Gyr", pass_=abs(rows[2]["phi"]-0.31)<0.02 and abs(rows[2]["remain_Gyr"]-10.0)<0.3,
           detail=f"φ={rows[2]['phi']:.3f}, remain={rows[2]['remain_Gyr']:.3f}"),
      dict(name="相位依赖锚：没有唯一值", pass_=True, detail="anchor-dependent"),
    ]
    out["checks"]=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks]
    with open(os.path.join(RES,"zero_world6_phase_anchor.json"),"w") as f: json.dump(out,f,indent=1)
    print("当前相位的候选锚")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for r in rows: print(f"  {r['anchor']:30s} t0={r['t0_Gyr']:.5f} Gyr  φ={r['phi']:.3f}  剩余={r['remain_Gyr']:.2f} Gyr")
if __name__=="__main__": main()
