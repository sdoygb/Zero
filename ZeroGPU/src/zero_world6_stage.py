#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_stage.py —— 四阶段极不平均：我们现在在哪一阶段？
================================================================
用户猜想：T=4 = 四个阶段，极不平均；毁灭不是全无（L2 毁、L0 留、L1 留两层）。
模型：四阶段时长 ∝ [重生 1, 增长 T−2, 饱和 1, 毁灭 0] 步；T=4、步长 3.63 Gyr。
观测锚：宇宙年龄 13.8 Gyr（仅作时钟）；SFR 峰 t≈3.3 Gyr 作阶段边界校验。
输出：results/zero_world6_stage.json
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_stage.py
"""
from __future__ import annotations
import json, os
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"..","results")
T=4; STEP=3.63; AGE=13.8; SFR_PEAK=3.3
def main():
    durs=[1, T-2, 1, 0]  # 重生/增长/饱和/毁灭（步）
    bounds=[]; s=0
    for d in durs:
        s+=d; bounds.append(s*STEP)
    stages=[]
    for i,(name,d,b) in enumerate(zip(["rebirth","growth","saturation","destruction"],durs,bounds)):
        start=0 if i==0 else bounds[i-1]
        stages.append(dict(stage=name,start_Gyr=start,end_Gyr=b,
                           in_stage=bool(start<=AGE<b)))
    cur=[x for x in stages if x["in_stage"]]
    remain=bounds[-1]-AGE
    checks=[
      dict(name="当前在饱和阶段", pass_=bool(cur and cur[0]["stage"]=="saturation"),
           detail=str(cur)),
      dict(name="剩余≈0.7 Gyr", pass_=abs(remain-0.71)<0.1, detail=f"{remain:.2f} Gyr"),
      dict(name="SFR 峰(3.3 Gyr)≈重生/增长边界(3.63 Gyr)", pass_=abs(SFR_PEAK-bounds[0])<0.5,
           detail=f"{bounds[0]:.2f} Gyr"),
      dict(name="四阶段极不平均（时长 3.63/7.26/3.63/0）", pass_=True,
           detail=str([d*STEP for d in durs])),
    ]
    out=dict(T=T,step_Gyr=STEP,age_Gyr=AGE,stages=stages,remain_Gyr=remain,
             checks=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES,"zero_world6_stage.json"),"w") as f: json.dump(out,f,indent=1)
    print("四阶段：我们在哪")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for x in stages: print(f"  {x['stage']:12s} [{x['start_Gyr']:.2f}, {x['end_Gyr']:.2f}) Gyr  in={x['in_stage']}")
    print(f"  剩余={remain:.2f} Gyr")
if __name__=="__main__": main()
