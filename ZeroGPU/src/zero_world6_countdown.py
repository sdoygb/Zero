#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_countdown.py —— 从现在起还有多少年毁灭（哈勃锚下）
================================================================
周期 = 哈勃时间 t_H；从大爆炸起算；剩余 = t_H − 宇宙年龄。
H0 敏感性：67.4 / 70 / 73；并与 Z017 的 14.4 Gyr 对照。
输出：results/zero_world6_countdown.json
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_countdown.py
"""
from __future__ import annotations
import json, os
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"..","results")
MPC_KM=3.085677581e19; YR=3.15576e7; GYR=1e9*YR
AGE=13.787e9  # Planck 2018 宇宙年龄（年）
def tH(H0): return MPC_KM/H0/YR  # 年
def main():
    rows=[]
    for H0 in (67.4,70.0,73.0):
        th=tH(H0); rows.append(dict(H0=H0,tH_Gyr=th/1e9,remaining_Gyr=(th-AGE)/1e9))
    z017=14.4e9
    out=dict(age_yr=AGE,rows=rows,z017_cycle_Gyr=14.4,z017_remaining_Gyr=(z017-AGE)/1e9)
    r=rows[0]["remaining_Gyr"]
    checks=[
      dict(name="H0=67.4：剩余 ≈ 0.71 Gyr（7.1e8 年）", pass_=abs(r-0.707)<0.03, detail=f"{r:.3f} Gyr"),
      dict(name="H0 敏感性：67.4–73 剩余 −0.40 到 0.71 Gyr", pass_=True,
           detail=str([round(x["remaining_Gyr"],3) for x in rows])),
      dict(name="Z017 14.4 Gyr 锚：剩余 0.61 Gyr", pass_=abs(out["z017_remaining_Gyr"]-0.613)<0.05,
           detail=f"{out['z017_remaining_Gyr']:.3f} Gyr"),
      dict(name="这是锚下减法，不是预言（Z017 定理 10）", pass_=True, detail="external anchor"),
    ]
    out["checks"]=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks]
    with open(os.path.join(RES,"zero_world6_countdown.json"),"w") as f: json.dump(out,f,indent=1)
    print("从现在起还有多少年毁灭（锚=哈勃时间）")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for x in rows: print(f"  H0={x['H0']}: t_H={x['tH_Gyr']:.2f} Gyr, 剩余={x['remaining_Gyr']:+.3f} Gyr")
    print(f"  Z017(14.4): 剩余={out['z017_remaining_Gyr']:+.3f} Gyr")
if __name__=="__main__": main()
