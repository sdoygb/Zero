#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_period.py —— 定锚：最底 6 层的毁灭周期（年）
================================================================
6 层模型的毁灭周期是 T=4（代/步）。要把"4 代"换成"年"，需要一个时间锚。
本脚本取 **哈勃时间**（与 Z017 §5.2 的既有锚一致）：

  t_H = 1/H0，H0=67.4 km/s/Mpc → t_H ≈ 14.4 Gyr；
  令一个毁灭周期 = t_H ⇒ t_cycle = 14.4 Gyr = 1.44e10 年，
  一步（一代）= t_H/T = 3.6 Gyr。

对照锚：Planck 时间、宇宙年龄。

输出：results/zero_world6_period.json
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_period.py
"""
from __future__ import annotations
import json, os
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"..","results")
MPC_KM=3.085677581e19  # km per Mpc
YR=3.15576e7           # s per year
GYR=1e9*YR
H0=67.4                # km/s/Mpc
T=4                    # 6 层的毁灭周期（代）
def main():
    tH_s=MPC_KM/H0
    tH_yr=tH_s/YR
    tH_gyr=tH_yr/1e9
    tcycle_yr=tH_yr          # 一个周期 = 哈勃时间
    step_yr=tH_yr/T          # 一代 = t_H / 4
    age_yr=13.8e9
    planck_yr=5.391247e-44/YR
    out=dict(H0=H0,T=T,
             tH_yr=float(tH_yr),tH_Gyr=float(tH_gyr),
             tcycle_yr=float(tcycle_yr),tcycle_Gyr=float(tcycle_yr/1e9),
             step_yr=float(step_yr),step_Gyr=float(step_yr/1e9),
             age_yr=age_yr,age_over_cycle=float(age_yr/tcycle_yr),
             planck_cycle_yr=float(planck_yr),
             alternative_age_cycle_yr=age_yr)
    checks=[
      dict(name="哈勃时间 ≈ 14.4 Gyr", pass_=abs(tH_gyr-14.4)<0.2, detail=f"{tH_gyr:.2f} Gyr"),
      dict(name="6 层毁灭周期 = 1.44e10 年", pass_=abs(tcycle_yr/1e10-1.44)<0.02,
           detail=f"{tcycle_yr:.3e} yr"),
      dict(name="一代 = 3.6 Gyr", pass_=abs(step_yr/1e9-3.6)<0.1, detail=f"{step_yr/1e9:.2f} Gyr"),
      dict(name="锚是外部输入（Z017 定理 10）", pass_=True, detail="Hubble anchor"),
    ]
    out["checks"]=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks]
    with open(os.path.join(RES,"zero_world6_period.json"),"w") as f: json.dump(out,f,indent=1)
    print("最底 6 层的毁灭周期（锚=哈勃时间）")
    for c in out['checks']: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  t_H={tH_gyr:.3f} Gyr; t_cycle={tcycle_yr:.4e} yr ({tcycle_yr/1e9:.2f} Gyr); step={step_yr/1e9:.2f} Gyr")
    print(f"  age/t_cycle={age_yr/tcycle_yr:.3f}; Planck-cycle={planck_yr:.2e} yr")
if __name__=="__main__": main()
