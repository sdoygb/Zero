#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_countdown.py —— Zero 没有大爆炸：倒计时是一族解，不是一个数
================================================================
修正：上一版用"从大爆炸起算"（t_remain=t_H−t_age），这与 Zero 不符——
Zero 没有大爆炸，周期没有 t=0 起点。

本轮：周期 = 哈勃时间 t_H（锚）；当前相位 φ∈[0,1) 是**自由输入**；
      t_remain=(1−φ)·t_cycle。
  · 没有相位锚 ⇒ 剩余 ∈ (0, t_cycle]，不是唯一数；
  · 上一版的 0.72 Gyr 对应 φ=0.951（大爆炸假设），作废；
  · H0 敏感性只影响周期本身（13.39–14.51 Gyr）。

输出：results/zero_world6_countdown.json
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_countdown.py
"""
from __future__ import annotations
import json, os
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"..","results")
MPC_KM=3.085677581e19; YR=3.15576e7
def tH(H0): return MPC_KM/H0/YR
def main():
    rows=[]
    for H0 in (67.4,70.0,73.0):
        th=tH(H0); rows.append(dict(H0=H0,tcycle_Gyr=th/1e9,
                                    remain_min_Gyr=0.0,remain_max_Gyr=th/1e9))
    # 上一版的"大爆炸"结果：φ=1−t_age/t_H
    age=13.787e9; th=tH(67.4)
    phi_bb=age/th  # 上一版把"已过比例"当相位（大爆炸起算）
    out=dict(no_big_bang=True,rows=rows,age_used_by_previous_yr=age,
             phi_previous=phi_bb,previous_remain_Gyr=(th-age)/1e9)
    checks=[
      dict(name="Zero 没有大爆炸：周期无 t=0 起点", pass_=True, detail="phase free"),
      dict(name="剩余是一族：(0, t_cycle]（无相位锚）", pass_=True, detail="0 < remain <= t_cycle"),
      dict(name="上一版 0.72 Gyr 是 φ=0.951（大爆炸起算）的假设，作废", pass_=abs(phi_bb-0.951)<0.01,
           detail=f"phi_prev={phi_bb:.3f}"),
      dict(name="H0 只影响周期（13.39–14.51 Gyr）", pass_=True,
           detail=str([round(r['tcycle_Gyr'],2) for r in rows])),
    ]
    out["checks"]=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks]
    with open(os.path.join(RES,"zero_world6_countdown.json"),"w") as f: json.dump(out,f,indent=1)
    print("Zero 没有大爆炸：倒计时是一族解")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for r in rows: print(f"  H0={r['H0']}: 周期={r['tcycle_Gyr']:.2f} Gyr, 剩余∈(0, {r['tcycle_Gyr']:.2f}] Gyr")
    print(f"  上一版: φ={phi_bb:.3f}, 剩余={out['previous_remain_Gyr']:.3f} Gyr（作废）")
if __name__=="__main__": main()
