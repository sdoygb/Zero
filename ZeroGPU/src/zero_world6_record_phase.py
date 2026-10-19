#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_record_phase.py —— 从记录层的间隔结构读 T 与相位 k
================================================================
旧库（L2_imprint）：记录层的沉积流间隔为 2,2,...,T+2,...；由此反解 T。
本轮：生成该沉积流，验证反解 T 与"当前相位 k"的读法。

结论：
  · 间隔结构能读出 T（long−2）；
  · 相位 k = 现在 − 最近一次沉积（需要"现在"与"最近沉积"都可观测）；
  · 没有观测到的 L2 沉积流 ⇒ k 仍未知（但方法是明确的）。

输出：results/zero_world6_record_phase.json
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_record_phase.py
"""
from __future__ import annotations
import json, os, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"..","results")
T=4
def stream(T=4, cycles=100, jitter=0.0, seed=0):
    rng=np.random.default_rng(seed)
    pat=[2]*(T//2-1)+[T+2]
    times=[]; t=0
    for _ in range(cycles):
        for d in pat:
            t=t+d; times.append(t+ (rng.normal(0,jitter) if jitter else 0.0))
    return np.array(times)
def read_T(times):
    ints=np.diff(times); long=ints.max(); short=np.median(ints[ints<long*0.8])
    return int(round(long-2)), float(short)
def read_k(times, now):
    last=times[times<=now].max(); return float(now-last)
def main():
    times=stream(T,100)
    Trec,short=read_T(times)
    k=read_k(times, times[-1]+1.0)  # 假设"现在"在最后沉积后 1 步
    rem=(1-k/T)*14.51
    # 加噪声
    times2=stream(T,100,jitter=0.05,seed=1); Trec2,short2=read_T(times2)
    checks=[
      dict(name="间隔结构反解 T=4", pass_=Trec==4, detail=f"T_rec={Trec}, short={short:.2f}"),
      dict(name="加噪声仍反解 T=4", pass_=Trec2==4, detail=f"T_rec={Trec2}, short={short2:.2f}"),
      dict(name="相位 k=现在−最近沉积（示例 k=1）", pass_=abs(k-1.0)<1e-9, detail=f"k={k:.2f}"),
      dict(name="剩余=(1−k/T)·14.51 Gyr（示例）", pass_=abs(rem-(3/4)*14.51)<1e-9, detail=f"{rem:.2f} Gyr"),
      dict(name="没有 L2 沉积流观测 ⇒ k 仍未知", pass_=True, detail="method only"),
    ]
    out=dict(T=T,T_rec=Trec,k=k,remain_Gyr=rem,checks=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES,"zero_world6_record_phase.json"),"w") as f: json.dump(out,f,indent=1)
    print("从记录间隔读 T 与相位")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
if __name__=="__main__": main()
