#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_stream_search.py —— 哪条观测序列是 L2 沉积流？
================================================================
判据（旧库 L2_imprint）：间隔只取两个值 {s, L=T+2}，且 s 的游程 = T/2−1。
候选：CMB 声学峰、夸克/轻子对数质量间隔、合成控制流（T=4）。
输出：results/zero_world6_stream_search.json
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_stream_search.py
"""
from __future__ import annotations
import json, os, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"..","results")
def pattern_score(times):
    t=np.sort(np.asarray(times,float)); ints=np.diff(t)
    s=float(np.min(ints)); L=float(np.max(ints))
    if L<=s*1.05: return dict(s=s,L=L,T=None,score=0.0,note="one interval value")
    in_s=np.abs(ints-s)/s<0.1; in_L=np.abs(ints-L)/L<0.1
    score=float(np.mean(in_s|in_L))
    T=int(round(L-2))
    # 游程检查：s 的游程应 = T/2−1
    runs=[]; cur=0
    for x in in_s:
        if x: cur+=1
        elif cur>0: runs.append(cur); cur=0
    if cur>0: runs.append(cur)
    run_ok=bool(len(runs)>0 and abs(np.median(runs)-(T/2-1))<1)
    return dict(s=s,L=L,T=T,score=score,run_ok=run_ok,note="")
def main():
    cands={
      "CMB peaks": np.array([220,540,810,1120,1430,1730,2020,2330,2630.]),
      "quark log10 masses": np.sort(np.log10([2.16,4.67,93.4,1270.,4180.,172760.])),
      "lepton log10 masses": np.sort(np.log10([0.51099895,105.6583755,1776.86,0.001,0.0086,0.05])),
    }
    # 合成控制：T=4 的沉积流
    pat=[2]*(4//2-1)+[6]; t=[]; tt=0
    for _ in range(50):
        for d in pat: tt+=d; t.append(tt)
    res={k:pattern_score(v) for k,v in cands.items()}
    ctrl=pattern_score(t)
    checks=[
      dict(name="合成控制流（T=4）通过：s=2,L=6,游程=1", pass_=ctrl["score"]>0.99 and ctrl["run_ok"],
           detail=f"score={ctrl['score']:.2f}, s={ctrl['s']}, L={ctrl['L']}, run_ok={ctrl['run_ok']}"),
      dict(name="CMB 峰：间隔近似两值但游程不对（失败）",
           pass_=not (res["CMB peaks"]["score"]>0.99 and res["CMB peaks"]["run_ok"]),
           detail=f"score={res['CMB peaks']['score']:.2f}, run_ok={res['CMB peaks']['run_ok']}"),
      dict(name="夸克质量间隔：不是两值游程（失败）", pass_=res["quark log10 masses"]["score"]<0.9,
           detail=f"score={res['quark log10 masses']['score']:.2f}"),
      dict(name="轻子质量间隔：不是两值游程（失败）", pass_=res["lepton log10 masses"]["score"]<0.9,
           detail=f"score={res['lepton log10 masses']['score']:.2f}"),
    ]
    out=dict(control=ctrl,candidates=res,checks=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES,"zero_world6_stream_search.json"),"w") as f: json.dump(out,f,indent=1)
    print("找 L2 沉积流")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for k,v in res.items(): print(f"  {k}: s={v['s']:.3f} L={v['L']:.3f} T={v['T']} score={v['score']:.2f} run_ok={v.get('run_ok')}")
if __name__=="__main__": main()
