#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_closure_events.py —— 闭合事件：从零走到零
================================================================
Zero 的闭合：词和 0 → 走 ±1 → 回到 0 即闭合。
本脚本算：
  · 第 2n 步回到 0 的路径数 R_n = C(2n,n)（全部闭合）；
  · 首次在第 2n 步回到 0 的路径数 F_n = 2·Catalan(n−1)（n≥1）；
  · 更新恒等式 R_n = Σ_{k=1..n} F_k R_{n−k}；
  · 返回概率 R_n/2^{2n}；
  · T=4（n=2）的闭合事件数。
输出：results/zero_closure_events.json
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_closure_events.py
"""
from __future__ import annotations
import json, os
from math import comb
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"..","results")
def catalan(n): return comb(2*n,n)//(n+1)
def main():
    N=8
    rows=[]
    for n in range(0,N+1):
        R=comb(2*n,n)
        F=(2*catalan(n-1) if n>=1 else 1)
        rows.append(dict(n=n,step=2*n,R=R,F=F,prob=R/2**(2*n)))
    # 更新恒等式
    ok=True
    for n in range(1,N+1):
        s=sum(rows[k]["F"]*rows[n-k]["R"] for k in range(1,n+1))
        # F_0=1, R_0=1 的项：F_0*R_n 不算首次；用 n>=1 的和
        if s!=rows[n]["R"]: ok=False
    T=4; nT=T//2
    R_tot=sum(rows[n]["R"] for n in range(0,nT+1))
    F_tot=sum(rows[n]["F"] for n in range(1,nT+1))
    checks=[
      dict(name="第 2n 步返回数 R_n=C(2n,n)", pass_=rows[1]["R"]==2 and rows[2]["R"]==6,
           detail=str([r["R"] for r in rows[:5]])),
      dict(name="首次返回 F_n=2·Catalan(n−1)", pass_=rows[1]["F"]==2 and rows[2]["F"]==2,
           detail=str([r["F"] for r in rows[:5]])),
      dict(name="更新恒等式 R_n=Σ F_k R_{n−k}", pass_=ok, detail="holds n=1..8"),
      dict(name="T=4：9 次返回（含平凡 n=0）、4 次首次返回", pass_=(R_tot==9 and F_tot==4),
           detail=f"R_tot={R_tot}, F_tot={F_tot}"),
    ]
    out=dict(rows=rows,T=T,R_tot=R_tot,F_tot=F_tot,
             checks=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES,"zero_closure_events.json"),"w") as f: json.dump(out,f,indent=1)
    print("闭合事件：从零走到零")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print("  n  step  R_n   F_n   prob")
    for r in rows: print(f"  {r['n']}  {r['step']:3d}  {r['R']:4d}  {r['F']:3d}  {r['prob']:.4f}")
    print(f"  T=4: R_tot={R_tot}, F_tot={F_tot}")
if __name__=="__main__": main()
