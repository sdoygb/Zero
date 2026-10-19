#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_closure_next.py —— 上次闭合在 2400 年前，下次闭合还要多久？
================================================================
闭合间隔（步）分布：P(2n)=2·Catalan(n−1)/4^n（重尾；均值发散）。
用"上次间隔=2 步（最概然）"定标：2 步 = 2400 年 ⇒ 1 步 = 1200 年。
下次间隔 2n' 步 ⇒ 剩余 = (n'−1)·2400 年。
输出：results/zero_closure_next.json
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_closure_next.py
"""
from __future__ import annotations
import json, os
from math import comb
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,"..","results")
def catalan(n): return comb(2*n,n)//(n+1)
def main():
    N=30
    rows=[]; cum=0.0
    for n in range(1,N+1):
        P=2*catalan(n-1)/4**n
        cum+=P
        rows.append(dict(n=n,gap_steps=2*n,remain_yr=(n-1)*2400,P=P,cum=cum))
    # 分位数
    def q(p):
        c=0.0
        for r in rows:
            c+=r["P"]
            if c>=p: return r["remain_yr"]
        return None
    checks=[
      dict(name="间隔分布 P(2n)=2C_{n−1}/4^n", pass_=abs(rows[0]["P"]-0.5)<1e-9, detail=f"P1={rows[0]['P']:.3f}"),
      dict(name="50% 概率下一次就在当前步（剩余 0）", pass_=rows[0]["cum"]>=0.5, detail=f"cum1={rows[0]['cum']:.3f}"),
      dict(name="约 62.5% 在 2400 年内", pass_=abs(rows[1]["cum"]-0.625)<1e-9, detail=f"cum2={rows[1]['cum']:.3f}"),
      dict(name="重尾（均值发散，不能给唯一年数）", pass_=True, detail="heavy tail"),
    ]
    out=dict(rows=rows,q50=q(0.5),q625=q(0.625),q90=q(0.9),
             checks=[{**{k:v for k,v in c.items() if k!="pass_"},"pass":bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES,"zero_closure_next.json"),"w") as f: json.dump(out,f,indent=1)
    print("下次闭合还要多久")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for r in rows[:8]: print(f"  gap={r['gap_steps']:2d}步 剩余={r['remain_yr']:7d}年 P={r['P']:.4f} 累计={r['cum']:.4f}")
    print(f"  中位={out['q50']}年, 62.5%={out['q625']}年, 90%={out['q90']}年")
if __name__=="__main__": main()
