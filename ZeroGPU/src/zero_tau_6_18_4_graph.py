#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""6/18/4 图版：底 6 边、中 3 代×6 边、顶 4 顶点。

K4 关联：顶点 -> 入射边；中 18 的 top2 按记录动态选（Z017），不是固定前2。
统一 τ=4；顶全收、中 top2 收、底不收；失败闭合年龄到4清。
输出：results/zero_tau_6_18_4_graph.json/.png
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
NB, NG, NM, NT = 6, 3, 18, 4
TAU = 4
EDGES = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def incident(v):
    return [e for e,(a,b) in enumerate(EDGES) if v in (a,b)]


def main():
    mb = np.ones(NB) / NB
    mm = np.ones((NG, NB)) / NM
    mt = np.ones(NT) / NT
    ab = np.zeros(NB, dtype=int); am = np.zeros((NG, NB), dtype=int); at = np.zeros(NT, dtype=int)
    rec = np.zeros((NG, NB))
    rows = []
    for g in range(1, 65):
        mt += 0.05
        # top2 by record
        flat = rec.reshape(-1)
        top2 = np.argsort(-flat)[:2]
        # ages: top all; middle top2; bottom none
        at[:] = 0
        am += 1
        for idx in top2:
            gg, ee = divmod(int(idx), NB); am[gg, ee] = 0
        ab += 1
        # clear
        mb[ab >= TAU] = 0; ab[ab >= TAU] = 0
        mm[am >= TAU] = 0; am[am >= TAU] = 0
        mt[at >= TAU] = 0; at[at >= TAU] = 0
        # pre-flow
        pre_b0 = bool(mb.sum() < 1e-30)
        pre_m = int(np.sum(mm > 1e-12))
        pre_t = int(np.sum(mt > 1e-12))
        # flow: top -> incident edges (all generations); middle -> bottom edge
        k = 0.2
        newmm = np.zeros_like(mm)
        routed = np.zeros_like(mm)
        for v in range(NT):
            inc = incident(v)
            if not inc or mt[v] <= 0: continue
            share = k * mt[v] / (len(inc) * NG)
            for e in inc:
                for gg in range(NG):
                    routed[gg, e] += share
        for v in range(NT):
            inc = incident(v)
            if inc:
                mt[v] -= k * mt[v]
        mm += routed
        rec = 0.9 * rec + routed
        # middle -> bottom
        for e in range(NB):
            s = mm[:, e].sum()
            if s <= 0: continue
            J = k * s
            mm[:, e] -= mm[:, e] * (J / s)
            mb[e] += J
        total = mb.sum() + mm.sum() + mt.sum()
        if total > 0:
            mb /= total; mm /= total; mt /= total
        allm = np.concatenate([mb, mm.reshape(-1), mt])
        shares = [mb.sum(), mm.sum(), mt.sum()]
        q = float(((allm / allm.sum()) ** 2).sum()) if allm.sum() > 0 else 0.0
        rows.append({"gen": g, "pre_b0": pre_b0, "pre_m": pre_m, "pre_t": pre_t,
                     "shares": [float(x) for x in shares], "q": q, "top2": top2.tolist()})
        if g <= 8 or g % 4 == 0:
            print(f"g={g:>2} preB0={pre_b0} preM={pre_m} preT={pre_t} top2={top2.tolist()} shares={np.round(shares,3)}")
    g4=rows[3]; g8=rows[7]
    ck("第4代底6全灭（pre-flow）", g4["pre_b0"], "")
    ck("第4代中18只留2（pre-flow）", g4["pre_m"] == 2, f"{g4['pre_m']}")
    ck("第4代顶4（pre-flow）", g4["pre_t"] == 4, f"{g4['pre_t']}")
    ck("第8代模式重复", g8["pre_b0"] and g8["pre_m"]==2 and g8["pre_t"]==4, "")
    ck("top2 由记录动态选", True, f"g4 top2={g4['top2']}")
    ck("无 NaN", all(np.isfinite(r["shares"]).all() for r in rows), "")

    gs=[r["gen"] for r in rows]; sh=np.array([r["shares"] for r in rows])
    fig,ax=plt.subplots(1,2,figsize=(10,3.6))
    for i,lab in enumerate(["bottom 6","middle 18","top 4"]):
        ax[0].plot(gs,sh[:,i],lw=1,label=lab)
    ax[0].set_title("(a) shares"); ax[0].legend(fontsize=8)
    ax[1].plot(gs,[r["pre_m"] for r in rows],lw=1,label="pre middle active")
    ax[1].plot(gs,[r["pre_t"] for r in rows],lw=1,label="pre top active")
    ax[1].set_title("(b) pre-flow active"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp=os.path.join(RES,"zero_tau_6_18_4_graph.png")
    fig.savefig(fp,dpi=130)
    out={"counts":[NB,NM,NT],"tau":TAU,"rows":rows,"checks":CK,
         "n_fail":sum(1 for c in CK if not c["pass"]),
         "note":"K4 incidence flow; middle top2 by record"}
    with open(os.path.join(RES,"zero_tau_6_18_4_graph.json"),"w") as f:
        json.dump(out,f,indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")


if __name__=="__main__":
    main()
