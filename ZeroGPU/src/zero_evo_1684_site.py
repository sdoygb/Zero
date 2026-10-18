#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""16/8/4 + cycle16：site-resolved 映射敏感性与 geometry/wave 谱。

五种 class->顶点双射；各 2000 代；比较：
  - 聚合量 q16/content q/份额（预期不变）
  - site-resolved：n 向量、top1/top2（预期随映射变）
  - geometry/wave：h_energy/psi_coh/fiedler 及其周期 4
输出：results/zero_evo_1684_site.json/.png
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402
from zero_evo import Params  # noqa: E402
from zero_evo_1684_cycle16 import ZeroEvo1684Cycle16, MAPS  # noqa: E402

RES = os.path.join(HERE, "..", "results")
T = 4
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def run(mapping, steps=2000):
    e = ZeroEvo1684Cycle16(zs.cycle(16), Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05, alpha=0.1),
                           seed=np.ones(16))
    e.perm = MAPS[mapping]
    rows = []
    for g in range(1, steps+1):
        e.step()
        s = e.snapshot()
        rows.append({"gen": g, "q": s["q"], "q16": e.q16(),
                     "h_E": s["h_energy"], "psi_coh": s["psi_coh"],
                     "fiedler": s["fiedler"], "top1": s["top1"], "top2": s["top2"],
                     "low": s["low"], "high": s["high"]})
    return {"rows": rows, "n": e.n.copy(), "h": e.h.copy(), "psi": e.psi.copy()}


def p4(a):
    a = np.asarray(a) - np.mean(a); f = np.abs(np.fft.rfft(a)); return float(f[len(a)//T])


def main():
    names = list(MAPS)
    res = {m: run(m) for m in names}
    print("map        q       q16     h_E       psi_coh   top1 top2")
    for m in names:
        r = res[m]["rows"][-100:]
        print(f"{m:10s} {np.mean([x['q'] for x in r]):.4f}  {np.mean([x['q16'] for x in r]):.4f}  "
              f"{np.mean([x['h_E'] for x in r]):.3e}  {np.mean([x['psi_coh'] for x in r]):.4f}  "
              f"{r[-1]['top1']:>4} {r[-1]['top2']:>4}")
    q16s = [np.mean([x["q16"] for x in res[m]["rows"][-100:]]) for m in names]
    qs = [np.mean([x["q"] for x in res[m]["rows"][-100:]]) for m in names]
    ck("q16 聚合不随映射变", max(q16s)-min(q16s) < 1e-9, f"spread={max(q16s)-min(q16s):.1e}")
    ck("content q 聚合不随映射变", max(qs)-min(qs) < 1e-6, f"spread={max(qs)-min(qs):.1e}")
    # site-resolved: n 向量是否不同
    n0 = res["seq"]["n"]
    dn = [float(np.linalg.norm(res[m]["n"] - n0)) for m in names]
    ck("n 向量随映射改变（site-resolved 不不变）", max(dn) > 1e-6, f"max Δ={max(dn):.3e}")
    # top labels
    tops = [(res[m]["rows"][-1]["top1"], res[m]["rows"][-1]["top2"]) for m in names]
    print("top labels", dict(zip(names, tops)))

    # 谱：h_E / psi_coh / q16 的周期4
    print("p4: h_E/psi_coh/q16")
    for m in names:
        r = res[m]["rows"]
        print(f"{m:10s} {p4([x['h_E'] for x in r]):.3e}  {p4([x['psi_coh'] for x in r]):.3e}  "
              f"{p4([x['q16'] for x in r]):.3f}")
    ck("q16 周期4仍在", all(p4([x['q16'] for x in res[m]['rows']]) > 1e-6 for m in names), "")
    # 画 q16 与 h_E/psi_coh（seq）
    g = np.arange(1, 2001)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].plot(g, [x["q16"] for x in res["seq"]["rows"]], lw=1, label="q16")
    ax[0].plot(g, [x["q"] for x in res["seq"]["rows"]], lw=1, label="content q")
    ax[0].set_title("(a) q"); ax[0].legend(fontsize=8)
    ax[1].plot(g, [x["h_E"] for x in res["seq"]["rows"]], lw=1, label="h_E")
    ax[1].plot(g, [x["psi_coh"] for x in res["seq"]["rows"]], lw=1, label="psi_coh")
    ax[1].set_title("(b) geometry/wave"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_evo_1684_site.png")
    fig.savefig(fp, dpi=130)
    out = {"q16": q16s, "q": qs, "dn": dn, "tops": dict(zip(names, tops)),
           "p4": {m: [p4([x['h_E'] for x in res[m]['rows']]),
                      p4([x['psi_coh'] for x in res[m]['rows']]),
                      p4([x['q16'] for x in res[m]['rows']])] for m in names},
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "site-resolved and geometry/wave spectra; 16/8/4 fixed"}
    with open(os.path.join(RES, "zero_evo_1684_site.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_evo_1684_site.json"), "和", fp)


if __name__ == "__main__":
    main()
