#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""16/8/4 + cycle16：beta 扫描（注意：beta 是 content n -> phase psi，不是 h -> psi）。

固定映射 seq；g_src=1, alpha=0.1；扫 beta={0,0.01,0.1,1.0}。
看 psi_coh/psi 周期4/h_E 是否稳定。
输出：results/zero_evo_1684_beta.json/.png
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


def run(beta, steps=2000):
    e = ZeroEvo1684Cycle16(zs.cycle(16),
                           Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05,
                                  alpha=0.1, g_src=1.0, beta=beta),
                           seed=np.ones(16))
    e.perm = MAPS["seq"]
    rows = []
    for g in range(1, steps+1):
        e.step()
        s = e.snapshot()
        rows.append({"gen": g, "q": s["q"], "q16": e.q16(), "h_E": s["h_energy"],
                     "psi_coh": s["psi_coh"], "n_std": s["n_std"]})
    return rows


def p4(a):
    a = np.asarray(a) - np.mean(a); f = np.abs(np.fft.rfft(a)); return float(f[len(a)//T])


def main():
    bs = [0.0, 0.01, 0.1, 1.0]
    res = {b: run(b) for b in bs}
    print("beta   q       h_E       psi_coh  n_std   p4psi   p4h")
    for b in bs:
        r = res[b][-200:]
        print(f"{b:.3f}  {np.mean([x['q'] for x in r]):.4f}  {np.mean([x['h_E'] for x in r]):.3e}  "
              f"{np.mean([x['psi_coh'] for x in r]):.4f}  {np.mean([x['n_std'] for x in r]):.4f}  "
              f"{p4([x['psi_coh'] for x in res[b]]):.3e}  {p4([x['h_E'] for x in res[b]]):.3e}")
    ck("无 NaN", all(np.all(np.isfinite([x['psi_coh'] for x in res[b]])) for b in bs), "")
    ck("psi_coh 有界（<1）", all(np.max([x['psi_coh'] for x in res[b]]) <= 1.0 + 1e-9 for b in bs), "")
    ck("beta>0 时 psi 有周期4", any(p4([x['psi_coh'] for x in res[b]]) > 1e-6 for b in bs if b > 0),
       str([round(p4([x['psi_coh'] for x in res[b]]),4) for b in bs]))
    ck("beta>0 时 h_E 仍在", all(np.mean([x['h_E'] for x in res[b][-200:]]) > 0 for b in bs if b > 0), "")

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    for b in bs:
        ax[0].plot([x['psi_coh'] for x in res[b]], lw=1, label=f"β={b}")
        ax[1].plot([x['h_E'] for x in res[b]], lw=1, label=f"β={b}")
    ax[0].set_title("(a) psi_coh"); ax[0].legend(fontsize=8)
    ax[1].set_title("(b) h_E"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_evo_1684_beta.png")
    fig.savefig(fp, dpi=130)
    out = {"results": {str(b): {"tail": {k: float(np.mean([x[k] for x in res[b][-200:]]))
                                       for k in ["q","q16","h_E","psi_coh","n_std"]},
                              "p4_psi": p4([x['psi_coh'] for x in res[b]]),
                              "p4_h": p4([x['h_E'] for x in res[b]]),
                              "rows": res[b]} for b in bs},
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "note": "beta couples content n -> psi, not h -> psi"}
    with open(os.path.join(RES, "zero_evo_1684_beta.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_evo_1684_beta.json"), "和", fp)


if __name__ == "__main__":
    main()
