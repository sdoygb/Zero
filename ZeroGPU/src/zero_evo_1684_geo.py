#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""16/8/4 + cycle16：固定物理映射 seq，开几何源 g_src。

映射由 16 槽载体的字典序 index=4*site+fiber 固定（seq），不再当自由置换。
扫 g_src={0,0.1,1.0}，看 h_E/psi 是否响应、是否出现周期 4、是否稳定。
输出：results/zero_evo_1684_geo.json/.png
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


def run(g_src, steps=2000):
    # 固定映射：16 槽载体 index=4*site+fiber
    e = ZeroEvo1684Cycle16(zs.cycle(16),
                           Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05,
                                  alpha=0.1, g_src=g_src, Gamma_h=0.0),
                           seed=np.ones(16))
    e.perm = MAPS["seq"]
    rows = []
    for g in range(1, steps+1):
        e.step()
        s = e.snapshot()
        rows.append({"gen": g, "q": s["q"], "q16": e.q16(),
                     "h_E": s["h_energy"], "h_max": s["h_max"],
                     "psi_coh": s["psi_coh"], "n_std": s["n_std"],
                     "fiedler": s["fiedler"]})
    return rows


def p4(a):
    a = np.asarray(a) - np.mean(a); f = np.abs(np.fft.rfft(a)); return float(f[len(a)//T])


def main():
    gs = [0.0, 0.1, 1.0]
    res = {g: run(g) for g in gs}
    print("g_src  q       q16     h_E(mean)  h_max   psi_coh  n_std   p4h     p4psi")
    for g in gs:
        r = res[g][-200:]
        print(f"{g:.2f}   {np.mean([x['q'] for x in r]):.4f}  {np.mean([x['q16'] for x in r]):.4f}  "
              f"{np.mean([x['h_E'] for x in r]):.3e}  {np.max([x['h_max'] for x in r]):.3e}  "
              f"{np.mean([x['psi_coh'] for x in r]):.4f}  {np.mean([x['n_std'] for x in r]):.4f}  "
              f"{p4([x['h_E'] for x in res[g]]):.3e}  {p4([x['psi_coh'] for x in res[g]]):.3e}")
    ck("g_src=0 时 h_E=0（对照）", abs(np.mean([x['h_E'] for x in res[0.0][-200:]])) < 1e-15, "")
    ck("g_src>0 时 h_E>0（几何响应）", all(np.mean([x['h_E'] for x in res[g][-200:]]) > 0 for g in gs[1:]),
       str([round(float(np.mean([x['h_E'] for x in res[g][-200:]])),6) for g in gs[1:]]))
    ck("所有 g_src 有界（h_max<1e6）", all(np.max([x['h_max'] for x in res[g]]) < 1e6 for g in gs),
       str([round(float(np.max([x['h_max'] for x in res[g]])),6) for g in gs]))
    ck("所有 g_src 无 NaN", all(np.all(np.isfinite([x['h_E'] for x in res[g]])) for g in gs), "")
    ck("g_src>0 时 h_E 或 psi 有周期4", any(p4([x['h_E'] for x in res[g]]) > 1e-6 or
                                         p4([x['psi_coh'] for x in res[g]]) > 1e-6 for g in gs[1:]),
       str([round(p4([x['h_E'] for x in res[g]]),3) for g in gs[1:]]))

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    for g in gs:
        ax[0].plot([x['h_E'] for x in res[g]], lw=1, label=f"g_src={g}")
        ax[1].plot([x['psi_coh'] for x in res[g]], lw=1, label=f"g_src={g}")
    ax[0].set_title("(a) h_E"); ax[0].legend(fontsize=8)
    ax[1].set_title("(b) psi_coh"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_evo_1684_geo.png")
    fig.savefig(fp, dpi=130)
    out = {"results": {str(g): {"tail": {k: float(np.mean([x[k] for x in res[g][-200:]]))
                                       for k in ["q","q16","h_E","h_max","psi_coh","n_std","fiedler"]},
                              "p4_h": p4([x['h_E'] for x in res[g]]),
                              "p4_psi": p4([x['psi_coh'] for x in res[g]]),
                              "rows": res[g]} for g in gs},
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "mapping fixed seq; g_src scan; no dimensionful anchor"}
    with open(os.path.join(RES, "zero_evo_1684_geo.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_evo_1684_geo.json"), "和", fp)


if __name__ == "__main__":
    main()
