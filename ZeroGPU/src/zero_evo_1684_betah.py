#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""新增 h->psi 项：beta_h 扫描。

16/8/4 + cycle16，固定 seq；beta=0（content->psi 关），
只开 beta_h (h-hbar)->psi；扫 beta_h={0,0.01,0.1,1.0}。
这是新的具名耦合，探索性。
输出：results/zero_evo_1684_betah.json/.png
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


def run(beta_h, steps=2000):
    e = ZeroEvo1684Cycle16(zs.cycle(16),
                           Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05,
                                  alpha=0.1, g_src=1.0, beta=0.0, beta_h=beta_h),
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
    bs = [0.1, 0.3, 1.0, 3.0, 10.0]
    res = {b: run(b, 2000) for b in bs}
    print("beta_h  psi_coh  h_E       p4psi    q16     n_std")
    for b in bs:
        r = res[b][-200:]
        print(f"{b:5.1f}  {np.mean([x['psi_coh'] for x in r]):.4f}  "
              f"{np.mean([x['h_E'] for x in r]):.3e}  "
              f"{p4([x['psi_coh'] for x in res[b]]):.4f}  "
              f"{np.mean([x['q16'] for x in r]):.4f}  {np.mean([x['n_std'] for x in r]):.4f}")
    ck("所有 beta_h 无 NaN", all(np.all(np.isfinite([x['psi_coh'] for x in res[b]])) for b in bs), "")
    ck("所有 beta_h psi_coh 有界（<=1）", all(np.max([x['psi_coh'] for x in res[b]]) <= 1.0 + 1e-9 for b in bs), "")
    ck("beta_h 增大时 psi 周期4 不消失", all(p4([x['psi_coh'] for x in res[b]]) > 1e-6 for b in bs),
       str([round(p4([x['psi_coh'] for x in res[b]]),3) for b in bs]))
    ck("所有 beta_h h_E 仍在", all(np.mean([x['h_E'] for x in res[b][-200:]]) > 0 for b in bs), "")

    print("长跑 beta_h=1.0, 10000 代 ...")
    rl = run(1.0, 10000)
    psi = np.array([x['psi_coh'] for x in rl]); h = np.array([x['h_E'] for x in rl])
    q16 = np.array([x['q16'] for x in rl])
    print("long tail psi_coh", round(float(psi[-200:].mean()),4), "p4psi", round(p4(psi),4),
          "p4h", round(p4(h),3), "q16", round(float(q16[-200:].mean()),4))
    ck("长跑 psi 有周期4", p4(psi) > 1e-6, f"{p4(psi):.4f}")
    ck("长跑 psi_coh 有界", float(np.max(psi)) <= 1.0 + 1e-9, f"max={np.max(psi):.4f}")
    ck("长跑无 NaN", np.all(np.isfinite(psi)), "")

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].semilogx(bs, [p4([x['psi_coh'] for x in res[b]]) for b in bs], marker="o", label="p4 psi")
    ax[0].semilogx(bs, [np.mean([x['psi_coh'] for x in res[b][-200:]]) for b in bs], marker="o", label="psi_coh")
    ax[0].set_title("(a) vs beta_h"); ax[0].legend(fontsize=8)
    ax[1].plot(psi, lw=1, label="psi_coh")
    ax[1].set_title("(b) long run beta_h=1"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_evo_1684_betah.png")
    fig.savefig(fp, dpi=130)
    out = {"scan": {str(b): {"tail": {k: float(np.mean([x[k] for x in res[b][-200:]]))
                                       for k in ["q","q16","h_E","psi_coh","n_std"]},
                              "p4_psi": p4([x['psi_coh'] for x in res[b]])} for b in bs},
           "long": {"psi_coh": float(psi[-200:].mean()), "p4_psi": p4(psi),
                    "p4_h": p4(h), "q16": float(q16[-200:].mean()),
                    "rows": rl},
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "note": "new exploratory h->psi coupling beta_h; beta=0"}
    with open(os.path.join(RES, "zero_evo_1684_betah.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_evo_1684_betah.json"), "和", fp)


if __name__ == "__main__":
    main()
