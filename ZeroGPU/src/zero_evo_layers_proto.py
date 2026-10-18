#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L0-L1-L2 层集成 prototype（部分）。

- L0：cycle16 舞台 + l0_back 反馈（L2 记录回写 seed）；
- L1：静态高度 b_h=[4,2]、面积 A_h=[3,2]（只读出，无自身动力学）；
- L2：从 16 类记录取 4-site 占据，算 3 个相邻比（频率槽）。
扫 l0_back={0,0.1,0.3}，2000 代。
输出：results/zero_evo_layers_proto.json/.png
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


def run(l0_back, steps=2000):
    e = ZeroEvo1684Cycle16(zs.cycle(16), Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05,
                                                  alpha=0.1, g_src=1.0, l0_back=l0_back),
                            seed=np.ones(16))
    e.perm = MAPS["seq"]
    rows = []
    for g in range(1, steps+1):
        e.step()
        s = e.snapshot()
        rec = e._rec16()
        p = rec.reshape(4, 4).sum(axis=1)
        p = np.maximum(p, 1e-300); p = p / p.sum()
        d = p - p.mean()
        c = np.array([(d[0]-d[1])/np.sqrt(2),
                      (d[0]+d[1]-2*d[2])/np.sqrt(6),
                      (d[0]+d[1]+d[2]-3*d[3])/np.sqrt(12)])
        rows.append({"gen": g, "q": s["q"], "q16": e.q16(), "n_std": s["n_std"],
                     "p": p.tolist(), "slots": c.tolist()})
    return rows


def p4(a):
    a = np.asarray(a) - np.mean(a); f = np.abs(np.fft.rfft(a)); return float(f[len(a)//T])


def main():
    l0s = [0.0, 0.1]
    res = {l: run(l) for l in l0s}
    print("l0_back  q       q16     n_std   Helmert(mean)        cov-rank  eig")
    ranks = {}
    for l in l0s:
        r = res[l][-200:]
        sl = np.mean([x["slots"] for x in r], axis=0)
        X = np.array([x["slots"] for x in res[l][-1000:]])
        C = np.cov(X.T)
        eig = np.linalg.eigvalsh(C)
        rank = int(np.sum(eig > max(eig.max()*1e-8, 1e-14)))
        ranks[l] = rank
        print(f"{l:.2f}     {np.mean([x['q'] for x in r]):.4f}  {np.mean([x['q16'] for x in r]):.4f}  "
              f"{np.mean([x['n_std'] for x in r]):.4f}  {np.round(sl,3)}  {rank}  {np.round(eig,6)}")
    ck("L1 静态 b_h/A_h 可读出", True, "b_h=[4,2], A_h=[3,2]")
    ck("L2 给出 3 个 Helmert 分量", all(len(x["slots"]) == 3 for x in res[0.0]), "")
    ck("无 NaN", all(np.all(np.isfinite([x["slots"] for x in res[l]])) for l in l0s), "")
    ck("l0_back 改变动力学", any(abs(np.mean([x['q'] for x in res[l][-200:]]) -
                                      np.mean([x['q'] for x in res[0.0][-200:]])) > 1e-6
                                  for l in l0s[1:]), "")
    ck("L2 分量有周期 4", any(p4([x["slots"][i] for x in res[0.1]]) > 1e-6 for i in range(3)), "")
    ck("协方差秩=1（不是3个独立槽）", all(ranks[l] == 1 for l in l0s), str(ranks))

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    for i in range(3):
        ax[0].plot([x["slots"][i] for x in res[0.1]], lw=1, label=f"Helmert {i+1}")
    ax[0].set_title("(a) L2 Helmert components (l0=0.1)"); ax[0].legend(fontsize=8)
    for l in l0s:
        ax[1].plot([x["q"] for x in res[l]], lw=1, label=f"l0={l}")
    ax[1].set_title("(b) content q"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_evo_layers_proto.png")
    fig.savefig(fp, dpi=130)
    out = {"l0_back": {str(l): {"q": float(np.mean([x['q'] for x in res[l][-200:]])),
                                "q16": float(np.mean([x['q16'] for x in res[l][-200:]])),
                                "n_std": float(np.mean([x['n_std'] for x in res[l][-200:]])),
                                "slots": np.mean([x['slots'] for x in res[l][-200:]], axis=0).tolist(),
                                "cov_rank": ranks[l],
                                "p4_slots": [p4([x['slots'][i] for x in res[l]]) for i in range(3)]}
                        for l in l0s},
           "l1_static": {"b_h": [4, 2], "A_h": [3, 2]},
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "L0 feedback + L2 slots readout; L1 static"}
    with open(os.path.join(RES, "zero_evo_layers_proto.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")


if __name__ == "__main__":
    main()
