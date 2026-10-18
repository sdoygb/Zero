#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""6/18/4 count prototype：底 6、中 18、顶 4。

6=K4 边；18=3×6（3 个 L2 槽 × 6 边）；4=K4 顶点。
统一 τ=4；顶全收闭合、中只收 top2、底不收；失败闭合年龄到 4 清。
这是 count-level 测试，不是图引擎。
输出：results/zero_tau_6_18_4.json/.png
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
NB, NM, NT = 6, 18, 4
TAU = 4
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def main():
    mb = np.ones(NB) / NB
    mm = np.ones(NM) / NM
    mt = np.ones(NT) / NT
    ab = np.zeros(NB, dtype=int); am = np.zeros(NM, dtype=int); at = np.zeros(NT, dtype=int)
    rows = []
    for g in range(1, 65):
        # seed to top
        mt += 0.05
        # ages: top all closure; middle top2 closure; bottom none
        at[:] = 0
        am += 1; am[:2] = 0   # top2 = first two classes (input)
        ab += 1
        # clear
        mb[ab >= TAU] = 0.0; ab[ab >= TAU] = 0
        mm[am >= TAU] = 0.0; am[am >= TAU] = 0
        mt[at >= TAU] = 0.0; at[at >= TAU] = 0
        # pre-flow observation
        pre_b0 = bool(np.abs(mb).sum() < 1e-30)
        pre_m = int(np.sum(mm > 1e-12))
        pre_t = int(np.sum(mt > 1e-12))
        # flow top->middle->bottom, uniform
        k = 0.2
        J = k * mt.sum()
        if mt.sum() > 0:
            mt -= mt * (J / mt.sum())
            mm += J / NM
        J2 = k * mm.sum()
        if mm.sum() > 0:
            mm -= mm * (J2 / mm.sum())
            mb += J2 / NB
        total = mb.sum() + mm.sum() + mt.sum()
        if total > 0:
            mb /= total; mm /= total; mt /= total
        # shares/q
        allm = np.concatenate([mb, mm, mt])
        shares = [mb.sum(), mm.sum(), mt.sum()]
        q = float(((allm / allm.sum()) ** 2).sum()) if allm.sum() > 0 else 0.0
        rows.append({"gen": g, "pre_b0": pre_b0, "pre_m": pre_m, "pre_t": pre_t,
                     "shares": [float(x) for x in shares], "q": q})
        if g <= 8 or g % 4 == 0:
            print(f"g={g:>2} preB0={pre_b0} preM={pre_m} preT={pre_t} shares={np.round(shares,3)} q={q:.3f}")
    g4 = rows[3]; g8 = rows[7]
    ck("第4代底6全灭（pre-flow）", g4["pre_b0"], "")
    ck("第4代中18只留2（pre-flow）", g4["pre_m"] == 2, f"{g4['pre_m']}")
    ck("第4代顶4（pre-flow）", g4["pre_t"] == 4, f"{g4['pre_t']}")
    ck("第8代模式重复", g8["pre_b0"] and g8["pre_m"] == 2 and g8["pre_t"] == 4, "")
    ck("无 NaN", all(np.isfinite(r["shares"]).all() for r in rows), "")

    gs = [r["gen"] for r in rows]; sh = np.array([r["shares"] for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    for i, lab in enumerate(["bottom 6", "middle 18", "top 4"]):
        ax[0].plot(gs, sh[:, i], lw=1, label=lab)
    ax[0].set_title("(a) shares"); ax[0].legend(fontsize=8)
    ax[1].plot(gs, [r["pre_m"] for r in rows], lw=1, label="pre middle active")
    ax[1].plot(gs, [r["pre_t"] for r in rows], lw=1, label="pre top active")
    ax[1].set_title("(b) pre-flow active"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_tau_6_18_4.png")
    fig.savefig(fp, dpi=130)
    out = {"counts": [NB, NM, NT], "tau": TAU, "rows": rows, "checks": CK,
           "n_fail": sum(1 for c in CK if not c["pass"]),
           "note": "count-level; middle top2 = first two classes (input)"}
    with open(os.path.join(RES, "zero_tau_6_18_4.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")


if __name__ == "__main__":
    main()
