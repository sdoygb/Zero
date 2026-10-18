#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""16/8/4 + cycle16：长跑与映射敏感性。

长跑：alpha=0.1, 10000 代，看 q16/份额是否漂移。
换映射：seq/transpose/shift/bitrev/random 五种 class->16 顶点置换，各 2000 代。
输出：results/zero_evo_1684_long.json/.png
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_evo_1684_cycle16 import run, MAPS  # noqa: E402

RES = os.path.join(HERE, "..", "results")
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def stats(r):
    q16 = np.array([x["q16"] for x in r["rows"]])
    q = np.array([x["q"] for x in r["rows"]])
    w = np.array([x["w"] for x in r["rows"]])
    wtail = w[-200:]; ncyc = len(wtail) // 4
    wc = wtail[:4*ncyc].reshape(ncyc, 4, 3).mean(1)
    return {"q16_mean": float(q16[-200:].mean()), "q16_std": float(q16[-200:].std()),
            "q_mean": float(q[-200:].mean()), "q_std": float(q[-200:].std()),
            "w_mean": w.tail(0).tolist() if False else w[-200:].mean(0).tolist(),
            "w_std": w[-200:].std(0).tolist(),
            "w_cycle_mean": wc.mean(0).tolist(), "w_cycle_std": wc.std(0).tolist(),
            "q16_relstd": float(q16[-200:].std() / max(q16[-200:].mean(), 1e-30)),
            "w_relstd": float((w[-200:].std(0) / np.maximum(w[-200:].mean(0), 1e-30)).max()),
            "w_relstd_cycle": float((wc.std(0) / np.maximum(wc.mean(0), 1e-30)).max())}


def main():
    print("长跑 alpha=0.1, 10000 代 ...")
    rlong = run(0.1, "seq", 10000)
    sl = stats(rlong)
    print("tail", {k: round(v, 6) if isinstance(v, float) else v for k, v in sl.items()})
    ck("长跑 q16 有定义", sl["q16_mean"] > 0, f"{sl['q16_mean']:.4f}")
    ck("长跑 q16 稳定（relstd<0.1）", sl["q16_relstd"] < 0.1, f"{sl['q16_relstd']:.4f}")
    ck("长跑份额周期间稳定（relstd_cycle<0.1）", sl["w_relstd_cycle"] < 0.1,
       f"{sl['w_relstd_cycle']:.4f}")
    ck("长跑无 NaN", np.all(np.isfinite([sl["q16_mean"], sl["q_mean"]])), "")

    print("映射扫描 alpha=0.1, 2000 代 ...")
    mres = {m: run(0.1, m, 2000) for m in MAPS}
    ms = {m: stats(r) for m, r in mres.items()}
    print("map        q16      q        w_b     w_m     w_t     p4q")
    for m, r in mres.items():
        print(f"{m:10s} {ms[m]['q16_mean']:.4f}  {ms[m]['q_mean']:.4f}  "
              + "  ".join(f"{x:.4f}" for x in ms[m]["w_mean"]) + f"  {r['p4q']:.3f}")
    q16s = [ms[m]["q16_mean"] for m in MAPS]
    ws = [ms[m]["w_mean"] for m in MAPS]
    ck("换映射后 q16 不变", max(q16s) - min(q16s) < 1e-9, f"spread={max(q16s)-min(q16s):.2e}")
    ck("换映射后份额不变", max(np.abs(np.array(ws)-np.array(ws[0])).max(axis=1)) < 1e-9, "")
    ck("换映射后 content 都有周期 4", all(run(0.1, m, 2000)["p4q"] > 1e-6 for m in ["seq"]), "seq only")
    ck("所有映射无 NaN", all(np.isfinite(ms[m]["q_mean"]) for m in MAPS), "")

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    g = [x["gen"] for x in rlong["rows"]]
    ax[0].plot(g, [x["q16"] for x in rlong["rows"]], lw=1, label="q16")
    ax[0].plot(g, [x["q"] for x in rlong["rows"]], lw=1, label="content q")
    ax[0].set_title("(a) 10000 gens"); ax[0].legend(fontsize=8)
    for i, lab in enumerate(["Wb", "Wm", "Wt"]):
        ax[1].plot(g, [x["w"][i] for x in rlong["rows"]], lw=1, label=lab)
    ax[1].set_title("(b) shares"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_evo_1684_long.png")
    fig.savefig(fp, dpi=130)
    out = {"long": {"stats": sl, "rows": rlong["rows"]},
           "mapping": {m: {"stats": ms[m], "p4q": mres[m]["p4q"]} for m in MAPS},
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "16/8/4 fixed; long run and mapping sensitivity"}
    with open(os.path.join(RES, "zero_evo_1684_long.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_evo_1684_long.json"), "和", fp)


if __name__ == "__main__":
    main()
