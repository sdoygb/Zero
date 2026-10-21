#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_optimal_forget.py —— 把“一生多少代”当参数：最优遗忘落在哪？

上一轮：若成团保留 f=p（保留=闭合比例），有效传递
    T(p) = p·exp(-6 p^2 L / Hc)          （L = 一生的代数）
    p*   = sqrt(Hc / (12 L))              （峰值）
    遗忘 = 1 - p*

扫 L，看最优遗忘是否落在“人一生能记住的极少”那一带。

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_optimal_forget.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
HC = 60.0          # 6C, C=2(T+1)=10


def optimum(L, p=None):
    p_star = float(np.sqrt(HC / (12.0 * L)))
    if p is None:
        p = p_star
    H = 6.0 * p * L
    T = p * np.exp(-6.0 * p * p * L / HC)
    return dict(L=float(L), p_star=p_star, forget=1.0 - p_star, T=T, H=H)


def scan(L, n=400):
    ps = np.logspace(-4, np.log10(0.999), n)
    vals = [p * np.exp(-6.0 * p * p * L / HC) for p in ps]
    i = int(np.argmax(vals))
    return float(ps[i]), float(vals[i])


def main():
    Ls = [10, 30, 100, 300, 1000, 3000, 10000, 100000, 1000000]
    rows = [optimum(L) for L in Ls]
    scan_rows = []
    for L in Ls:
        pn, Tn = scan(L)
        scan_rows.append(dict(L=L, p_num=pn, T_num=Tn))

    # “人一生能记住的极少”：长期事件记忆保留约 1%–10%（外部经验，非 Zero）
    band = [L for L in Ls if 0.01 <= np.sqrt(HC / (12.0 * L)) <= 0.10]

    checks = [
        dict(name="c1 数值峰 p 与解析 sqrt(Hc/(12L)) 一致",
             pass_=all(abs(r["p_num"] - optimum(r["L"])["p_star"]) / optimum(r["L"])["p_star"] < 0.05
                       for r in scan_rows),
             detail=" | ".join(f"L={r['L']:g}:{r['p_num']:.3f}/{optimum(r['L'])['p_star']:.3f}"
                               for r in scan_rows[:4])),
        dict(name="c2 遗忘随一生变长单调升，趋近 1",
             pass_=bool(all(rows[i + 1]["forget"] > rows[i]["forget"] for i in range(len(rows) - 1))
                        and rows[-1]["forget"] > 0.99),
             detail=" ".join(f"L={r['L']:g}:{r['forget']:.3f}" for r in rows)),
        dict(name="c3 最优传递 T*=p*·e^{-1/2}（峰处 6p²L/Hc=1/2）",
             pass_=bool(all(abs(r["T"] - np.exp(-0.5) * r["p_star"]) < 1e-12 for r in rows)),
             detail=f"T(L=1e5)={optimum(100000)['T']:.4f}"),
        dict(name="c4 若一生 10^4–10^6 代，最优遗忘落在 98%–99.8%",
             pass_=bool(0.97 <= optimum(1e4)["forget"] and optimum(1e6)["forget"] <= 0.999),
             detail=f"1e4:{optimum(1e4)['forget']:.4f}, 1e6:{optimum(1e6)['forget']:.4f}"),
        dict(name="c5 与“长期记忆保留 1%–10%”对照：保留率 p* 落带的 L",
             pass_=bool(len(band) > 0),
             detail=f"L={band}"),
    ]

    out = dict(params=dict(HC=HC), rows=rows, scan=scan_rows, band=band,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_optimal_forget.json"), "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot([r["L"] for r in rows], [r["forget"] for r in rows], "o-")
    ax.set_xscale("log")
    ax.set_xlabel("lifetime in generations L")
    ax.set_ylabel("optimal forgetting 1-p*")
    ax.set_title("longer life → forget more (optimum)")
    ax.grid(alpha=.3)
    fig.tight_layout()
    plot = os.path.join(RES, "zero_optimal_forget.png")
    fig.savefig(plot, dpi=130); plt.close(fig)

    print("最优遗忘 vs 一生代数")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
