#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_fit.py —— 4（顶）→ 18（中）→ 6（底=我们的世界）：先拟合
=====================================================================
按用户方向：4 最顶层，18 中间层（不是我们的世界），最底 6 层是我们的世界。
先拟合，拟合通过再证明。

两个模型：
  A（约束）：4 源权重 w → 6 条边 (w_i+w_j)/2 → 质量 m = C e^p（5 参数）；
  B（中18介导）：3 代 × 6 边，每代一个幂 p_g；底 6 = Σ_g e^{p_g}；m = C·b
    （8 参数：w 4 + p 3 + C 1）。

拟合目标（6 个观测，二选一）：
  · 6 个夸克质量（18 = 3 色 × 6 味 的自然读法）；
  · 6 个轻子质量（e,μ,τ + 3 个中微子）。

结论（先说）：
  · A 拟合失败（夸克 RMS 3.9 dex、轻子 1.6 dex）——4→6 的约束产不出观测层级；
  · B 拟合能过（夸克 0.50 dex、轻子 0.69 dex），但 **8 参数对 6 数据（过拟合）**，
    且参数极端（|p|~6–18、log w 跨 3 个量级）——这不是"证明"，是先"拟合通过"。

输出：results/zero_world6_fit.json、zero_world6_fit.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_fit.py
"""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
E = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
QUARKS = np.array([2.16, 4.67, 93.4, 1270.0, 4180.0, 172760.0])       # MeV
LEPTONS = np.array([0.51099895, 105.6583755, 1776.86,
                    0.001, 0.0086, 0.05])                              # MeV (nu ~ eV)


def edges(w):
    return np.array([(w[u] + w[v]) / 2.0 for (u, v) in E])


def fit_A(target, restarts=30, seed=0):
    lm = np.log(np.asarray(target, float))
    order = np.argsort(lm)
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(restarts):
        x0 = np.concatenate([rng.normal(0, 1, 4), [0.0, 1.0]])
        def loss(x):
            w = np.exp(np.clip(x[:4], -50, 50)); w /= w.sum()
            e = edges(w)
            if np.any(e <= 0) or not np.all(np.isfinite(e)):
                return 1e9
            pred = x[4] + x[5] * np.sort(np.log(e))
            return float(np.sum((pred - lm[order]) ** 2))
        res = minimize(loss, x0, method="Nelder-Mead",
                       options={"maxiter": 20000, "xatol": 1e-10, "fatol": 1e-12})
        if best is None or res.fun < best.fun:
            best = res
    w = np.exp(np.clip(best.x[:4], -50, 50)); w /= w.sum()
    pred = np.exp(best.x[4] + best.x[5] * np.sort(np.log(edges(w))))
    return dict(rms=float(np.sqrt(best.fun / 6)), w=w.tolist(),
                p=float(best.x[5]), C=float(np.exp(best.x[4])),
                pred=pred.tolist())


def fit_B(target, restarts=40, seed=1):
    lm = np.log(np.asarray(target, float))
    order = np.argsort(lm)
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(restarts):
        x0 = np.concatenate([rng.normal(0, 1, 4), [0.0, 0.0, 0.0], [0.0]])
        def loss(x):
            w = np.exp(np.clip(x[:4], -50, 50)); w /= w.sum()
            e = edges(w)
            if np.any(e <= 0) or not np.all(np.isfinite(e)):
                return 1e9
            p = x[4:7]
            b = np.sum(np.exp(np.clip(p[None, :] * np.log(e)[:, None], -700, 700)), axis=1)
            if np.any(b <= 0) or not np.all(np.isfinite(b)):
                return 1e9
            pred = x[7] + np.log(b)
            return float(np.sum((np.sort(pred) - lm[order]) ** 2))
        res = minimize(loss, x0, method="Nelder-Mead",
                       options={"maxiter": 40000, "xatol": 1e-10, "fatol": 1e-12})
        if best is None or res.fun < best.fun:
            best = res
    w = np.exp(np.clip(best.x[:4], -50, 50)); w /= w.sum()
    e = edges(w)
    b = np.sum(np.exp(np.clip(np.array(best.x[4:7])[None, :] * np.log(e)[:, None],
                              -700, 700)), axis=1)
    pred = np.exp(best.x[7] + np.log(b))
    return dict(rms=float(np.sqrt(best.fun / 6)), w=w.tolist(),
                p=[float(x) for x in best.x[4:7]], C=float(best.x[7]),
                pred=pred.tolist())


def main() -> None:
    rows = {}
    for name, target in (("quarks", QUARKS), ("leptons", LEPTONS)):
        rows[name] = dict(A=fit_A(target), B=fit_B(target), target=target.tolist())

    extreme = any(abs(p) > 5 for k in ("quarks", "leptons") for p in rows[k]["B"]["p"])
    checks = [
        dict(name="A（4→6 约束）拟合失败（RMS > 1 dex）",
             pass_=all(rows[k]["A"]["rms"] > 1.0 for k in rows),
             detail={k: round(rows[k]["A"]["rms"], 3) for k in rows}),
        dict(name="B（18 介导）拟合通过（RMS < 0.7 dex）",
             pass_=all(rows[k]["B"]["rms"] < 0.7 for k in rows),
             detail={k: round(rows[k]["B"]["rms"], 3) for k in rows}),
        dict(name="B 是过拟合（8 参数 > 6 数据）",
             pass_=True, detail="8 params vs 6 data"),
        dict(name="B 的参数极端（|p| > 5 或 w 跨量级）",
             pass_=bool(extreme),
             detail={k: [round(p, 2) for p in rows[k]["B"]["p"]] for k in rows}),
        dict(name="约束模型 A 的 w 退化（夸克：全 1/4，边全等）",
             pass_=bool(abs(max(rows["quarks"]["A"]["w"]) - min(rows["quarks"]["A"]["w"])) < 1e-3),
             detail={k: [round(x, 3) for x in rows[k]["A"]["w"]] for k in rows}),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, name in zip(axes, ("quarks", "leptons")):
        t = np.sort(np.array(rows[name]["target"]))
        ax.loglog(t, t, "k--", lw=0.8, label="obs=pred")
        for lab, key, mk in (("A (4→6)", "A", "s"), ("B (18-mediated)", "B", "o")):
            r = rows[name][key]
            ax.loglog(t, np.sort(r["pred"]), mk, ms=5,
                      label=f"{lab}, RMS={r['rms']:.2f} dex")
        ax.set_xlabel("observed mass (MeV)")
        ax.set_ylabel("predicted (MeV)")
        ax.set_title(f"{name}: A={rows[name]['A']['rms']:.2f}, "
                     f"B={rows[name]['B']['rms']:.2f} dex")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    plot = os.path.join(RES, "zero_world6_fit.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    out = dict(rows=rows,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_world6_fit.json"), "w") as f:
        json.dump(out, f, indent=1)

    print("4→18→6 先拟合")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for k in rows:
        print(f"  {k}: A rms={rows[k]['A']['rms']:.3f} dex, "
              f"B rms={rows[k]['B']['rms']:.3f} dex, B p={[round(x,2) for x in rows[k]['B']['p']]}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
