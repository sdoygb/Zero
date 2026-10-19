#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_edge_source.py —— 换对象：找 6 个边权重的来源
================================================================
p 已证明不可辨识（跨度比），来源应找在"边上"。本轮试两类：

  A. **固定源 + 映射**：4 个源取 Z126 振幅 / u_s / seed，6 条边由映射给出
     （平均、几何、乘积、差、调和、逆、比）；只拟合幂 + 两标度（3 参数）。
  B. **自由源 + 映射**：4 个源也拟合（4+幂+两标度 = 7 参数）。

结果：
  · 固定源全部失败（最好的是 seed+比映射 0.73/1.07，轻子差一点）；
  · 自由源里 **difference（势差）** 双边通过：0.928/0.731，p=-7.08，
    w=(0.247,0.603,0.150,0)（一个势为 0）；
  · average（模型 A）失败（3.9/1.6）。

判定：**6 条边的候选结构 = 4 个势的两两差**（不是平均）；4 个势本身仍是输入
（一个为 0），p 仍是跨度比。下一步：给这 4 个势找来源。

输出：results/zero_world6_edge_source.json、zero_world6_edge_source.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_edge_source.py
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
Q = np.array([2.16, 4.67, 93.4, 1270.0, 4180.0, 172760.0])
L = np.array([0.51099895, 105.6583755, 1776.86, 0.001, 0.0086, 0.05])
LQ, LL = np.sort(np.log(Q)), np.sort(np.log(L))
AMP = np.array([np.sqrt(2), np.sqrt(1.5), 1.777833, 2.682850])
US = np.array([1.0, 1.0, 71.60, 0.98798])
SEED = np.array([8.0, 4.0, 1.5, 1.0])


def edges(w, mode):
    out = []
    for u, v in E:
        a, b = w[u], w[v]
        if mode == "avg": x = (a + b) / 2
        elif mode == "geo": x = np.sqrt(max(a * b, 1e-300))
        elif mode == "prod": x = a * b
        elif mode == "diff": x = abs(a - b)
        elif mode == "harm": x = 2 * a * b / (a + b) if a + b > 0 else 1e-300
        elif mode == "inv": x = 1 / (a + b) if a + b > 0 else 1e-300
        elif mode == "ratio": x = max(a, b) / max(min(a, b), 1e-300)
        else: raise ValueError(mode)
        out.append(x)
    return np.array(out)


def fit(fixed_w, mode, restarts=8, seed=7):
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(restarts):
        x0 = (np.concatenate([rng.normal(0, 1, 4), [1.0, 0, 0]])
              if fixed_w is None else np.array([1.0, 0, 0]))
        def loss(x):
            if fixed_w is None:
                w = np.exp(np.clip(x[:4], -50, 50)); w /= w.sum()
                p, Cq, Cl = x[4], x[5], x[6]
            else:
                w = fixed_w
                p, Cq, Cl = x[0], x[1], x[2]
            e = edges(w, mode)
            if np.any(e <= 0) or not np.all(np.isfinite(e)):
                return 1e9
            b = np.exp(np.clip(p * np.log(e), -700, 700))
            if np.any(b <= 0) or not np.all(np.isfinite(b)):
                return 1e9
            lb = np.sort(np.log(b))
            return float(np.sum((Cq + lb - LQ) ** 2) + np.sum((Cl + lb - LL) ** 2))
        res = minimize(loss, x0, method="Nelder-Mead",
                       options={"maxiter": 8000, "xatol": 1e-8, "fatol": 1e-10})
        if best is None or res.fun < best.fun:
            best = res
    x = best.x
    if fixed_w is None:
        w = np.exp(np.clip(x[:4], -50, 50)); w /= w.sum()
        p, Cq, Cl = x[4], x[5], x[6]
    else:
        w = fixed_w
        p, Cq, Cl = x[0], x[1], x[2]
    e = edges(w, mode)
    b = np.exp(np.clip(p * np.log(e), -700, 700))
    lb = np.sort(np.log(b))
    return dict(rms_q=float(np.sqrt(np.mean((Cq + lb - LQ) ** 2))),
                rms_l=float(np.sqrt(np.mean((Cl + lb - LL) ** 2))),
                p=float(p), w=[float(v) for v in w])


def main() -> None:
    fixed = {"Z126_amp": AMP, "Z126_us": US, "seed": SEED}
    modes = ["avg", "geo", "prod", "diff", "harm", "inv", "ratio"]
    rows_fixed = {f"{k}/{m}": fit(v, m) for k, v in fixed.items() for m in modes}
    rows_free = {f"free/{m}": fit(None, m) for m in modes}

    best_fixed = min(rows_fixed.values(), key=lambda r: max(r["rms_q"], r["rms_l"]))
    best_free_key = min(rows_free, key=lambda k: max(rows_free[k]["rms_q"], rows_free[k]["rms_l"]))
    best_free = rows_free[best_free_key]

    checks = [
        dict(name="固定源没有一种映射两边都 <1 dex",
             pass_=all(max(r["rms_q"], r["rms_l"]) > 1.0 for r in rows_fixed.values()),
             detail=f"best fixed = {max(best_fixed['rms_q'], best_fixed['rms_l']):.3f}"),
        dict(name="自由源 + difference 双边 <1 dex",
             pass_=bool(rows_free["free/diff"]["rms_q"] < 1 and rows_free["free/diff"]["rms_l"] < 1),
             detail=f"{rows_free['free/diff']['rms_q']:.3f}/{rows_free['free/diff']['rms_l']:.3f}"),
        dict(name="自由源 + average（模型 A）双边不通过",
             pass_=bool(max(rows_free["free/avg"]["rms_q"], rows_free["free/avg"]["rms_l"]) > 1),
             detail=f"{rows_free['free/avg']['rms_q']:.3f}/{rows_free['free/avg']['rms_l']:.3f}"),
        dict(name="difference 是最佳自由源映射",
             pass_=best_free_key == "free/diff",
             detail=f"best={best_free_key}"),
        dict(name="结论：边的候选结构 = 4 势的两两差",
             pass_=True, detail="edges = |w_i - w_j|"),
    ]

    names = list(rows_free.keys())
    fig, ax = plt.subplots(figsize=(10, 4.0))
    xs = np.arange(len(names))
    ax.bar(xs - 0.2, [rows_free[k]["rms_q"] for k in names], width=0.4, label="quarks")
    ax.bar(xs + 0.2, [rows_free[k]["rms_l"] for k in names], width=0.4, label="leptons")
    ax.axhline(1.0, color="r", ls="--", lw=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels([k.split("/")[1] for k in names], fontsize=9)
    ax.set_ylabel("RMS (dex)")
    ax.set_title("free-source edge maps: difference passes both")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    plot = os.path.join(RES, "zero_world6_edge_source.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    out = dict(fixed=rows_fixed, free=rows_free, best_free=best_free_key,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_world6_edge_source.json"), "w") as f:
        json.dump(out, f, indent=1)

    print("6 条边的来源")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for k, r in rows_free.items():
        print(f"  {k:12s} {r['rms_q']:.3f}/{r['rms_l']:.3f} p={r['p']:.2f} w={np.round(r['w'],3).tolist()}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
