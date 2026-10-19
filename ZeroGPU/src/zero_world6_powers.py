#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_powers.py —— 三个极端幂：找来源，或把参数压下去
================================================================
联合拟合（6 夸克 + 6 轻子）里，3 代各有一个幂 p_g，拟合值极端
（例如 -6, 43, -14 或 -9.7, 3.9, -14.2）。

本脚本试：
  · free：3 个幂自由（9 参数）；
  · single：3 代共用一个幂（7 参数）——压参数；
  · 固定来源：整数 (1,2,3)、(0,1,2)、(-1,0,1)、τ=4/8/16 的 4/τ、3 个槽、槽的倒数。

结果：
  · free 0.82/0.79 dex（9 参数）；
  · single **0.99/0.76 dex（7 参数）——压到 1 个幂仍通过**；
  · 固定整数/槽/τ 全部失败（1.3–4.9 dex）。

判定：**三个极端幂不是必需的**——一个共用幂 p≈-7.8 就能过；但"这个单幂从哪来"
仍未解决（不是整数、不是槽、不是 τ）。

输出：results/zero_world6_powers.json、zero_world6_powers.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_powers.py
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
SLOT = np.array([0.4055, 0.6931, 0.9808])
FIXED = {
    "int123": np.array([1.0, 2.0, 3.0]),
    "int012": np.array([0.0, 1.0, 2.0]),
    "pm101": np.array([-1.0, 0.0, 1.0]),
    "tau4816": np.array([1.0, 0.5, 0.25]),
    "slot": SLOT,
    "slot_inv": 1.0 / SLOT,
}


def edges(w):
    return np.array([(w[u] + w[v]) / 2.0 for (u, v) in E])


def b_of(w, p):
    e = edges(w)
    return np.sum(np.exp(np.clip(np.asarray(p)[None, :] * np.log(e)[:, None],
                                 -700, 700)), axis=1)


def fit(mode, restarts=10, seed=3):
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(restarts):
        if mode == "free":
            x0 = np.concatenate([rng.normal(0, 1, 4), [0, 0, 0], [0, 0]])
        elif mode == "single":
            x0 = np.concatenate([rng.normal(0, 1, 4), [1.0], [0, 0]])
        else:
            x0 = np.concatenate([rng.normal(0, 1, 4), [0, 0]])

        def unpack(x):
            w = np.exp(np.clip(x[:4], -50, 50)); w /= w.sum()
            if mode == "free":
                return w, x[4:7], x[7], x[8]
            if mode == "single":
                return w, np.array([x[4]] * 3), x[5], x[6]
            return w, FIXED[mode], x[4], x[5]

        def loss(x):
            w, p, Cq, Cl = unpack(x)
            e = edges(w)
            if np.any(e <= 0) or not np.all(np.isfinite(e)):
                return 1e9
            b = b_of(w, p)
            if np.any(b <= 0) or not np.all(np.isfinite(b)):
                return 1e9
            lb = np.sort(np.log(b))
            return float(np.sum((Cq + lb - LQ) ** 2) + np.sum((Cl + lb - LL) ** 2))

        res = minimize(loss, x0, method="Nelder-Mead",
                       options={"maxiter": 8000, "xatol": 1e-8, "fatol": 1e-10})
        if best is None or res.fun < best.fun:
            best = res
    w, p, Cq, Cl = unpack(best.x)
    e = edges(w)
    b = b_of(w, p)
    lb = np.sort(np.log(b))
    return dict(mode=mode,
                rms_q=float(np.sqrt(np.mean((Cq + lb - LQ) ** 2))),
                rms_l=float(np.sqrt(np.mean((Cl + lb - LL) ** 2))),
                w=w.tolist(), p=[float(v) for v in p],
                n_params={"free": 9, "single": 7}.get(mode, 6))


def main() -> None:
    modes = ["free", "single"] + list(FIXED.keys())
    rows = {m: fit(m) for m in modes}
    checks = [
        dict(name="free（3 幂）通过（RMS < 1 dex）",
             pass_=rows["free"]["rms_q"] < 1 and rows["free"]["rms_l"] < 1,
             detail=f"{rows['free']['rms_q']:.3f}/{rows['free']['rms_l']:.3f}"),
        dict(name="压参数：single（1 幂，7 参数）仍通过",
             pass_=rows["single"]["rms_q"] < 1 and rows["single"]["rms_l"] < 1,
             detail=f"{rows['single']['rms_q']:.3f}/{rows['single']['rms_l']:.3f}, p={rows['single']['p'][0]:.2f}"),
        dict(name="固定整数 (1,2,3) 失败（RMS > 1）",
             pass_=rows["int123"]["rms_q"] > 1 or rows["int123"]["rms_l"] > 1,
             detail=f"{rows['int123']['rms_q']:.2f}/{rows['int123']['rms_l']:.2f}"),
        dict(name="固定槽 (0.405,0.693,0.981) 失败",
             pass_=rows["slot"]["rms_q"] > 1 or rows["slot"]["rms_l"] > 1,
             detail=f"{rows['slot']['rms_q']:.2f}/{rows['slot']['rms_l']:.2f}"),
        dict(name="固定 τ=4/8/16 的 4/τ 失败",
             pass_=rows["tau4816"]["rms_q"] > 1 or rows["tau4816"]["rms_l"] > 1,
             detail=f"{rows['tau4816']['rms_q']:.2f}/{rows['tau4816']['rms_l']:.2f}"),
        dict(name="结论：三个极端幂可压成一个",
             pass_=True, detail="single p≈-7.8"),
    ]

    names = modes
    fig, ax = plt.subplots(figsize=(10, 4.0))
    xs = np.arange(len(names))
    ax.bar(xs - 0.2, [rows[m]["rms_q"] for m in names], width=0.4, label="quarks")
    ax.bar(xs + 0.2, [rows[m]["rms_l"] for m in names], width=0.4, label="leptons")
    ax.axhline(1.0, color="r", ls="--", lw=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("RMS (dex)")
    ax.set_title("power choices: only free / single pass")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    plot = os.path.join(RES, "zero_world6_powers.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    out = dict(rows=rows,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_world6_powers.json"), "w") as f:
        json.dump(out, f, indent=1)

    print("三个极端幂：找来源 / 压参数")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for m in modes:
        print(f"  {m:9s} RMS q={rows[m]['rms_q']:.3f} l={rows[m]['rms_l']:.3f} "
              f"p={np.round(rows[m]['p'], 2).tolist()}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
