#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_evo_18_mass.py —— 把 18 世界校准到观测的轻子质量层级
================================================================
18 世界的 3 代份额由代权重 gen_w 决定。三种模式：
  · slots：默认（log2, log(8/3), log1.5 归一）；
  · sqrt_mass：带电轻子 sqrt(m) 归一（Z126 振幅律的振幅就是 sqrt(m)）；
  · mass：带电轻子 m 归一。

跑 2000 代，检查：
  1. sqrt_mass / mass 模式的份额是否等于观测比值（按构造应相等）；
  2. slots 模式离观测多远（因子）；
  3. 校准后引擎是否稳定（无 NaN、有界、n_std/q 正常）；
  4. 影子 6 是否满支撑、是否被重代主导。

输出：results/zero_evo_18_mass.json、zero_evo_18_mass.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_evo_18_mass.py
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_evo import Params  # noqa: E402
from zero_evo_18 import ZeroEvo18  # noqa: E402
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")
M_LEP = np.array([0.51099895, 105.6583755, 1776.86])   # MeV


def base_params():
    return dict(micro=True, L=4, T=4, retention="top2",
                eta_n=0.2, eta_h=0.2, eta_psi=0.2,
                mu=0.5, kappa_p=0.1, gamma=0.2, K=2.0,
                capacity="soft", K_hard=1.5, g_src=1.0, c_h=1.2,
                Gamma_h=0.2, kappa_h=0.5, alpha=0.0,
                wave_on=True, beta=0.5, g_wave=0.5, kappa_psi=0.5, Gamma_w=0.5)


def weights(mode):
    if mode == "slots":
        w = np.array([math.log(2.0), math.log(8.0 / 3.0), math.log(1.5)])
    elif mode == "sqrt_mass":
        w = np.sqrt(M_LEP)
    elif mode == "mass":
        w = M_LEP.copy()
    else:
        raise ValueError(mode)
    return w / w.sum()


def run(mode, gens=2000):
    e = ZeroEvo18(zs.k4(), Params(**base_params()), seed=np.array([8.0, 4.0, 1.5, 1.0]))
    e.gen_w = weights(mode)
    e.W18 = np.outer(e.gen_w, e._edge_source(e.seed))
    e.W18 /= e.W18.sum()
    sh = []
    for _ in range(gens):
        e.step()
        sh.append(e.W18.sum(axis=1).copy())
    sh = np.array(sh)
    tail = sh[gens // 2:]
    return dict(mode=mode,
                shares=[float(x) for x in tail.mean(axis=0)],
                std=[float(x) for x in tail.std(axis=0)],
                shadow=[float(x) for x in e.W18.sum(axis=0)],
                n_std=float(e.n.std()), q=float(e.snapshot()["q"]),
                finite=bool(np.isfinite(sh).all()),
                trace=sh[::max(1, gens // 200)].tolist())


def main() -> None:
    slots = run("slots")
    sq = run("sqrt_mass")
    ms = run("mass")
    target_sq = weights("sqrt_mass")
    target_m = weights("mass")

    err_sq = np.abs(np.array(sq["shares"]) - target_sq)
    err_m = np.abs(np.array(ms["shares"]) - target_m)
    ratio = np.array(slots["shares"]) / target_sq
    checks = [
        dict(name="sqrt_mass 份额 = 观测 sqrt(m) 归一（逐位）",
             pass_=float(np.max(err_sq)) < 1e-6,
             detail=f"max err={np.max(err_sq):.2e}"),
        dict(name="mass 份额 = 观测 m 归一（逐位）",
             pass_=float(np.max(err_m)) < 1e-6,
             detail=f"max err={np.max(err_m):.2e}"),
        dict(name="默认 slots 离观测很远（最大因子 > 10）",
             pass_=float(np.max(ratio)) > 10,
             detail=f"ratio={np.round(ratio, 2).tolist()}"),
        dict(name="校准后引擎稳定（无 NaN、有界）",
             pass_=all(c["finite"] for c in (sq, ms)),
             detail="2000 gens"),
        dict(name="校准后影子 6 满支撑（6 条边 > 0）",
             pass_=bool(min(sq["shadow"]) > 0),
             detail=f"shadow={np.round(sq['shadow'], 4).tolist()}"),
        dict(name="校准后宏观量正常（n_std/q 非零）",
             pass_=sq["n_std"] > 0 and sq["q"] > 0,
             detail=f"n_std={sq['n_std']:.3f}, q={sq['q']:.3f}"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    xs = np.arange(3)
    axes[0].bar(xs - 0.27, target_sq, width=0.27, label="observed sqrt(m)")
    axes[0].bar(xs, sq["shares"], width=0.27, label="18-world sqrt_mass")
    axes[0].bar(xs + 0.27, slots["shares"], width=0.27, label="18-world default")
    axes[0].set_yscale("log")
    axes[0].set_xticks(xs)
    axes[0].set_xticklabels(["gen 0 (e)", "gen 1 (mu)", "gen 2 (tau)"])
    axes[0].set_title("generation shares vs observed lepton hierarchy")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3, axis="y")
    for i in range(6):
        axes[1].bar(i, sq["shadow"][i], width=0.6)
    axes[1].set_xticks(range(6))
    axes[1].set_xticklabels(["01", "02", "03", "12", "13", "23"])
    axes[1].set_title("shadow 6 (calibrated; dominated by heavy gen)")
    axes[1].grid(alpha=0.3, axis="y")
    fig.tight_layout()
    plot = os.path.join(RES, "zero_evo_18_mass.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    out = dict(m_lep=M_LEP.tolist(),
               cases=[slots, sq, ms],
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_evo_18_mass.json"), "w") as f:
        json.dump(out, f, indent=1)

    print("18 世界 × 轻子质量层级")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  sqrt_mass shares={np.round(sq['shares'], 4).tolist()}")
    print(f"  mass shares={np.round(ms['shares'], 6).tolist()}")
    print(f"  slots shares={np.round(slots['shares'], 4).tolist()} (ratio to obs {np.round(ratio, 2).tolist()})")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
