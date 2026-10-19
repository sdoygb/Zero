#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_evo_cores.py —— 18 世界 vs 16/8/4：同一宏观下的并排比较
================================================================
(a) 代依赖宏观耦合：ZeroEvo18 的 gamma_gen（用**边梯度**区分 3 个匹配；
    顶点平均对 3 个完美匹配相同，无效）。
(b) 并排：18（gamma=0/2）与 16/8/4 在同一套 12 通道宏观参数下跑 5000 代，
    比较 n_std / q / h_E / psi_coh / 层份额。

输出：results/zero_evo_cores.json、zero_evo_cores.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_evo_cores.py
"""
from __future__ import annotations

import json
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
from zero_evo_1684 import ZeroEvo1684  # noqa: E402
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")
SEED = np.array([8.0, 4.0, 1.5, 1.0])


def base_params():
    return dict(micro=True, L=4, T=4, retention="top2",
                eta_n=0.2, eta_h=0.2, eta_psi=0.2,
                mu=0.5, kappa_p=0.1, gamma=0.2, K=2.0,
                capacity="soft", K_hard=1.5, g_src=1.0, c_h=1.2,
                Gamma_h=0.2, kappa_h=0.5, alpha=0.0,
                wave_on=True, beta=0.5, g_wave=0.5, kappa_psi=0.5, Gamma_w=0.5)


def run(cls, gamma, gens):
    e = cls(zs.k4(), Params(**base_params()), seed=SEED.copy())
    if hasattr(e, "gamma_gen"):
        e.gamma_gen = gamma
    sh = []
    ns, qs, he, pc = [], [], [], []
    for _ in range(gens):
        e.step()
        if isinstance(e, ZeroEvo18):
            sh.append(e.W18.sum(axis=1).copy())
        else:
            s = np.array([e.Wb.sum(), e.Wm.sum(), e.Wt.sum()])
            sh.append(s / s.sum())
        snap = e.snapshot()
        ns.append(snap["n_std"]); qs.append(snap["q"])
        he.append(snap["h_energy"]); pc.append(snap["psi_coh"])
    sh = np.array(sh)
    tail = slice(gens // 2, gens)
    return dict(name=f"{cls.__name__}_g{gamma:g}",
                shares=[float(x) for x in sh[tail].mean(axis=0)],
                shares_std=[float(x) for x in sh[tail].std(axis=0)],
                n_std=float(np.mean(ns[tail])), q=float(np.mean(qs[tail])),
                h_energy=float(np.mean(he[tail])), psi_coh=float(np.mean(pc[tail])),
                finite=bool(np.isfinite(sh).all()),
                trace=sh[::max(1, gens // 200)].tolist())


def main() -> None:
    cases = [run(ZeroEvo18, 0.0, 5000),
             run(ZeroEvo18, 2.0, 5000),
             run(ZeroEvo1684, 0.0, 5000)]

    c18, c18g, c16 = cases
    same_ns = abs(c18["n_std"] - c16["n_std"]) < 1e-3
    same_q = abs(c18["q"] - c16["q"]) < 1e-3
    checks = [
        dict(name="三个 core 都跑完 5000 代（无 NaN）",
             pass_=all(c["finite"] for c in cases),
             detail="all finite"),
        dict(name="18(gamma=0) 份额固定（std=0）",
             pass_=max(c18["shares_std"]) < 1e-9,
             detail=f"shares={np.round(c18['shares'], 3).tolist()}"),
        dict(name="18(gamma=2) 份额被宏观改变（std>0.01 且均值差>0.05）",
             pass_=(max(c18g["shares_std"]) > 0.01
                    and max(np.abs(np.array(c18g["shares"]) - np.array(c18["shares"]))) > 0.05),
             detail=f"shares={np.round(c18g['shares'], 3).tolist()}, std={np.round(c18g['shares_std'], 3).tolist()}"),
        dict(name="16/8/4 份额归一、正",
             pass_=bool(min(c16["shares"]) > 0 and abs(sum(c16["shares"]) - 1) < 1e-6),
             detail=f"shares={np.round(c16['shares'], 3).tolist()}"),
        dict(name="宏观 n_std 对 core 不敏感（18 vs 16/8/4 差<1e-3）",
             pass_=bool(same_ns),
             detail=f"{c18['n_std']:.4f} vs {c16['n_std']:.4f}"),
        dict(name="宏观 q 对 core 不敏感（差<1e-3）",
             pass_=bool(same_q),
             detail=f"{c18['q']:.4f} vs {c16['q']:.4f}"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))
    for c in cases:
        t = np.array(c["trace"])
        for i in range(3):
            axes[0].plot(t[:, i], lw=1.0, label=f"{c['name']} {i}")
    axes[0].set_title("layer/generation shares")
    axes[0].legend(fontsize=7, ncol=2)
    axes[0].grid(alpha=0.3)
    names = [c["name"] for c in cases]
    axes[1].bar([i - 0.2 for i in range(3)], [c["n_std"] for c in cases],
                width=0.4, label="n_std")
    axes[1].bar([i + 0.2 for i in range(3)], [c["q"] for c in cases],
                width=0.4, label="q")
    axes[1].set_xticks(range(3))
    axes[1].set_xticklabels(names, fontsize=7, rotation=10)
    axes[1].set_title("macro observables: insensitive to micro core")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3, axis="y")
    fig.tight_layout()
    plot = os.path.join(RES, "zero_evo_cores.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    out = dict(cases=cases,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_evo_cores.json"), "w") as f:
        json.dump(out, f, indent=1)

    print("18 vs 16/8/4")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for c in cases:
        print(f"  {c['name']}: shares={np.round(c['shares'], 3).tolist()} "
              f"n_std={c['n_std']:.4f} q={c['q']:.4f} h_E={c['h_energy']:.3f} psi_coh={c['psi_coh']:.3f}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
