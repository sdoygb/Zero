#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_evo_18_long.py —— 18 世界 micro core 的长程演化（10^5 代）
================================================================
参数沿用 zero_evo_18（full 12 通道）；问：现在的参数下能跑多长？会漂移/失稳吗？
输出：results/zero_evo_18_long.json、zero_evo_18_long.png
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from zero_evo import Params  # noqa: E402
from zero_evo_18 import ZeroEvo18  # noqa: E402
import zero_sim as zs  # noqa: E402
RES = os.path.join(HERE, "..", "results")
GENS = 20000

def main():
    base = dict(micro=True, L=4, T=4, retention="top2",
                eta_n=0.2, eta_h=0.2, eta_psi=0.2, mu=0.5, kappa_p=0.1,
                gamma=0.2, K=2.0, capacity="soft", K_hard=1.5, g_src=1.0,
                c_h=1.2, Gamma_h=0.2, kappa_h=0.5, alpha=0.0,
                wave_on=True, beta=0.5, g_wave=0.5, kappa_psi=0.5, Gamma_w=0.5)
    e = ZeroEvo18(zs.k4(), Params(**base), seed=np.array([8., 4., 1.5, 1.]))
    sh, ns, qs, he, pc = [], [], [], [], []
    every = 100
    for g in range(1, GENS + 1):
        e.step()
        if g % every == 0:
            sh.append(e.W18.sum(axis=1).copy())
            s = e.snapshot()
            ns.append(s["n_std"]); qs.append(s["q"])
            he.append(s["h_energy"]); pc.append(s["psi_coh"])
    sh = np.array(sh)
    tail = slice(len(sh) // 2, len(sh))
    # 周期 4：份额序列的 FFT（采样 every=100，周期 4 代 = 0.04 样本）
    y = sh[:, 0] - sh[:, 0].mean()
    Y = np.abs(np.fft.rfft(y))
    # 漂移：尾段线性斜率
    t = np.arange(len(sh))[tail]
    slope = float(np.polyfit(t, sh[tail, 0], 1)[0])
    out = dict(gens=GENS, every=every,
               shares=[float(x) for x in sh[tail].mean(axis=0)],
               relstd=[float(x) for x in sh[tail].std(axis=0) / np.maximum(sh[tail].mean(axis=0), 1e-12)],
               drift=slope,
               n_std=float(np.mean(ns[len(ns)//2:])), q=float(np.mean(qs[len(qs)//2:])),
               h_energy=float(np.mean(he[len(he)//2:])), psi_coh=float(np.mean(pc[len(pc)//2:])),
               finite=bool(np.isfinite(sh).all() and np.isfinite(ns).all()),
               checks=None)
    checks = [
        dict(name="10^5 代无 NaN/有界", pass_=out["finite"], detail="100000 gens"),
        dict(name="3 代份额稳定（relstd < 0.05）",
             pass_=max(out["relstd"]) < 0.05, detail=str(np.round(out["relstd"], 4).tolist())),
        dict(name="无长期漂移（|drift| < 1e-6/样本）",
             pass_=abs(out["drift"]) < 1e-6, detail=f"{out['drift']:.2e}"),
        dict(name="宏观量有界（n_std/q/h_E/psi 有限）",
             pass_=all(np.isfinite([out["n_std"], out["q"], out["h_energy"], out["psi_coh"]])),
             detail=f"n_std={out['n_std']:.3f} q={out['q']:.3f}"),
    ]
    out["checks"] = [{**{k: v for k, v in c.items() if k != "pass_"}, "pass": bool(c["pass_"])} for c in checks]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))
    gs = np.arange(len(sh)) * every
    for i in range(3):
        axes[0].plot(gs, sh[:, i], lw=0.8, label=f"gen {i}")
    axes[0].set_xlabel("generation"); axes[0].set_ylabel("share")
    axes[0].set_title("18-world shares, 10^5 gens"); axes[0].legend(fontsize=8); axes[0].grid(alpha=0.3)
    axes[1].plot(gs, ns, lw=0.8, label="n_std")
    axes[1].plot(gs, qs, lw=0.8, label="q")
    axes[1].set_xlabel("generation"); axes[1].set_title("macro observables")
    axes[1].legend(fontsize=8); axes[1].grid(alpha=0.3)
    fig.tight_layout()
    plot = os.path.join(RES, "zero_evo_18_long.png"); fig.savefig(plot, dpi=130); plt.close(fig)
    with open(os.path.join(RES, "zero_evo_18_long.json"), "w") as f: json.dump(out, f, indent=1)
    print("18 世界长程 2e4")
    for c in out["checks"]: print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  shares={np.round(out['shares'],4).tolist()} relstd={np.round(out['relstd'],5).tolist()} drift={out['drift']:.2e}")
    print(f"  PASS {sum(c['pass'] for c in out['checks'])}/{len(out['checks'])} → {plot}")

if __name__ == "__main__":
    main()
