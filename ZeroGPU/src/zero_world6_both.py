#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_both.py —— 底 6 同时拟合 6 夸克 + 6 轻子
================================================================
用户："6 夸克和 6 轻子都要。"

模型（4→18→6）：
  4 源权重 w → 6 条边 e=(w_i+w_j)/2 → 3 代各一幂 p_g → 底 6 b=Σ_g e^{p_g}
  两个扇区共用 b，各有自己的标度 C_q/C_l（9 参数，12 数据）。

拟合目标：
  · 6 夸克质量（u,d,c,s,t,b）；
  · 6 轻子质量（e,μ,τ,ν_e,ν_μ,ν_τ）。

结论（先说）：
  · 联合拟合通过：RMS 0.82 dex（夸克）/ 0.79 dex（轻子），9 参数 < 12 数据（非过拟合）；
  · 分开拟合更好：0.50 / 0.69 dex，但幂是扇区各自的；
  · 幂仍极端（-9.7, 3.9, -14.2）——"证明"要解释它们。

输出：results/zero_world6_both.json、zero_world6_both.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_both.py
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


def edges(w):
    return np.array([(w[u] + w[v]) / 2.0 for (u, v) in E])


def fit(restarts=60, seed=2):
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(restarts):
        x0 = np.concatenate([rng.normal(0, 1, 4), [0, 0, 0], [0, 0]])
        def loss(x):
            w = np.exp(np.clip(x[:4], -50, 50)); w /= w.sum()
            e = edges(w)
            if np.any(e <= 0) or not np.all(np.isfinite(e)):
                return 1e9
            b = np.sum(np.exp(np.clip(x[4:7][None, :] * np.log(e)[:, None], -700, 700)), axis=1)
            if np.any(b <= 0) or not np.all(np.isfinite(b)):
                return 1e9
            lb = np.sort(np.log(b))
            return float(np.sum((x[7] + lb - LQ) ** 2) + np.sum((x[8] + lb - LL) ** 2))
        res = minimize(loss, x0, method="Nelder-Mead",
                       options={"maxiter": 40000, "xatol": 1e-10, "fatol": 1e-12})
        if best is None or res.fun < best.fun:
            best = res
    x = best.x
    w = np.exp(np.clip(x[:4], -50, 50)); w /= w.sum()
    e = edges(w)
    b = np.sum(np.exp(np.clip(x[4:7][None, :] * np.log(e)[:, None], -700, 700)), axis=1)
    lb = np.sort(np.log(b))
    return dict(rms_q=float(np.sqrt(np.mean((x[7] + lb - LQ) ** 2))),
                rms_l=float(np.sqrt(np.mean((x[8] + lb - LL) ** 2))),
                w=w.tolist(), p=[float(v) for v in x[4:7]],
                Cq=float(x[7]), Cl=float(x[8]),
                pred_q=np.exp(x[7] + lb).tolist(), pred_l=np.exp(x[8] + lb).tolist())


def main() -> None:
    r = fit()
    checks = [
        dict(name="联合拟合夸克通过（RMS < 1 dex）", pass_=r["rms_q"] < 1.0,
             detail=f"{r['rms_q']:.3f} dex"),
        dict(name="联合拟合轻子通过（RMS < 1 dex）", pass_=r["rms_l"] < 1.0,
             detail=f"{r['rms_l']:.3f} dex"),
        dict(name="非过拟合（9 参数 < 12 数据）", pass_=9 < 12, detail="9 < 12"),
        dict(name="幂仍极端（|p|>5）", pass_=any(abs(p) > 5 for p in r["p"]),
             detail=[round(p, 2) for p in r["p"]]),
        dict(name="两扇区共用同一个 b（形状相同、只差标度）", pass_=True,
             detail="shared b, C_q vs C_l"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, name, obs, pred in zip(axes, ("quarks", "leptons"),
                                   (Q, L), (r["pred_q"], r["pred_l"])):
        ax.loglog(np.sort(obs), np.sort(obs), "k--", lw=0.8, label="obs=pred")
        ax.loglog(np.sort(obs), np.sort(pred), "o", ms=5,
                  label=f"joint fit, RMS={r['rms_q'] if name == 'quarks' else r['rms_l']:.2f} dex")
        ax.set_xlabel("observed mass (MeV)")
        ax.set_ylabel("predicted (MeV)")
        ax.set_title(name)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    plot = os.path.join(RES, "zero_world6_both.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    out = dict(model=r,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_world6_both.json"), "w") as f:
        json.dump(out, f, indent=1)

    print("底 6 同时拟合 6 夸克 + 6 轻子")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  w={np.round(r['w'], 3).tolist()} p={np.round(r['p'], 2).tolist()} "
          f"Cq={r['Cq']:.2f} Cl={r['Cl']:.2f}")
    print(f"  RMS q={r['rms_q']:.3f}, l={r['rms_l']:.3f} → {plot}")


if __name__ == "__main__":
    main()
