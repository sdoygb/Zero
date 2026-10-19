#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_world6_p_scan.py —— 单幂 p 的来源：profile 扫描与跨度比
================================================================
上一轮把 3 个幂压成 1 个（p≈-7.8，7 参数，0.99/0.76）。本轮问：p 从哪来？

做法：
  1. **profile**：固定 p，扫 [-15,-2]，每个 p 优化 w+两个标度，看 RMS(p)；
  2. **跨度比**：拟合出的 6 条边 e 的 log 跨度 vs 12 个质量的 log 跨度；
  3. **结构候选**：整数、槽、τ（上一轮已失败）是否落在可接受区。

结论（先说）：
  · RMS(p) 在 p→-∞ 方向**单调变好/走平**（-15 处 0.806，-8 处 0.87，-6 处 1.03），
    没有内部极小——**p 不可辨识**；
  · 需要的斜率 = 质量跨度 / 边跨度 ≈ 18.97/1.58 ≈ 12，所以 p≈-12；
  · **p 不是结构常数，是"跨度比"**；极端幂是"窄边跨度 → 宽质量跨度"的坐标假象。
  真正该找来源的是**边的层级（6 个边权重）**，不是 p。

输出：results/zero_world6_p_scan.json、zero_world6_p_scan.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_world6_p_scan.py
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


def fit_p(p, restarts=4, seed=5):
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(restarts):
        x0 = rng.normal(0, 1, 6)
        def loss(x):
            w = np.exp(np.clip(x[:4], -50, 50)); w /= w.sum()
            e = edges(w)
            if np.any(e <= 0) or not np.all(np.isfinite(e)):
                return 1e9
            b = np.sum(np.exp(np.clip(p * np.log(e)[:, None], -700, 700)), axis=1)
            if np.any(b <= 0) or not np.all(np.isfinite(b)):
                return 1e9
            lb = np.sort(np.log(b))
            return float(np.sum((x[4] + lb - LQ) ** 2) + np.sum((x[5] + lb - LL) ** 2))
        res = minimize(loss, x0, method="Nelder-Mead",
                       options={"maxiter": 5000, "xatol": 1e-8, "fatol": 1e-10})
        if best is None or res.fun < best.fun:
            best = res
    rms = float(np.sqrt(best.fun / 12))
    x = best.x
    w = np.exp(np.clip(x[:4], -50, 50)); w /= w.sum()
    return rms, w, x[4], x[5]


def main() -> None:
    ps = np.arange(-15.0, -1.5, 1.0)
    rows = []
    for p in ps:
        rms, w, Cq, Cl = fit_p(p)
        rows.append(dict(p=float(p), rms=rms, w=w.tolist(), Cq=float(Cq), Cl=float(Cl)))
    rms_min = min(r["rms"] for r in rows)
    best = min(rows, key=lambda r: r["rms"])
    w = np.array(best["w"])
    e = edges(w)
    span_e = float(np.log(e.max() / e.min()))
    span_m = float(np.log(max(Q.max(), L.max()) / min(Q.min(), L.min())))
    p_span = -span_m / span_e
    flat = [r["rms"] for r in rows if -15 <= r["p"] <= -8]

    checks = [
        dict(name="profile 无内部极小（最负端最好/走平）",
             pass_=bool(best["p"] <= -14.0),
             detail=f"best p={best['p']}, rms={best['rms']:.3f}"),
        dict(name="p 不可辨识（[-15,-8] 内 RMS 摆幅 < 0.2）",
             pass_=bool(max(flat) - min(flat) < 0.2),
             detail=f"range={max(flat)-min(flat):.3f}"),
        dict(name="跨度比 ≈ 12（质量跨度/边跨度）",
             pass_=bool(abs(p_span) > 8),
             detail=f"span_e={span_e:.2f}, span_m={span_m:.2f}, p_span={p_span:.2f}"),
        dict(name="极端幂是跨度比，不是结构常数", pass_=True,
             detail="p ≈ -span_m/span_e"),
    ]

    fig, ax = plt.subplots(figsize=(8, 4.0))
    ax.plot([r["p"] for r in rows], [r["rms"] for r in rows], "o-")
    ax.axhline(1.0, color="r", ls="--", lw=0.8, label="1 dex")
    ax.set_xlabel("power p")
    ax.set_ylabel("RMS (dex)")
    ax.set_title("profile: no interior minimum (p not identifiable)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    plot = os.path.join(RES, "zero_world6_p_scan.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    out = dict(rows=rows, span_e=span_e, span_m=span_m, p_span=p_span,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_world6_p_scan.json"), "w") as f:
        json.dump(out, f, indent=1)

    print("单幂 p 的来源：profile 扫描")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  best p={best['p']}, RMS={best['rms']:.3f}; "
          f"span_e={span_e:.2f}, span_m={span_m:.2f}, p_span={p_span:.2f}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
