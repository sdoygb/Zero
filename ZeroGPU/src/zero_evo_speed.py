#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_evo_speed.py —— 高 k 带的行波速度 vs 几何响应速度 c_h
==========================================================
背景：`zero_evo_phase.py` 用 n_std 的时间涨落分类；匀速平移的行波
      n_std 恒定，会被误标成“冻结”。本脚本直接测相位斜率，把“带”
      和“是否平移”分开：

  · k=7 复振幅 z_7(t) = Σ_j n_j(t) e^{-2πi·7j/16}/16
  · 相位斜率 ω = d arg z_7/dt；速度 v = -ω/(2π·7/16)（sites/time）
  · 振幅 |z_7| 与相对涨落（判断是否匀速）

输出：results/zero_evo_speed.json、results/zero_evo_speed.png

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_evo_speed.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_evo import ZeroEvo  # noqa: E402
from zero_evo_phase import base_params, init_random  # noqa: E402
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")
V = 16


def measure(c_h, seed=1, warm=60000, tail=4000):
    adj = zs.cycle(V)
    e = ZeroEvo(adj, base_params(c_h=c_h), n0=1.0 + init_random(V, seed))
    e.run(warm, every=500)
    xs = []
    for _ in range(tail):
        e.step()
        xs.append(e.n.copy())
    xs = np.asarray(xs)
    ts = np.arange(len(xs)) * e.p.dt
    z = (xs @ np.exp(-2j * np.pi * 7 * np.arange(V) / V)) / V
    amp = np.abs(z)
    rel = float(amp.std() / max(1e-12, amp.mean()))
    ph = np.unwrap(np.angle(z))
    slope = float(np.polyfit(ts, ph, 1)[0])
    v = float(-slope / (2 * np.pi * 7 / V))
    # 低 k 对照
    z1 = (xs @ np.exp(-2j * np.pi * 1 * np.arange(V) / V)) / V
    return {"c_h": c_h, "seed": seed,
            "k7_amp": float(amp.mean()), "k7_relstd": rel,
            "k7_phase_slope": slope, "k7_speed": v,
            "k1_amp": float(np.abs(z1).mean()),
            "n_std_mean": float(np.mean([np.std(x) for x in xs]))}


def main():
    os.makedirs(RES, exist_ok=True)
    rows = []
    for c in [0.8, 1.0, 1.2, 1.35, 1.45, 1.5, 1.52, 1.54, 1.55, 1.58, 1.6]:
        r = measure(c, seed=1)
        rows.append(r)
        print(f"c_h={c:4.2f}  k7 amp={r['k7_amp']:.3f} rel={r['k7_relstd']:.3f} "
              f"slope={r['k7_phase_slope']:+.3f} v={r['k7_speed']:+.3f} "
              f"k1 amp={r['k1_amp']:.3f}", flush=True)
    # 方向：换初值
    dirs = []
    for c in (1.2, 1.55):
        for s in (1, 4):
            r = measure(c, seed=s)
            dirs.append(r)
            print(f"direction c_h={c} seed={s}: v={r['k7_speed']:+.3f}", flush=True)

    # 渐近点：1.55 的交接态在 150k 后分到纯 k7 行波或 FID；这里取 seed 1 的行波支
    asym = []
    for c in (1.52, 1.55):
        r = measure(c, seed=1, warm=150000)
        asym.append(r)
        print(f"asymptotic c_h={c}: k7 amp={r['k7_amp']:.3f} rel={r['k7_relstd']:.3f} "
              f"v={r['k7_speed']:+.3f}", flush=True)
    out = {"speed": rows, "direction": dirs, "asymptotic": asym,
           "note": "n_std 恒定的“冻结”标签会把匀速平移误判；本表用相位斜率识别。"}
    fp = os.path.join(RES, "zero_evo_speed.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
    cs = [r["c_h"] for r in rows]
    ax = axes[0]
    ax.plot(cs, [abs(r["k7_speed"]) for r in rows], "o-", label="60k warm-up")
    ax.plot([r["c_h"] for r in asym], [abs(r["k7_speed"]) for r in asym], "*",
            ms=14, label="150k warm-up (asymptotic branch)")
    ax.set_xlabel("$c_h$")
    ax.set_ylabel("$|v|$ (sites/time)")
    ax.set_title("(a) traveling speed of the $k=7$ carrier")
    ax.legend(fontsize=8)
    ax = axes[1]
    ax.plot(cs, [r["k7_amp"] for r in rows], "o-", label="$|z_7|$")
    ax.plot(cs, [r["k1_amp"] for r in rows], "s--", label="$|z_1|$")
    ax.set_xlabel("$c_h$")
    ax.set_ylabel("Fourier amplitude")
    ax.set_title("(b) carrier weight: high-$k$ to low-$k$ handover")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "zero_evo_speed.png"), dpi=150)
    print(f"写入 {fp} 与 zero_evo_speed.png")


if __name__ == "__main__":
    main()
