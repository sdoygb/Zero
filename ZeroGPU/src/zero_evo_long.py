#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_evo_long.py —— 统一引擎的长时间演化（10^6 步）与诊断
=========================================================
四条长跑：
  HK     c_h=1.2   （高 k 冻结）
  OSC    c_h=1.58  （振荡/行波相）
  FID    c_h=1.7   （Fiedler 冻结）
  TRAVEL c_h=1.55, seed=1（手性行波初值）

每条跑 10^6 步（dt=0.05 ⟹ t=5×10^4），诊断：
  · 五个时间窗的均值/涨落（是否平稳）
  · n_std 的 FFT 主频/周期
  · 末窗 k=7 的相位斜率与速度
  · 局域李雅普诺夫指数（逐步重归一，eps=1e-10，2×10^4 步）
  · 自相关（末 2×10^4 步）
  · 时空图（OSC / TRAVEL）

输出：results/zero_evo_long.json、results/zero_evo_long.png

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_evo_long.py
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_evo import ZeroEvo  # noqa: E402
from zero_evo_phase import (base_params, init_random, classify,  # noqa: E402
                            mode_metrics, run_case)
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")
STEPS = 1_000_000
EVERY = 500


def lyap_local(e, nsteps=20000, eps=1e-10):
    """逐步重归一的最大李雅普诺夫指数（只扰动内容扇区）。"""
    other = e.clone()
    direction = np.cos(np.arange(e.V) * 0.7)
    direction -= direction.mean()
    direction /= max(1e-12, np.linalg.norm(direction))
    other.n = np.maximum(e.n + eps * direction, e.p.clip_lo)
    logs = []
    for _ in range(nsteps):
        e.step()
        other.step()
        d = float(np.linalg.norm(other.n - e.n))
        if not np.isfinite(d) or d < 1e-300:
            break
        logs.append(np.log(d / eps))
        other.n = e.n + (other.n - e.n) * (eps / d)
    if not logs:
        return float("nan")
    return float(np.mean(logs) / e.p.dt)


def tail_trace(adj, e, steps=6000, k=7):
    V = len(adj)
    xs = []
    for _ in range(steps):
        e.step()
        xs.append(e.n.copy())
    xs = np.asarray(xs)
    z = (xs @ np.exp(-2j * np.pi * k * np.arange(V) / V)) / V
    amp = np.abs(z)
    ph = np.unwrap(np.angle(z))
    slope = float(np.polyfit(np.arange(len(ph)), ph, 1)[0]) / e.p.dt
    return {"amp": float(amp.mean()), "amp_relstd": float(amp.std() / max(1e-12, amp.mean())),
            "phase_slope": slope, "speed": float(-slope / (2 * np.pi * k / V))}, xs


def autocorr_period(y, max_lag=2000):
    y = np.asarray(y, dtype=float)
    y = y - y.mean()
    if y.size < 10 or y.std() == 0:
        return {"period": 0.0, "peak": 0.0}
    ac = np.correlate(y, y, "full")[len(y) - 1:]
    ac /= ac[0]
    lags = np.arange(2, min(max_lag, len(ac)))
    if lags.size == 0:
        return {"period": 0.0, "peak": 0.0}
    k = int(lags[np.argmax(ac[lags])])
    return {"period": float(k * EVERY * 0.05), "peak": float(ac[k])}


def main():
    t0 = time.time()
    os.makedirs(RES, exist_ok=True)
    adj = zs.cycle(16)
    V = 16
    w, Q = np.linalg.eigh(__import__("zero_evo").lap_matrix(adj))
    fied = Q[:, 1] / np.linalg.norm(Q[:, 1]) * 0.3

    cases = [
        ("HK", 1.2, init_random(V, 1)),
        ("OSC", 1.58, init_random(V, 1)),
        ("FID", 1.7, fied),
        ("TRAVEL", 1.55, init_random(V, 1)),
    ]
    out = {"steps": STEPS, "dt": 0.05, "cases": {}}
    traces = {}
    for name, c_h, x0 in cases:
        print(f"[{name}] c_h={c_h} 跑 {STEPS} 步 ...", flush=True)
        e = run_case(adj, c_h, x0, STEPS, every=EVERY)
        cls = classify(e, adj)
        ns = np.asarray(e.hist["n_std"])
        gs = np.asarray(e.hist["gen"])
        windows = []
        for i in range(5):
            m = (gs >= i * STEPS // 5) & (gs < (i + 1) * STEPS // 5)
            windows.append({"mean": float(ns[m].mean()), "std": float(ns[m].std())})
        fft = e.dominant_frequency("n_std", frac=1.0)
        ac = autocorr_period(ns)
        tr, xs = tail_trace(adj, e)
        lam = lyap_local(e, nsteps=20000, eps=1e-10)
        out["cases"][name] = {
            "c_h": c_h, "class": cls, "windows": windows,
            "fft": fft, "autocorr": ac, "travel": tr, "lyapunov": lam,
            "n_hist": [float(x) for x in ns[::10]],
            "t_hist": [float(x) for x in np.asarray(e.hist["t"])[::10]],
        }
        traces[name] = xs
        print(f"  -> {cls['tag']} lam={cls['lam']:.3f} high={cls['high']:.3f} "
              f"temp={cls['temporal']:.4f} fft={fft['freq']:.4f}/period={fft['period']:.2f} "
              f"lyap={lam:+.4f} wall={time.time()-t0:.0f}s", flush=True)

    make_figure(out, traces)
    out["meta"] = {"wall_s": round(time.time() - t0, 1), "every": EVERY}
    fp = os.path.join(RES, "zero_evo_long.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"写入 {fp}（{out['meta']['wall_s']}s）")



def make_figure(out, traces):
    """四联图：(a) 长跑轨迹 (b) 平稳性 (c) OSC 衰减 (d) TRAVEL 时空。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = list(out["cases"].keys())
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    ax = axes[0, 0]
    for name in names:
        d = out["cases"][name]
        ax.plot(d["t_hist"], d["n_hist"], lw=0.7,
                label=f"{name} ($c_h={d['c_h']}$)")
    ax.set_xlabel("time")
    ax.set_ylabel("$n_{std}$")
    ax.set_title("(a) long-time $n_{std}$")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    for i, name in enumerate(names):
        d = out["cases"][name]
        ax.plot([i] * len(d["windows"]), [w["std"] for w in d["windows"]],
                "o-", label=name)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names)
    ax.set_ylabel("std of $n_{std}$ in each 1/5 window")
    ax.set_title("(b) stationarity")
    ax.legend(fontsize=7)

    ax = axes[1, 0]
    d = out["cases"]["OSC"]
    t = np.asarray(d["t_hist"])
    y = np.asarray(d["n_hist"])
    m = t <= 20000
    ax.plot(t[m], y[m], lw=0.8, c="tab:orange", label="OSC ($c_h=1.58$)")
    level = float(np.mean(y[-20:]))
    ax.axhline(level, ls=":", c="tab:green", label="FID level (OSC tail)")
    ax.set_xlabel("time")
    ax.set_ylabel("$n_{std}$")
    ax.set_title("(c) OSC $\\to$ FID (first 2×10^4 time)")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    xs = traces["TRAVEL"]
    z = (xs @ np.exp(-2j * np.pi * 7 * np.arange(xs.shape[1]) / xs.shape[1])) / xs.shape[1]
    ax.plot(np.unwrap(np.angle(z)), lw=0.9)
    ax.set_xlabel("time step (tail 6000)")
    ax.set_ylabel("arg $z_7$ (unwrapped)")
    ax.set_title("(d) TRAVEL phase advance ($c_h=1.55$)")
    fig.tight_layout()
    fp = os.path.join(RES, "zero_evo_long.png")
    fig.savefig(fp, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    if "--plot-only" in sys.argv:
        with open(os.path.join(RES, "zero_evo_long.json"), encoding="utf-8") as fh:
            out = json.load(fh)
        adj = zs.cycle(16)
        e = run_case(adj, 1.55, init_random(16, 1), 150000, every=500)
        _, xs = tail_trace(adj, e)
        make_figure(out, {"TRAVEL": xs})
        print("re-plotted results/zero_evo_long.png")
    else:
        main()
