#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_evo_phase.py —— 统一引擎的相图：几何回写如何反转图案选择
=============================================================
跑四组可复跑实验，落 `results/zero_evo_phase.json` 与两张图：

  P1 机制对照     cycle16：无几何 / 只给源 / 只回写 / 双向
  P2 c_h 相图     c_h = 0.8 … 2.0，3 个初值；末态分类（HK / OSC / FID）
  P3 迟滞         c_h = 1.52 / 1.55 / 1.58，三种初值（random / fiedler / highest）
  P4 行波         c_h = 1.55（高 k）、1.58（低 k），10 / 3 个初值；测 k=7 相位斜率
  P5 换图         cycle16 / path16 / torus4x4 / K4 / Heawood，c_h = 1.2 / 1.6
  P6 层间回写     α 扫描（micro 扇区接入强度），看宏观图案被拉到哪里

口径：
  · "末态" = 最后 10% 步的瞬时 n；模式功率 = |Q^T (n/n.sum - 1/V)|^2 归一化
  · HK = 高 k 冻结（λ_peak > 2.5 且时间涨落 < 0.02）
  · FID = Fiedler 冻结（λ_peak < 1 且时间涨落 < 0.02）
  · OSC = 时间涨落 ≥ 0.02（极限环/行波/振荡）
  · 图：results/zero_evo_phase.png, results/zero_evo_travel.png

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_evo_phase.py
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_evo import ZeroEvo, Params, lap_matrix  # noqa: E402
import zero_sim as zs  # noqa: E402
from z0_engine import torus_adj  # noqa: E402

RES = os.path.join(HERE, "..", "results")
OUT: dict = {}


# ────────────────────────────── 工具 ──────────────────────────────
def mode_metrics(adj, vec):
    """归一化向量的模式统计（向量先除以自身和，避免 1e180 量级的平方溢出）。"""
    V = len(adj)
    w, Q = np.linalg.eigh(lap_matrix(adj))
    vec = np.asarray(vec, dtype=float)
    s = vec.sum()
    if s <= 0 or not np.isfinite(s):
        return {"lam": float("nan"), "high": float("nan"), "low": float("nan"),
                "pr": float("nan")}
    p = vec / s
    x = p - 1.0 / V
    c = Q.T @ x
    P = c * c
    if P.sum() <= 0:
        return {"lam": 0.0, "high": 0.0, "low": 0.0, "pr": float("nan")}
    P /= P.sum()
    pr = float(1.0 / max(1e-300, (p * p).sum()))
    return {"lam": float(P @ w), "high": float(P[w > 2.5].sum()),
            "low": float(P[w < 1.0].sum()), "pr": pr}


def classify(e, adj, tail=50):
    """末态分类 + 瞬时模式统计。"""
    ns = np.asarray(e.hist["n_std"][-tail:])
    temporal = float(ns.std())
    m = mode_metrics(adj, e.n)
    if temporal >= 0.02:
        tag = "OSC"
    elif m["high"] > 0.5:
        tag = "HK"
    elif m["low"] > 0.5:
        tag = "FID"
    else:
        tag = "MIX"
    return {"tag": tag, "temporal": temporal, "lam": m["lam"],
            "high": m["high"], "low": m["low"], "pr": m["pr"],
            "n_mean": float(e.n.mean()), "n_std_mean": float(ns.mean())}


def base_params(**kw):
    p = dict(micro=False, mu=0.5, kappa_p=0.1, gamma=0.2, K=2.0,
             capacity="soft", K_hard=1.5, g_src=1.0, c_h=1.2,
             Gamma_h=0.2, kappa_h=0.5, wave_on=False, poisson_on=False)
    p.update(kw)
    return Params(**p)


def init_random(V, seed, amp=0.05):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=V)
    x -= x.mean()
    x *= amp / max(1e-12, np.abs(x).max())
    return x


def run_case(adj, c_h, x0, steps, every=200, **kw):
    e = ZeroEvo(adj, base_params(c_h=c_h, **kw), n0=1.0 + x0)
    e.run(steps, every=every)
    return e


# ────────────────────────────── P1 机制 ──────────────────────────────
def p1_controls(adj):
    V = len(adj)
    x0 = init_random(V, 1)
    cases = {
        "no_geom": dict(g_src=0.0, kappa_h=0.0),
        "source_only": dict(g_src=1.0, kappa_h=0.0),
        "backreact_only": dict(g_src=0.0, kappa_h=0.5),
        "both": dict(g_src=1.0, kappa_h=0.5),
    }
    out = {}
    for name, kw in cases.items():
        e = run_case(adj, 1.2, x0, 60000, **kw)
        out[name] = classify(e, adj)
        print(f"  P1 {name:14s} -> {out[name]['tag']:4s} lam={out[name]['lam']:.3f} "
              f"high={out[name]['high']:.3f} temp={out[name]['temporal']:.4f}")
    return out


# ────────────────────────────── P2 相图 ──────────────────────────────
def p2_phase(adj, c_list, seeds=(0, 1), steps=100000):
    rows = []
    for c in c_list:
        for s in seeds:
            x0 = init_random(len(adj), s)
            e = run_case(adj, c, x0, steps)
            r = classify(e, adj)
            r.update({"c_h": c, "seed": s})
            rows.append(r)
            print(f"  P2 c_h={c:4.2f} s={s} -> {r['tag']:4s} lam={r['lam']:.3f} "
                  f"high={r['high']:.3f} temp={r['temporal']:.4f}")
    return rows


# ────────────────────────────── P3 迟滞 ──────────────────────────────
def p3_hysteresis(adj, c_list=(1.52, 1.55, 1.58), steps=100000):
    V = len(adj)
    w, Q = np.linalg.eigh(lap_matrix(adj))
    fied = Q[:, 1]
    fied = fied / np.linalg.norm(fied) * 0.3
    high = Q[:, -1]
    high = high / np.linalg.norm(high) * 0.3
    inits = {"random": init_random(V, 1), "fiedler": fied, "highest": high}
    rows = []
    for c in c_list:
        for name, x0 in inits.items():
            e = run_case(adj, c, x0, steps)
            r = classify(e, adj)
            r.update({"c_h": c, "init": name})
            rows.append(r)
            print(f"  P3 c_h={c:4.2f} init={name:9s} -> {r['tag']:4s} "
                  f"lam={r['lam']:.3f} high={r['high']:.3f} temp={r['temporal']:.4f}")
    return rows


# ────────────────────────────── P4 行波 ──────────────────────────────
def _tail_travel(adj, e, k=7, steps=3000):
    V = len(adj)
    xs, ts = [], []
    for _ in range(steps):
        e.step()
        xs.append(e.n.copy())
        ts.append(e.t)
    xs = np.asarray(xs)
    ts = np.asarray(ts)
    z = (xs @ np.exp(-2j * np.pi * k * np.arange(V) / V)) / V
    amp = np.abs(z)
    rel = float(amp.std() / max(1e-12, amp.mean()))
    ph = np.unwrap(np.angle(z))
    slope = float(np.polyfit(ts, ph, 1)[0])
    v = float(-slope / (2 * np.pi * k / V))
    return {"k": k, "amp": float(amp.mean()), "amp_relstd": rel,
            "phase_slope": slope, "speed": v}


def p4_travel(adj, c_h, seeds, steps=150000):
    rows = []
    for s in seeds:
        x0 = init_random(len(adj), s)
        e = run_case(adj, c_h, x0, steps)
        r = classify(e, adj)
        r.update({"c_h": c_h, "seed": s})
        if r["tag"] == "OSC":
            tr = _tail_travel(adj, e)
            r["travel"] = tr
            print(f"  P4 c_h={c_h:.2f} s={s} -> OSC amp={tr['amp']:.3f} "
                  f"rel={tr['amp_relstd']:.3f} v={tr['speed']:+.3f}")
        else:
            print(f"  P4 c_h={c_h:.2f} s={s} -> {r['tag']:4s} lam={r['lam']:.3f}")
        rows.append(r)
    return rows


# ────────────────────────────── P5 换图 ──────────────────────────────
def p5_graphs(graphs, c_list=(1.2, 1.6), steps=40000):
    rows = []
    for name, adj in graphs:
        for c in c_list:
            for s in (1, 2):
                x0 = init_random(len(adj), s)
                e = run_case(adj, c, x0, steps)
                r = classify(e, adj)
                r.update({"graph": name, "c_h": c, "seed": s})
                rows.append(r)
                print(f"  P5 {name:9s} c_h={c:.2f} s={s} -> {r['tag']:4s} "
                      f"lam={r['lam']:.3f} high={r['high']:.3f} temp={r['temporal']:.4f}")
    return rows


# ────────────────────────────── P6 α 扫描 ──────────────────────────────
def p6_alpha(adj, alphas=(0.0, 0.001, 0.005, 0.02, 0.1, 0.5, 1.0), steps=12000):
    V = len(adj)
    x0 = init_random(V, 1)
    rows = []
    for a in alphas:
        e = ZeroEvo(adj, base_params(c_h=1.55, micro=True, alpha=a, eta_n=0.2),
                    seed=np.array([2.0 ** (V - 1 - i) for i in range(V)]),
                    n0=1.0 + x0)
        e.run(steps, every=10)
        nm = mode_metrics(adj, e.n)
        rm = mode_metrics(adj, e.record)
        rows.append({"alpha": a, "n_lam": nm["lam"], "n_high": nm["high"],
                     "n_pr": nm["pr"], "rec_lam": rm["lam"], "rec_high": rm["high"],
                     "rec_pr": rm["pr"], "lock": float(e.hist["lock_cn"][-1])})
        print(f"  P6 alpha={a:5.3f} n(lam,high,PR)=({nm['lam']:.2f},{nm['high']:.2f},"
              f"{nm['pr']:.2f}) rec(lam,high,PR)=({rm['lam']:.2f},{rm['high']:.2f},"
              f"{rm['pr']:.2f}) lock={rows[-1]['lock']:.2f}")
    return rows


# ────────────────────────────── 画图 ──────────────────────────────
def make_plots(phase_rows, travel_rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cs = sorted(set(r["c_h"] for r in phase_rows))
    lam = [np.mean([r["lam"] for r in phase_rows if r["c_h"] == c]) for c in cs]
    lam_sd = [np.std([r["lam"] for r in phase_rows if r["c_h"] == c]) for c in cs]
    temp = [np.mean([r["temporal"] for r in phase_rows if r["c_h"] == c]) for c in cs]

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5))
    ax = axes[0, 0]
    ax.errorbar(cs, lam, yerr=lam_sd, marker="o", capsize=3)
    ax.axhline(0.27, ls=":", c="tab:green", label="Fiedler (0.27)")
    ax.axhline(3.85, ls=":", c="tab:red", label="high-k (3.85)")
    ax.set_xlabel("geometry speed $c_h$")
    ax.set_ylabel("selected $\\lambda_{peak}$")
    ax.set_title("(a) selection vs geometry response speed")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.semilogy(cs, np.maximum(temp, 1e-8), marker="s", c="tab:purple")
    ax.set_xlabel("$c_h$")
    ax.set_ylabel("temporal std of $n_{std}$ (tail)")
    ax.set_title("(b) oscillation / traveling window")

    ax = axes[1, 0]
    for c, col in [(1.2, "tab:red"), (1.58, "tab:orange"), (1.7, "tab:green")]:
        rr = [r for r in phase_rows if r["c_h"] == c]
        if rr:
            ax.scatter([r["lam"] for r in rr], [r["high"] for r in rr], label=f"$c_h={c}$")
    ax.set_xlabel("$\\lambda_{peak}$")
    ax.set_ylabel("high-k fraction")
    ax.set_title("(c) final states (HK / OSC / FID)")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    # 一个 HK 末态与一个 FID 末态的模式谱
    for c, col in [(1.2, "tab:red"), (1.58, "tab:orange"), (1.7, "tab:green")]:
        rr = [r for r in phase_rows if r["c_h"] == c]
        if rr:
            ax.bar([c], [np.mean([r["high"] for r in rr])], color=col, width=0.4)
    ax.set_xlabel("$c_h$")
    ax.set_ylabel("mean high-k fraction")
    ax.set_title("(d) high-k fraction vs $c_h$")
    fig.tight_layout()
    fp = os.path.join(RES, "zero_evo_phase.png")
    fig.savefig(fp, dpi=150)
    plt.close(fig)

    # 行波图：取一个 OSC 的时空图
    osc = [r for r in travel_rows if r.get("travel") and "amp" in r["travel"]]
    if osc:
        r = max(osc, key=lambda z: z["travel"]["amp"])
        e = run_case(zs.cycle(16), r["c_h"], init_random(16, r["seed"]), 200000)
        # 取最后 600 步时空
        xs = []
        for _ in range(600):
            e.step()
            xs.append(e.n.copy())
        xs = np.array(xs).T
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
        ax = axes[0]
        im = ax.imshow(xs, aspect="auto", origin="lower", cmap="RdBu_r")
        ax.set_xlabel("time step (last 600)")
        ax.set_ylabel("site")
        ax.set_title(f"traveling wave, $c_h={r['c_h']}$")
        fig.colorbar(im, ax=ax, fraction=0.046)
        ax = axes[1]
        z = (xs.T @ np.exp(-2j * np.pi * 7 * np.arange(16) / 16)) / 16
        ax.plot(np.unwrap(np.angle(z)), lw=0.8)
        ax.set_xlabel("time step")
        ax.set_ylabel("arg $z_7$ (unwrapped)")
        ax.set_title("phase advance of the $k=7$ mode")
        fig.tight_layout()
        fp = os.path.join(RES, "zero_evo_travel.png")
        fig.savefig(fp, dpi=150)
        plt.close(fig)


# ────────────────────────────── main ──────────────────────────────
def main():
    t0 = time.time()
    os.makedirs(RES, exist_ok=True)
    cycle = zs.cycle(16)

    print("=" * 78)
    print("P1 机制对照（cycle16, c_h=1.2, 60k 步）")
    print("=" * 78)
    OUT["P1_controls"] = p1_controls(cycle)

    print("=" * 78)
    print("P2 相图（cycle16, c_h=0.8..2.0, 2 初值, 100k 步）")
    print("=" * 78)
    OUT["P2_phase"] = p2_phase(cycle, [0.8, 1.0, 1.2, 1.35, 1.45, 1.5, 1.54, 1.58, 1.6, 1.8])

    print("=" * 78)
    print("P3 迟滞（三种初值）")
    print("=" * 78)
    OUT["P3_hysteresis"] = p3_hysteresis(cycle)

    print("=" * 78)
    print("P4 行波（c_h=1.55, 10 初值；c_h=1.58, 3 初值）")
    print("=" * 78)
    OUT["P4_travel_155"] = p4_travel(cycle, 1.55, range(10))
    OUT["P4_travel_158"] = p4_travel(cycle, 1.58, range(3))

    print("=" * 78)
    print("P5 换图（c_h=1.2 / 1.6）")
    print("=" * 78)
    graphs = [("cycle16", zs.cycle(16)), ("path16", zs.path(16)),
              ("torus4x4", torus_adj(2, 4)), ("K4", zs.k4()),
              ("Heawood", zs.heawood())]
    OUT["P5_graphs"] = p5_graphs(graphs)

    print("=" * 78)
    print("P6 层间回写 α 扫描（c_h=1.55, micro=True）")
    print("=" * 78)
    OUT["P6_alpha"] = p6_alpha(cycle)

    make_plots(OUT["P2_phase"], OUT["P4_travel_155"])
    OUT["meta"] = {"wall_s": round(time.time() - t0, 1),
                   "engine": "zero_evo.py",
                   "note": "exploratory couplings: logistic saturation, soft capacity, scalar geometry"}
    fp = os.path.join(RES, "zero_evo_phase.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(OUT, f, ensure_ascii=False, indent=1)
    print(f"\n写入 {fp}（{OUT['meta']['wall_s']}s）")
    print("图：results/zero_evo_phase.png / results/zero_evo_travel.png")


if __name__ == "__main__":
    main()
