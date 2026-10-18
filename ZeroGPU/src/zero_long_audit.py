#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Zero 忠实引擎的长时间演化审计：卡住/冻结/爆掉/下溢。

用法：
    /usr/bin/python3 zero_long_audit.py --case full --steps 10000000
    /usr/bin/python3 zero_long_audit.py --case layer --steps 10000000
    /usr/bin/python3 zero_long_audit.py --case both --steps 2000000

判据（每 block 采样）：
  · 非有限值 → blowup
  · record/seed 为 0（micro 开）→ underflow
  · 所有监控量在连续 5 个 block 的相对变化 < 1e-12 → frozen
  · block 墙钟 > 5× 中位数 → slow
  · 否则 ok（若 n_std 有周期，另报周期）
"""
from __future__ import annotations
import argparse, json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_evo import ZeroEvo  # noqa: E402
from zero_evo_phase import base_params, init_random  # noqa: E402
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")


def make_case(case: str):
    adj = zs.cycle(16)
    V = len(adj)
    x0 = init_random(V, 1)
    if case == "full":
        p = base_params(micro=True, wave_on=True, poisson_on=True,
                        c_h=1.7, g_src=1.0, kappa_h=0.5, Gamma_h=0.2,
                        kappa_p=0.1, gamma=0.2, mu=0.5, K=2.0)
        e = ZeroEvo(adj, p, n0=1.0 + x0)
        return e, "cycle16/full/12ch/c_h=1.7"
    if case == "layer":
        p = base_params(micro=True, wave_on=False, poisson_on=False,
                        c_h=0.0, g_src=0.0, kappa_h=0.0)
        p.micro_semantics = "persistent"
        e = ZeroEvo(adj, p, seed=np.array([2.0 ** (V - 1 - i) for i in range(V)]))
        return e, "cycle16/micro-only/persistent"
    raise SystemExit(f"unknown case {case}")


def audit(case: str, steps: int, block: int):
    e, label = make_case(case)
    t0 = time.perf_counter()
    prev = None
    frozen_run = 0
    blocks = []
    times = []
    status = "ok"
    last_gen = 0
    for g in range(1, steps + 1):
        tb = time.perf_counter()
        e.step()
        if g % block != 0:
            continue
        s = e.snapshot()
        rec = float(e.record.sum())
        seed = float(e.seed.sum())
        met = {
            "gen": g, "wall_s": round(time.perf_counter() - t0, 3),
            "n_std": s["n_std"], "q": s["q"], "h_energy": s["h_energy"],
            "psi_norm": s["psi_norm"], "psi_coh": s["psi_coh"],
            "record_sum": rec, "seed_sum": seed,
            "q_live": s["q_live"], "q_hist": s["q_hist"],
            "log_scale": float(e.log_scale),
        }
        vals = np.array([met[k] for k in ("n_std", "q", "h_energy", "psi_norm",
                                          "psi_coh", "record_sum", "q_live", "q_hist")],
                        dtype=float)
        if not np.all(np.isfinite(vals)):
            status = "blowup"
        if e.p.micro and (rec <= 0.0 or seed <= 0.0):
            status = "underflow"
        if prev is not None and status == "ok":
            rel = [abs(met[k] - prev[k]) / max(1.0, abs(prev[k])) for k in
                   ("n_std", "q", "h_energy", "psi_norm", "psi_coh", "q_live", "q_hist")]
            if max(rel) < 1e-12:
                frozen_run += 1
            else:
                frozen_run = 0
            if frozen_run >= 5:
                status = "frozen"
        blocks.append(met)
        last_gen = g
        times.append(time.perf_counter() - tb)
        prev = met
        if len(blocks) <= 5 or len(blocks) % 10 == 0:
            print(f"[{case}] gen={g:>10} wall={met['wall_s']:>8.1f}s "
                  f"n_std={met['n_std']:.6g} q={met['q']:.6g} "
                  f"h_E={met['h_energy']:.6g} psi={met['psi_norm']:.6g} "
                  f"rec={rec:.3g} seed={seed:.3g} status={status}", flush=True)
        if status != "ok":
            break
    # summary
    ns = np.array([b["n_std"] for b in blocks], dtype=float)
    qs = np.array([b["q"] for b in blocks], dtype=float)
    hs = np.array([b["h_energy"] for b in blocks], dtype=float)
    tail = max(1, len(blocks) // 5)
    def stat(a):
        b = a[-tail:]
        return {"mean": float(b.mean()), "std": float(b.std()),
                "min": float(b.min()), "max": float(b.max())}
    # crude period from zero crossings of n_std tail
    period_blocks = 0.0
    y = ns[-tail:] - ns[-tail:].mean()
    if len(y) > 10 and y.std() > 0:
        zc = np.where(np.sign(y[:-1]) * np.sign(y[1:]) < 0)[0]
        if len(zc) >= 2:
            period_blocks = float(2 * np.median(np.diff(zc)))
    out = {
        "case": case, "label": label, "requested_steps": steps,
        "steps_run": last_gen, "block": block,
        "blocks": blocks, "status": status,
        "wall_s": round(time.perf_counter() - t0, 3),
        "steps_per_s": round(last_gen / max(1e-9, time.perf_counter() - t0), 1),
        "tail_n_std": stat(ns), "tail_q": stat(qs), "tail_h_energy": stat(hs),
        "period_blocks": period_blocks,
        "period_steps": period_blocks * block,
        "block_time_median_s": float(np.median(times)) if times else 0.0,
        "block_time_max_s": float(np.max(times)) if times else 0.0,
    }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="full", choices=["full", "layer", "both"])
    ap.add_argument("--steps", type=int, default=10_000_000)
    ap.add_argument("--block", type=int, default=10_000)
    args = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    cases = ["full", "layer"] if args.case == "both" else [args.case]
    out = {"steps": args.steps, "block": args.block, "cases": {}}
    for c in cases:
        print(f"=== {c} {args.steps} steps, block={args.block} ===", flush=True)
        out["cases"][c] = audit(c, args.steps, args.block)
        r = out["cases"][c]
        print(f"=== {c}: status={r['status']} wall={r['wall_s']}s "
              f"speed={r['steps_per_s']} steps/s period={r['period_steps']:.1f} ===", flush=True)
        if r["status"] != "ok":
            break
    fp = os.path.join(RES, "zero_long_audit.json" if args.case == "both"
                      else f"zero_long_audit_{args.case}.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print("写入", fp, flush=True)


if __name__ == "__main__":
    main()
