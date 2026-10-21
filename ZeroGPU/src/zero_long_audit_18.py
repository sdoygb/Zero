#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_long_audit_18.py —— 18 世界 micro core 的长时间演化审计

背景：2026-10-18 的 `zero_long_audit.py` 用的是 16/8/4 micro core，结论是
「形状冻结」（~2e4 步后 n_std/q/h_E 尾段 std=0，只有 psi 相位还在转）。
此后用户改判当前亚层模型为 **18 世界**（4→18→6，3 代 × 6 边；见
`zero_l2_partitions.md` §8.2、`zero_evo_18.md`），但 18 core 的长跑只做到
1e4 代、且当时没有冻结判据。本脚本把同一套判据搬到 18 core 上。

监视量分三组：
  · 宏观形状：n_std, q, h_E, psi_coh, q_live, q_hist
  · 层份额（18 世界的三个代）：share_g = W18.sum(axis=1)（新，16/8/4 时没有）
  · 影子 6：W18 对代求和后的边分布（新）

判据（每 block 采样）：
  · 非有限值 → blowup
  · record/seed 为 0（micro 开）→ underflow
  · 全部监视量连续 5 个 block 相对变化 < 1e-12 → frozen
  · 否则 ok

用法：
    /usr/bin/python3 zero_long_audit_18.py --case full  --steps 1000000
    /usr/bin/python3 zero_long_audit_18.py --case micro --steps 1000000
    /usr/bin/python3 zero_long_audit_18.py --case both  --steps 1000000
"""
from __future__ import annotations
import argparse, json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_evo import ZeroEvo, Params  # noqa: E402
from zero_evo_18 import ZeroEvo18     # noqa: E402
import zero_sim as zs                 # noqa: E402

RES = os.path.join(HERE, "..", "results")


def make_case(case: str):
    adj = zs.cycle(16)
    V = len(adj)
    if case == "full":
        # 与 zero_long_audit.py 的 full 同参（12 通道全开，c_h=1.7）
        p = Params(micro=True, micro_semantics="persistent", T=4, L=4,
                   c_h=1.7, g_src=1.0, kappa_h=0.5, Gamma_h=0.2,
                   kappa_p=0.1, gamma=0.2, mu=0.5, K=2.0,
                   wave_on=True, poisson_on=True)
        e = ZeroEvo18(adj, p)
        return e, "cycle16/18core/full/12ch/c_h=1.7"
    if case == "coupled":
        # 同上，但把 micro->宏观 的边打开：alpha>0（边③ 闭合→内容）+ eta_*
        # 这是唯一让 18 core 与 16/8/4 core 产生区别的配置（见审计正文 §3）
        p = Params(micro=True, micro_semantics="persistent", T=4, L=4,
                   c_h=1.7, g_src=1.0, kappa_h=0.5, Gamma_h=0.2,
                   kappa_p=0.1, gamma=0.2, mu=0.5, K=2.0,
                   wave_on=True, poisson_on=True,
                   alpha=0.3, eta_n=0.2, eta_h=0.2, eta_psi=0.2)
        e = ZeroEvo18(adj, p)
        return e, "cycle16/18core/coupled/12ch+alpha=0.3+eta=0.2"
    if case == "micro":
        p = Params(micro=True, micro_semantics="persistent", T=4, L=4,
                   c_h=0.0, g_src=0.0, kappa_h=0.0,
                   wave_on=False, poisson_on=False)
        e = ZeroEvo18(adj, p, seed=np.array([2.0 ** (V - 1 - i) for i in range(V)]))
        return e, "cycle16/18core/micro-only/persistent"
    raise SystemExit(f"unknown case {case}")


def layer_stats(e):
    """18 世界的层份额与影子 6。"""
    tot = float(e.W18.sum())
    if tot > 0:
        share = (e.W18.sum(axis=1) / tot)
    else:
        share = np.zeros(3)
    sh = e.W18.sum(axis=0)
    sh = sh / sh.sum() if sh.sum() > 0 else sh
    return share, sh


def audit(case: str, steps: int, block: int, seed_tag: str = "cycle16"):
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
        sd = float(e.seed.sum())
        share, sh = layer_stats(e)
        met = {
            "gen": g, "wall_s": round(time.perf_counter() - t0, 3),
            "n_std": s["n_std"], "q": s["q"], "h_energy": s["h_energy"],
            "psi_norm": s["psi_norm"], "psi_coh": s["psi_coh"],
            "record_sum": rec, "seed_sum": sd,
            "q_live": s["q_live"], "q_hist": s["q_hist"],
            "log_scale": float(e.log_scale),
            "share0": float(share[0]), "share1": float(share[1]), "share2": float(share[2]),
            "shadow": [float(v) for v in sh],
        }
        keys = ("n_std", "q", "h_energy", "psi_norm", "psi_coh", "q_live", "q_hist",
                "share0", "share1", "share2")
        vals = np.array([met[k] for k in keys], dtype=float)
        if not np.all(np.isfinite(vals)):
            status = "blowup"
        if e.p.micro and (rec <= 0.0 or sd <= 0.0):
            status = "underflow"
        if prev is not None and status == "ok":
            rel = [abs(met[k] - prev[k]) / max(1.0, abs(prev[k])) for k in keys]
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
                  f"share=({met['share0']:.4f},{met['share1']:.4f},{met['share2']:.4f}) "
                  f"rec={rec:.3g} seed={sd:.3g} status={status}", flush=True)
        if status != "ok":
            break

    def stat(a):
        a = np.asarray(a, dtype=float)
        tail = max(1, len(a) // 5)
        b = a[-tail:]
        return {"mean": float(b.mean()), "std": float(b.std()),
                "min": float(b.min()), "max": float(b.max())}

    ns = [b["n_std"] for b in blocks]
    qs = [b["q"] for b in blocks]
    hs = [b["h_energy"] for b in blocks]
    tail = max(1, len(blocks) // 5)
    period_blocks = 0.0
    y = np.asarray(ns[-tail:], float)
    if len(y) > 10:
        y = y - y.mean()
        if y.std() > 0:
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
        "tail_share0": stat([b["share0"] for b in blocks]),
        "tail_share1": stat([b["share1"] for b in blocks]),
        "tail_share2": stat([b["share2"] for b in blocks]),
        "period_blocks": period_blocks,
        "period_steps": period_blocks * block,
        "block_time_median_s": float(np.median(times)) if times else 0.0,
        "block_time_max_s": float(np.max(times)) if times else 0.0,
    }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="full",
                    choices=["full", "micro", "coupled", "both", "all"])
    ap.add_argument("--steps", type=int, default=1_000_000)
    ap.add_argument("--block", type=int, default=10_000)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    cases = {"full": ["full"], "micro": ["micro"], "coupled": ["coupled"],
             "both": ["full", "micro"], "all": ["full", "micro", "coupled"]}[args.case]
    out = {"steps": args.steps, "block": args.block, "core": "18",
           "cases": {}}
    for c in cases:
        print(f"=== {c} {args.steps} steps, block={args.block} ===", flush=True)
        out["cases"][c] = audit(c, args.steps, args.block)
        r = out["cases"][c]
        print(f"=== {c}: status={r['status']} wall={r['wall_s']}s "
              f"speed={r['steps_per_s']} steps/s period={r['period_steps']:.1f} ===",
              flush=True)
    suffix = f"_{args.tag}" if args.tag else ""
    if len(cases) > 1:
        fp = os.path.join(RES, f"zero_long_audit_18{suffix}.json")
    else:
        fp = os.path.join(RES, f"zero_long_audit_18_{cases[0]}{suffix}.json")
    with open(fp, "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"→ {fp}")


if __name__ == "__main__":
    main()