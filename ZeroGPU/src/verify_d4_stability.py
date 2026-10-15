# -*- coding: utf-8 -*-
"""
ZeroGPU · D4 窗中等长度验证（数据保留 + 分析）

对扫描命中的进窗配置跑中等长度演化（4000 代），验证三层稳态是否稳定保持：
  q_hist / q_live 是否全程落在 D4 窗 (0.5, 0.6)；
  是否漂移出窗、进出窗次数；
  稳态形态（structured / single / None）；
  record 累积、Gini、位宽。

数据保留：每 samp_every 代采样一次，存入 ZeroGPU/results/ 的 JSON。

用法：
    /usr/bin/python3 verify_d4_stability.py        # 跑并分析
"""
from __future__ import annotations
import sys, os, json, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import ZeroCore, path_adj, torus_adj, star_adj
from layer_closure import LayerClosure, QW_MIN, QW_MAX

RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def make_seed(V, kind):
    if kind == "hub2":
        s = [1] * V
        s[0] = 100; s[1] = 100
        return s
    if kind == "hub3":
        s = [1] * V
        s[0] = s[1] = s[2] = 100
        return s
    if kind == "geom":
        return [2 ** (V - 1 - v) for v in range(V)]
    if kind == "all1":
        return [1] * V
    raise ValueError(kind)

# 验证配置（扫描命中者）
CONFIGS = [
    dict(name="path6_hub2_last_L4", adj="path6", seed="hub2", ret="last", L=4),
    dict(name="path6_hub2_all_L6", adj="path6", seed="hub2", ret="all", L=6),
    dict(name="path8_hub2_last_L4", adj="path8", seed="hub2", ret="last", L=4),
    dict(name="path8_hub2_top2_L6", adj="path8", seed="hub2", ret="top2", L=6),
    dict(name="path10_hub2_all_L6", adj="path10", seed="hub2", ret="all", L=6),
    dict(name="path10_hub2_last_L4", adj="path10", seed="hub2", ret="last", L=4),
]

GRAPHS = {
    "path6": path_adj(6), "path8": path_adj(8), "path10": path_adj(10),
}

GENS = 4000
SAMP = 50          # 每 50 代采样一次


def run_one(cfg):
    adj = GRAPHS[cfg["adj"]]
    V = len(adj)
    lc = LayerClosure(adj, seed=make_seed(V, cfg["seed"]),
                      retention=cfg["ret"], L=cfg["L"])
    samp = []
    t0 = time.time()
    for g in range(GENS):
        lc.core.run_generation()
        lc.gens += 1
        if (g + 1) % SAMP == 0:
            qh, ql = lc.q_hist(), lc.q_live()
            samp.append({
                "gen": g + 1,
                "q_hist": qh, "q_live": ql,
                "in_hist": QW_MIN < qh < QW_MAX if qh else False,
                "in_live": QW_MIN < ql < QW_MAX if ql else False,
            })
    s = lc.summary()
    rec = {
        "config": cfg,
        "generations": GENS,
        "sample_every": SAMP,
        "n_samples": len(samp),
        "elapsed_s": round(time.time() - t0, 1),
        "q_hist_final": s["q_hist(L1→L0)"],
        "q_live_final": s["q_live(L2→L0)"],
        "d4_hist_final": s["d4_hist"],
        "d4_live_final": s["d4_live"],
        "hist_form": s["hist_form"],
        "live_form": s["live_form"],
        "record_total": s["record_total"],
        "gini": s["gini"],
        "bits": s["bits"],
        "series": samp,
    }
    # 分析：全程窗内？进出窗次数？
    n_in_hist = sum(1 for p in samp if p["in_hist"])
    n_in_live = sum(1 for p in samp if p["in_live"])
    rec["n_in_hist"] = n_in_hist
    rec["n_in_live"] = n_in_live
    rec["frac_in_hist"] = round(n_in_hist / len(samp), 3)
    rec["frac_in_live"] = round(n_in_live / len(samp), 3)
    # 漂移：q_hist 首窗 vs 末窗
    qh_first = samp[0]["q_hist"]; qh_last = samp[-1]["q_hist"]
    rec["q_hist_first"] = qh_first; rec["q_hist_last"] = qh_last
    rec["q_hist_drift"] = round(qh_last - qh_first, 5)
    return rec


def main():
    os.makedirs(RES, exist_ok=True)
    out = []
    print("=" * 80)
    print(f"中等长度验证：{GENS} 代，每 {SAMP} 代采样，数据保留到 {RES}")
    print("=" * 80)
    for cfg in CONFIGS:
        print(f"\n▶ {cfg['name']}")
        rec = run_one(cfg)
        out.append(rec)
        fname = f"d4_window_{cfg['name']}.json"
        with open(os.path.join(RES, fname), "w", encoding="utf-8") as f:
            json.dump(rec, f, ensure_ascii=False, indent=1)
        print(f"  q_hist 终值={rec['q_hist_final']:.4f}  q_live={rec['q_live_final']:.4f}")
        print(f"  进窗占比 hist={rec['frac_in_hist']:.1%} live={rec['frac_in_live']:.1%}")
        print(f"  q_hist 漂移={rec['q_hist_drift']:+.5f}")
        print(f"  稳态形态 hist={rec['hist_form']} live={rec['live_form']}  Gini={rec['gini']} bits={rec['bits']}")
        print(f"  已存 {RES}/{fname}")

    # 汇总
    print("\n" + "=" * 80)
    print("汇总表（配置 × 进窗占比 × 漂移 × 形态）")
    print("=" * 80)
    print(f"{'配置':<20}{'q_hist终':>8}{'q_live终':>8}{'进窗%hist':>9}{'进窗%live':>9}{'漂移':>9}{'形态':>10}")
    for r in out:
        print(f"{r['config']['name']:<20}{r['q_hist_final']:>8.4f}{r['q_live_final']:>8.4f}"
              f"{r['frac_in_hist']:>9.1%}{r['frac_in_live']:>9.1%}{r['q_hist_drift']:>+9.5f}"
              f"{r['hist_form'] or '-':>10}")
    # 汇总文件
    with open(os.path.join(RES, "d4_window_verify_summary.json"), "w", encoding="utf-8") as f:
        json.dump({"generations": GENS, "configs": out}, f, ensure_ascii=False, indent=1)
    print(f"\n汇总已存 {RES}/d4_window_verify_summary.json")


if __name__ == "__main__":
    main()
