# -*- coding: utf-8 -*-
"""
ZeroGPU · D4 窗更长程确认（path6·last，20000 代）

目标：确认三层自洽不动点在更长时间尺度上成立——
  ① q_hist 是否全程落在 D4 窗、漂移≈0；
  ② ω 分布（record/Σrecord 相对份额）本身是否收敛到稳定形态（不只是 q）；
  ③ 稳态是有结构（多顶点承载）还是坍塌为单点；
  ④ 数据保留（results/ JSON），可复现。

配置：path6, hub2 seed, retention=last, L=4（扫描命中最优，q≈0.5992）。

用法：
    /usr/bin/python3 verify_d4_longrun.py
"""
from __future__ import annotations
import sys, os, json, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import ZeroCore, path_adj
from layer_closure import LayerClosure, QW_MIN, QW_MAX

RES = os.path.normpath(os.path.join(HERE, "..", "results"))

GENS = 20000
SAMP = 200          # 每 200 代采样
EPS = 1e-6

def make_hub2(V):
    s = [1] * V
    s[0] = 100; s[1] = 100
    return s


def main():
    os.makedirs(RES, exist_ok=True)
    adj = path_adj(6)
    V = len(adj)
    lc = LayerClosure(adj, seed=make_hub2(V), retention="last", L=4)

    print("=" * 80)
    print(f"更长程确认：path6·hub2·last·L4 — {GENS} 代，每 {SAMP} 代采样")
    print("=" * 80)

    t0 = time.time()
    samples = []
    for g in range(GENS):
        lc.core.run_generation()
        lc.gens += 1
        if (g + 1) % SAMP == 0:
            w = lc.hist_partition()          # ω = record/Σrecord（6 顶点相对份额）
            qh, ql = lc.q_hist(), lc.q_live()
            samples.append({
                "gen": g + 1,
                "q_hist": qh, "q_live": ql,
                "omega": [round(x, 10) for x in w],
                "in_hist": QW_MIN < qh < QW_MAX if qh else False,
                "record_total": sum(lc.core.record),
            })

    s = lc.summary()
    rec = {
        "config": {"adj": "path6", "seed": "hub2", "retention": "last", "L": 4},
        "generations": GENS, "sample_every": SAMP, "n_samples": len(samples),
        "elapsed_s": round(time.time() - t0, 1),
        "q_hist_final": s["q_hist(L1→L0)"], "q_live_final": s["q_live(L2→L0)"],
        "d4_hist_final": s["d4_hist"], "d4_live": s["d4_live"],
        "record_total": s["record_total"], "gini": s["gini"], "bits": s["bits"],
        "series": samples,
    }

    # ── 分析 ①：q 全程窗内 / 漂移 ──
    qh_series = [p["q_hist"] for p in samples]
    qh_first, qh_last = qh_series[0], qh_series[-1]
    n_in = sum(1 for p in samples if p["in_hist"])
    rec["frac_in_window"] = round(n_in / len(samples), 4)
    rec["q_hist_drift"] = round(qh_last - qh_first, 8)
    rec["q_hist_first"] = qh_first; rec["q_hist_last"] = qh_last
    rec["q_hist_range"] = round(max(qh_series) - min(qh_series), 8)

    # ── 分析 ②：ω 分布收敛（前一半 vs 后一半的 L∞ 距离）──
    half = len(samples) // 2
    w_first_half = samples[half // 2]["omega"]     # 前四分之一末
    w_mid = samples[half]["omega"]                  # 中
    w_last = samples[-1]["omega"]                   # 末
    def linf(a, b):
        return max(abs(x - y) for x, y in zip(a, b))
    rec["omega_Linf_mid_vs_last"] = round(linf(w_mid, w_last), 8)
    rec["omega_Linf_firstQ_vs_last"] = round(linf(w_first_half, w_last), 8)

    # ── 分析 ③：稳态形态（ω 结构）──
    top = sorted(enumerate(w_last), key=lambda x: -x[1])
    rec["omega_final"] = w_last
    rec["omega_top_vertices"] = [(v, round(share, 6)) for v, share in top]
    rec["steady_form"] = s["hist_form"]

    # 保存
    fname = "d4_longrun_path6_last_L4.json"
    with open(os.path.join(RES, fname), "w", encoding="utf-8") as f:
        json.dump(rec, f, ensure_ascii=False, indent=1)

    # ── 打印 ──
    print(f"耗时 {rec['elapsed_s']}s，{len(samples)} 个采样点")
    print(f"① q_hist: 首={qh_first:.6f} 末={qh_last:.6f} 极差={rec['q_hist_range']:.2e} 漂移={rec['q_hist_drift']:+.2e}")
    print(f"   进窗占比 {rec['frac_in_window']:.1%}  D4窗 终判={rec['d4_hist_final']}")
    print(f"② ω 收敛: L∞(中,末)={rec['omega_Linf_mid_vs_last']:.2e}  L∞(前¼,末)={rec['omega_Linf_firstQ_vs_last']:.2e}")
    print(f"③ ω 末态分布: {rec['omega_final']}")
    print(f"   top 顶点: {rec['omega_top_vertices']}")
    print(f"   Gini={rec['gini']:.4f}  bits={rec['bits']}")
    print(f"④ 已存 {RES}/{fname}")


if __name__ == "__main__":
    main()
