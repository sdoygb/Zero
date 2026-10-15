# -*- coding: utf-8 -*-
"""
ZeroGPU · K₄ 多频率长期稳定共存验证

背景：K₄ 无反馈 record 精确保持 seed 的 3 个相邻频率（k4_multi_seed）。本脚本
长程验证：在忠实层间动力学（top2 记忆 + T=4 毁灭）下，K₄ 的 3 频率能否长期
稳定并存（不衰减到单频率）？并扫毁灭周期 T 与记忆规则 retention。

频率定义：record 相邻比对数（K₄ 3 邻比→3 候选频率）。

用法：/usr/local/bin/python3 k4_stability.py [gens]
"""
from __future__ import annotations
import sys, os, math, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from loyalty_engine import LoyaltyEngine

RES = os.path.normpath(os.path.join(HERE, "..", "results"))
K4 = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]


def freqs_of(arr):
    fs = [abs(math.log(arr[i] / arr[i + 1])) for i in range(len(arr) - 1) if arr[i] > 0 and arr[i + 1] > 0]
    return sorted(round(f, 4) for f in fs)


def branches(arr, signif=1e-4):
    w = arr / arr.sum()
    mx = w.max()
    sig = [(i, float(w[i])) for i in range(len(w)) if w[i] > signif * mx]
    sig.sort(key=lambda x: -x[1])
    return len(sig), sig


def run(gens=60000):
    seed_trip = [8.0, 4.0, 1.5, 1.0]   # 3 频率：0.405, 0.693, 0.981
    seed_quad = [16.0, 8.0, 3.0, 1.0]  # 3 频率：0.693, 0.981, 1.099
    configs = [
        ("trip_top2_T4",  seed_trip, "top2", 4),
        ("trip_top2_T2",  seed_trip, "top2", 2),
        ("trip_top2_T8",  seed_trip, "top2", 8),
        ("trip_last_T4",  seed_trip, "last", 4),
        ("trip_all_T4",   seed_trip, "all",  4),
        ("quad_top2_T4",  seed_quad, "top2", 4),
    ]
    print("=" * 86)
    print(f"K₄ 多频率长期稳定共存：{gens} 代（每 5000 代采样 record 频率）")
    print("=" * 86)
    out = {}
    for name, seed, ret, T in configs:
        eng = LoyaltyEngine(K4, seed=seed, T=T, retention=ret)
        SAMP = 5000
        samples = []
        for g in range(gens):
            eng.run_generation()
            if (g + 1) % SAMP == 0:
                w = eng.record / eng.record.sum()
                samples.append({"gen": g + 1, "freqs": freqs_of(w),
                                "branches": [round(float(x), 4) for x in w],
                                "q": float(w @ w)})
        w0 = eng.record / eng.record.sum()
        last_freqs = freqs_of(w0)
        # 稳定性：各频率在全程是否保持（相对首采样漂移）
        first = samples[0]["freqs"]
        stable = len(last_freqs) == len(first) and all(
            abs(a - b) < 1e-3 for a, b in zip(first, last_freqs))
        out[name] = {"seed_freqs": freqs_of(seed), "retention": ret, "T": T,
                     "gens": gens, "final_freqs": last_freqs,
                     "stable": stable, "final_q": float(w0 @ w0),
                     "final_branches": [round(float(x), 4) for x in w0],
                     "series": samples}
        print(f"  [{name:<13}] 末代频率={last_freqs}  稳定={stable}  "
              f"q={out[name]['final_q']:.4f}")
        print(f"             record={out[name]['final_branches']}")
        if stable:
            print(f"             ✓ 多频率并存：首={first} 全程保持", flush=True)
        else:
            print(f"             ✗ 频率漂移/衰减：首={first} → 末={last_freqs}", flush=True)
    print()
    os.makedirs(RES, exist_ok=True)
    fp = os.path.join(RES, "k4_stability.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"已存 {fp}")


if __name__ == "__main__":
    g = int(sys.argv[1]) if len(sys.argv) > 1 else 60000
    run(gens=g)
