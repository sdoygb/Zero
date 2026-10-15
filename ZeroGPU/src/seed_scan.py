# -*- coding: utf-8 -*-
"""
ZeroGPU · 换 seed：检验 log2 是普适还是 geom seed 特例

扫多种 seed 类型（path16/24），看稳态的支结构（支数、比、q、内禀频率）是否变化。
若都收敛 2:1 → log2 普适（seed 只定位置不定比）；若不同 seed 给不同比/支 → seed 决定谱。

用法：/usr/local/bin/python3 seed_scan.py [path16] [path24]
"""
from __future__ import annotations
import sys, os, time, json, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import path_adj
from z0_vectorized import ZeroCoreVec

RES = os.path.normpath(os.path.join(HERE, "..", "results"))
SIGNIF = 1e-4


def seed_types(V):
    rng = np.random.default_rng(12345)
    return {
        "geom":     [float(2 ** (V - 1 - v)) for v in range(V)],
        "geom2":    [float(1.5 ** (V - 1 - v)) for v in range(V)],   # 不同底
        "uniform":  [1.0] * V,
        "linear_d": [float(V - v) for v in range(V)],
        "linear_i": [float(v + 1) for v in range(V)],
        "quad_d":   [float((V - v) ** 2) for v in range(V)],
        "rand":     list(rng.uniform(0.1, 1.0, V)),
        "end_peaks":[0.0] * V,  # 两端高中间低
    }


def scan(n, lock_gens=80000):
    adj = path_adj(n); V = len(adj)
    seeds = seed_types(V)
    seeds["end_peaks"][0] = 10.0; seeds["end_peaks"][V - 1] = 10.0
    seeds["end_peaks"] = [max(x, 0.001) for x in seeds["end_peaks"]]
    out = {}
    for name, seed in seeds.items():
        v = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
        for _ in range(lock_gens):
            v.run_generation()
        w = v.record / v.record.sum()
        q = float((w * w).sum())
        mx = w.max()
        sig = [(i, float(w[i])) for i in range(V) if w[i] > SIGNIF * mx]
        sig.sort(key=lambda x: -x[1])
        if len(sig) >= 2:
            ratio = sig[0][1] / sig[1][1]
            freq = abs(math.log(ratio))
        else:
            ratio = freq = None
        out[name] = {"seed": "…", "n_branch": len(sig),
                     "branches": [(i, round(x, 5)) for i, x in sig],
                     "q": q, "in_win": 0.5 < q < 0.6,
                     "top2_ratio": (round(ratio, 4) if ratio else None),
                     "freq": (round(freq, 4) if freq else None)}
        print(f"path{n} seed={name:<10}: 支={len(sig):>2} q={q:.5f} "
              f"{'入窗' if 0.5<q<0.6 else '出窗'}  支={[(i,round(x,3)) for i,x in sig][:5]}  "
              f"比={ratio and round(ratio,3)}  频率={freq and round(freq,4)}", flush=True)
    return out


def main():
    names = sys.argv[1:] or ["path16", "path24"]
    full = {}
    for name in names:
        n = int(name.replace("path", ""))
        print(f"\n==== {name}（锁 80000 代）====", flush=True)
        full[name] = scan(n)
    os.makedirs(RES, exist_ok=True)
    path = os.path.join(RES, "seed_scan.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(full, f, ensure_ascii=False, indent=1)
    print(f"\n已存 {path}")


if __name__ == "__main__":
    main()
