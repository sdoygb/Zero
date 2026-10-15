# -*- coding: utf-8 -*-
"""
ZeroGPU · 多 seed 分量：K₄ 上多频率能否并存

背景：无反馈 K₄ record 精确继承 seed 比（多支保持）。若 seed 给不同邻对比
（非几何），record 是否同时承载多个频率（多粒子谱）？反馈（Hebbian）下
多频能否并存，还是单点化？

频率定义：record 相邻比的对数（K₄ 有 3 个邻对比 → 3 个候选频率）。

用法：/usr/local/bin/python3 k4_multi_seed.py [cycles] [gens_per]
"""
from __future__ import annotations
import sys, os, math, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from loyalty_engine import LoyaltyEngine
from selforg_k4 import SelfOrgK4

RES = os.path.normpath(os.path.join(HERE, "..", "results"))


def freqs_of(arr):
    """record 的相邻频率集（按顶点顺序的相邻比对数，非降序）。"""
    fs = [abs(math.log(arr[i] / arr[i + 1])) for i in range(len(arr) - 1) if arr[i] > 0 and arr[i + 1] > 0]
    return sorted(round(f, 4) for f in fs)


def run(cycles=60, gens_per=200):
    K4 = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]
    seeds = {
        "geom(2:2:2)":   [8.0, 4.0, 2.0, 1.0],       # 单频率 log2
        "trip(2:8/3:3/2)": [8.0, 4.0, 1.5, 1.0],      # 三频率不同
        "quad(2:8/3:3:2)":[16.0, 8.0, 3.0, 1.0],       # 2, 8/3, 3
        "rand_multi":    [10.0, 4.0, 1.7, 1.0],        # 随机多频
        "slow_trip":     [9.0, 5.0, 2.0, 1.0],         # 9/5, 5/2, 2
    }
    print("=" * 86)
    print("多 seed 分量：K₄ 上多频率能否并存（无反馈 vs Hebbian 反馈）")
    print("=" * 86)
    out = {}
    for name, seed in seeds.items():
        # 无反馈
        eng = LoyaltyEngine(K4, seed=seed, T=4, retention="top2")
        for _ in range(cycles * gens_per):
            eng.run_generation()
        w = eng.record / eng.record.sum()
        nofb = {"branches": [round(float(x), 4) for x in w], "freqs": freqs_of(w),
                "q": float(w @ w)}
        # Hebbian 反馈
        fb = SelfOrgK4(K4, seed=seed, T=4, retention="top2", beta=1.0, lam=0.05)
        for _ in range(cycles * gens_per):
            fb.run_generation()
        w2 = fb.record / fb.record.sum()
        has = {"branches": [round(float(x), 4) for x in w2], "freqs": freqs_of(w2),
               "q": float(w2 @ w2)}
        out[name] = {"seed": seed, "seed_freqs": freqs_of(seed),
                     "nofeedback": nofb, "feedback": has}
        print(f"  [{name:<16}] seed={seed}")
        print(f"      seed 频率       : {freqs_of(seed)}")
        print(f"      无反馈 record   : {nofb['branches']}  频率={nofb['freqs']}  q={nofb['q']:.4f}")
        print(f"      反馈   record   : {has['branches']}  频率={has['freqs']}  q={has['q']:.4f}", flush=True)
    print()
    for name in out:
        s = out[name]
        print(f"  {name:<16}: seed多频={len(s['seed_freqs'])} | 无反馈多频={len(s['nofeedback']['freqs'])} "
              f"| 反馈多频={len(s['feedback']['freqs'])}")
    os.makedirs(RES, exist_ok=True)
    fp = os.path.join(RES, "k4_multi_seed.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"\n已存 {fp}")


if __name__ == "__main__":
    c = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    g = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    run(cycles=c, gens_per=g)
