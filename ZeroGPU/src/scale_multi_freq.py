# -*- coding: utf-8 -*-
"""
ZeroGPU · 更大图完全收敛：能否出现多分支（多频率）

path8/12/16 收敛到二分支 2:1（单频率 log2，Z087）。path24/32 有泄漏分支
（中间两主支两侧的小值）。本脚本跑超长程，判断泄漏分支是完全收敛到 0
（q→5/9，纯二分支）还是保持有限（多支稳态 → 多内禀频率）。

用法：/usr/local/bin/python3 scale_multi_freq.py [path24] [path32] [path64]
"""
from __future__ import annotations
import sys, os, time, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import path_adj
from z0_vectorized import ZeroCoreVec

RES = os.path.normpath(os.path.join(HERE, "..", "results"))
SIGNIF = 1e-4   # 显著支阈值（相对 max）


def significant(w, rel=SIGNIF):
    """显著支：权重 > rel × max 的顶点。"""
    mx = w.max()
    return [int(i) for i in range(len(w)) if w[i] > rel * mx]


def run(n, lock_gens, run_gens, sample_every):
    adj = path_adj(n); V = len(adj)
    seed = [float(2 ** (V - 1 - v)) for v in range(V)]
    v = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
    for _ in range(lock_gens):
        v.run_generation()
    trace = []
    t0 = time.time()
    for g in range(run_gens):
        v.run_generation()
        if (g + 1) % sample_every == 0:
            w = v.record / v.record.sum()
            sig = significant(w)
            q = float((w * w).sum())
            wt = [(i, round(float(w[i]), 6)) for i in sig]
            trace.append({"gen": lock_gens + g + 1, "n_branch": len(sig),
                          "branches": wt, "q": q})
    return {"graph": f"path{n}", "V": V, "lock_gens": lock_gens,
            "run_gens": run_gens, "sample_every": sample_every,
            "wall_s": round(time.time() - t0, 1), "trace": trace}


def main():
    names = sys.argv[1:] or ["path24", "path32", "path64"]
    LOCK = {"path24": 80000, "path32": 80000, "path64": 80000}
    RUN = {"path24": 5_000_000, "path32": 5_000_000, "path64": 3_000_000}
    EVERY = 200_000
    out = {}
    for name in names:
        n = int(name.replace("path", ""))
        print(f"==== {name}: 锁 {LOCK[name]} 代 + 跑 {RUN[name]} 代（每 {EVERY} 采样）====", flush=True)
        r = run(n, LOCK[name], RUN[name], EVERY)
        out[name] = r
        for t in r["trace"][::2]:  # 隔行打印
            print(f"  gen {t['gen']:>7}: 支={t['n_branch']}  q={t['q']:.6f}  支权重={t['branches']}")
        print(f"  耗时 {r['wall_s']}s\n", flush=True)
    os.makedirs(RES, exist_ok=True)
    path = os.path.join(RES, "scale_multi_freq.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"已存 {path}")


if __name__ == "__main__":
    main()
