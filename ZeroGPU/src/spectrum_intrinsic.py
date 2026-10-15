# -*- coding: utf-8 -*-
"""
ZeroGPU · 内禀频率谱（粒子候选）

Z027 命题 3 的应用：稳态 ω_v 若为"块内两支不等权"→ 自发非中心态
→ 模流非平凡 → 内禀时间（频率 log(p/q)）。粒子质量候选 = 内禀频率。

块 = 稳态非零顶点集；两支 = 前两大非零顶点的占比。
频率 = |log(p/q)|，来自各"块"两支的不等权比。

数据保留到 results/。
用法：/usr/local/bin/python3 spectrum_intrinsic.py [path8] [path12] ...
"""
from __future__ import annotations
import sys, os, time, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import path_adj
from z0_vectorized import ZeroCoreVec

RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOW = 1e-9   # 非零阈值


def intrinsic(n, n_steady=60000):
    adj = path_adj(n); V = len(adj)
    seed = [float(2 ** (V - 1 - v)) for v in range(V)]
    vec = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
    for _ in range(n_steady):
        vec.run_generation()
    w = vec.record / vec.record.sum()
    q = float((w * w).sum())
    # 非零顶点集
    nz = [(i, float(w[i])) for i in range(V) if w[i] > LOW]
    nz.sort(key=lambda x: -x[1])
    # 两支：前两大（含零填充到2）
    top = nz[:2]
    p = top[0][1] if len(top) > 0 else 0.0
    qq = top[1][1] if len(top) > 1 else 0.0
    freq = abs(np.log(p / qq)) if (p > 0 and qq > 0) else None
    return {"graph": f"path{n}", "V": V, "n_steady": n_steady, "q": q,
            "q_59": 5 / 9, "branches": [(i, round(x, 6)) for i, x in nz],
            "n_branches": len(nz),
            "top2": [round(top[0][1], 9)] if top else [],
            "ratio_pq": (round(p / qq, 9) if qq > 0 else None),
            "intrinsic_freq_log_pq": (round(freq, 9) if freq is not None else None)}


def main():
    names = sys.argv[1:] or ["path8", "path12", "path16", "path24", "path32"]
    out = {}
    for name in names:
        n = int(name.replace("path", ""))
        t0 = time.time()
        r = intrinsic(n)
        r["wall_s"] = round(time.time() - t0, 1)
        out[name] = r
        print(f"path{n:>2}: q={r['q']:.6f}  非零分支={r['n_branches']}  "
              f"top2={r['top2']}  比={r['ratio_pq']}  "
              f"频率 log(p/q)={r['intrinsic_freq_log_pq']}")
    os.makedirs(RES, exist_ok=True)
    path = os.path.join(RES, "intrinsic_spectrum.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"\n已存 {path}")


if __name__ == "__main__":
    main()
