# -*- coding: utf-8 -*-
"""
ZeroGPU · 多峰 seed：能否突破"恒为 2 支"产生多频率

此前所有 seed 稳态都收敛到中间两顶点（2 支），比 = seed 中间两顶点比。
本脚本测多峰/多分离峰 seed，看能否让稳态保留 >2 支（→ 多频率）。

用法：/usr/local/bin/python3 seed_multipeak.py [path16] [path24]
"""
from __future__ import annotations
import sys, os, json, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import path_adj
from z0_vectorized import ZeroCoreVec

RES = os.path.normpath(os.path.join(HERE, "..", "results"))
SIGNIF = 1e-4


def multipeak_seeds(V):
    c = V // 2 - 1          # path 中心偏左顶点（path16→7, path24→11）
    def mk(base):
        s = [1.0] * V
        for i, val in base:
            s[i] = val
        return s
    return {
        # 中间 3 峰（相邻）
        "3peak_adj":  mk([(c, 100.0), (c + 1, 50.0), (c + 2, 25.0)]),
        # 中间 4 峰（相邻）
        "4peak_adj":  mk([(c, 100.0), (c + 1, 50.0), (c + 2, 25.0), (c + 3, 12.5)]),
        # 中间 3 峰居中
        "3peak_ctr":  mk([(c - 1, 100.0), (c, 50.0), (c + 1, 25.0)]),
        # 平台（中间 4 顶点等高）
        "plateau4":   mk([(c, 10.0), (c + 1, 10.0), (c + 2, 10.0), (c + 3, 10.0)]),
        # 3 个分离峰（不相邻）
        "3peak_sep":  mk([(c - 3, 100.0), (c, 50.0), (c + 3, 25.0)]),
        # 2 个分离峰
        "2peak_sep":  mk([(c - 2, 100.0), (c + 2, 100.0)]),
        # 中间 5 峰（跨度更大）
        "5peak_adj":  mk([(c - 1, 100.0), (c, 50.0), (c + 1, 25.0),
                          (c + 2, 12.5), (c + 3, 6.25)]),
        # 全带偏（中间半段指数，其余低）
        "mid_band":   mk([(i, 2 ** (c + 3 - i)) for i in range(max(0, c - 2), min(V, c + 5))]),
    }


def scan(n, lock_gens=80000):
    adj = path_adj(n); V = len(adj)
    out = {}
    for name, seed in multipeak_seeds(V).items():
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
        out[name] = {"n_branch": len(sig), "branches": [(i, round(x, 4)) for i, x in sig],
                     "q": q, "in_win": 0.5 < q < 0.6,
                     "top2_ratio": (round(ratio, 3) if ratio else None),
                     "freq": (round(freq, 4) if freq else None)}
        print(f"{name:<12}: 支={len(sig):>2} q={q:.5f} {'入窗' if 0.5<q<0.6 else '出窗'}  "
              f"支={[(i,round(x,3)) for i,x in sig][:6]}  频率={freq and round(freq,4)}", flush=True)
    return out


def main():
    names = sys.argv[1:] or ["path16", "path24"]
    full = {}
    for name in names:
        n = int(name.replace("path", ""))
        print(f"\n==== {name}（锁 80000 代）====", flush=True)
        full[name] = scan(n)
    os.makedirs(RES, exist_ok=True)
    path = os.path.join(RES, "seed_multipeak.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(full, f, ensure_ascii=False, indent=1)
    print(f"\n已存 {path}")


if __name__ == "__main__":
    main()
