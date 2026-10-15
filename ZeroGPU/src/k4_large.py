# -*- coding: utf-8 -*-
"""
ZeroGPU · 更大 K₄ 构造：能否承载多频率（多粒子谱）

背景：K₄ 正则不消灭多支、结构可从 seed 稳定保持（k4_growth 已证）。本脚本
构造更大的 K₄ 衍生图，检验能否承载多个独立频率：
  * K₄ 双覆盖（8顶点、度3、二部、girth4）
  * Cartesian 积 K₄□K₄（16顶点、度6、正则）
  * 张量积 K₄×K₄（16顶点、度9、正则）

seed 设计：分块 seed——不同顶点域给不同 seed 比（不同频率）。
问：record 是同时保持多个频率（多粒子谱），还是收敛到单一频率？

用法：/usr/local/bin/python3 k4_large.py [cycles] [gens_per]
"""
from __future__ import annotations
import sys, os, math, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from loyalty_engine import LoyaltyEngine
from z0_vectorized import _adj_matrix

RES = os.path.normpath(os.path.join(HERE, "..", "results"))


def K4():
    return [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]


def double_cover(k4):
    """K₄ 双覆盖：8 顶点、度3、二部、girth4。"""
    from itertools import combinations
    n = 8
    adj = [[] for _ in range(n)]
    def vid(a, s): return 2 * a + s
    for a, b in combinations(range(4), 2):
        adj[vid(a, 0)].append(vid(b, 1)); adj[vid(b, 1)].append(vid(a, 0))
        adj[vid(a, 1)].append(vid(b, 0)); adj[vid(b, 0)].append(vid(a, 1))
    return [sorted(set(a)) for a in adj]


def cartesian(k4):
    """K₄□K₄：16 顶点、度6。顶点(u,v)。"""
    V = len(k4)
    adj = [[] for _ in range(V * V)]
    def idx(u, v): return u * V + v
    for u in range(V):
        for v in range(V):
            i = idx(u, v)
            for un in k4[u]:
                adj[i].append(idx(un, v))
            for vn in k4[v]:
                adj[i].append(idx(u, vn))
    return [sorted(set(a)) for a in adj]


def tensor(k4):
    """K₄×K₄（张量积）：16 顶点、度9。顶点(u,v)。"""
    V = len(k4)
    adj = [[] for _ in range(V * V)]
    def idx(u, v): return u * V + v
    for u in range(V):
        for v in range(V):
            i = idx(u, v)
            for un in k4[u]:
                for vn in k4[v]:
                    adj[i].append(idx(un, vn))
    return [sorted(set(a)) for a in adj]


def gini(arr):
    x = sorted(float(v) for v in arr if v)
    n = len(x)
    if n == 0:
        return 0.0
    s = sum(x)
    num = sum((2 * (i + 1) - n - 1) * x[i] for i in range(n))
    return num / (n * s)


def run(cycles=60, gens_per=200):
    k4 = K4()
    dc = double_cover(k4)     # 8 顶点
    cart = cartesian(k4)      # 16 顶点
    ten = tensor(k4)          # 16 顶点

    # 分块 seed：把顶点按"域"分组，不同域给不同 seed 比 → 不同频率
    def geom_in_block(indices, vals):
        s = [1.0] * (max(indices) + 1) if indices else [1.0]
        # 用图的真实规模
        return None  # 占位，下面按图规模造

    def block_seed(V, blocks):
        """blocks: list of (indices, ratios)；ratio 给块内顶点。"""
        s = np.ones(V, dtype=np.float64)
        for idxs, ratios in blocks:
            for k, i in enumerate(idxs):
                if k < len(ratios):
                    s[i] = ratios[k]
        return list(s)

    seeds = {
        # 单频率 geom（整体几何衰减）
        "K4_geom":        ([8.0, 4.0, 2.0, 1.0],),
        # 双覆盖：两半各给不同频率 (8:4:2:1 和 4:3:2:1)
        "dc_2freq":       None,  # 下方构建
    }

    configs = [
        ("K4",        k4,   [8.0, 4.0, 2.0, 1.0]),
        # 双覆盖：上半 geom，下半反 geom → 两块不同频率
        ("dc_2freq",  dc,   block_seed(len(dc), [([0, 2, 4, 6], [8, 4, 2, 1]),
                                                 ([1, 3, 5, 7], [1, 2, 4, 8])])),
        # Cartesian：4 个 K₄ 层，各层不同比
        ("cart_4freq", cart, block_seed(len(cart), [
            ([0 + 4 * v for v in range(4)], [8, 4, 2, 1]) if False else
            ([i for i in range(len(cart)) if i % 4 == v], [8, 4, 2, 1][:4]) for v in range(4)])),
        # 张量积：geom 分块
        ("tensor_geom", ten, None),
    ]
    # 修正 seed
    configs[3] = ("tensor_geom", ten, block_seed(len(ten), [([i for i in range(len(ten)) if i % 4 == v],
                                                             [16, 8, 4, 2][v % 4] and [4, 2, 1, 0.5][:4]) for v in range(4)]))
    configs[3] = ("tensor_geom", ten, block_seed(len(ten), [
        ([i for i in range(len(ten)) if i // 4 == u], [8, 4, 2, 1]) for u in range(4)]))

    print("=" * 86)
    print(f"更大 K₄ 构造：多频率承载测试，{cycles}周期×{gens_per}代")
    print("=" * 86)
    out = {}
    for name, adj, seed in configs:
        V = len(adj)
        degs = [len(a) for a in adj]
        eng = LoyaltyEngine(adj, seed=seed, T=4, retention="top2")
        last = {}
        for cyc in range(cycles):
            for _ in range(gens_per):
                eng.run_generation()
            nb, br, q, ratio, freq = eng.branches()
            g = gini(eng.record)
            last = {"n_branch": nb, "branches": br, "q": q,
                    "top2_ratio": (round(ratio, 3) if ratio else None),
                    "freq": (round(freq, 4) if freq else None),
                    "gini": g}
            if (cyc + 1) % 15 == 0:
                print(f"  [{name:<11} V={V} 度{set(degs)}] gen{eng.gen:>5}: 支={nb:>2} "
                      f"q={q:.4f} Gini={g:.4f} 主支={br[:6]} 频率={freq and round(freq,4)}", flush=True)
        out[name] = {"V": V, "deg": sorted(set(degs)), "final": last}
    print()
    for name in out:
        f = out[name]["final"]
        print(f"  {name:<11}: 支={f['n_branch']:>2} q={f['q']:.4f} Gini={f['gini']:.4f} "
              f"主支={f['branches'][:8]} 频率={f['freq']}")
    os.makedirs(RES, exist_ok=True)
    fp = os.path.join(RES, "k4_large.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"\n已存 {fp}")


if __name__ == "__main__":
    c = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    g = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    run(cycles=c, gens_per=g)
