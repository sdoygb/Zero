# -*- coding: utf-8 -*-
"""
ZeroGPU · 从 K₄ 长出结构：正则图上能否自发产生结构

背景：gamma_origin 导出 Γ 局部 = K₄（4顶点、度3、正则）；gamma_gap 判定
"结构只能从度不均匀来，度均匀⟹Gini=0"。本脚本在正则 K₄ 上跑忠实引擎
（top2 记忆 + T=4 毁灭 + 可选容量），检验：
  (a) K₄ 无权重：演化回馈能否自发长出 Gini>0 / 多支？
  (b) K₄ 带权重（耦合不均匀）：有效度不均匀能否给结构？
  (c) 对比路径图（非正则，已知有结构）。

用法：/usr/local/bin/python3 k4_growth.py [T] [cycles] [gens_per]
"""
from __future__ import annotations
import sys, os, math, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from loyalty_engine import LoyaltyEngine
from z0_engine import path_adj
from z0_vectorized import NORM_THRESHOLD, _adj_matrix

RES = os.path.normpath(os.path.join(HERE, "..", "results"))


def gini(arr):
    x = sorted(float(v) for v in arr if v)
    n = len(x)
    if n == 0:
        return 0.0
    s = sum(x)
    num = sum((2 * (i + 1) - n - 1) * x[i] for i in range(n))
    return num / (n * s)


class K4Engine(LoyaltyEngine):
    """在 K₄（带可选边权重）上跑忠实引擎。权重通过 A 矩阵注入 → 有效度可不均匀。"""

    def __init__(self, adj, weights=None, L=4, seed=None, T=4, retention="top2", cap=None):
        super().__init__(adj, L=L, seed=seed, T=T, retention=retention, cap=cap)
        if weights is None:
            self.A = _adj_matrix(adj)
        else:
            A = np.zeros((len(adj), len(adj)), dtype=np.float64)
            for u, nb in enumerate(adj):
                for w_ in nb:
                    A[u][w_] = weights[u].get(w_, 1.0)
            self.A = A

    def run_generation(self):
        V, L = self.V, self.L
        W = self.W.copy()
        W[np.arange(V), np.arange(V)] += self.seed
        for _ in range(L):
            W = W @ self.A
            if self.cap is not None:
                np.minimum(W, float(self.cap), out=W)
        diag = W[np.arange(V), np.arange(V)]
        self.gen_closures = diag.copy()
        self.record = self.record + diag
        W2 = W.copy()
        W2[np.arange(V), np.arange(V)] = 0
        self.W = W2
        self.seed = self._retained()
        if float(self.seed.sum()) == 0:
            self.seed = np.ones(self.V, dtype=np.float64)
        self.gen += 1
        if self.gen % self.T == 0:
            self.destructions.append((self.gen, self.branches()))
            self.W = np.zeros((self.V, self.V), dtype=np.float64)
        self._normalize()
        return self


def run(T=4, cycles=60, gens_per=200):
    K4 = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]  # 完全图 K₄
    P16 = path_adj(16)

    configs = {
        # (名称, adj, weights, seed)
        "K4_uniform":  (K4, None, [1.0] * 4),                      # 正则无权，seed 均匀
        "K4_geom":     (K4, None, [8.0, 4.0, 2.0, 1.0]),           # 正则无权，seed 不均匀
        "K4_weighted": (K4,
                        # 边权重不均匀 → 有效度不均匀（顶点0强耦合、顶点3弱）
                        [{1: 3.0, 2: 2.0, 3: 1.0},
                         {0: 3.0, 2: 2.0, 3: 1.0},
                         {0: 2.0, 1: 2.0, 3: 1.0},
                         {0: 1.0, 1: 1.0, 2: 1.0}],
                        [1.0] * 4),
        "path16_geom": (P16, None, [float(2 ** (15 - v)) for v in range(16)]),  # 对照：非正则
    }

    print("=" * 86)
    print(f"从 K₄ 长出结构：忠实引擎 top2 记忆 + T={T} 毁灭，{cycles}周期×{gens_per}代")
    print("=" * 86)
    out = {}
    for name, (adj, w, seed) in configs.items():
        V = len(adj)
        eng = K4Engine(adj, weights=w, seed=seed, T=T, retention="top2")
        last = {}
        gini_traj = []
        for cyc in range(cycles):
            for _ in range(gens_per):
                eng.run_generation()
            nb, br, q, ratio, freq = eng.branches()
            g = gini(eng.record)
            gini_traj.append((eng.gen, g, q, nb))
            last = {"n_branch": nb, "branches": br, "q": q,
                    "top2_ratio": (round(ratio, 3) if ratio else None),
                    "freq": (round(freq, 4) if freq else None),
                    "gini_record": g}
            if (cyc + 1) % 15 == 0:
                print(f"  [{name:<13}] gen{eng.gen:>5}: 支={nb:>2} q={q:.4f} "
                      f"Gini={g:.4f} 主支={br[:4]} 频率={freq and round(freq,4)}", flush=True)
        out[name] = {"config": {"V": V, "weights": w is not None},
                     "final": last, "gini_traj": gini_traj}
    print()
    for name in out:
        f = out[name]["final"]
        print(f"  {name:<13}: 支={f['n_branch']} q={f['q']:.4f} Gini={f['gini_record']:.4f} "
              f"主支={f['branches'][:5]} 频率={f['freq']}")
    os.makedirs(RES, exist_ok=True)
    fp = os.path.join(RES, f"k4_growth_T{T}.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"\n已存 {fp}")


if __name__ == "__main__":
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    run(T=T)
