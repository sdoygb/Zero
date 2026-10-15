# -*- coding: utf-8 -*-
"""
ZeroGPU · 自组织 K₄：闭合回馈改写边权 → 有效度自发不均匀

背景：正则 K₄ 只继承 seed 形状，不产生结构（gamma_gap：度均匀⟹无自发结构）。
本脚本让边权由闭合回馈**自发演化**（Hebbian 型：共同闭合增强边权），
使有效度（加权度）不均匀 → 演化才能"长出"结构。

测两种关键情形：
  (a) 均匀 seed + 边权反馈：能否**对称破缺**（无 seed 偏置下自发产生结构）？
  (b) geom seed + 边权反馈：放大 seed 不均匀，还是产生新结构？

用法：/usr/local/bin/python3 selforg_k4.py [cycles] [gens_per]
"""
from __future__ import annotations
import sys, os, math, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from loyalty_engine import LoyaltyEngine
from z0_vectorized import NORM_THRESHOLD

RES = os.path.normpath(os.path.join(HERE, "..", "results"))


def gini(arr):
    x = sorted(float(v) for v in arr if v)
    n = len(x)
    if n == 0:
        return 0.0
    s = sum(x)
    num = sum((2 * (i + 1) - n - 1) * x[i] for i in range(n))
    return num / (n * s)


class SelfOrgK4(LoyaltyEngine):
    """边权由闭合回馈演化的 K₄。A 不再固定，每代按 Hebbian 更新。"""

    def __init__(self, adj, L=4, seed=None, T=4, retention="top2", beta=1.0,
                 lam=0.05, cap=None):
        super().__init__(adj, L=L, seed=seed, T=T, retention=retention, cap=cap)
        V = len(adj)
        self.A = np.ones((V, V), dtype=np.float64) * 0.0
        for u, nb in enumerate(adj):
            for w in nb:
                self.A[u][w] = 1.0
        self.beta = beta          # 反馈强度
        self.lam = lam            # 平滑系数（0<λ<1，防突变）
        self.edge_history = []

    def _update_weights(self):
        """闭合回馈改写边权：边(u,v)权重随两端闭合**比例**共同增长（Hebbian）。
        用归一化闭合份额（有界）而非绝对量，避免 K₄ 高度连通下的数值爆炸。"""
        c = self.gen_closures
        cs = float(c.sum())
        cn = c / cs if cs > 0 else np.ones(self.V) / self.V
        V = self.V
        target = np.ones_like(self.A)
        for u in range(V):
            for w in range(V):
                if self.A[u][w] > 0:
                    target[u][w] = 1.0 + self.beta * cn[u] * cn[w]
        # 行归一化（保总流尺度），平滑更新
        row_sum = target.sum(axis=1, keepdims=True)
        target = target / np.maximum(row_sum, 1e-300) * V
        self.A = (1 - self.lam) * self.A + self.lam * target

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
        # 闭合回馈改写边权（L0 内容累积改变结构条件）
        self._update_weights()
        self.gen += 1
        if self.gen % self.T == 0:
            self.W = np.zeros((self.V, self.V), dtype=np.float64)
        self._normalize()
        return self

    def eff_deg(self):
        """有效度（加权度）分布。"""
        return self.A.sum(axis=1)


def run(cycles=80, gens_per=200):
    K4 = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]
    configs = [
        # (名称, seed, beta, lam)
        ("uniform_b1",   [1.0, 1.0, 1.0, 1.0], 1.0, 0.05),
        ("uniform_b10",  [1.0, 1.0, 1.0, 1.0], 10.0, 0.05),
        ("geom_b1",      [8.0, 4.0, 2.0, 1.0], 1.0, 0.05),
        ("geom_b10",     [8.0, 4.0, 2.0, 1.0], 10.0, 0.05),
    ]
    print("=" * 86)
    print("自组织 K₄：闭合回馈改写边权 → 有效度能否自发不均匀")
    print("=" * 86)
    out = {}
    for name, seed, beta, lam in configs:
        eng = SelfOrgK4(K4, seed=seed, T=4, retention="top2", beta=beta, lam=lam)
        last = {}
        for cyc in range(cycles):
            for _ in range(gens_per):
                eng.run_generation()
            nb, br, q, ratio, freq = eng.branches()
            ed = eng.eff_deg()
            g_rec = gini(eng.record)
            g_deg = gini(ed)
            last = {"n_branch": nb, "branches": br, "q": q,
                    "freq": (round(freq, 4) if freq else None),
                    "gini_record": g_rec, "gini_effdeg": g_deg,
                    "eff_deg": [round(float(x), 3) for x in ed],
                    "weight_mat": [[round(float(x), 2) for x in r] for r in eng.A]}
            if (cyc + 1) % 20 == 0:
                print(f"  [{name:<11}] gen{eng.gen:>5}: 支={nb:>2} q={q:.4f} "
                      f"Gini(record)={g_rec:.4f} Gini(有效度)={g_deg:.4f} "
                      f"有效度={[round(float(x),2) for x in ed]} 频率={freq and round(freq,4)}", flush=True)
        out[name] = {"beta": beta, "lam": lam, "final": last}
    print()
    for name in out:
        f = out[name]["final"]
        print(f"  {name:<11}: 支={f['n_branch']:>2} q={f['q']:.4f} "
              f"Gini(record)={f['gini_record']:.4f} Gini(有效度)={f['gini_effdeg']:.4f} "
              f"有效度={f['eff_deg']} 频率={f['freq']}")
        print(f"           权重矩阵 = {f['weight_mat']}")
    os.makedirs(RES, exist_ok=True)
    fp = os.path.join(RES, "selforg_k4.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"\n已存 {fp}")


if __name__ == "__main__":
    c = int(sys.argv[1]) if len(sys.argv) > 1 else 80
    g = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    run(cycles=c, gens_per=g)
