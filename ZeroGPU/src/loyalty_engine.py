# -*- coding: utf-8 -*-
"""
ZeroGPU · 忠实层间动力学引擎（LoyaltyEngine）

按 Zero 正式篇参数忠实实现：
  * 记忆 = 保留最高两层（Z017 定理 7）→ retention='top2'（seed = 前两层闭合 + 本代）
  * 毁灭周期 T = L = 4（Z022 定理 4）→ 每 T 代清空 L2 活动层（W=0），L1 保留最高两层重播种，L0 保留规则
  * 闭合 = 回到起点（closure_verdict，读法①）——已忠实

观测：每个毁灭周期内的稳态支结构（支数、主两支比、内禀频率）。
问：毁灭重生后是反复长出同一频率（单稳态），还是长出不同频率/新支（级联多粒子）？

用法：/usr/local/bin/python3 loyalty_engine.py [path16] [T=4] [top2]
"""
from __future__ import annotations
import sys, os, time, math, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import path_adj
from z0_vectorized import ZeroCoreVec, NORM_THRESHOLD, _adj_matrix

RES = os.path.normpath(os.path.join(HERE, "..", "results"))
SIGNIF = 1e-4


class LoyaltyEngine:
    """忠实三层闭环：top2 记忆 + 周期毁灭 T + 闭合回灌监测。"""

    def __init__(self, adj, L=4, seed=None, T=4, retention="top2", cap=None):
        self.adj = [tuple(a) for a in adj]
        self.V = len(adj)
        self.A = _adj_matrix(adj)
        self.L = int(L)
        self.T = int(T)             # 毁灭周期（Z022 定理4 → 4）
        self.retention = retention
        self.cap = cap
        self.seed = np.array(seed if seed is not None else [1.0] * self.V, dtype=np.float64)
        self.W = np.zeros((self.V, self.V), dtype=np.float64)
        self.record = np.zeros(self.V, dtype=np.float64)   # L1 历史（L0 累积）
        self.term = np.zeros(self.V, dtype=np.float64)
        self.gen = 0
        self.gen_closures = np.zeros(self.V, dtype=np.float64)
        self.layers = []            # 每代闭合（L1 亚层栈）
        self.norm_events = 0
        self.destructions = []      # (gen, 毁灭前支结构)

    def _normalize(self):
        for a in (self.record, self.W, self.seed):
            t = a.sum()
            if t > NORM_THRESHOLD:
                a /= t
                self.norm_events += 1

    def _retained(self):
        # 保留最高两层（Z017 定理7）：前一层 + 本代闭合
        self.layers.append(self.gen_closures.copy())
        if self.retention == "top2":
            if len(self.layers) >= 2:
                return self.layers[-2] + self.gen_closures
            return self.gen_closures.copy()
        if self.retention == "last":
            return self.gen_closures.copy()
        if self.retention == "all":
            return self.record.copy()
        raise ValueError(self.retention)

    def branches(self):
        w = self.record / self.record.sum()
        mx = w.max()
        sig = [(i, float(w[i])) for i in range(self.V) if w[i] > SIGNIF * mx]
        sig.sort(key=lambda x: -x[1])
        q = float(w @ w)
        ratio = sig[0][1] / sig[1][1] if len(sig) >= 2 else None
        freq = abs(math.log(ratio)) if ratio else None
        return len(sig), [(i, round(x, 4)) for i, x in sig], q, ratio, freq

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
        # 重播种（L1 → L2，保留最高两层）
        self.seed = self._retained()
        if float(self.seed.sum()) == 0:
            self.seed = np.ones(self.V, dtype=np.float64)
        # 周期毁灭（Z022 定理4：T=4）：L2 全毁，L1 最高两层保留作 seed，L0 规则保留
        self.gen += 1
        if self.gen % self.T == 0:
            self.destructions.append((self.gen, self.branches()))
            self.W = np.zeros((self.V, self.V), dtype=np.float64)
        self._normalize()
        return self


def run(name, T=4, retention="top2", cycles=200, gens_per=200):
    n = int(name.replace("path", ""))
    adj = path_adj(n); V = len(adj)
    seed = [float(2 ** (V - 1 - v)) for v in range(V)]
    eng = LoyaltyEngine(adj, seed=seed, T=T, retention=retention)
    print("=" * 76)
    print(f"忠实引擎：path{n} geom seed  {retention}记忆  T={T}毁灭  {cycles}周期×{gens_per}代")
    print("=" * 76)
    per_cycle = []
    for cyc in range(cycles):
        for _ in range(gens_per):
            eng.run_generation()
        nb, br, q, ratio, freq = eng.branches()
        per_cycle.append({"cycle": cyc + 1, "gen": eng.gen, "n_branch": nb,
                          "branches": br, "q": q,
                          "top2_ratio": (round(ratio, 3) if ratio else None),
                          "freq": (round(freq, 4) if freq else None)})
        if (cyc + 1) % 20 == 0:
            print(f"  周期{cyc+1:>3} gen{eng.gen:>5}: 支={nb:>2} q={q:.5f} "
                  f"主支={br[:3]} 频率={freq and round(freq,4)}", flush=True)
    # 频率是否稳定
    freqs = [c["freq"] for c in per_cycle if c["freq"] is not None]
    n_branches = set(c["n_branch"] for c in per_cycle)
    uniq_freq = sorted(set(round(f, 4) for f in freqs))
    print(f"\n跨周期支数集={sorted(n_branches)}  不同频率数={len(uniq_freq)}  "
          f"频率集={uniq_freq[:12]}")
    rec = {"config": {"adj": name, "seed": "geom", "retention": retention,
                      "T": T, "L": 4},
           "cycles": cycles, "gens_per_cycle": gens_per, "total_gen": eng.gen,
           "uniq_freq": uniq_freq, "n_branch_set": sorted(n_branches),
           "per_cycle": per_cycle}
    os.makedirs(RES, exist_ok=True)
    fp = os.path.join(RES, f"loyalty_{name}_T{T}_{retention}.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(rec, f, ensure_ascii=False, indent=1)
    print(f"已存 {fp}")
    return rec


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "path16"
    T = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    retention = sys.argv[3] if len(sys.argv) > 3 else "top2"
    run(name, T=T, retention=retention)
