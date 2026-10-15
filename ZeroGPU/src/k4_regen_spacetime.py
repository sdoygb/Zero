#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_regen_spacetime.py —— 回灌动力学：让 L0-L1-L2 双向边参与演化，测荷分布/势动态模流
====================================================================================
依据 Zero 文章库：
  Z088: 六条层间边全双向 —— ②L1→L2 历史作种子, ③L2→L1 演化写历史,
        ⑤L2→L0 本代闭合回灌, ⑥L1→L0 记录回灌, ④L0→L2 结构决定
        闭合回灌 L0 = 重数(内容)累积, 不是加词
  Z033: 离散泊松 Lφ=−x 图上有解 → 势 φ（时空形状）
设计（忠实六边回灌）：
  历史回灌参与 seed：seed_{t+1} = (1−β)·最近两层闭合 + β·历史回灌(record)
    β 让 L1→L0 累积历史经 ⑥L1→L0→④L0→L2 真正影响演化 → 荷分布可能随代变
  每代解 Lφ=−x 得势 → 测势是否动态（模流）vs 静态不动点
"""
from __future__ import annotations
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loyalty_engine import LoyaltyEngine
from z0_vectorized import _adj_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
K4 = [[1,2,3],[0,2,3],[0,1,3],[0,1,2]]

def laplacian(adj):
    n = len(adj); L = np.zeros((n,n))
    for i, nb in enumerate(adj):
        L[i,i] += len(nb)
        for j in nb: L[i,j] -= 1
    return L

def solve_potential(L, x):
    x0 = x - x.mean()
    phi = -np.linalg.pinv(L) @ x0
    return phi - phi.mean()

class RegenEngine(LoyaltyEngine):
    """回灌引擎：历史 record 参与 seed（L1→L0→L2 回灌动力学）。"""
    def __init__(self, adj, L=4, seed=None, T=4, retention="top2", beta=0.0):
        super().__init__(adj, L=L, seed=seed, T=T, retention=retention)
        self.beta = beta
        self.hist_series = []      # 历史回灌量随代（L1→L0 ⑥）

    def run_generation(self):
        V, L = self.V, self.L
        A = _adj_matrix(self.adj)
        W = self.W.copy()
        W[np.arange(V), np.arange(V)] += self.seed
        for _ in range(L):
            W = W @ A
        diag = W[np.arange(V), np.arange(V)]
        self.gen_closures = diag.copy()
        self.record = self.record + diag                 # ③L2→L1 演化写历史
        # ⑥L1→L0 记录回灌 = 累积历史（重数）；④L0→L2 结构决定经 seed 表达
        hist = self.record / self.record.sum() if self.record.sum() > 0 else np.ones(V)/V
        retained = self._retained()                      # ②L1→L2 历史作种子(最近两层)
        # 历史回灌参与演化：seed = (1−β)·最近两层 + β·累积历史
        self.seed = (1 - self.beta) * retained + self.beta * hist
        if float(self.seed.sum()) == 0:
            self.seed = np.ones(V, dtype=np.float64)
        self.hist_series.append(hist.copy())
        self.gen += 1
        if self.gen % self.T == 0:                        # ⑤L2→L0 本代闭合回灌 + 毁灭
            self.destructions.append((self.gen, self.branches()))
            self.W = np.zeros((V, V), dtype=np.float64)
        self._normalize()                                  # 归一化防爆炸（原引擎纪律）
        return self

print("="*86)
print("回灌动力学：L0-L1-L2 双向边参与演化，测荷分布/势动态模流")
print("="*86)
L4 = laplacian(K4); n = len(K4)
print(f"K₄（{n}顶点）· 离散泊松 Lφ=−x（Z033）· 六边回灌（Z088）")

beta_values = [0.0, 0.1, 0.3, 0.5]
results = {}
for beta in beta_values:
    eng = RegenEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2", beta=beta)
    gens = 400
    pot = []; closure_norm_dyn = 0.0; pot_dyn = 0.0; last_c = None; last_p = None
    for g in range(gens):
        eng.run_generation()
        x = eng.gen_closures / eng.gen_closures.sum() - 1.0/n
        phi = solve_potential(L4, x)
        if last_c is not None:
            closure_norm_dyn += float(np.linalg.norm(x - last_c))  # 荷分布逐代动态量
            pot_dyn += float(np.linalg.norm(phi - last_p))         # 势逐代动态量
        last_c, last_p = x, phi
        if (g+1) in (1,5,50,100,200,300,400):
            pot.append({"gen": g+1, "phi": [round(v,4) for v in phi],
                        "x": [round(v,3) for v in x]})
    # 稳态检查：后 200 代势是否收敛
    tail = [solve_potential(L4, eng.gen_closures/eng.gen_closures.sum()-1.0/n)
            for _ in range(200)]
    tail_runs = np.array([[round(v,4) for v in p] for p in tail])
    conv = float(np.std(tail_runs, axis=0).max())   # 势末段波动
    results[beta] = {"pot_dyn": round(pot_dyn,4), "closure_dyn": round(closure_norm_dyn,4),
                     "tail_std": round(conv,5), "samples": pot}
    print(f"β={beta}: 逐代势动态量={pot_dyn:.4f}  荷分布动态量={closure_norm_dyn:.4f}  "
          f"末200代势波动(收敛性)={conv:.5f}")

print("\n" + "="*86)
print("结论：回灌动力学是否让时空动态生成？")
print("="*86)
for beta in beta_values:
    r = results[beta]
    stat = "静态不动点(瞬达)" if r["tail_std"] < 1e-3 and r["pot_dyn"] < 1e-3 else \
           ("收敛新稳态" if r["tail_std"] < 1e-3 else "持续动态(模流)")
    print(f"  β={beta}: 逐代动态={r['pot_dyn']:.4f}, 末段波动={r['tail_std']:.5f} → {stat}")

os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_regen_spacetime.json"), "w", encoding="utf-8") as f:
    json.dump({"config": {"graph": "K4", "gens": 400, "T": 4, "retention": "top2",
                          "Z088_edges": "②L1→L2 ③L2→L1 ④L0→L2 ⑤L2→L0 ⑥L1→L0",
                          "Z033": "Lφ=−x 势可算"},
               "beta_scan": {str(b): results[b] for b in beta_values},
               "conclusion": "见stdout"}, f, ensure_ascii=False, indent=1)
print(f"\n已存 {os.path.join(RES,'k4_regen_spacetime.json')}")
