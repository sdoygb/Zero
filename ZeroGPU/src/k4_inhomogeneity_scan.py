#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_inhomogeneity_scan.py —— 不均匀度扫描：物质分布不均匀性 → 时空动态活跃度
====================================================================================
验证 Z101/打破边界的深层结论：时空动态活跃度依赖物质分布的不均匀性。
  质量分布 m_i ∝ w_i^p (w=record权重, p=不均匀度参数)
    p→0: 均匀(平坦时空)   p大: 不均匀(弯曲/动态)
  双向耦合(λ=1.0, feedback=0.5) → 测模流活跃度(末段波动/逐代动态/斜率)
预期: p→0 末段波动→0(静态); p大 末段波动大(动态) → 验证"不均匀→动态"定量关系
"""
from __future__ import annotations
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from k4_particle_spacetime import ParticleSpacetimeEngine, laplacian, solve_potential
from loyalty_engine import LoyaltyEngine
from z0_vectorized import _adj_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
K4 = [[1,2,3],[0,2,3],[0,1,3],[0,1,2]]
L4 = laplacian(K4); n = 4

class InhomogEngine(LoyaltyEngine):
    """粒子源不均匀度 p 控制的耦合引擎。"""
    def __init__(self, adj, L=4, seed=None, T=4, retention="top2", p=1.0, beta=0.2, feedback=0.5):
        super().__init__(adj, L=L, seed=seed, T=T, retention=retention)
        self.p = p; self.beta = beta; self.feedback = feedback

    def mass_dist(self, record):
        w = record / record.sum() if record.sum() > 0 else np.ones(self.V)/self.V
        m = w ** self.p                      # p=0 均匀, p大 不均匀
        if m.sum() > 0:
            m = m / m.sum()
        return m

    def run_generation(self):
        V, L = self.V, self.L; A = _adj_matrix(self.adj)
        W = self.W.copy(); W[np.arange(V), np.arange(V)] += self.seed
        for _ in range(L): W = W @ A
        diag = W[np.arange(V), np.arange(V)]
        self.gen_closures = diag.copy(); self.record = self.record + diag
        m = self.mass_dist(self.record)
        hist = self.record / self.record.sum() if self.record.sum() > 0 else np.ones(V)/V
        retained = self._retained()
        if self.feedback > 0:
            xs = m - 1.0/V                      # 粒子源(零和), 纯粒子弯曲时空
            phi = solve_potential(laplacian(self.adj), xs)
            self.seed = (1-self.feedback)*((1-self.beta)*retained + self.beta*hist) + self.feedback*(phi - phi.min() + 1e-9)
        else:
            self.seed = (1-self.beta)*retained + self.beta*hist
        if float(self.seed.sum()) == 0: self.seed = np.ones(V)
        self.gen += 1
        if self.gen % self.T == 0:
            self.destructions.append((self.gen, self.branches())); self.W = np.zeros((V,V))
        self._normalize()
        return self

print("="*90)
print("不均匀度扫描：物质分布不均匀性 → 时空动态活跃度")
print("="*90)
print("质量分布 m_i ∝ w_i^p  (p=不均匀度: 0均匀 → 大不均匀)")

ps = [0.0, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
N = 2000; results = {}
for p in ps:
    eng = InhomogEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2", p=p)
    phis = []
    for _ in range(N):
        eng.run_generation()
        m = eng.mass_dist(eng.record)
        xs = m - 1.0/n
        phis.append(solve_potential(L4, xs))
    phis = np.array(phis)
    tail_std = float(np.std(phis[-600:], axis=0).max())     # 末段波动(动态活跃度)
    dyn_total = float(np.sum(np.linalg.norm(np.diff(phis, axis=0), axis=1)))  # 全程逐代动态
    sl = [round(np.polyfit(np.arange(N), phis[:,v], 1)[0], 7) for v in range(n)]
    # 不均匀度指标: 粒子源的质量集中度(有效均匀度)
    w = eng.record / eng.record.sum(); m = w ** p; m = m / m.sum()
    gini_like = float(np.var(m))                            # 源方差=不均匀度
    results[p] = {"tail_std": tail_std, "dyn_total": round(dyn_total,4),
                  "slopes": sl, "src_var": round(gini_like,6)}
    print(f"  p={p:>4}: 源方差(不均匀度)={gini_like:.5f}  末段波动={tail_std:.6f}  "
          f"全程动态={dyn_total:.4f}  斜率={sl}")

print("\n" + "="*90)
print("定量关系：不均匀度(源方差) → 动态活跃度(末段波动)")
print("="*90)
for p in ps:
    r = results[p]
    print(f"  p={p:>4}: 源方差={r['src_var']:.5f} → 末段波动={r['tail_std']:.6f}  "
          f"({'动态' if r['tail_std']>5e-4 else ('微动' if r['tail_std']>1e-4 else '静态/收敛')})")

os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_inhomogeneity_scan.json"), "w", encoding="utf-8") as f:
    json.dump({"config": {"graph": "K4", "N": N, "coupling": "lambda=1.0, feedback=0.5",
                          "mass": "m_i ∝ w_i^p"},
               "results": {str(p): results[p] for p in ps},
               "thesis": "物质分布不均匀性 → 时空动态活跃度(GR: 物质分布决定时空弯曲)"},
              f, ensure_ascii=False, indent=1)
print(f"\n已存 {os.path.join(RES,'k4_inhomogeneity_scan.json')}")
