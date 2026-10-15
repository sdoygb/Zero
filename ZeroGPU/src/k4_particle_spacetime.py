#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_particle_spacetime.py —— 粒子对时空的动力学：物质弯曲时空，时空势随粒子动态
====================================================================================
用户洞察：时空刚开始静态的可能，但有了粒子就不可能是静态了——粒子对时空有反作用。
  = 广义相对论核心：物质(能量-动量)弯曲时空，时空反过来引导物质 → 双向耦合。
Zero 依据：
  Z033: 离散泊松 Lφ=−x 图上有解 → 势 φ（时空形状）
  Z089/Z091: 粒子谱 = 稳态 + 3 频率槽（质量归一 [0.413,0.707,1.0]）
  Z088: L0-L1-L2 层间双向动力学
设计（粒子↔时空双向耦合）：
  粒子质量分布(频率槽读出) 作为时空源 → 物质弯曲时空
  时空源 x_t = (1−λ)·荷分布 + λ·粒子质量分布（粒子参与弯曲时空）
  粒子随时间(代)演化(毁灭-重播种+历史回灌) → 时空源随代变 → 势 φ 动态模流
  势反馈演化(时空引导粒子) → 双向耦合
"""
from __future__ import annotations
import sys, os, json, math
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

def freq_mass(record, V):
    """从 record 读 3 频率槽质量（Z089/Z091）：相邻比值对数 → 质量归一。"""
    w = np.sort(record / record.sum())[::-1]
    masses = np.zeros(V)
    for i in range(min(3, V-1)):
        r = w[i] / (w[i+1] + 1e-12)
        masses[i] = abs(math.log(r)) if r != 0 else 0.0
    if masses.max() > 0:
        masses = masses / masses.max()
    return masses

class ParticleSpacetimeEngine(LoyaltyEngine):
    """粒子↔时空双向耦合：粒子弯曲时空(源)，时空势引导演化。"""
    def __init__(self, adj, L=4, seed=None, T=4, retention="top2", lam=0.0, beta=0.2, feedback=0.0):
        super().__init__(adj, L=L, seed=seed, T=T, retention=retention)
        self.lam = lam          # 粒子参与时空源权重
        self.beta = beta        # 历史回灌参与 seed
        self.feedback = feedback  # 时空势反馈到演化(时空引导粒子)权重

    def run_generation(self):
        V, L = self.V, self.L
        A = _adj_matrix(self.adj)
        W = self.W.copy()
        W[np.arange(V), np.arange(V)] += self.seed
        for _ in range(L):
            W = W @ A
        diag = W[np.arange(V), np.arange(V)]
        self.gen_closures = diag.copy()
        self.record = self.record + diag                       # L1 历史
        # 粒子质量分布（3 频率槽）作为时空源 → 物质弯曲时空
        masses = freq_mass(self.record, V)
        hist = self.record / self.record.sum() if self.record.sum() > 0 else np.ones(V)/V
        retained = self._retained()
        # 时空源 x = (1−λ)·荷分布 + λ·粒子质量（粒子参与弯曲时空）
        # 势反馈到 seed：时空引导粒子
        if self.feedback > 0:
            # 用上一代势(由当前荷+粒子源解得)反馈 seed（双向耦合）
            x_src = (1-self.lam) * (diag/diag.sum()) + self.lam * (masses/masses.sum() if masses.sum()>0 else np.ones(V)/V)
            x_src -= x_src.mean()
            phi = solve_potential(laplacian(self.adj), x_src)
            seed_fb = (1 - self.feedback) * ((1-self.beta)*retained + self.beta*hist) + self.feedback * (phi - phi.min() + 1e-9)
            self.seed = seed_fb
        else:
            self.seed = (1-self.beta) * retained + self.beta * hist
        if float(self.seed.sum()) == 0:
            self.seed = np.ones(V, dtype=np.float64)
        self.gen += 1
        if self.gen % self.T == 0:
            self.destructions.append((self.gen, self.branches()))
            self.W = np.zeros((V, V), dtype=np.float64)
        self._normalize()
        return self

print("="*88)
print("粒子对时空的动力学：物质弯曲时空，时空势随粒子动态")
print("="*88)
L4 = laplacian(K4); n = len(K4)
print(f"K₄ · Lφ=−x(Z033) · 粒子谱3频率槽质量(Z089/Z091) · 粒子↔时空双向耦合")

scenarios = [
    ("无粒子(λ=0,基线)", dict(lam=0.0, beta=0.2, feedback=0.0)),
    ("粒子弯曲时空(λ=0.5)", dict(lam=0.5, beta=0.2, feedback=0.0)),
    ("粒子弯曲时空(λ=1.0)", dict(lam=1.0, beta=0.2, feedback=0.0)),
    ("双向耦合(λ=0.5+反馈)", dict(lam=0.5, beta=0.2, feedback=0.5)),
    ("双向耦合(λ=1.0+反馈)", dict(lam=1.0, beta=0.2, feedback=0.5)),
]
out = {}
for label, cfg in scenarios:
    eng = ParticleSpacetimeEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2", **cfg)
    gens = 400
    pot_dyn = 0.0; closure_dyn = 0.0; last_p = None; last_c = None
    traj = []
    for g in range(gens):
        eng.run_generation()
        x = eng.gen_closures / eng.gen_closures.sum() - 1.0/n
        masses = freq_mass(eng.record, n)
        x_src = (1-cfg["lam"]) * (x + 1.0/n) + cfg["lam"] * (masses/masses.sum() if masses.sum()>0 else np.ones(n)/n)
        x_src -= x_src.mean()
        phi = solve_potential(L4, x_src)
        if last_p is not None:
            pot_dyn += float(np.linalg.norm(phi - last_p))
            closure_dyn += float(np.linalg.norm(x - last_c))
        last_p, last_c = phi, x
        if (g+1) in (1,5,50,200,400):
            traj.append({"gen": g+1, "phi": [round(v,4) for v in phi],
                         "mass": [round(v,3) for v in masses]})
    # 末段波动（收敛性）
    eng2 = ParticleSpacetimeEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2", **cfg)
    for _ in range(gens): eng2.run_generation()
    tail_phi = []
    for _ in range(200):
        eng2.run_generation()
        m2 = freq_mass(eng2.record, n)
        x2 = eng2.gen_closures/eng2.gen_closures.sum()-1.0/n
        xs = (1-cfg["lam"])*(x2+1.0/n) + cfg["lam"]*(m2/m2.sum() if m2.sum()>0 else np.ones(n)/n)
        xs -= xs.mean()
        tail_phi.append(solve_potential(L4, xs))
    tail_std = float(np.std(np.array(tail_phi), axis=0).max())
    out[label] = {"pot_dyn": round(pot_dyn,4), "closure_dyn": round(closure_dyn,4),
                  "tail_std": round(tail_std,6), "traj": traj}
    stat = "静态不动点" if tail_std < 1e-4 and pot_dyn < 1e-3 else ("持续动态(模流)" if tail_std >= 1e-3 else "收敛新稳态")
    print(f"  {label}: 逐代势动态={pot_dyn:.4f}  荷动态={closure_dyn:.4f}  末段波动={tail_std:.6f} → {stat}")

print("\n" + "="*88)
print("结论：粒子对时空的动力学是否让时空动态？")
print("="*88)
for label, cfg in scenarios:
    r = out[label]
    stat = "静态不动点" if r["tail_std"] < 1e-4 and r["pot_dyn"] < 1e-3 else ("持续动态(模流)" if r["tail_std"] >= 1e-3 else "收敛新稳态")
    print(f"  {label}: 势动态={r['pot_dyn']:.4f} 末段波动={r['tail_std']:.6f} → {stat}")

os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_particle_spacetime.json"), "w", encoding="utf-8") as f:
    json.dump({"config": {"graph": "K4", "gens": 400, "T": 4,
                          "GR_analog": "物质弯曲时空, 时空引导物质(双向耦合)",
                          "Z033": "Lφ=−x 势可算", "Z089": "粒子谱3频率槽质量"},
               "scenarios": out}, f, ensure_ascii=False, indent=1)
print(f"\n已存 {os.path.join(RES,'k4_particle_spacetime.json')}")
