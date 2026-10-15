#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
z_integrate.py —— 场与力集成进 Zero 模拟器：多代演化看是否聚合
================================================================
把 Z093 的"场（S₃/Z₃ 方向结构）+ 力（约束回复力 D）"集成进忠实引擎，
让场和力作为粒子的演化参与者（回写调制 L2），一代代演化，观测是否聚合。

集成语义（忠实 Zero）：
  * 粒子 = L2 稳态 record 分布 p = record/Σ（K₄，3 频率槽 = 三代）
  * 场   = 粒子分布的"方向结构"：主分支数 / 频率槽（Z093：S₃⊃Z₃ 三代）
  * 力   = 约束回复力 D_v = 1/V − p_v（闭合约束把粒子权重"拉向"对称平衡的回复强度）
  * 回写 = seed 重播种时叠加力项 seed += γ·D（L1→L2 边②，力驱动粒子演化）
  * 聚合观测：record 是否收敛到稳定分布(吸引子)、频率槽数、均匀度 q=Σw²、香农熵

γ（力强度）扫描，看强/弱力回写下粒子谱是聚合稳定、破坏、还是出现新模式。
诚实边界：力 D 是标量约束回复力（Z093），非矢量力；回写是否真实对应
L0-L1-L2 双向动力学，另行核对（Z088 六边、Z092 回写负结论）。
"""
import sys, os, math, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loyalty_engine import LoyaltyEngine
from z0_engine import path_adj

RES = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
SIGNIF = 1e-4

# K₄（Z089：3 频率槽承载），seed[8,4,1.5,1]（稳态 record = seed 归一化）
K4 = [[1,2,3],[0,2,3],[0,1,3],[0,1,2]]
SEED4 = [8.0, 4.0, 1.5, 1.0]


class IntegEngine(LoyaltyEngine):
    """场+力集成引擎：seed 重播种叠加力回写，观测聚合。"""

    def __init__(self, adj, seed, T=4, retention="top2", gamma=0.0, anchor="attractor"):
        super().__init__(adj, seed=seed, T=T, retention=retention)
        self.gamma = float(gamma)          # 力强度
        self.V = len(adj)
        self.anchor = anchor               # 力锚点: 'attractor'(seed比稳态) | 'uniform'(均匀1/V)
        self.attractor = np.array(seed, dtype=np.float64)
        self.attractor = self.attractor / self.attractor.sum()  # 闭合稳态分布
        self.force_hist = []               # 每代力向量
        self.field_hist = []               # 每代场(主分支)结构

    def _force(self, p):
        # 约束回复力：回复到锚点。attractor=粒子谱稳态(seed比), uniform=对称平衡(1/V)
        if self.anchor == "uniform":
            return (1.0 / self.V - p)
        return (self.attractor - p)        # 回复到闭合稳态(粒子谱吸引子)

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

        # 场/力：从当前 record 分布读
        tot = self.record.sum()
        p = self.record / tot if tot > 0 else np.ones(V)/V
        D = self._force(p)
        self.force_hist.append(D.copy())

        # 重播种：top2 + 力回写（γ·D 调制 L2 种子）
        self.seed = self._retained()
        if self.gamma != 0.0:
            self.seed = self.seed + self.gamma * D
        if float(self.seed.sum()) == 0:
            self.seed = np.ones(V, dtype=np.float64)

        self.gen += 1
        if self.gen % self.T == 0:
            self.destructions.append((self.gen, self.branches()))
            self.W = np.zeros((V, V), dtype=np.float64)
        self._normalize()
        return self

    def field_structure(self):
        """场的方向结构：主分支(频率槽)数、q(均匀度)、香农熵。"""
        w = self.record / self.record.sum()
        mx = w.max()
        nb = int((w > SIGNIF * mx).sum())
        q = float(w @ w)
        ent = float(-np.sum(w * np.log(w + 1e-12)))
        # 频率槽：取显著分支的相邻比 → log 频率
        sig = sorted([(i, float(wi)) for i, wi in enumerate(w) if wi > SIGNIF * mx],
                     key=lambda x: -x[1])
        ratio = sig[0][1] / sig[1][1] if len(sig) >= 2 else None
        freq = abs(math.log(ratio)) if ratio else None
        return {"n_field": nb, "q": q, "entropy": ent,
                "top": [round(wi, 4) for _, wi in sig], "ratio": ratio,
                "freq": round(freq, 4) if freq else None}


def run_integration(adj, seed, gamma, cycles, gens_per, anchor="attractor"):
    eng = IntegEngine(adj, seed=seed, T=4, retention="top2", gamma=gamma, anchor=anchor)
    V = len(adj)
    per = []
    for cyc in range(cycles):
        for _ in range(gens_per):
            eng.run_generation()
        st = eng.field_structure()
        per.append({"cycle": cyc+1, "gen": eng.gen,
                    "n_field": st["n_field"], "q": round(st["q"],6),
                    "entropy": round(st["entropy"],6),
                    "top": st["top"], "freq": st["freq"]})
    # 末态
    w = eng.record / eng.record.sum()
    return eng, per, {"final_record": [round(float(x),6) for x in eng.record],
                      "final_p": [round(float(x),6) for x in w],
                      "n_field_final": st["n_field"], "q_final": round(st["q"],6)}


def main():
    gammas = [0.0, 0.1, 1.0, 10.0]
    anchors = ["attractor", "uniform"]
    cycles, gens_per = 200, 200     # 200 周期 × 200 代 = 40000 代
    print("="*80)
    print("场+力集成进 Zero 模拟器：多代演化聚合检验  (K₄, T=4, top2)")
    print("力回写 seed += γ·D ; 锚点: attractor=粒子谱稳态(seed比) | uniform=均匀1/V")
    print("每配置 %d×%d=%d 代" % (cycles, gens_per, cycles*gens_per))
    print("="*80)
    summary = {}
    for anchor in anchors:
        for g in gammas:
            eng, per, fin = run_integration(K4, SEED4, g, cycles, gens_per, anchor=anchor)
            nfs = set(c["n_field"] for c in per)
            freqs = sorted(set(c["freq"] for c in per if c["freq"] is not None))
            qs = [c["q"] for c in per]
            tag = "基线" if g == 0.0 else ""
            print(f"\n[{anchor}] γ={g:<5} 末p={fin['final_p']}  "
                  f"场分支={fin['n_field_final']}  q末={fin['q_final']:.5f}  "
                  f"q范围=[{min(qs):.5f},{max(qs):.5f}]  频率集={freqs[:5]}")
            key = f"{anchor}_γ{g}"
            summary[key] = {"anchor": anchor, "gamma": g,
                            "final_record": fin["final_record"],
                            "final_p": fin["final_p"],
                            "n_field_final": fin["n_field_final"],
                            "q_final": fin["q_final"],
                            "q_range": [round(min(qs),5), round(max(qs),5)],
                            "n_field_set": sorted(nfs),
                            "freq_set": freqs[:8]}
    fp = os.path.join(RES, "z_integrate.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump({"anchor_scan": summary, "config": {"adj": "K4",
                  "seed": SEED4, "T": 4, "retention": "top2",
                  "cycles": cycles, "gens_per": gens_per}}, f,
                  ensure_ascii=False, indent=1)
    print(f"\n已存 {fp}")


if __name__ == "__main__":
    main()
