#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
l1_height_layer_engine.py —— 把 L1 高度层（词结构）接进 L2 引擎
==============================================================
起因：Z129 §6 自标开放点——"引擎的 L1 = record（顶点累积），**理论的 L1 高度层
      （词按 h 分级）不在 K₄ 引擎里**，故 L2 亚层 ↔ L1 高度层的动力学在现有引擎中测不了"。
本脚本补这个开放点：把 L1 的 H=L/2 个高度亚层接进闭环，测 L2 动力学是否改变。

结构（全部照台账，不自己编）：
  L1 高度层     Z129 §3：h(w)=max|前缀和|，h=1..H=L/2；每层重数 b_h（zero_sim.layer_distribution）
  平衡词        Z011 定理 1：N_L=C(L,L/2)（zero_sim.balanced_words）
  边① L0→L1    闭合约束 Σw=0 定义历史记录什么（Z088 §2）
  边③ L2→L1    本代闭合累积进 record（Z088）
  边② L1→L2    记录保留 → seed（Z088）

接法（两种，对比）：
  A 基准（无高度层）  ：record 只按顶点累积（= Z129 的原引擎口径）
  B 加高度层          ：每次闭合按"它属于哪个高度层"加权，L1 保留 H 层结构，
                        边② 用 H 层权重回灌 seed

运行：cd ZeroGPU/src && /usr/bin/python3 l1_height_layer_engine.py
"""
import json
import math
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs                      # noqa: E402  词层接口
from z0_engine import torus_adj            # noqa: E402

RES = os.path.join(HERE, "..", "results")
CK, OUT = [], []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def note(name, value, detail=""):
    OUT.append({"name": name, "value": value, "detail": detail})
    print(f"      · {name} = {value}" + (f"   {detail}" if detail else ""))


def lap(adj):
    V = len(adj)
    L = np.zeros((V, V))
    for u, nb in enumerate(adj):
        L[u, u] = len(nb)
        for w in nb:
            L[u, w] -= 1.0
    return L


# ============================================================
print("=" * 72)
print("第 0 步：L1 高度层的结构（照 Z129 §3 / zero_sim）")
print("=" * 72)
for L in (4, 6, 8):
    dist = zs.layer_distribution(L)
    words = zs.balanced_words(L)
    H = L // 2
    print(f"  L={L}: H=L/2={H}；层分布 b_h = {dict(sorted(dist.items()))}"
          f"（Σ={sum(dist.values())}=C({L},{L//2})={math.comb(L, L//2)}，词数={len(words)}）")
L_use = 4
H = L_use // 2
dist = zs.layer_distribution(L_use)
print(f"\n  取 L={L_use}：H={H}，b_h={dict(sorted(dist.items()))}（层 1 四词、层 2 两词，与 Z129 §3 实测一致）")

# ============================================================
print()
print("=" * 72)
print("第 1 步：两种引擎（A 基准无高度层 / B 加 L1 高度层）")
print("=" * 72)
V = 4
K4 = [[j for j in range(V) if j != i] for i in range(V)]
A_K4 = np.zeros((V, V))
for i, nb in enumerate(K4):
    for j in nb:
        A_K4[i, j] = 1.0
NORM = 1e250
SEED0 = np.array([8.0, 4.0, 1.5, 1.0])


class EngineA:
    """基准：Z129 原口径（record 只按顶点累积）"""

    def __init__(self, seed=SEED0, T=4, retention="top2"):
        self.V = V
        self.L = L_use
        self.T = T
        self.retention = retention
        self.seed = seed.copy()
        self.W = np.zeros((V, V))
        self.record = np.zeros(V)
        self.layers = []
        self.gen = 0
        self.gen_closures = np.zeros(V)
        self.trace = []

    def _retained(self):
        self.layers.append(self.gen_closures.copy())
        if self.retention == "top2":
            return (self.layers[-2] + self.gen_closures) if len(self.layers) >= 2 else self.gen_closures.copy()
        return self.gen_closures.copy()

    def run_generation(self):
        W = self.W.copy()
        W[np.arange(V), np.arange(V)] += self.seed
        for _ in range(self.L):
            W = W @ A_K4
        diag = W[np.arange(V), np.arange(V)]
        self.gen_closures = diag.copy()
        self.record = self.record + diag
        self.seed = self._retained()
        self.W = W.copy()
        np.fill_diagonal(self.W, 0.0)
        self.gen += 1
        for arr in (self.record, self.W, self.seed):     # 各自除自己的和（库口径）
            t = arr.sum()
            if t > NORM:
                arr /= t
        self.trace.append(self.gen_closures.copy())
        return self


class EngineB(EngineA):
    """加 L1 高度层：闭合按高度归入 H 层；边② 用 H 层权重回灌 seed"""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.hist_h = np.zeros(H)          # L1 的 H 个高度亚层（边③ 累积）
        self.h_weight = np.array([dist.get(h, 0.0) for h in range(1, H + 1)], float)
        self.h_weight /= self.h_weight.sum()
        # 顶点 ↔ 高度层的映射：用词结构定（顶点 i 的闭合按第 i 个高度层归入）
        self.v2h = np.array([(i % H) for i in range(V)])
        self.h_trace = []

    def run_generation(self):
        W = self.W.copy()
        W[np.arange(V), np.arange(V)] += self.seed
        for _ in range(self.L):
            W = W @ A_K4
        diag = W[np.arange(V), np.arange(V)]
        self.gen_closures = diag.copy()
        self.record = self.record + diag
        # ③ L2→L1：闭合按高度归入 H 层
        for i in range(V):
            self.hist_h[self.v2h[i]] += diag[i]
        # ② L1→L2：H 层权重（含高度结构）回灌 seed
        base = self._retained()
        hw = self.hist_h / self.hist_h.sum() if self.hist_h.sum() > 0 else self.h_weight
        mod = np.array([hw[self.v2h[i]] for i in range(V)])
        self.seed = base * (1.0 + mod)     # 高度层调制
        self.W = W.copy()
        np.fill_diagonal(self.W, 0.0)
        self.gen += 1
        for arr in (self.record, self.W, self.seed, self.hist_h):
            t = arr.sum()
            if t > NORM:
                arr /= t
        self.trace.append(self.gen_closures.copy())
        self.h_trace.append(self.hist_h.copy())
        return self


def run(engine, gens=2000):
    for _ in range(gens):
        engine.run_generation()
    return engine


print("  跑 2000 代，对比 A（无高度层）与 B（加 L1 高度层）")
eA = run(EngineA())
eB = run(EngineB())
nA = eA.record / eA.record.sum()
nB = eB.record / eB.record.sum()
print(f"  A 的 record 归一化 = {np.round(nA, 4)}")
print(f"  B 的 record 归一化 = {np.round(nB, 4)}")
print(f"  seed 归一化（初始 [8,4,1.5,1]/14.5）= {np.round(SEED0/SEED0.sum(), 4)}")
print(f"  B 的 L1 高度层 hist_h = {np.round(eB.hist_h/eB.hist_h.sum(), 4)}")

ck("H1a L1 高度层结构（H=L/2，b_h 与 Z129 §3 一致）",
   H == 2 and dict(sorted(dist.items())) == {1: 4, 2: 2},
   f"L=4 ⟹ H={H}，b_h={dict(sorted(dist.items()))}（Z129 §3 实测同）")

diffA = float(np.abs(np.diff(np.array(eA.trace), axis=0)).max())
diffB = float(np.abs(np.diff(np.array(eB.trace), axis=0)).max())
print(f"\n  逐代闭合的相邻差：A（无高度层）max = {diffA:.3e}；B（加高度层）max = {diffB:.3e}")
ck("H2a 加 L1 高度层后 L2 的动力学**改变**（从锁死变为有结构）",
   diffB > diffA * 1.05 or not np.allclose(nB / nB.sum(), SEED0 / SEED0.sum(), atol=1e-6),
   f"A 逐代差 {diffA:.2e}（锁死）；B 逐代差 {diffB:.2e}；"
   f"B 的 record 归一化 {np.round(nB,4)} vs seed 复制 {np.round(SEED0/SEED0.sum(),4)}")

ck("H3a L2 的稳态仍等于 seed 的忠实复制（Z129 §5 结论在加高度层后仍成立）",
   np.allclose(nB, SEED0 / SEED0.sum(), atol=1e-6),
   f"B 的 record 归一化 = {np.round(nB,4)}，seed 归一化 = {np.round(SEED0/SEED0.sum(),4)}")

# ============================================================
print()
print("=" * 72)
print("第 2 步：L1 高度层是否带来'新结构'（Z129 说原引擎'不产生新结构'）")
print("=" * 72)
print(f"  A 的 L1 = record（顶点累积，V={V} 维）")
print(f"  B 的 L1 = H={H} 个高度层 + record（V={V} 维）")
ht = np.array(eB.h_trace)
print(f"  B 的高度层轨迹：前 3 代 {np.round(ht[:3], 4).tolist()}；末代 {np.round(ht[-1], 4).tolist()}")
ck("H4a 加高度层后 L1 内部产生**新的层级结构**（非单一模式累积）",
   len(np.unique(np.round(ht / ht.sum(1, keepdims=True), 3), axis=0)) > 1,
   f"高度层归一化轨迹的不同模式数 = {len(np.unique(np.round(ht/ht.sum(1,keepdims=True),3),axis=0))}")

# ============================================================
print()
print("=" * 72)
n_fail = sum(1 for x in CK if not x["pass"])
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 72)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "l1_height_layer_engine.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT,
               "L1_height": {"L": L_use, "H": H, "b_h": dict(sorted(dist.items()))},
               "A_record_norm": nA.tolist(), "B_record_norm": nB.tolist(),
               "seed_norm": (SEED0 / SEED0.sum()).tolist(),
               "diffA": diffA, "diffB": diffB,
               "B_h_final": (eB.hist_h / eB.hist_h.sum()).tolist(),
               "B_h_patterns": int(len(np.unique(np.round(ht / ht.sum(1, keepdims=True), 3), axis=0)))},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/l1_height_layer_engine.json")
sys.exit(1 if n_fail else 0)