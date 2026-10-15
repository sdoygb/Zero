# -*- coding: utf-8 -*-
"""
ZeroGPU · 批量并行 MPS 引擎（真 GPU 落地）

把 NumPy 向量化升到 PyTorch MPS（Metal）。GPU 的价值在**吞吐**：
把 B 个独立配置放进 (B,V,V) 批量张量同时演化（torch.bmm）。

  ① 单配置（B=1）用于与 NumPy 版一致性验证；
  ② 批量（B 大）展示 MPS 并行加速 —— 为大规模多配置扫描铺路。

与 ZeroCoreVec 语义等价（float32 + 归一化，比例严格保持）。

用法：
    /usr/local/bin/python3 z0_mps.py --verify      # 与 NumPy 版一致
    /usr/local/bin/python3 z0_mps.py --bench       # MPS vs CPU 测速
"""
from __future__ import annotations
import sys, os, time, json

import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import path_adj, torus_adj, star_adj
from z0_vectorized import ZeroCoreVec, QW_MIN, QW_MAX

DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"
NORM_THRESHOLD = 1e30   # float32 上限 ~3.4e38，必须远低于它触发归一化防溢出


def _adj_tensor(adj_batch, device):
    B = len(adj_batch); V = len(adj_batch[0])
    A = torch.zeros(B, V, V, device=device)
    for b, adj in enumerate(adj_batch):
        for u, nb in enumerate(adj):
            for w in nb:
                A[b, u, w] = 1.0
    return A


class ZeroBatchMPS:
    """批量并行 MPS 引擎（B 个配置同时演化）。"""

    def __init__(self, adj_batch, seeds_batch, L=4, retention="last",
                 device=DEVICE):
        self.B = len(adj_batch)
        self.V = len(adj_batch[0])
        self.A = _adj_tensor(adj_batch, device)
        self.L = int(L)
        self.retention = retention
        self.device = device
        self.seed = torch.tensor(seeds_batch, dtype=torch.float32, device=device)
        self.W = torch.zeros(self.B, self.V, self.V, device=device)
        self.record = torch.zeros(self.B, self.V, device=device)
        self.term = torch.zeros(self.B, self.V, device=device)
        self.running = torch.zeros(self.B, self.V, device=device)
        self.layers = []
        self.gen = 0
        self.gen_closures = torch.zeros(self.B, self.V, device=device)
        self.idx = torch.tensor(list(range(self.V)), dtype=torch.long, device=device)

    def run_generation(self):
        W = self.W.clone()
        W[:, self.idx, self.idx] += self.seed
        for _ in range(self.L):
            W = torch.bmm(W, self.A)
        diag = W[:, self.idx, self.idx]
        self.gen_closures = diag.clone()
        self.record = self.record + diag
        W2 = W.clone()
        W2[:, self.idx, self.idx] = 0
        self.W = W2
        self.term = self.term + (W2.sum(dim=1) - W2.sum(dim=2))
        self._reseed()
        self._normalize()
        self.gen += 1
        return self

    def _retained(self):
        r = self.retention
        if r == "last":
            return self.gen_closures.clone()
        if r == "all":
            self.running = self.running + self.gen_closures
            return self.running.clone()
        if r == "top2":
            if len(self.layers) >= 2:
                return self.layers[-2] + self.gen_closures
            return self.gen_closures.clone()
        raise ValueError(f"未知保留规则 {r}")

    def _reseed(self):
        w = self._retained()
        zero = w.sum(dim=1) == 0
        if zero.any():
            w[zero] = 1.0
        self.layers.append(self.gen_closures.clone())
        self.seed = w

    def _normalize(self):
        for a in (self.record, self.W, self.term, self.running, self.seed):
            t = a.sum()
            if float(t) > NORM_THRESHOLD:
                a /= t

    def q_weighted(self, arr):
        tot = arr.sum(dim=1, keepdim=True)
        w = arr / tot
        return (w * w).sum(dim=1)

    def q_hist(self):
        return self.q_weighted(self.record)

    def q_live(self):
        return self.q_weighted(self.gen_closures)


# ────────────────────────── 一致性验证 ──────────────────────────
def verify(N=200):
    ok = [0]; bad = [0]
    def ck(name, cond, detail=""):
        (ok if cond else bad)[0] += 1
        print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    print("=" * 76); print(f"MPS 一致性：path16 geom last L4 — {N} 代 vs NumPy"); print("=" * 76)
    adj = path_adj(16); V = len(adj)
    seed = [float(2 ** (V - 1 - v)) for v in range(V)]

    ref = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
    mps = ZeroBatchMPS([adj], [seed], L=4, retention="last")
    max_dq = 0.0
    for g in range(N):
        ref.run_generation()
        mps.run_generation()
        max_dq = max(max_dq, abs(float(ref.q_hist()) - float(mps.q_hist()[0])))
    ck("q_hist 与 NumPy 一致（float32 精度）", max_dq < 1e-5, f"最大差 {max_dq:.2e}")
    ck("MPS 在 mps device", str(mps.device) == "mps", str(mps.device))
    ck("q_hist 进窗", QW_MIN < float(mps.q_hist()[0]) < QW_MAX,
       f"q={float(mps.q_hist()[0]):.6f}")
    print(); print(f"断言：通过 {ok[0]} / 不符 {bad[0]}")
    return 0 if bad[0] == 0 else 1


# ────────────────────────── 测速 ──────────────────────────
def bench(gens=3000):
    print("=" * 76); print("MPS vs CPU NumPy 测速"); print("=" * 76)
    adj = path_adj(16); V = len(adj)
    seed = [float(2 ** (V - 1 - v)) for v in range(V)]

    # CPU NumPy 单配置
    t0 = time.time()
    ref = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
    for _ in range(gens): ref.run_generation()
    cpu_t = time.time() - t0
    print(f"CPU NumPy 单配置: {gens} 代 {cpu_t:.2f}s  {gens/cpu_t:.0f} 代/秒")

    # MPS 单配置 B=1
    mps1 = ZeroBatchMPS([adj], [seed], L=4, retention="last")
    t0 = time.time()
    for _ in range(gens): mps1.run_generation()
    mps1_t = time.time() - t0
    print(f"MPS 单配置 B=1:   {gens} 代 {mps1_t:.2f}s  {gens/mps1_t:.0f} 代/秒")

    # MPS 批量 B=64（多配置并行）
    for B in (16, 64):
        adjb = [adj] * B
        seedb = [seed] * B
        mb = ZeroBatchMPS(adjb, seedb, L=4, retention="last")
        t0 = time.time()
        for _ in range(gens): mb.run_generation()
        bt = time.time() - t0
        print(f"MPS 批量 B={B}:   {gens} 代 {bt:.2f}s  {B*gens/bt:.0f} 配置·代/秒"
              f"  (每个 {gens/bt:.0f} 代/秒)")

    # CPU NumPy 循环 B=64（对比）
    B = 64
    t0 = time.time()
    for _ in range(B):
        r = ZeroCoreVec(adj, seed=seed, retention="last", L=4)
        for _g in range(gens): r.run_generation()
    cpuB_t = time.time() - t0
    print(f"CPU NumPy ×{B}:   {gens} 代 ×{B} {cpuB_t:.2f}s  {B*gens/cpuB_t:.0f} 配置·代/秒")

    print("=" * 76)


if __name__ == "__main__":
    if "--verify" in sys.argv:
        sys.exit(verify())
    if "--bench" in sys.argv:
        bench()
        sys.exit(0)
    sys.exit(verify())
