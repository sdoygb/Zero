#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_scalar_wave.py —— 标量引力波检查(稳健版): 局域源脉冲能否产生跨顶点波前传播
========================================================================
标量引力波 = 标量势扰动传播。Zero 势 φ=−L⁺x 全局瞬时(单步无传播),
传播只能体现于跨代动力学。在 path8 稳态模流顶点0制造局域脉冲,
测各顶点峰值响应时间是否随距离延迟(波前=标量波候选) 或全局同步(无波)。
"""
from __future__ import annotations
import sys, os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def path_laplacian(n):
    adj = [[i-1] for i in range(n)]
    for i in range(n):
        adj[i] = [j for j in (i-1, i+1) if 0 <= j < n]
    L = np.zeros((n, n))
    for i in range(n):
        L[i, i] = len(adj[i])
        for j in adj[i]:
            L[i, j] = -1
    return L, adj

def robust_solve(L, x):
    """Lφ=−x, 零空间投影(零和)。特征分解, 排除零本征。"""
    n = L.shape[0]
    x0 = x - x.mean()
    evals, evecs = np.linalg.eigh(L)
    mask = evals > 1e-8
    lam = evals[mask]
    U = evecs[:, mask]
    proj = U @ (U.T @ x0 / lam)          # U (Uᵀ x / λ)
    phi = -proj
    return phi - phi.mean()

def evolve(n, steps, seed0, T=4, retention=2, pulse=None):
    L, adj = path_laplacian(n)
    record = np.zeros(n)
    seed = np.array(seed0, dtype=float)
    hist = []
    for t in range(steps):
        W = np.zeros(n)
        for i in range(n):
            W[i] = sum(seed[j] for j in adj[i])
        diag = W.copy()
        if pulse is not None and t == pulse[0]:
            diag[pulse[1]] += pulse[2]
        record = record + diag
        # 归一化量级(避免累积溢出)
        s = record.sum()
        if s > 0: record = record / s * n
        idx = np.argsort(seed)[::-1][:retention]
        seed = np.zeros(n)
        seed[idx] = record[idx]
        if (t+1) % T == 0:
            pass  # 毁灭由 seed 清零体现(下次 W 邻居贡献)
        x = record - record.mean()
        phi = robust_solve(L, x)
        hist.append(phi.copy())
    return np.array(hist)

def main():
    n = 8
    L, _ = path_laplacian(n)
    out = {"config": {"graph": "path8", "test": "local pulse vertex0, cross-vertex wavefront check", "solve": "eig projection"}}

    hist0 = evolve(n, 1500, [8.0,4.0,1.5,1.0,0.8,0.5,0.3,0.2])
    base = hist0[-100:].mean(axis=0)
    if not np.all(np.isfinite(base)):
        base = np.zeros(n)
    out["baseline_steady_phi"] = [round(float(v),5) for v in base]
    print("稳态势:", [round(float(v),5) for v in base])

    pulse_step = 1000
    hist = evolve(n, 1500, [8.0,4.0,1.5,1.0,0.8,0.5,0.3,0.2],
                  pulse=(pulse_step, 0, 50.0))
    pre = hist[pulse_step-10].copy()
    dev = hist - pre[None, :]
    # 忽略 NaN/含 inf 的代
    peak_t = np.zeros(n); peak_a = np.zeros(n)
    for v in range(n):
        a = np.abs(dev[:, v])
        a[~np.isfinite(a)] = 0.0
        pk = int(np.argmax(a[pulse_step:]))
        peak_t[v] = pk + pulse_step
        peak_a[v] = float(a.max())
    out["pulse"] = {"step": pulse_step, "vertex": 0, "amp": 50.0}
    out["peak_time"] = {int(v): int(peak_t[v]) for v in range(n)}
    out["peak_amp"] = {int(v): round(peak_a[v],5) for v in range(n)}
    print("\n局域脉冲(顶点0): 各顶点峰值响应")
    for v in range(n):
        print(f"  顶点{v}: 峰值 t={int(peak_t[v])}  幅度={peak_a[v]:.5f}")

    delay = {int(v): int(peak_t[v]-pulse_step) for v in range(n)}
    out["delay"] = delay
    # 传播判据: 远端延迟 > 近端延迟 → 波前; 全等 → 全局同步无波
    delays = [delay[v] for v in range(n)]
    spread = max(delays) - min(delays)
    wavefront = (delay.get(7,0) - delay.get(1,0)) > 1
    verdict = ("wavefront" if wavefront else
               "global_sync_no_wave" if spread <= 1 else "diffusive")
    out["verdict"] = verdict
    print(f"\n传播判据: 延迟 {delay}")
    print(f"  延迟跨度={spread}  远端-近端={delay.get(7,0)-delay.get(1,0)}")
    print(f"  判定: {verdict}")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_scalar_wave.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_scalar_wave.json")

if __name__ == "__main__":
    main()
