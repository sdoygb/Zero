#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_L1_wave.py —— 势的历史入 L1：L 驱动势自身演化, 验证 ω²∝λ_k 色散 + 波前
====================================================================================
用户洞察: 势的自身演化(∂ₜ²φ 二阶时间惯性)藏在 L1 历史记录层。
实现: L1 记录势过去 φ(t−1)(惯性) → 势演化 φ(t+1)=2φ(t)−φ(t−1)+c²(Lφ(t)+x(t)) 蛙跳。
验证1(核心): 自由波 + 局域脉冲初值 → 波前传播 + 投影 L 本征基测 ω²∝λ_k 色散。
验证2(Zero整合): L1 势历史 + 模流时变源 → 测本征模式激发。
"""
from __future__ import annotations
import sys, os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def path_laplacian(n):
    adj = [[] for _ in range(n)]
    for i in range(n):
        if i > 0: adj[i].append(i-1)
        if i < n-1: adj[i].append(i+1)
    L = np.zeros((n, n))
    for i in range(n):
        L[i, i] = len(adj[i])
        for j in adj[i]:
            L[i, j] = -1
    return L, adj

def leapfrog(L, c2, steps, phi0, phi1, x_src=None):
    """∂ₜ²φ=−c²Lφ+x 蛙跳(物理波动: Δ负定=−L, 双曲稳定)。L1=prev(势历史/惯性)。
       x_src: (steps, n) 或 None。"""
    n = L.shape[0]
    phis = [phi1.copy()]
    prev, cur = phi0.copy(), phi1.copy()
    for t in range(1, steps):
        src = np.zeros(n) if x_src is None else x_src[t]
        nxt = 2*cur - prev - c2*(L @ cur) + c2*src
        prev, cur = cur, nxt
        phis.append(cur.copy())
    return np.array(phis)   # (steps, n)

def moduli_source(n, steps, seed0, T=4, retention=2):
    """L2 模流时变源 x(t): seed 演化 → record(物质历史) → 时变源。"""
    _, adj = path_laplacian(n)
    record = np.zeros(n); seed = np.array(seed0, dtype=float)
    srcs = []
    for _ in range(steps):
        W = np.zeros(n)
        for i in range(n):
            W[i] = sum(seed[j] for j in adj[i])
        record = record + W
        s0 = record.sum()
        if s0 > 0: record = record / s0 * n
        idx = np.argsort(seed)[::-1][:retention]
        seed = np.zeros(n); seed[idx] = record[idx]
        if (0) % T == 0: pass
        srcs.append(record - record.mean())
    return np.array(srcs)

def dominant_freq(series, dt=1.0):
    x = series - series.mean()
    N = len(x)
    w = np.hanning(N); x = x * w
    X = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(N, d=dt) * 2 * np.pi
    power = np.abs(X)**2
    lo = max(1, int(0.001/(2*np.pi*dt)*N))
    seg = power[lo:]
    if seg.sum() <= 0: return 0.0, 0.0
    pk = np.argmax(seg) + lo
    return float(freqs[pk]), float(power[pk]/power.sum())

def main():
    n = 16
    L, adj = path_laplacian(n)
    evals, evecs = np.linalg.eigh(L)
    w_theory = np.sqrt(np.clip(evals, 0, None))
    c2 = 0.8   # 稳定: c²·λ_max(≈3.92)<4
    out = {"config": {"graph": "path16", "n": n, "c2": c2,
                      "method": "L1 记录势历史 φ(t-1) → ∂ₜ²φ=c²Lφ+x 蛙跳 → ω²∝λ_k 色散"}}

    print("="*90)
    print("势的历史入 L1：L 驱动势自身演化, 验证色散 + 波前")
    print("="*90)

    # === 验证1: 自由波 + 局域脉冲初值 → 波前传播 + 色散 ===
    print("\n【验证1】自由波 ∂ₜ²φ=c²Lφ, 顶点0 初值脉冲")
    steps = 2048
    phi0 = np.zeros(n); phi1 = np.zeros(n); phi1[0] = 1.0
    phis = leapfrog(L, c2, steps, phi0, phi1)
    # 波前: 各顶点|φ|峰值时间(应随距离延迟)
    peak_t = {}
    for v in range(n):
        a = np.abs(phis[:, v])
        peak_t[v] = int(np.argmax(a[:steps//2]))   # 只查前半,避免边界反射
    out["test1_wavefront_peak_time"] = {int(v): int(peak_t[v]) for v in range(n)}
    print("  顶点峰值时间(代):", {v: peak_t[v] for v in range(n)})
    # 波前速度: 每代传播几个顶点
    # 色散: 投影 L 本征基 a_k(t), 测各模式主导频率 vs √λ_k
    A = phis @ evecs
    print("  模式k  √λ_k(理论)  实测主导频率  (ω≈√λ ⇒ 色散)")
    test1_modes = {}
    for k in range(1, n):
        fk, pw = dominant_freq(A[:, k])
        test1_modes[int(k)] = {"sqrt_lambda": round(float(w_theory[k]),4),
                               "dom_freq": round(fk,4), "power": round(pw,4)}
        if k <= 6:
            print(f"    {k:3d}   {w_theory[k]:7.3f}   {fk:9.4f}")
    out["test1_modes"] = test1_modes
    # 色散相关: 低频模式(k=1..6, √λ<Nyquist可测) 实测频率 vs √λ
    ks = [k for k in range(1, n) if w_theory[k] < np.pi]   # 可测模式
    xs = [w_theory[k] for k in ks]; ys = [test1_modes[k]["dom_freq"] for k in ks]
    if len(ks) >= 2 and np.std(xs) > 0:
        corr = float(np.corrcoef(xs, ys)[0, 1])
    else:
        corr = 0.0
    out["test1_dispersion_corr"] = round(corr, 4)
    print(f"  色散相关 ω vs √λ (k=1..{ks[-1]}) = {corr:.4f}  (≈1 ⇒ ω²∝λ 色散成立)")

    # === 验证2: L1 势历史 + 模流时变源 → 激发 ===
    print("\n【验证2】L1 势历史(惯性) + 模流时变源 → ∂ₜ²φ=c²Lφ+x(t)")
    steps2 = 2048
    srcs = moduli_source(n, steps2, list(np.linspace(8.0, 0.2, n)))
    phi0 = np.zeros(n); phi1 = np.zeros(n)
    phis2 = leapfrog(L, c2, steps2, phi0, phi1, x_src=srcs)
    A2 = phis2 @ evecs
    print("  模式k  √λ_k(理论)  实测主导频率  (模流源是否激发本征)")
    test2_modes = {}
    for k in range(1, n):
        fk, pw = dominant_freq(A2[:, k])
        test2_modes[int(k)] = {"sqrt_lambda": round(float(w_theory[k]),4),
                               "dom_freq": round(fk,4), "power": round(pw,4)}
        if k <= 6:
            print(f"    {k:3d}   {w_theory[k]:7.3f}   {fk:9.4f}")
    out["test2_modes"] = test2_modes
    xs2 = [w_theory[k] for k in ks]; ys2 = [test2_modes[k]["dom_freq"] for k in ks]
    corr2 = float(np.corrcoef(xs2, ys2)[0, 1]) if len(ks) >= 2 and np.std(xs2) > 0 else 0.0
    out["test2_dispersion_corr"] = round(corr2, 4)
    print(f"  色散相关 ω vs √λ = {corr2:.4f}")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_L1_wave.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_L1_wave.json")

if __name__ == "__main__":
    main()
