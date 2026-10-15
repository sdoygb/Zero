#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_grid_wave.py —— 更大图色散普适性(逐模式激发) + 模流源辐射功率
====================================================================================
验证1(普适性): 对每个可测本征模式 k, 初值=该本征向量 φ_k 单独激发, 蛙跳自由波,
                 测该模式振荡频率 f_k vs c√λ_k (色散 ω²=c²λ_k 逐模式严格)。
验证2(辐射功率): 模流时变源 x(t) 驱动 ∂ₜ²φ=−c²Lφ+x(t), 稳态波振幅 → 辐射功率
                  P=Σω_k²⟨|a_k|²⟩, 扫描源强度 scale → 幂律(标量辐射功率类比)。
"""
from __future__ import annotations
import sys, os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def grid_graph(nr, nc):
    adj = {v: [] for v in range(nr*nc)}
    for i in range(nr):
        for j in range(nc):
            v = i*nc + j
            for di, dj in [(1,0),(-1,0),(0,1),(0,-1)]:
                ni, nj = i+di, j+dj
                if 0 <= ni < nr and 0 <= nj < nc:
                    adj[v].append(ni*nc + nj)
    n = nr*nc
    L = np.zeros((n, n))
    for v in range(n):
        L[v, v] = len(adj[v])
        for j in adj[v]:
            L[v, j] = -1
    return L, adj

def leapfrog(L, c2, steps, phi0, phi1, x_src=None):
    n = L.shape[0]
    phis = [phi1.copy()]
    prev, cur = phi0.copy(), phi1.copy()
    for t in range(1, steps):
        src = np.zeros(n) if x_src is None else x_src[t]
        nxt = 2*cur - prev - c2*(L @ cur) + c2*src
        prev, cur = cur, nxt
        phis.append(cur.copy())
    return np.array(phis)

def dominant_freq(series, dt=1.0):
    x = series - series.mean()
    N = len(x)
    w = np.hanning(N); x = x * w
    X = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(N, d=dt) * 2 * np.pi
    power = np.abs(X)**2
    lo = max(1, int(0.001/(2*np.pi*dt)*N))
    seg = power[lo:]
    if seg.sum() <= 1e-30: return 0.0, 0.0
    pk = np.argmax(seg) + lo
    return float(freqs[pk]), float(power[pk]/power.sum())

def moduli_source(n, adj, steps, seed0, T=4, retention=2, scale=1.0):
    record = np.zeros(n); seed = np.array(seed0, dtype=float)
    srcs = []
    for _ in range(steps):
        W = np.zeros(n)
        for v in range(n):
            W[v] = sum(seed[j] for j in adj[v])
        record = record + W
        s0 = record.sum()
        if s0 > 0: record = record / s0 * n
        idx = np.argsort(seed)[::-1][:retention]
        seed = np.zeros(n); seed[idx] = record[idx]
        srcs.append((record - record.mean()) * scale)
    return np.array(srcs)

def main():
    out = {"config": {"method": "per-mode eigen excitation dispersion + moduli-source radiation power"}}

    # === 验证1: 逐模式激发测色散普适性 (grid8x8 与 grid5x20) ===
    for nr, nc in [(8, 8), (5, 20)]:
        L, adj = grid_graph(nr, nc)
        n = nr*nc
        evals, evecs = np.linalg.eigh(L)
        w_th = np.sqrt(np.clip(evals, 0, None))
        lam_max = float(evals[-1])
        c2 = min(0.8, 3.5 / lam_max)    # CFL: c²λ_max<4, 自适应保稳定
        c = np.sqrt(c2)
        steps = 2048
        modes = {}
        for k in range(1, n):
            if w_th[k] >= np.pi*0.95: continue      # 欠采样跳过
            phi0 = np.zeros(n); phi1 = evecs[:, k].copy()   # 单独激发模式k
            phis = leapfrog(L, c2, steps, phi0, phi1)
            A = phis @ evecs
            fk, _ = dominant_freq(A[:, k])          # 模式k振幅振荡频率
            modes[int(k)] = {"sqrt_lambda": round(float(w_th[k]),4),
                             "c_sqrt": round(float(c*w_th[k]),4),
                             "dom_freq": round(fk,4)}
        ks = [k for k in modes if modes[k]["dom_freq"] > 0.01]
        xs = [modes[k]["sqrt_lambda"] for k in ks]
        ys = [modes[k]["dom_freq"] for k in ks]
        corr = float(np.corrcoef(xs, ys)[0,1]) if len(ks) >= 3 and np.std(xs) > 0 and np.std(ys) > 0 else 0.0
        dev = np.mean(np.abs(np.array(ys) - c*np.array(xs)) / np.maximum(c*np.array(xs), 1e-9)) * 100
        key = f"grid{nr}x{nc}"
        out[key] = {"n": n, "n_excited": len(ks), "dispersion_corr": round(corr,4),
                    "mean_dev_pct": round(float(dev),1)}
        print(f"[验证1] grid{nr}x{nc} (n={n}): 激发模式={len(ks)}  色散相关={corr:.4f}  平均偏差={dev:.1f}%")

    # === 验证2: 模流源辐射功率 ===
    print("\n[验证2] 模流时变源 → 辐射功率 P=Σω_k²⟨|a_k|²⟩, 扫描源强度")
    L, adj = grid_graph(8, 8)
    n = 64; evals, evecs = np.linalg.eigh(L)
    w_th = np.sqrt(np.clip(evals, 0, None))
    lam_max = float(evals[-1])
    c2 = min(0.8, 3.5 / lam_max)
    steps = 4096
    rad_powers = {}
    for scale in [0.1, 0.5, 1.0, 2.0, 5.0]:
        srcs = moduli_source(n, adj, steps, list(np.linspace(8.0, 0.01, n)), scale=scale)
        phi0 = np.zeros(n); phi1 = np.zeros(n)
        phis = leapfrog(L, c2, steps, phi0, phi1, x_src=srcs)
        A = phis @ evecs
        tail = A[-steps//3:]
        P = 0.0
        for k in range(1, n):
            P += float(w_th[k]**2 * np.mean(tail[:, k]**2))
        rad_powers[scale] = round(P, 6)
        print(f"  源强度={scale:.1f}: 辐射功率 P={P:.6f}")
    xs = list(rad_powers.keys()); ys = list(rad_powers.values())
    slope = float(np.polyfit(np.log(xs), np.log(np.maximum(ys,1e-30)), 1)[0]) if min(ys) > 0 else 0.0
    out["rad_power"] = {"P_vs_scale": {str(k): v for k, v in rad_powers.items()},
                        "powerlaw_slope": round(slope, 3)}
    print(f"  幂律: P ∝ scale^{slope:.3f}")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_grid_wave.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_grid_wave.json")

if __name__ == "__main__":
    main()
