#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_eigen_wave.py —— Zero 内生验证：时空动态投影到 L 本征基, 测 ω²∝λ_k 色散
====================================================================================
Zero 的 GR 理论里引力波答案 = L 的谱 λ_k 是本征频率平方。
不外部加波动方程, 纯用 Zero 已有结构: 更大图(path16 非简并谱) + 双向耦合(fb)
+ 自参考(s) 产生时空动态 → 投影到 L 本征基 φ=Σa_k φ_k → 测各模式主导频率 f_k
→ 对比 f_k vs √λ_k (ω²∝λ_k ⇒ 色散 ⇒ Zero 内生时空波)。
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

def robust_solve(L, x):
    x0 = x - x.mean()
    evals, evecs = np.linalg.eigh(L)
    mask = evals > 1e-8
    U = evecs[:, mask]; lam = evals[mask]
    phi = -(U @ (U.T @ x0 / lam))
    return phi - phi.mean()

def run(n, steps, fb, s, seed0, T=4, retention=2, warmup=100):
    L, adj = path_laplacian(n)
    evals, evecs = np.linalg.eigh(L)          # λ_k, φ_k (本征基)
    record = np.zeros(n)
    seed = np.array(seed0, dtype=float)
    phis = []
    for t in range(warmup + steps):
        W = np.zeros(n)
        for i in range(n):
            W[i] = sum(seed[j] for j in adj[i])
        record = record + W
        s0 = record.sum()
        if s0 > 0: record = record / s0 * n
        phi = robust_solve(L, record - record.mean())
        # 自参考: 时空自身能量(曲率)参与源 → 放大拉伸
        if s > 0:
            kappa = -L @ phi; sel = kappa**2
            sel = sel / sel.sum() if sel.sum() > 0 else np.ones(n)/n
            xs = (1-s) * (record - record.mean()) + s * (sel - sel.mean())
            phi = robust_solve(L, xs)
        # 双向耦合: seed = (1-fb)*保留 + fb*势反馈(归一化量级, 避免爆炸)
        idx = np.argsort(seed)[::-1][:retention]
        keep = np.zeros(n); keep[idx] = record[idx]
        fb_feed = phi - phi.min() + 1e-9
        kmean = keep.mean()
        fmean = fb_feed.mean()
        if fmean > 1e-12:
            fb_feed = fb_feed * (kmean / fmean)   # 势反馈量级匹配保留 seed
        seed = (1-fb) * keep + fb * fb_feed
        if (t+1) % T == 0:
            keep = np.zeros(n); seed = (1-fb)*keep + fb*fb_feed
        if t >= warmup:
            phis.append(phi.copy())
    return np.array(phis), evals, evecs, L

def dominant_freq(series, dt=1.0):
    """去均值 + 谱分析, 返回主导(角)频率与相对功率。"""
    x = series - series.mean()
    N = len(x)
    w = np.hanning(N); x = x * w
    X = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(N, d=dt) * 2 * np.pi   # 角频率
    power = np.abs(X)**2
    # 忽略 DC(索引0) 与极低频噪声(带通: 只查 0.005~Nyquist)
    lo = int(0.005 / (2*np.pi*dt) * N) if (0.005/(2*np.pi*dt))*N < N//2 else 1
    if lo < 1: lo = 1
    seg = power[lo:]
    if seg.sum() <= 0: return 0.0, 0.0
    pk = np.argmax(seg) + lo
    return float(freqs[pk]), float(power[pk]/power.sum())

def main():
    n = 16
    L, adj = path_laplacian(n)
    evals, evecs = np.linalg.eigh(L)
    # 理论本征角频率 ω=√λ (图上波动色散)
    w_theory = np.sqrt(np.clip(evals, 0, None))
    out = {"config": {"graph": "path16", "n": n, "solve": "eig projection",
                      "method": "space-time dynamic → project onto L eigenbasis → mode dominant freq vs sqrt(lambda)"}}

    print("="*90)
    print("Zero 内生验证：时空动态投影到 L 本征基, 测 ω²∝λ_k 色散")
    print("="*90)

    out["scans"] = {}
    # 扫描: 双向耦合(fb) × 自参考(s) → 时空动态 → 本征模式激发
    for fb, s in [(0.8, 0.0), (0.8, 0.6), (0.6, 0.6), (0.4, 0.0)]:
        phis, ev, evc, _ = run(n, steps=4096, fb=fb, s=s,
                               seed0=list(np.linspace(8.0, 0.2, n)))
        phis = np.array(phis)
        dyn = float(np.std(phis[-500:], axis=0).max())
        # 投影到本征基: A[:,k] = φ_k·φ(t)
        A = phis @ evc                      # (steps, n) 每列 = a_k(t)
        key = f"fb{fb}_s{s}"
        out["scans"][key] = {"phi_dyn": round(dyn,6)}
        print(f"\n--- fb={fb:.1f} s={s:.1f}: 时空动态={dyn:.6f} ---")
        print("   模式k  λ_k    √λ_k(理论)  实测主导频率  功率占比")
        mode_freqs = {}
        for k in range(n):
            fk, pw = dominant_freq(A[:, k])
            mode_freqs[int(k)] = {"lambda": round(float(ev[k]),4),
                                  "sqrt_lambda": round(float(w_theory[k]),4),
                                  "dom_freq": round(fk,4), "power": round(pw,4)}
            if k < 8:   # 只打印前8模式(高频>Nyquist欠采样)
                print(f"   {k:3d}   {ev[k]:6.3f}   {w_theory[k]:7.3f}   {fk:9.4f}   {pw:.3f}")
        out["scans"][key]["modes"] = mode_freqs

    # 色散判定: 对低频模式(k≤3, √λ<0.5 可测), 实测主导频率是否≈√λ
    out["verdict_note"] = ("模流耦合时空动态若激发本征模式, 低频模式(k≤3)主导频率应≈√λ_k; "
                           "若全≈模流慢频(0.001~0.01)则模流与L谱脱节, 无本征色散")
    print("\n" + "="*90)
    print("色散判定: 低频模式(√λ<0.5)主导频率是否≈√λ_k")
    print("="*90)
    print(out["verdict_note"])

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_eigen_wave.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_eigen_wave.json")

if __name__ == "__main__":
    main()
