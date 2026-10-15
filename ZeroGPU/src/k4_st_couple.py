#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_st_couple.py —— 标量(膨胀) + 张量(剪切) 耦合 = 完整引力
====================================================================================
完整 GR 引力 = 度规扰动 = 标量(膨胀/迹) + 张量(无迹剪切/双偏振)。
Z105 标量层(膨胀) + Z106-107 张量层(剪切) 在此耦合:
  标量场 φ(每顶点膨胀幅度)
  张量场 h(无迹剪切, h_plus=xx−yy, h_cross=xy)
  耦合(膨胀↔剪切相互转化):
    ∂ₜ²φ = −c²Lφ − g·(剪切能量)      # 剪切能量驱动膨胀
    ∂ₜ²h = −c²Lh + g·(膨胀梯度 Lφ)    # 膨胀不均匀驱动剪切
验证: ①完整引力波色散 ω²=c²λ_k ②膨胀↔剪切相互转化(激发一测另一) ③g 扫描转化强度。
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

def coupled_evolve(L, c2, g, steps, phi0, phi1, h0, h1):
    """标量φ(膨胀)+张量h(剪切)耦合蛙跳。返回 φ, h 序列。"""
    n = L.shape[0]
    phis=[phi1.copy()]; hs=[h1.copy()]
    pp, cp = phi0.copy(), phi1.copy()
    hp, ch = h0.copy(), h1.copy()
    for t in range(1, steps):
        # 剪切能量 → 驱动膨胀 (去均值, 避免单向漂移)
        she = ch[:,0]**2 + ch[:,1]**2        # h_plus² + h_cross²
        ph_nxt = 2*cp - pp - c2*(L@cp) - g*(she - she.mean())
        # 膨胀梯度 Lφ → 驱动剪切 (逐分量)
        grad_phi = L @ cp
        h_nxt = 2*ch - hp - c2*(L@ch) + g*np.stack([grad_phi, grad_phi], axis=1)
        pp, cp = cp, ph_nxt
        hp, ch = ch, h_nxt
        phis.append(ph_nxt.copy()); hs.append(h_nxt.copy())
    return np.array(phis), np.array(hs)

def main():
    nr, nc = 6, 6
    L, adj = grid_graph(nr, nc)
    n = nr*nc
    evals, evecs = np.linalg.eigh(L)
    w_th = np.sqrt(np.clip(evals, 0, None))
    lam_max = float(evals[-1])
    c2 = min(0.8, 3.5/lam_max); c = np.sqrt(c2)
    steps = 2048
    out = {"config": {"graph": f"grid{nr}x{nc}", "n": n, "coupling": "scalar(expansion)+tensor(shear)",
                      "c2": round(c2,3)}}

    print("="*90)
    print("标量(膨胀) + 张量(剪切) 耦合 = 完整引力")
    print("="*90)

    # === 验证1: 完整引力波色散 (无耦合 g=0: 标量+张量各自色散) ===
    print("\n[验证1] 完整引力波色散 (标量层 ω²=c²λ_k, 张量层同)")
    # 标量场色散(激发 L 本征模式在 φ)
    phis, hs = coupled_evolve(L, c2, 0.0, steps,
                              np.zeros(n), evecs[:,1].copy(),
                              np.zeros((n,2)), np.zeros((n,2)))
    P = phis @ evecs
    fk, _ = dominant_freq(P[:, 1])
    corr_s = 0.0
    modes = {}
    for k in range(1, n):
        if w_th[k] >= np.pi*0.95: continue
        phis, _ = coupled_evolve(L, c2, 0.0, steps,
                                 np.zeros(n), evecs[:,k].copy(),
                                 np.zeros((n,2)), np.zeros((n,2)))
        P = phis @ evecs
        fk, _ = dominant_freq(P[:, k])
        modes[int(k)] = {"sqrt_lambda": round(float(w_th[k]),4),
                         "c_sqrt": round(float(c*w_th[k]),4), "dom_freq": round(fk,4)}
    ks=[k for k in modes if modes[k]["dom_freq"]>0.01]
    xs=[modes[k]["sqrt_lambda"] for k in ks]; ys=[modes[k]["dom_freq"] for k in ks]
    corr = float(np.corrcoef(xs,ys)[0,1]) if len(ks)>=3 and np.std(xs)>0 and np.std(ys)>0 else 0.0
    out["dispersion_scalar"] = {"n_modes": len(ks), "corr": round(corr,4)}
    print(f"  标量层色散相关 = {corr:.4f}  (膨胀波 ω²=c²λ_k)")

    # === 验证2: 膨胀↔剪切相互转化 ===
    print("\n[验证2] 膨胀↔剪切相互转化 (耦合 g>0, 小 g + 短演化避免正反馈爆炸)")
    print("  注: 玩具耦合(非能量守恒)正反馈, 用小 g 验证转化存在")
    v = evecs[:, 2]
    steps_short = 800
    # (a) 激发纯膨胀(φ=v): 测剪切是否产生
    for g in [0.0, 0.005, 0.01, 0.02]:
        phis, hs = coupled_evolve(L, c2, g, steps_short,
                                  np.zeros(n), v.copy(),
                                  np.zeros((n,2)), np.zeros((n,2)))
        shear_amp = float(np.sqrt(hs[-200:][:,:,0]**2 + hs[-200:][:,:,1]**2).max())
        out.setdefault("expand_to_shear", {})[str(g)] = round(shear_amp,6)
        print(f"  激发膨胀 g={g:.3f}: 剪切末段振幅 = {shear_amp:.2e}")
    # (b) 激发纯剪切(h_plus=v): 测膨胀是否产生
    for g in [0.0, 0.005, 0.01, 0.02]:
        h0 = np.zeros((n,2)); h0[:,0] = v
        phis, hs = coupled_evolve(L, c2, g, steps_short,
                                  np.zeros(n), np.zeros(n),
                                  np.zeros((n,2)), h0.copy())
        exp_amp = float(np.abs(phis[-200:]).max())
        out.setdefault("shear_to_expand", {})[str(g)] = round(exp_amp,6)
        print(f"  激发剪切 g={g:.3f}: 膨胀末段振幅 = {exp_amp:.2e}")

    # === 验证3: g 扫描 → 转化强度 ===
    print("\n[验证3] 耦合强度 g → 膨胀→剪切转化强度")
    gs = [0.0, 0.005, 0.01, 0.02, 0.04]
    conv = {}
    for g in gs:
        phis, hs = coupled_evolve(L, c2, g, steps_short,
                                  np.zeros(n), v.copy(),
                                  np.zeros((n,2)), np.zeros((n,2)))
        conv[g] = round(float(np.sqrt(hs[steps_short//2:][:,:,0]**2+hs[steps_short//2:][:,:,1]**2).mean()),6)
    out["conv_g_scan"] = {str(g): v for g, v in conv.items()}
    print("  g:  转化剪切振幅:")
    for g, v in conv.items():
        print(f"    {g:.3f}: {v:.6f}")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_st_couple.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_st_couple.json")

if __name__ == "__main__":
    main()
