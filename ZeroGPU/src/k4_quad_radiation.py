#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_quad_radiation.py —— 四极辐射: 守恒耦合时变源 → P∝(d³I/dt³)²
====================================================================================
GR 核心预言: 时变质量四极矩辐射引力波, 辐射功率 P∝(d³I/dt³)²。
Zero 守恒耦合引擎(Z109) + 时变四极矩源(双质量块振荡):
  源: 两质量块(顶点v1,v2)振幅 A 振荡 sin(ω_s t) → 时变四极矩 I(t)
  张量波: ∂ₜ²h=−c²Lh−g·Lφ+源  (守恒耦合+时变源)
  辐射功率: P=Σ_k ω_k²⟨|a_k|²⟩ (本征基投影, 稳态段)
  四极矩:  I(t)=Σ_v ρ_v(t)·x_v² ; d³I/dt³=三阶导数
验证: P ∝ (d³I/dt³)² —— 扫描源振幅 A, 拟合 logP vs log(d³I), slope=2
"""
from __future__ import annotations
import sys, os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def grid_graph(nr, nc):
    adj = {v: [] for v in range(nr*nc)}
    coords = {}
    for i in range(nr):
        for j in range(nc):
            v = i*nc + j; coords[v] = (i, j)
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
    return L, adj, coords

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

def evolve_src(L, c2, g, m2, steps, src_fn, phi0, phi1, h0, h1):
    """守恒耦合 + 时变源(驱动张量). 返回 φ,h 历史(采样), 能量历史."""
    n = L.shape[0]
    pp, cp = phi0.copy(), phi1.copy()
    hp_, ch = h0.copy(), h1.copy()
    E_hist=[]; hsave=[]
    for t in range(1, steps):
        S = src_fn(t)                       # (n,2) 时变源(剪切)
        Lcp = L @ cp
        Lch0 = L @ ch[:,0]; Lch1 = L @ ch[:,1]
        pn = 2*cp - pp - c2*Lcp - m2*cp - g*(Lch0 + Lch1)
        hn = 2*ch - hp_ - c2*(L @ ch) - g*np.stack([Lcp, Lcp], axis=1) + S
        E_hist.append(0.5*np.sum((pn-cp)**2) + 0.5*np.sum((hn-ch)**2)
                     + 0.5*c2*(pn@L@pn + np.sum(hn*(L@hn)))
                     + 0.5*m2*np.sum(pn**2) + g*np.sum((L@pn)*np.sum(hn,axis=1)))
        if t % 8 == 0: hsave.append(ch.copy())
        pp, cp = cp, pn; hp_, ch = ch, hn
    return np.array(E_hist), np.array(hsave)

def rad_power(hsave, L, evecs):
    """辐射功率: P=Σ ω_k²⟨|a_k|²⟩ 稳态段(后1/3). 投影 h_plus 到本征基."""
    n = L.shape[0]
    evals = np.diag(evecs.T @ L @ evecs); evals = np.clip(evals, 0, None)
    Hp = hsave[:, :, 0]                     # h_plus 历史 (T, n)
    A = Hp @ evecs
    T = len(hsave); seg = slice(int(T*2/3), T)
    power = 0.0
    for k in range(n):
        a = A[seg, k]
        w2 = evals[k]
        power += w2 * np.mean(a**2)
    return float(power)

def quadrupole_third(coords, rho_history, vset):
    """四极矩三阶导 |d³I/dt³|: I(t)=Σ_v ρ_v(t)·(x_v²+y_v²), 数值三阶导. rho_history:(T,n)"""
    T = rho_history.shape[0]
    w = np.array([coords[v][0]**2 + coords[v][1]**2 for v in vset])
    I = np.array([np.dot(rho_history[t][vset], w) for t in range(T)])
    return float(np.max(np.abs(np.diff(I, n=3, prepend=[0,0,0]))))

def main():
    nr, nc = 8, 8
    L, adj, coords = grid_graph(nr, nc)
    n = nr*nc
    evals, evecs = np.linalg.eigh(L)
    lam_max = float(evals[-1])
    c2 = min(0.8, 1.2/lam_max); c = np.sqrt(c2)
    m2 = 0.3
    g = 0.1
    # 双质量块(四极配置): 对称两顶点
    v1, v2 = 0, (nr-1)*nc + (nc-1)          # 对角
    ws = 0.5                                 # 源振荡频率(须低于色散, 落在图模式)
    steps = 6000
    out = {"config": {"graph": f"grid{nr}x{nc}", "n": n, "source": "dual-mass quadrupole",
                      "ws": ws, "g": g, "c2": round(c2,3), "m2": m2}}

    print("="*90)
    print("四极辐射: 守恒耦合时变源 → P∝(d³I/dt³)²")
    print("="*90)

    # 对每个源振幅 A, 测辐射功率 P 和四极矩三阶导 d³I
    amps = [0.1, 0.3, 0.5, 1.0, 1.5, 2.0]
    P_list=[]; dI_list=[]
    rho_static = np.zeros(n); rho_static[v1]=1.0; rho_static[v2]=1.0
    for A in amps:
        rho_hist = np.stack([rho_static*(A*(1+0.8*np.sin(ws*t))) for t in range(0,steps,8)])  # 振荡∝A
        def src_fn(t):
            S = np.zeros((n,2))
            S[v1,0] = A*np.sin(ws*t)
            S[v2,0] = A*np.sin(ws*t)
            return S
        phi0=np.zeros(n); phi1=np.zeros(n); h0=np.zeros((n,2)); h1=np.zeros((n,2))
        E, hsave = evolve_src(L, c2, g, m2, steps, src_fn, phi0, phi1, h0, h1)
        P = rad_power(hsave, L, evecs)
        dI = quadrupole_third(coords, rho_hist, [v1, v2])
        P_list.append(P); dI_list.append(dI)
        print(f"  A={A:.1f}: P={P:.2e}  |d³I/dt³|={dI:.2e}")
    # 拟合 logP vs log(d³I): slope 应为 2 (P∝(d³I)²)
    x = np.log(np.array(dI_list)); y = np.log(np.array(P_list))
    slope = float(np.polyfit(x, y, 1)[0])
    corr = float(np.corrcoef(x, y)[0,1])
    out["radiation"] = {"amps": amps, "P": [round(p,6) for p in P_list],
                        "d3I": [round(d,6) for d in dI_list],
                        "slope_logP_logdI": round(slope,3), "corr": round(corr,3)}
    print(f"\n  logP vs log(d³I) slope = {slope:.3f}  (GR 四极: 2.0)  corr={corr:.3f}")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_quad_radiation.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_quad_radiation.json")

if __name__ == "__main__":
    main()
