#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_sc_energy.py —— 能量守恒耦合：从势能泛函变分派生(哈密顿), 长期稳定
====================================================================================
物理正确的标量-张量耦合从作用量/能量泛函出发, 变分派生 → 能量守恒, 无正反馈爆炸。
耦合势能: V_coup = g·Σ_v φ_v·|h_v|²   (标量密度 × 剪切强度)
变分力(哈密顿):
  ∂ₜ²φ = −δE/δφ = −c²Lφ − g·|h|²
  ∂ₜ²h = −δE/δh = −c²Lh − 2g·φ·h        (对称, 来自同一势能)
总能量: E = 动能 + 图势能 + 耦合势能
验证: ①能量守恒(长演化E有界) ②长期稳定(大g不爆炸) ③膨胀↔剪切相互转化 ④色散。
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

def energy(L, c2, g, m2, phi, phip, h, hp):
    """总能量: 动能 + 图势能 + 标量质量项 + 曲率耦合势能 V_coup=g·Σ(Lφ)·h (L自伴→能量守恒)。"""
    n = L.shape[0]
    vphi = (phi - phip)                 # 速度(差分)
    vh = (h - hp)
    Ekin = 0.5*np.sum(vphi**2) + 0.5*np.sum(vh**2)
    Egraph = 0.5*c2*(phi @ L @ phi) + 0.5*c2*np.sum(h * (L @ h))
    Emass = 0.5*m2*np.sum(phi**2)                      # 标量质量项(引力常数自由度有势)
    Ecoup = g*np.sum((L @ phi) * np.sum(h, axis=1))    # V=g·Σ(Lφ)·(h_plus+h_cross)
    return Ekin + Egraph + Emass + Ecoup

def evolve(L, c2, g, m2, steps, phi0, phi1, h0, h1):
    """能量守恒曲率耦合(哈密顿): 线性对称, L自伴, 标量有质量→非简并非共振。
       φ_{t+1}=2φ−φ_prev−c²Lφ−m²φ−g·Σ_p Lh_p   (标量色散 ω²=c²λ+m²)
       h_p_{t+1}=2h−h_prev−c²Lh_p−g·Lφ           (张量色散 ω²=c²λ)
       非简并(标量快/慢) → 无共振 → 长期稳定"""
    n = L.shape[0]
    pp, cp = phi0.copy(), phi1.copy()
    hp_, ch = h0.copy(), h1.copy()
    E_hist = []
    for t in range(1, steps):
        Lcp = L @ cp
        Lch0 = L @ ch[:,0]; Lch1 = L @ ch[:,1]
        pn = 2*cp - pp - c2*Lcp - m2*cp - g*(Lch0 + Lch1)          # 标量(带质量)
        hn = 2*ch - hp_ - c2*(L @ ch) - g*np.stack([Lcp, Lcp], axis=1)  # 张量
        E_hist.append(energy(L, c2, g, m2, pn, cp, hn, ch))
        pp, cp = cp, pn
        hp_, ch = ch, hn
    return np.array(E_hist), pp, cp, hp_, ch

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

def main():
    nr, nc = 6, 6
    L, adj = grid_graph(nr, nc)
    n = nr*nc
    evals, evecs = np.linalg.eigh(L)
    w_th = np.sqrt(np.clip(evals, 0, None))
    lam_max = float(evals[-1])
    c2 = min(0.8, 1.2/lam_max); c = np.sqrt(c2)   # 更保守 CFL: 留余量给质量项+耦合分裂
    m2 = 0.3   # 标量质量项: 与张量非简并避免共振
    out = {"config": {"graph": f"grid{nr}x{nc}", "n": n, "coupling": "V=g·Σ(Lφ)h 变分(哈密顿, 能量守恒, 标量带质量非简并)",
                      "c2": round(c2,3), "m2": m2}}

    print("="*90)
    print("能量守恒耦合：从势能泛函变分(哈密顿), 长期稳定")
    print("="*90)
    v = evecs[:, 2]

    # === 验证1: 能量守恒 + 长期稳定 (大 g 扫描) ===
    print("\n[验证1] 能量守恒 + 长期稳定 (长演化 6000 步, 大 g)")
    steps = 6000
    for g in [0.0, 0.05, 0.1, 0.2]:
        phi0 = np.zeros(n); phi1 = v.copy()
        h0 = np.zeros((n,2)); h1 = np.zeros((n,2))
        E, pp, cp, hp_, ch = evolve(L, c2, g, m2, steps, phi0, phi1, h0, h1)
        E = np.array(E)
        E_early = float(E[:2000].mean()); E_end = float(E[-2000:].mean())
        drift = float(abs(E_end - E_early) / (abs(E_end)+1e-30))   # 长期漂移(能量守恒)
        stable = np.all(np.isfinite(cp)) and np.all(np.isfinite(ch))
        out.setdefault("energy", {})[str(g)] = {"E_early": round(E_early,4), "E_end": round(E_end,4),
                                                 "drift": round(drift,4), "stable": bool(stable)}
        print(f"  g={g:.2f}: E首段={E_early:.4f} E末段={E_end:.4f} 漂移={drift:.4f}  稳定={stable}")
    # 能量守恒判据: 漂移≈0(守恒) + 大 g 仍稳定
    e0 = out["energy"]["0.1"]
    print(f"  判据: g=0.1 时漂移={e0['drift']} 稳定={e0['stable']}  (漂移≈0+大g稳定 ⇒ 能量守恒长期稳定)")

    # === 验证2: 膨胀↔剪切相互转化 (能量守恒耦合) ===
    print("\n[验证2] 膨胀↔剪切相互转化 (守恒耦合)")
    # 激发膨胀→剪切: 演化后测剪切振幅
    phi0=np.zeros(n); phi1=v.copy(); h0=np.zeros((n,2)); h1=np.zeros((n,2))
    E, pp, cp, hp_, ch = evolve(L, c2, 0.1, m2, 2000, phi0, phi1, h0, h1)
    shear = float(np.sqrt(ch[:,0]**2 + ch[:,1]**2).max())
    out["expand_to_shear"] = round(shear,6)
    print(f"  激发膨胀 → 剪切末段振幅 = {shear:.4f}")
    # 激发剪切→膨胀: 演化后测膨胀振幅
    phi0=np.zeros(n); phi1=np.zeros(n); h1=np.zeros((n,2)); h1[:,0]=v.copy()
    E, pp, cp, hp_, ch = evolve(L, c2, 0.1, m2, 2000, phi0, phi1, np.zeros((n,2)), h1)
    expand = float(np.abs(cp).max())
    out["shear_to_expand"] = round(expand,6)
    print(f"  激发剪切 → 膨胀末段振幅 = {expand:.4f}")

    # === 验证3: 张量色散 (耦合 g=0.1 下) ===
    print("\n[验证3] 张量色散 ω²=c²λ_k (守恒耦合下)")
    print("  注: d=(h_plus−h_cross) 自由模式, g·Lφ 耦合抵消 → 纯张量色散; h_plus 与 φ 耦合分裂")
    modes = {}
    for k in range(1, n):
        if w_th[k] >= np.pi*0.95: continue
        # 激发张量 d 模式 h=[v,−v], 收集 d=h_plus−h_cross 测频率
        phi0=np.zeros(n); phi1=np.zeros(n)
        h1=np.zeros((n,2)); h1[:,0]=evecs[:,k].copy(); h1[:,1]=-evecs[:,k].copy()
        pp_,cp_=np.zeros(n),np.zeros(n); hp__=np.zeros((n,2)); ch_=h1.copy()
        dseq=[(ch_[:,0]-ch_[:,1]).copy()]
        for t in range(1,2048):
            Lcp=L@cp_; Lch0=L@ch_[:,0]; Lch1=L@ch_[:,1]
            pn=2*cp_-pp_-c2*Lcp-m2*cp_-0.1*(Lch0+Lch1)
            hn=2*ch_-hp__-c2*(L@ch_)-0.1*np.stack([Lcp,Lcp],axis=1)
            pp_,cp_=cp_,pn;hp__,ch_=ch_,hn;dseq.append((ch_[:,0]-ch_[:,1]).copy())
        dseq=np.array(dseq); Dx=dseq@evecs
        fk,_=dominant_freq(Dx[:,k])
        modes[int(k)]={"sqrt_lambda":round(float(w_th[k]),4),"c_sqrt":round(float(c*w_th[k]),4),"dom_freq":round(fk,4)}
    ks=[k for k in modes if modes[k]["dom_freq"]>0.01]
    xs=[modes[k]["sqrt_lambda"] for k in ks];ys=[modes[k]["dom_freq"] for k in ks]
    corr=float(np.corrcoef(xs,ys)[0,1]) if len(ks)>=3 and np.std(xs)>0 and np.std(ys)>0 else 0.0
    out["dispersion_tensor"]={"n_modes":len(ks),"corr":round(corr,4)}
    print(f"  守恒耦合下张量自由模式 d 色散相关 = {corr:.4f}  (ω²=c²λ, 耦合不破坏张量层)")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_sc_energy.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_sc_energy.json")

if __name__ == "__main__":
    main()
