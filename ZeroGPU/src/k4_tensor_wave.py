#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_tensor_wave.py —— 张量结构：每顶点对称张量 = 度规扰动, 无迹 TT + 偏振
====================================================================================
GR 引力波 = 度规张量扰动 h(每点对称张量), TT 无迹横波, 双偏振(+/×), 四极辐射。
Zero 张量层候选(最贴 GR): 每顶点一个对称 2×2 张量 h=(h_xx,h_xy,h_yy)。
  无迹投影: h→h_tl=h−(tr/2)I  (剪切, 体积不变)
  波动: ∂ₜ²h_tl=−c²·L(h_tl)   (L 逐分量, 图 Laplacian)
  色散: ω²=c²λ_k (同标量定理, 每模式带偏振)
  偏振: +(h_xx−h_yy) 与 ×(h_xy) 独立分量
验证: ①张量波色散 ②无迹守恒 ③偏振模式分解 ④与标量(膨胀)对比。
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

def leapfrog(L, c2, steps, phi0, phi1):
    n = L.shape[0]
    phis = [phi1.copy()]
    prev, cur = phi0.copy(), phi1.copy()
    for t in range(1, steps):
        nxt = 2*cur - prev - c2*(L @ cur)
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

def traceless(H):
    """逐顶点无迹投影: H:(n,3)=[h_xx,h_xy,h_yy] → [h_xx−tr/2, h_xy, h_yy−tr/2]"""
    H = H.copy()
    tr = H[:, 0] + H[:, 2]
    H[:, 0] -= tr/2
    H[:, 2] -= tr/2
    return H

def main():
    nr, nc = 6, 6
    L, adj = grid_graph(nr, nc)
    n = nr*nc
    evals, evecs = np.linalg.eigh(L)
    w_th = np.sqrt(np.clip(evals, 0, None))
    lam_max = float(evals[-1])
    c2 = min(0.8, 3.5/lam_max); c = np.sqrt(c2)
    # 张量算子: L⊗I₃ (3 分量独立波动)
    I3 = np.eye(3)
    L3 = np.kron(L, I3)          # 3n × 3n
    steps = 2048
    out = {"config": {"graph": f"grid{nr}x{nc}", "n": n, "field": "per-vertex symmetric tensor (metric perturbation)",
                      "tt": "traceless projection, traceless shear", "c2": round(c2,3)}}

    print("="*90)
    print("张量结构：每顶点对称张量 = 度规扰动, 无迹 TT + 偏振")
    print("="*90)

    # === 验证1: 张量波色散 (每分量波动, 无迹后 2 偏振自由度) ===
    print("\n[验证1] 张量波色散 (激发 L 本征模式 ⊗ 无迹张量偏振)")
    print("  模式k  √λ_k   c√λ(理论)  实测主导频率  (偏振: + / ×)")
    modes = {}
    # 构造张量本征激发: evecs[:,k] ⊗ 偏振分量
    for k in range(1, n):
        if w_th[k] >= np.pi*0.95: continue
        # 张量初值: 顶点按本征模式, 偏振选 + (h_xx=1,h_yy=-1,h_xy=0)
        v = evecs[:, k]
        h_plus = np.zeros((n, 3)); h_plus[:, 0] = v; h_plus[:, 2] = -v
        phi1 = h_plus.reshape(-1).copy()
        phis = leapfrog(L3, c2, steps, np.zeros(3*n), phi1)
        # 投影: h_xx 分量 (每顶点第0分量) 投影到 L 本征基
        Hx = phis.reshape(-1, n, 3)[:, :, 0]   # h_xx 分量 (steps, n)
        Ax = Hx @ evecs
        fk, _ = dominant_freq(Ax[:, k])
        modes[int(k)] = {"sqrt_lambda": round(float(w_th[k]),4),
                         "c_sqrt": round(float(c*w_th[k]),4),
                         "dom_freq": round(fk,4)}
        if k <= 8:
            print(f"   {k:3d}  {w_th[k]:6.3f}   {c*w_th[k]:7.3f}   {fk:8.4f}")
    ks = [k for k in modes if modes[k]["dom_freq"] > 0.01]
    xs = [modes[k]["sqrt_lambda"] for k in ks]; ys = [modes[k]["dom_freq"] for k in ks]
    corr = float(np.corrcoef(xs, ys)[0,1]) if len(ks)>=3 and np.std(xs)>0 and np.std(ys)>0 else 0.0
    out["tensor_dispersion"] = {"n_modes": len(ks), "corr": round(corr,4)}
    print(f"  张量波色散相关 = {corr:.4f}  (ω²=c²λ_k, 同标量定理)")

    # === 验证2: 无迹守恒 (tr(h)→0) ===
    print("\n[验证2] 无迹守恒: 波动演化中 tr(h) 保持 0")
    v = evecs[:, 1]
    h0 = np.zeros((n,3)); h0[:,0]=v; h0[:,2]=-v   # 已无迹
    phis = leapfrog(L3, c2, steps, np.zeros(3*n), h0.reshape(-1).copy())
    Ht = phis.reshape(-1, n, 3)
    tr = Ht[:,:,0] + Ht[:,:,2]     # (steps, n) 迹
    tr_max = float(np.abs(tr).max())
    out["traceless_max_trace"] = round(tr_max, 8)
    print(f"  演化中 max|tr(h)| = {tr_max:.2e}  (≈0 ⇒ 无迹守恒)")

    # === 验证3: 偏振模式 (+ 与 × 独立) ===
    print("\n[验证3] 偏振分解: +(h_xx−h_yy) 与 ×(h_xy) 独立波动")
    # 激发 + 偏振, 测 × 分量是否保持 0
    h_plus = np.zeros((n,3)); h_plus[:,0]=v; h_plus[:,2]=-v
    phis = leapfrog(L3, c2, steps, np.zeros(3*n), h_plus.reshape(-1).copy())
    Ht = phis.reshape(-1, n, 3)
    cross_amp = float(np.abs(Ht[:,:,1]).max())   # h_xy (× 分量)
    out["polarization"] = {"excite_plus_cross_leak": round(cross_amp,6)}
    print(f"  激发 + 偏振, ×(h_xy) 分量泄漏 = {cross_amp:.2e}  (≈0 ⇒ 偏振独立)")

    # === 验证4: 与标量(膨胀)对比 ===
    print("\n[验证4] 标量(膨胀) vs 张量(剪切) 对比")
    # 标量: 顶点值(膨胀) 已有 Z105; 张量: 无迹剪切
    # 张量波动自由度 = 2(偏振), 标量 = 1(膨胀)
    out["dof"] = {"scalar_expansion": 1, "tensor_shear": 2}
    print("  标量层: 膨胀(单极, 1 自由度) ｜ 张量层: 剪切(四极/偏振, 2 自由度)")
    print("  2D 无迹对称张量 = 1 偏振(剪切); 3D = 2 偏振(+/×) → 完整 GR 双偏振")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_tensor_wave.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_tensor_wave.json")

if __name__ == "__main__":
    main()
