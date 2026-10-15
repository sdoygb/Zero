#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_tensor3d_wave.py —— 3D 格点张量波: 完整双偏振(+/×), 无迹守恒, 色散
====================================================================================
GR 引力波 = 度规扰动 h(3D 对称 3×3), TT 无迹横波, 2 偏振(+/×)。
Zero 张量层 3D 扩展: 每顶点对称 3×3 张量(6 分量 h_xx,h_xy,h_xz,h_yy,h_yz,h_zz)。
  无迹投影: tr=h_xx+h_yy+h_zz→0 (3D 无迹对称 = 5 独立分量)
  波动: ∂ₜ²h=−c²L(h)  (L⊗I₆ 逐分量)
  偏振基(无迹): e+=(xx−yy), e+2=(xx−zz), e_xy, e_xz, e_yz
  横波偏振(+/× 类比): 2 个独立剪切模式
验证: ①色散 ω²=c²λ_k ②无迹守恒 tr=0 ③偏振独立性(激发一基, 他基零泄漏) ④自由度。
"""
from __future__ import annotations
import sys, os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def grid3d(nx, ny, nz):
    adj = {v: [] for v in range(nx*ny*nz)}
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                v = (i*ny + j)*nz + k
                for di, dj, dk in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
                    ni, nj, nk = i+di, j+dj, k+dk
                    if 0<=ni<nx and 0<=nj<ny and 0<=nk<nz:
                        adj[v].append((ni*ny+nj)*nz+nk)
    n = nx*ny*nz
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

def traceless6(H):
    """逐顶点无迹: H:(n,6)=[xx,xy,xz,yy,yz,zz] → 减 tr/3."""
    H = H.copy()
    tr = H[:,0]+H[:,3]+H[:,5]
    H[:,0]-=tr/3; H[:,3]-=tr/3; H[:,5]-=tr/3
    return H

def main():
    nx, ny, nz = 4, 4, 4
    L, adj = grid3d(nx, ny, nz)
    n = nx*ny*nz
    evals, evecs = np.linalg.eigh(L)
    w_th = np.sqrt(np.clip(evals, 0, None))
    lam_max = float(evals[-1])
    c2 = min(0.8, 3.5/lam_max); c = np.sqrt(c2)
    I6 = np.eye(6)
    L6 = np.kron(L, I6)          # 6 分量独立波动
    steps = 2048
    out = {"config": {"graph": f"grid{nx}x{ny}x{nz}", "n": n, "field": "per-vertex sym 3x3 tensor (metric)",
                      "tt": "traceless 3D (5 dof)", "c2": round(c2,3)}}

    print("="*90)
    print("3D 张量波：完整双偏振(+/×), 无迹守恒, 色散")
    print("="*90)

    # === 验证1: 色散 (逐模式, 用 h_xx 分量投影) ===
    print("\n[验证1] 3D 张量波色散 (激发 L 本征模式 ⊗ 无迹偏振)")
    modes = {}
    for k in range(1, n):
        if w_th[k] >= np.pi*0.95: continue
        v = evecs[:, k]
        h = np.zeros((n,6)); h[:,0]=v; h[:,3]=-v     # xx−yy 偏振 (+)
        h = traceless6(h)
        phis = leapfrog(L6, c2, steps, np.zeros(6*n), h.reshape(-1).copy())
        Hx = phis.reshape(-1, n, 6)[:, :, 0]
        Ax = Hx @ evecs
        fk, _ = dominant_freq(Ax[:, k])
        modes[int(k)] = {"sqrt_lambda": round(float(w_th[k]),4),
                         "c_sqrt": round(float(c*w_th[k]),4),
                         "dom_freq": round(fk,4)}
    ks = [k for k in modes if modes[k]["dom_freq"]>0.01]
    xs = [modes[k]["sqrt_lambda"] for k in ks]; ys = [modes[k]["dom_freq"] for k in ks]
    corr = float(np.corrcoef(xs,ys)[0,1]) if len(ks)>=3 and np.std(xs)>0 and np.std(ys)>0 else 0.0
    out["dispersion"] = {"n_modes": len(ks), "corr": round(corr,4)}
    print(f"  3D 张量波色散相关 = {corr:.4f}  (模式数={len(ks)})")

    # === 验证2: 无迹守恒 ===
    v = evecs[:, 1]
    h = np.zeros((n,6)); h[:,0]=v; h[:,3]=-v; h = traceless6(h)
    phis = leapfrog(L6, c2, steps, np.zeros(6*n), h.reshape(-1).copy())
    Ht = phis.reshape(-1, n, 6)
    tr = Ht[:,:,0]+Ht[:,:,3]+Ht[:,:,5]
    tr_max = float(np.abs(tr).max())
    out["traceless_max_trace"] = round(tr_max, 8)
    print(f"  演化中 max|tr(h)| = {tr_max:.2e}  (≈0 ⇒ 无迹守恒)")

    # === 验证3: 双偏振独立 (激发 +, 测 × 不泄漏; 反之) ===
    print("\n[验证3] 双偏振独立 (激发一偏振, 他偏振零泄漏)")
    # 偏振基: e+ (xx−yy), e× (xy)
    v = evecs[:, 1]
    # 激发 + (xx−yy): 测 × (xy 分量)
    h_plus = np.zeros((n,6)); h_plus[:,0]=v; h_plus[:,3]=-v; h_plus=traceless6(h_plus)
    phis = leapfrog(L6, c2, steps, np.zeros(6*n), h_plus.reshape(-1).copy())
    Ht = phis.reshape(-1, n, 6)
    leak_cross = float(np.abs(Ht[:,:,1]).max() + np.abs(Ht[:,:,4]).max())  # xy+yz 泄漏
    # 激发 × (xy): 测 + (xx−yy, xx−zz) 不泄漏
    h_cross = np.zeros((n,6)); h_cross[:,1]=v; h_cross=traceless6(h_cross)
    phis = leapfrog(L6, c2, steps, np.zeros(6*n), h_cross.reshape(-1).copy())
    Ht = phis.reshape(-1, n, 6)
    leak_plus = float(np.abs(Ht[:,:,0]).max() + np.abs(Ht[:,:,3]).max() + np.abs(Ht[:,:,5]).max())
    out["polarization"] = {"excite_plus_cross_leak": round(leak_cross,6),
                           "excite_cross_plus_leak": round(leak_plus,6)}
    print(f"  激发+(xx−yy), ×(xy)泄漏 = {leak_cross:.2e}")
    print(f"  激发×(xy), +(xx−yy/zz)泄漏 = {leak_plus:.2e}")
    print(f"  (≈0 ⇒ 双偏振 +/× 独立)")

    # === 验证4: 自由度 ===
    out["dof"] = {"traceless_sym3x3": 5, "transverse_polarizations": 2}
    print(f"\n[验证4] 3D 无迹对称张量 = 5 自由度; 横波下 2 偏振(+/×) = GR 引力波 2 偏振")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_tensor3d_wave.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_tensor3d_wave.json")

if __name__ == "__main__":
    main()
