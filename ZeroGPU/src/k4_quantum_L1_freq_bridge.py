#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_quantum_L1_freq_bridge.py —— B' 验证: L1 稳态(内禀频率 f)作 L2 波源 → 测波频率是否=f
======================================================================================
用户选定 B': 稳态的频率特征进波. L1 稳态 ω_v=record_v/Σ (中间两顶点 2:1, Z087)
编码内禀频率 f=log(p/q)=log2. 把 L1 稳态(2:1 签名)直接作 L2 波的时间演化源,
测图上波振荡频率 vs f.

关键: 源振幅按 L1 稳态的模流时钟调制 (f=log2, Z027 命题3 内禀时间),
     一阶量子行走 i∂tψ=c²Lψ, 测传播区振荡频率 → 对比 f.
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

def path_laplacian(n):
    L = np.zeros((n, n))
    for i in range(n):
        L[i, i] += (1 if i > 0 else 0) + (1 if i < n - 1 else 0)
        if i > 0: L[i, i-1] -= 1
        if i < n-1: L[i, i+1] -= 1
    return L

def main():
    P_("=" * 96)
    P_("B' 验证: L1 稳态(内禀频率 f)作 L2 波源 → 波频率是否=f (k4_quantum_L1_freq_bridge)")
    P_("=" * 96)
    out = {}

    # ---- L1 稳态: 中间两顶点 2:1, 内禀频率 f=log2 (Z087) ----
    f = np.log(2)
    lam_half = 0.5                              # 稳态两顶点分布 (2/3, 1/3)
    P_(f"L1 稳态签名: 中间两顶点 2:1, 内禀频率 f=log(p/q)=log2={f:.4f} (Z087 实测, path8-32 普适)")
    out["L1_steady"] = {"ratio": "2:1", "f_intrinsic": round(f, 5)}

    # ---- 图: path16 (Z104 引擎图), 源区=中间两顶点 ----
    n = 16
    src = [7, 8]
    c2 = 0.8                                    # Z104 值
    L = path_laplacian(n)
    lam, V = np.linalg.eigh(L)

    # ---- 源调制: 振幅按 L1 模流时钟 cos(f·t) 振荡, 2:1 分布注入 ----
    # 一阶量子行走 i∂tψ=c²Lψ, 蛙跳时间步进 + 源调制
    T = 4000
    dt = 0.1
    psi = np.zeros(n, complex)
    w = np.zeros(n, complex)                    # 一阶: 用 i∂tψ=Hψ, 蛙跳 w 为辅助
    # 本征演化 ψ(t)=Σc e^{-i c² λ t}, 源调制加在每步
    c_n = np.zeros(n, complex)
    # 源: 两顶点振幅 a1=√(2/3), a2=√(1/3) (稳态签名), 调制 cos(f t)
    a1, a2 = np.sqrt(2/3), np.sqrt(1/3)
    freqs = []
    probe = np.zeros(T, complex)
    for t in range(T):
        tt = t * dt
        # 源调制注入(每步加到两源顶点)
        amp = np.cos(f * tt)                    # L1 模流时钟调制
        # 用本征分解精确演化: 源作为非齐次项近似(蛙跳)
        psi[src[0]] += a1 * amp * dt
        psi[src[1]] += a2 * amp * dt
        # 演化半步 (i∂tψ=c²Lψ)
        psi = psi - 1j * c2 * (L @ psi) * dt
        probe[t] = psi[12]                      # 探测传播区顶点12
    # FFT 测主导频率
    from numpy.fft import rfft, rfftfreq
    fr = rfftfreq(T, dt)
    F = np.abs(rfft(probe.real))
    f_peak = fr[np.argmax(F[1:]) + 1]
    P_(f"源调制频率 = {f:.4f} (L1 内禀时钟 log2)")
    P_(f"一阶量子行走 T={T} 步, dt={dt}, 探测顶点12:")
    P_(f"  实测波主导频率 = {f_peak:.4f}")
    P_(f"  比值 f_peak/f = {f_peak/f:.4f}  (≈1 ⟹ L1 稳态频率进波 ✓)")
    out["bridge_B"] = {"source_modulation_f": round(f, 5), "measured_wave_freq": round(f_peak, 5),
                       "ratio": round(f_peak / f, 5), "engine": "i∂tψ=c²Lψ (path16, c²=0.8)"}

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_quantum_L1_freq_bridge.json"), "w", encoding="utf-8") as f_:
        json.dump(out, f_, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "k4_quantum_L1_freq_bridge.txt"), "w", encoding="utf-8") as f_:
        f_.write("\n".join(LOG))
    P_("\n已存 k4_quantum_L1_freq_bridge.json / .txt")

if __name__ == "__main__":
    main()
