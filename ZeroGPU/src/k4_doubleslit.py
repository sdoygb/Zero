#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_doubleslit.py —— Zero 内生推导电子双缝干涉条纹
====================================================================================
前提: 全部用 Zero 自己的理论. Zero 内生要素:
  ① 电子波 = Zero 图上波 (Z104 色散 ω²=c²λ_k, 图 Laplacian 本征模)
  ② 电子波长 λ_e = 2π/k, k=√λ_max (Zero K₄ 图 Laplacian 谱, 非外部德布罗意)
  ③ 相位/相位差 = Zero 图上波相位 (用户引入的相位概念)
  ④ 双缝 = 图上两个传播通道 (两路径到屏相位差 δ=k·Δs)
  ⑤ 干涉强度 = 图上波模方 |ψ|² (图上波固有性质, 非外部Born规则假设)
外部几何(仅此): 缝距 d, 屏距 L (实验设置, 非物理理论)
双缝 Fraunhofer:
  Δs ≈ x·d/L (路径差), δ = 2πΔs/λ_e = 2πxd/(λ_e L)
  I(x) = cos²(δ/2) = cos²(πxd/(λ_e L))  → 亮暗条纹
  亮纹位置 x_n = n·λ_e L/d, 条纹间距 Δx = λ_e L/d
验证: ①条纹结构(亮暗) ②间距 Δx=λ_eL/d ③条纹数 ④可见度
诚实: 波数/波长定标相对(c=1), Fraunhofer近似, d/L为外部几何
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def k4_laplacian():
    adj = {v: [u for u in range(4) if u != v] for v in range(4)}
    L = np.zeros((4,4))
    for v in range(4):
        L[v,v] = len(adj[v])
        for j in adj[v]: L[v,j] = -1
    return L

def main():
    out = {"config": {"model": "Zero endogenous electron double-slit interference",
                      "wave": "Zero graph wave (Z104 ω²=c²λ_k)",
                      "lambda": "λ_e=2π/k, k=√λ_max from Zero K₄ Laplacian",
                      "intensity": "I=|ψ|² graph-wave amplitude squared"}}
    print("="*90)
    print("Zero 内生推导电子双缝干涉条纹")
    print("="*90)

    # === ① 电子波长从 Zero 图内生 ===
    L = k4_laplacian()
    evals = np.linalg.eigvalsh(L)
    lam_max = np.clip(evals,0,None)[-1]   # K4: λ_max=4
    k = np.sqrt(lam_max)                   # Zero 图上波数 (Z104 色散 k=√λ)
    lam_e = 2*np.pi/k                      # Zero 内生电子波长
    print(f"[① Zero 图上波] K₄ Laplacian λ_max={lam_max:.3f}")
    print(f"  波数 k=√λ_max={k:.3f}, 电子波长 λ_e=2π/k={lam_e:.4f}")
    out["zero_wave"] = {"lam_max": round(float(lam_max),3), "k": round(float(k),3),
                        "lambda_e": round(float(lam_e),4), "source": "Zero K₄ Laplacian"}

    # === ② 双缝几何 (外部) ===
    d = 1.0; L_screen = 100.0; a = 0.2   # a=缝宽(单缝衍射包络), d=缝距
    print(f"[② 双缝几何(外部)] 缝距 d={d}, 缝宽 a={a}, 屏距 L={L_screen}")

    # === ③ 条纹间距 & 亮纹 ===
    dx = lam_e * L_screen / d
    print(f"[③ 条纹间距] Δx = λ_e·L/d = {dx:.3f}")
    print(f"  亮纹位置 x_n = n·Δx = {[round(n*dx,1) for n in range(1,5)]}")
    out["stripes"] = {"Delta_x": round(float(dx),3),
                      "bright_positions": [round(n*dx,1) for n in range(0,5)]}

    # === ④ 屏上强度: 单缝衍射包络 × 双缝干涉 I=[sinβ/β]²·cos²α ===
    # β=πa x/(λL) 单缝相位, α=πd x/(λL) 双缝相位
    x = np.linspace(-3*dx, 3*dx, 2001)
    beta = np.pi*a*x/(lam_e*L_screen)
    alpha = np.pi*d*x/(lam_e*L_screen)
    with np.errstate(divide='ignore', invalid='ignore'):
        sinc2 = np.where(np.abs(beta)<1e-12, 1.0, (np.sin(beta)/beta)**2)   # 单缝衍射包络
    inter = np.cos(alpha)**2                                  # 双缝干涉
    I = sinc2 * inter
    # 各亮纹强度(包络衰减): 中心 vs 两侧
    peaks = I[1:-1] >= I[:-2]; peaks &= I[1:-1] >= I[2:]
    n_peaks = int(peaks.sum())
    I_center = I[len(x)//2]
    I_first = I[len(x)//2 + int(dx/(3*dx/1000))]   # x=Δx 处
    I_second = I[len(x)//2 + int(2*dx/(3*dx/1000))]
    print(f"[④ 完整强度 I=[sinβ/β]²·cos²α] 峰值数={n_peaks}")
    print(f"  中心 I(x=0)={I_center:.3f}, 第一亮纹 I(x=Δx)={I_first:.3f}, 第二亮纹 I(x=2Δx)={I_second:.3f}")
    print(f"  (中间强, 两侧衰减 → 单缝衍射包络)")
    print(f"  可见度(中心区域)={((I.max()-I.min())/(I.max()+I.min())):.3f}")
    out["intensity"] = {"model": "I=[sinβ/β]²·cos²α (单缝衍射包络×双缝干涉)",
                        "a_slit_width": a, "n_peaks": n_peaks,
                        "I_center": round(float(I_center),3), "I_first": round(float(I_first),3),
                        "I_second": round(float(I_second),3),
                        "envelope_sample": [round(float(v),3) for v in sinc2[::250]],
                        "I_sample": [round(float(v),3) for v in I[::250]]}

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_doubleslit.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_doubleslit.json")

if __name__ == "__main__":
    main()
