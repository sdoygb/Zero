#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_superconductor.py —— 图上声子 → Debye 温度 → BCS Tc 结构退化
====================================================================================
超导 BCS 核心:  Tc ≈ 1.13·Θ_D·exp(−1/λ), Θ_D=Debye温度(晶格振动最高频率)。
Zero 图振动谱(Z104 色散 ω²=c²λ_k) = 图上声子谱 → 可内生算 Θ_D:
  声子频率:  ω_k = √(c²λ_k),  λ_k=图Laplacian特征值
  Debye温度: Θ_D ∝ ω_max = c√λ_max  (最高声子频率)
  BCS Tc:    Tc = 1.13·Θ_D·exp(−1/λ_ep), λ_ep=电子-声子耦合
Zero 增量: 不同图(晶格类型) 给不同声子谱 → 不同 Θ_D → Tc 结构差异预言
验证: ①各图声子谱→Θ_D ②BCS Tc(λ扫描) ③不同图Tc差异 ④图上声子密度→Debye模式数
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def make_graph(kind, n):
    if kind == "k4":
        # K4: 4顶点全连接
        adj = {v: [u for u in range(4) if u != v] for v in range(4)}
    elif kind.startswith("grid"):
        nr = nc = int(np.sqrt(n))
        adj = {}
        for i in range(nr):
            for j in range(nc):
                v = i*nc+j; adj[v]=[]
                for di,dj in [(1,0),(-1,0),(0,1),(0,-1)]:
                    ni,nj=i+di,j+dj
                    if 0<=ni<nr and 0<=nj<nc: adj[v].append(ni*nc+nj)
    elif kind == "path":
        adj = {v: [v-1] if v>0 else [] for v in range(n)}
        for v in range(n-1): adj[v].append(v+1)
    elif kind == "cycle":
        adj = {v: [] for v in range(n)}
        for v in range(n):
            adj[v].append((v-1)%n); adj[v].append((v+1)%n)
    nvert = len(adj)
    L = np.zeros((nvert, nvert))
    for v in range(nvert):
        L[v,v] = len(adj[v])
        for j in adj[v]: L[v,j] = -1
    return L, nvert

def main():
    graphs = {"K4": ("k4", 4), "Path16": ("path", 16), "Cycle16": ("cycle", 16),
              "Grid4x4": ("grid", 16), "Grid8x8": ("grid", 64)}
    out = {"config": {"model": "Zero graph phonon → Debye → BCS Tc",
                      "Tc": "Tc=1.13·Θ_D·exp(−1/λ_ep)"}}
    print("="*90)
    print("图上声子 → Debye 温度 → BCS Tc 结构退化")
    print("="*90)

    c = 1.0   # 图上波速(相对; 物理定标未决, 用相对 Θ_D)
    lam_ep_scan = [0.3, 0.5, 0.8, 1.0, 1.5]

    results = {}
    for name, (kind, n) in graphs.items():
        L, nv = make_graph(kind, n)
        evals = np.linalg.eigvalsh(L)
        evals = np.clip(evals, 0, None)
        lam_max = float(evals[-1])
        lam_1 = float(evals[1]) if nv > 1 else 0.0   # 最低非零声子
        # 声子频率谱
        omegas = np.sqrt(c**2 * evals)
        # Debye 温度 ∝ 最高声子频率
        Theta_D = c*np.sqrt(lam_max)
        # Debye 模式数(非零声子 = nv-1)
        n_modes = nv - 1
        # BCS Tc: Tc=1.13·Θ_D·exp(−1/λ)
        Tcs = {round(lam,2): round(1.13*Theta_D*np.exp(-1.0/lam), 4) for lam in lam_ep_scan}
        results[name] = {"n_vert": nv, "lam_max": round(lam_max,3), "lam_1": round(lam_1,3),
                         "n_phonon_modes": n_modes, "Theta_D_norm": round(Theta_D,3),
                         "Tc_vs_lambda": Tcs}
        print(f"\n[{name}] 顶点={nv}, 声子模式={n_modes}")
        print(f"  声子谱: λ_max={lam_max:.3f}, 最低非零 λ_1={lam_1:.3f}")
        print(f"  Debye 温度 Θ_D ∝ ω_max = {Theta_D:.3f}")
        print(f"  BCS Tc(λ_ep 扫描): ", {k: f"{v:.4f}" for k,v in Tcs.items()})

    # 对比: 不同图(晶格)的 Θ_D 与 Tc 差异
    print("\n[对比] 不同图(晶格类型) 的 Θ_D 与 Tc")
    for name in graphs:
        r = results[name]
        print(f"  {name}: Θ_D={r['Theta_D_norm']:.3f}, Tc(λ=0.8)={r['Tc_vs_lambda'][0.8]:.4f}")
    out["results"] = results
    # 比值: 各图 Tc 相对 K4
    base = results["K4"]["Tc_vs_lambda"][0.8]
    ratio = {name: round(results[name]["Tc_vs_lambda"][0.8]/base, 2) for name in graphs}
    out["Tc_ratio_to_K4"] = ratio
    print(f"\n  Tc 相对 K4: {ratio}")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_superconductor.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_superconductor.json")

if __name__ == "__main__":
    main()
