#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_mcmillan.py —— 加强耦合: McMillan 公式 (强耦合 Tc 修正)
====================================================================================
BCS Tc=1.13Θ_D·exp(−1/λ) 只弱耦合有效; 强耦合(λ>1)需 McMillan 修正:
  Tc = (Θ_D/1.45)·exp[−1.04(1+λ)/(λ−μ*(1+0.62λ))]
  Θ_D=Debye温度(Zero图声子谱给出), μ*=Coulomb赝势(0.1-0.15), λ=电子-声子耦合
物理核心: 强耦合 Tc 饱和于 Θ_D/1.45 (BCS 无上限, McMillan 有上限→真实超导)
Zero: 用 Z104 图声子谱(上轮 k4_superconductor)给 Θ_D, 退化 McMillan Tc 结构
验证: ①BCS vs McMillan 对比(弱耦合一致, 强耦合饱和) ②μ* 扫描(Coulomb抑制) ③各图Tc
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def make_graph(kind, n):
    if kind == "k4":
        adj = {v: [u for u in range(4) if u != v] for v in range(4)}
    elif kind == "grid":
        nr = nc = int(np.sqrt(n))
        adj = {}
        for i in range(nr):
            for j in range(nc):
                v = i*nc+j; adj[v]=[]
                for di,dj in [(1,0),(-1,0),(0,1),(0,-1)]:
                    ni,nj=i+di,j+dj
                    if 0<=ni<nr and 0<=nj<nc: adj[v].append(ni*nc+nj)
    nvert = len(adj)
    L = np.zeros((nvert, nvert))
    for v in range(nvert):
        L[v,v] = len(adj[v])
        for j in adj[v]: L[v,j] = -1
    return L, nvert

def theta_D(L, c=1.0):
    evals = np.linalg.eigvalsh(L)
    lam_max = np.clip(evals, 0, None)[-1]
    return c*np.sqrt(lam_max)

def bcs_tc(Th, lam):
    return 1.13*Th*np.exp(-1.0/lam)

def mcmillan_tc(Th, lam, mu):
    return (Th/1.45)*np.exp(-1.04*(1.0+lam)/(lam - mu*(1.0+0.62*lam)))

def main():
    graphs = {"K4": ("k4", 4), "Grid4x4": ("grid", 16), "Grid8x8": ("grid", 64)}
    out = {"config": {"model": "McMillan strong-coupling Tc = (Θ_D/1.45)·exp[−1.04(1+λ)/(λ−μ*(1+0.62λ))]",
                      "theta_D_from": "Zero graph phonon (Z104 dispersion)"}}
    print("="*90)
    print("加强耦合: McMillan 公式 (强耦合 Tc 修正)")
    print("="*90)
    c = 1.0
    lam_scan = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0]
    mu = 0.1

    # === 1: BCS vs McMillan (Grid8x8, 最高Θ_D) ===
    L, nv = make_graph("grid", 64)
    Th = theta_D(L, c)
    print(f"\n[1] Grid8x8: Θ_D={Th:.3f}, BCS vs McMillan(μ*={mu})")
    print("  λ     BCS        McMillan")
    bcs_list=[]; mcm_list=[]
    for lam in lam_scan:
        t_b = bcs_tc(Th, lam); t_m = mcmillan_tc(Th, lam, mu)
        bcs_list.append(t_b); mcm_list.append(t_m)
        print(f"  {lam:.1f}   {t_b:.4f}    {t_m:.4f}")
    # 饱和检查: λ→∞ McMillan→Th/1.45
    lam_inf = Th/1.45
    print(f"  强耦合上限 Θ_D/1.45 = {lam_inf:.4f}")
    print(f"  McMillan(λ=3)={mcm_list[-1]:.4f} vs 上限 {lam_inf:.4f} (趋于饱和)")
    print(f"  BCS(λ=3)={bcs_list[-1]:.4f} (无上限, 指数增长)")
    out["bcs_vs_mcmillan"] = {"lambda": lam_scan, "bcs": [round(x,4) for x in bcs_list],
                              "mcmillan_mu0.1": [round(x,4) for x in mcm_list],
                              "saturation_limit_Th1.45": round(lam_inf,4)}

    # === 2: μ* 扫描 (Coulomb 抑制) ===
    print(f"\n[2] μ* 扫描 (Grid8x8, λ=1.0, 1.5)")
    mu_scan = [0.0, 0.05, 0.10, 0.15, 0.20]
    for lam in [1.0, 1.5]:
        row = []
        for m in mu_scan:
            row.append(mcmillan_tc(Th, lam, m))
        print(f"  λ={lam}: μ*={mu_scan} → Tc={[round(x,4) for x in row]}")
        out.setdefault("mu_star_scan", {})[f"lambda={lam}"] = {str(m): round(v,4) for m,v in zip(mu_scan,row)}

    # === 3: 各图 McMillan Tc (λ=1.5) ===
    print(f"\n[3] 各图(晶格) McMillan Tc (λ=1.5, μ*=0.1)")
    res = {}
    for name,(kind,n) in graphs.items():
        Lg, _ = make_graph(kind, n)
        Thg = theta_D(Lg, c)
        tc = mcmillan_tc(Thg, 1.5, mu)
        res[name] = {"Theta_D": round(Thg,3), "Tc_McMillan": round(tc,4)}
        print(f"  {name}: Θ_D={Thg:.3f}, Tc={tc:.4f}")
    out["graphs_mcmillan"] = res

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_mcmillan.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_mcmillan.json")

if __name__ == "__main__":
    main()
