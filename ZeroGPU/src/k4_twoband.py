#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_twoband.py —— 图上双带 → 双能隙(各向异性 s), 对应 MgB₂
====================================================================================
MgB₂(39K): 双能隙各向异性 s-wave — σ带(硼层E₂g声子强耦合) Δ_σ≈7meV,
            π带(弱耦合) Δ_π≈1.7meV, 比≈4.1. 非d-wave节点, 全同号各向异性.
Zero: 图上分两带(σ强/π弱耦合顶点), 双带BCS能隙方程 → 双能隙
双带能隙方程:
  Δ_σ = Σ_σ V_σσΔ_σ/2E_σ + Σ_π V_σπΔ_π/2E_π
  Δ_π = Σ_σ V_πσΔ_σ/2E_σ + Σ_π V_ππΔ_π/2E_π
  带间 V_σπ 混合两带
验证: ①双能隙分离(σ大π小) ②比 Δ_σ/Δ_π  vs MgB₂≈4.1 ③各向异性s(同号) ④T_c
诚实: 图上声子介导Z114有符号变化; 双带建模用带内带间势, 检验能否给各向异性s
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def solve_twoband(Vss, Vpp, Vsp, Ns, Np, wD, tol=1e-10, maxit=20000):
    """双带自洽能隙 (两带各自费米能级 ξ 在 [-wD,wD])"""
    xis = np.linspace(-wD, wD, Ns); xip = np.linspace(-wD, wD, Np)
    Ds, Dp = 0.5, 0.3
    for it in range(maxit):
        Es = np.sqrt(xis**2 + Ds**2); Ep = np.sqrt(xip**2 + Dp**2)
        S_s = np.sum(Ds/Es); S_p = np.sum(Dp/Ep)
        Dsn = 0.5*(Vss*S_s + Vsp*S_p)
        Dpn = 0.5*(Vsp*S_s + Vpp*S_p)   # V_σπ=V_πσ
        # 阻尼
        Ds2 = 0.3*Ds + 0.7*Dsn; Dp2 = 0.3*Dp + 0.7*Dpn
        if abs(Ds2-Ds)+abs(Dp2-Dp) < tol:
            Ds, Dp = Ds2, Dp2; break
        Ds, Dp = Ds2, Dp2
    return Ds, Dp, it+1

def main():
    out = {"config": {"model": "Zero two-band → two-gap anisotropic s, vs MgB₂ σ/π",
                      "mgb2": "Δ_σ≈7meV, Δ_π≈1.7meV, ratio≈4.1, T_c=39K"}}
    print("="*90)
    print("图上双带 → 双能隙(各向异性 s), 对应 MgB₂")
    print("="*90)
    wD = 2.774  # Grid8x8 Θ_D (图上声子窗)
    Ns, Np = 32, 32   # σ带/π带 顶点数

    # 扫描耦合比 g_σ/g_π, 固定带间混合
    print(f"\n双带能隙 (wD={wD}, Ns={Ns}, Np={Np}):")
    print("  gσ:gπ   带间r    Δ_σ    Δ_π    Δ_σ/Δ_π")
    rows = []
    for (gs,gp), r in [((0.5,0.15),0.1), ((0.8,0.2),0.15), ((1.0,0.25),0.15), ((1.2,0.3),0.2),
                       ((1.5,0.4),0.2), ((1.5,0.6),0.2)]:   # (1.5,0.6)→Δσ/Δπ≈4.5 最接近MgB₂
        Vss, Vpp, Vsp = gs**2, gp**2, gs*gp*r
        Ds, Dp, its = solve_twoband(Vss, Vpp, Vsp, Ns, Np, wD)
        ratio = Ds/Dp if Dp>0 else 0
        sym = "各向异性s(同号)" if Ds>0 and Dp>0 else "?"
        rows.append({"g_sigma":gs, "g_pi":gp, "r_cross":r, "Delta_sigma":round(float(Ds),4),
                     "Delta_pi":round(float(Dp),4), "ratio":round(float(ratio),3)})
        print(f"  {gs}:{gp}   {r:.2f}   {Ds:.4f}  {Dp:.4f}   {ratio:.3f}")
    out["twoband_scan"] = rows

    # 参考: MgB₂ 比
    print(f"\n参考 MgB₂: Δ_σ/Δ_π ≈ 7/1.7 ≈ 4.12")
    out["mgb2_ratio"] = 7.0/1.7

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_twoband.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_twoband.json")

if __name__ == "__main__":
    main()
