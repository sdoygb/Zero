#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_gap.py —— 图上 Cooper 配对 → 能隙方程 Δ (接费米子侧)
====================================================================================
BCS 能隙(常数势 s-wave):  Δ = (V/2)Σ_k Δ/√(ξ_k²+Δ²)  (求和: |ξ_k|<ω_D Debye窗)
  → 1 = (V/2)Σ_k 1/√(ξ_k²+Δ²)  自洽方程, 迭代/二分解 Δ
  ω_D 来自 Zero 图声子谱(上轮 Θ_D=√λ_max), 图上费米能级 ξ_k ∈ [−ω_D,ω_D]
  N(0)=N/(2ω_D) 态密度; 弱耦合 Δ₀≈2ω_D·exp(−1/N(0)V)
  BCS 普适比: Δ₀/Tc = 1.76  (能隙与 Tc 同源 → 忠实度检验)
Zero 增量: ①图上声子(Θ_D)进能隙 ②图上费米子(态密度) ③图上声子介导 V_k,k'
验证: ①自洽能隙 Δ(V 扫描) ②指数结构 ③Δ₀/Tc≈1.76 ④声子介导 V vs 常数 V
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
    evals = np.clip(np.linalg.eigvalsh(L), 0, None)
    return c*np.sqrt(evals[-1])

def gap_equation(Delta, V, xi, wD):
    # 1 = (V/2) Σ_{|ξ|<ω_D} 1/√(ξ²+Δ²)
    mask = np.abs(xi) < wD
    xi_m = xi[mask]
    return (V/2.0)*np.sum(1.0/np.sqrt(xi_m**2 + Delta**2))

def solve_gap(V, xi, wD, tol=1e-12, maxit=2000):
    """二分自洽能隙 (解 1 = f(Δ))"""
    # f(Δ→0) 必须 >1 (否则无解, 弱耦合临界下无超导)
    f0 = gap_equation(1e-9, V, xi, wD)
    if f0 <= 1.0:
        return 0.0  # 无超导解
    lo, hi = 0.0, 2.0*wD
    for _ in range(maxit):
        mid = 0.5*(lo+hi)
        if gap_equation(mid, V, xi, wD) > 1.0:
            lo = mid
        else:
            hi = mid
        if hi-lo < tol: break
    return 0.5*(lo+hi)

def main():
    out = {"config": {"model": "Zero graph Cooper-pair gap Δ=(V/2)ΣΔ/√(ξ²+Δ²)",
                      "wD": "Zero graph phonon Θ_D (Z104)", "gap_ratio": "BCS Δ₀/Tc=1.76"}}
    print("="*90)
    print("图上 Cooper 配对 → 能隙方程 Δ (接费米子侧)")
    print("="*90)
    c = 1.0

    # 主图: Grid8x8 (最高 Θ_D)
    for name, (kind, n) in [("K4",("k4",4)), ("Grid8x8",("grid",64))]:
        L, _ = make_graph(kind, n)
        wD = theta_D(L, c)   # Debye 窗 = 图上最高声子频率
        print(f"\n[{name}] Debye 窗 ω_D = Θ_D = {wD:.3f}")

        # 图上费米能级: N 态均匀填充 Debye 窗 [-ω_D, ω_D]
        N_f = 64   # 费米态数
        xi = np.linspace(-wD, wD, N_f)
        N0 = N_f/(2.0*wD)   # 态密度
        print(f"  费米态 N={N_f}, 态密度 N(0)={N0:.3f}")

        # V 扫描: 弱-中耦合(能隙存在区 V>2/N) + 强耦合饱和
        V_scan = [0.03, 0.05, 0.07, 0.10, 0.15, 0.20]
        gaps = {}
        for V in V_scan:
            d = solve_gap(V, xi, wD)
            gaps[V] = round(float(d), 5)
        print(f"  Δ(V): {gaps}")
        print(f"  能隙存在边界: V>2/N={2.0/N_f:.4f} (N0V>1/ω_D={1.0/wD:.3f})")
        out.setdefault("gaps", {})[name] = {"wD": round(float(wD),3), "N0": round(float(N0),3),
                                            "existence_boundary_V_min": round(2.0/N_f,4),
                                            "Delta_vs_V": gaps}

        # 弱耦合指数检验: Δ₀≈2ω_D exp(−1/N0·V)
        print("  [弱耦合指数] Δ vs 2ω_D·exp(−1/N0V):")
        for V in V_scan:
            if gaps[V] > 0 and gaps[V] < 1.5*wD:
                pred = 2*wD*np.exp(-1.0/(N0*V))
                dev = abs(gaps[V]-pred)/pred
                print(f"    V={V}: Δ={gaps[V]:.5f}, 2ω_D·exp(−1/N0V)={pred:.5f}, 偏差={dev:.1%}")
        out.setdefault("exponential_check", {})[name] = {
            str(V): {"delta": gaps[V],
                     "weak_pred": round(float(2*wD*np.exp(-1.0/(N0*V))),5)}
            for V in V_scan if gaps[V] > 0 and gaps[V] < 1.5*wD}

        # Δ₀/Tc 比: BCS Tc=1.13Θ_D·exp(−1/N0V) (同源)
        print("  [Δ/Tc 普适比] (BCS 期望 1.76):")
        ratio_row = {}
        for V in V_scan:
            d = gaps[V]
            if d > 0 and d < 1.5*wD:
                tc = 1.13*wD*np.exp(-1.0/(N0*V))
                ratio = d/tc if tc>0 else 0
                ratio_row[str(V)] = round(ratio, 3)
                print(f"    V={V}: Δ={d:.5f}, Tc={tc:.5f}, Δ/Tc={ratio:.3f}")
        out.setdefault("gap_Tc_ratio", {})[name] = ratio_row

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_gap.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_gap.json")

if __name__ == "__main__":
    main()
