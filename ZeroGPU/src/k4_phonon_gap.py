#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_phonon_gap.py —— ①图上声子介导(s-wave) + ②图上磁介导/反铁磁(d-wave)
====================================================================================
配对介导机制决定配对对称性:
  ① 图上声子介导(Laplacian声子交换): V_kk'=g²Σφ_n(k)φ_n(k')/ω_n² (全正核)
     → s-wave 能隙 (全同号, 各向同性)  [BCS: Pb/Al/MgB₂ 声子超导]
  ② 图上磁介导(反铁磁自旋涨落, 二部图Néel序): V_kk'=J²·σ_k·σ_k'
     σ_k=±1 (二部Néel序) → 符号变化核
     → d-wave 能隙 (符号交替 + 节点过零)  [高温超导: 铜氧化物]
Zero 增量: 配对对称性由图上介导机制(声子/磁)内生决定, 非假设
验证: ①s-wave 全同号 ②d-wave 符号交替+节点 ③两种机制能隙对比
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def make_grid(n):
    nr = nc = int(np.sqrt(n))
    adj = {}
    coords = {}
    for i in range(nr):
        for j in range(nc):
            v = i*nc+j; adj[v]=[]; coords[v]=(i,j)
            for di,dj in [(1,0),(-1,0),(0,1),(0,-1)]:
                ni,nj=i+di,j+dj
                if 0<=ni<nr and 0<=nj<nc: adj[v].append(ni*nc+nj)
    nv = len(adj)
    L = np.zeros((nv, nv))
    for v in range(nv):
        L[v,v] = len(adj[v])
        for j in adj[v]: L[v,j] = -1
    return L, nv, coords

def solve_gap_vec(V, xi, tol=1e-8, maxit=8000):
    nv = V.shape[0]
    Delta = np.full(nv, 0.3)
    for it in range(maxit):
        E = np.sqrt(xi**2 + Delta**2)
        F = 0.5 * V @ (Delta / E)
        Dn = 0.4*Delta + 0.6*F
        if np.max(np.abs(Dn-Delta)) < tol:
            return Dn, it+1
        Delta = Dn
    return Delta, maxit

def main():
    out = {"config": {"model": "Pairing symmetry from Zero graph mediation mechanism",
                      "phonon": "V=g²Σφφ'/ω² → s-wave", "magnetic": "V=J²σ_kσ_k' Néel → d-wave"}}
    print("="*90)
    print("①图上声子介导(s-wave) + ②图上磁介导/反铁磁(d-wave)")
    print("="*90)
    for name, n in [("Grid4x4",16), ("Grid8x8",64)]:
        L, nv, coords = make_grid(n)
        evals, evecs = np.linalg.eigh(L)
        omegas = np.sqrt(np.clip(evals,0,None))
        wD = float(omegas[-1])
        xi = np.linspace(-wD, wD, nv)
        res = {"wD": round(wD,3)}

        # === ① 声子介导 (s-wave) ===
        inv = 1.0/(omegas[1:]**2 + 0.3)   # 去零模+正则化, 弱耦合
        Vp = (0.5**2) * evecs[:,1:] @ np.diag(inv) @ evecs[:,1:].T
        Dp, it1 = solve_gap_vec(Vp, xi)
        nneg_p = int((Dp < -1e-6).sum())
        sym_p = "s-wave(全同号)" if nneg_p==0 else "符号变化"
        res["phonon"] = {"gap_min": round(float(Dp.min()),5), "gap_max": round(float(Dp.max()),5),
                         "gap_mean": round(float(Dp.mean()),5), "n_neg": nneg_p, "sym": sym_p,
                         "V_offdiag_sign": "all>0" if Vp.min()>0 else "mixed"}
        print(f"\n[{name}] ω_D={wD:.3f}")
        print(f"  ①声子介导: 能隙min={Dp.min():.4f} max={Dp.max():.4f} mean={Dp.mean():.4f}, 负顶点={nneg_p} → {sym_p}")

        # === ② 磁介导 (反铁磁 Néel, d-wave) ===
        # σ_k = (-1)^(i+j) 二部Néel序; 磁介导势 V_m = -J²·σ_k·σ_k' (吸引在同子格, 排斥异子格)
        sigma = np.array([(-1)**(i+j) for (i,j) in coords.values()])
        J2 = 0.5
        Vm = -J2 * np.outer(sigma, sigma)
        Dm, it2 = solve_gap_vec(Vm, xi)
        nneg_m = int((Dm < -1e-6).sum())
        # 节点判定: 符号变化的边界(能隙近零)
        signed = np.sign(Dm)
        # 节点 = 相邻反号顶点间的连线数
        edges = [(a,b) for a in range(nv) for b in range(nv) if b>a and abs((coords[a][0]-coords[b][0])+(coords[a][1]-coords[b][1]))==1]
        n_nodes = sum(1 for a,b in edges if signed[a]!=signed[b])
        sym_m = "d-wave(符号交替)" if nneg_m>0 else "s-wave"
        res["magnetic"] = {"gap_min": round(float(Dm.min()),5), "gap_max": round(float(Dm.max()),5),
                           "gap_mean": round(float(Dm.mean()),5), "n_neg": nneg_m,
                           "n_Néel_edges_cross": n_nodes, "sym": sym_m,
                           "sigma_sample": [int(x) for x in sigma[:8]]}
        print(f"  ②磁介导: 能隙min={Dm.min():.4f} max={Dm.max():.4f} mean={Dm.mean():.4f}, 负顶点={nneg_m}")
        print(f"    跨Néel边节点数={n_nodes}, → {sym_m}")
        print(f"    Δ_k(前8)={np.round(Dm[:8],3).tolist()}")
        out.setdefault("results", {})[name] = res

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_phonon_gap.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_phonon_gap.json")

if __name__ == "__main__":
    main()
