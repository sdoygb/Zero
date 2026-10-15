#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_compete.py —— 两种介导竞争: 声子(s-wave) vs 磁(d-wave)
====================================================================================
总配对势 V = g²V_phonon + J²V_magnetic (声子吸引 + 反铁磁磁介导)
  V_phonon: 图上Laplacian声子介导 (近s, 各向同性)
  V_magnetic: 二部Néel序磁介导 (符号变化, d-wave)
竞争: 扫 g/J 相对强度 → 能隙对称性从 s-wave 转 d-wave
判定胜出: ①能隙符号结构(s全同号 / d符号交替) ②能隙幅度(胜出者给大Tc)
物理对应: 铜氧化物磁介导胜出(d-wave, 声子被抑制); BCS金属声子胜出(s-wave)
Zero增量: 两种介导竞争 → 配对对称性 + 主导介导由相对强度内生决定
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def make_grid(n):
    nr = nc = int(np.sqrt(n))
    adj = {}; coords = {}
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

def solve_gap(V, xi, tol=1e-8, maxit=8000):
    nv = V.shape[0]
    D = np.full(nv, 0.4)
    for it in range(maxit):
        E = np.sqrt(xi**2 + D**2)
        F = 0.5 * V @ (D / E)
        Dn = 0.4*D + 0.6*F
        if np.max(np.abs(Dn-D)) < tol: return Dn
        D = Dn
    return D

def main():
    out = {"config": {"model": "Competition: V=g²V_phonon(s-wave)+J²V_magnetic(d-wave)"}}
    print("="*90)
    print("两种介导竞争: 声子(s-wave) vs 磁(d-wave)")
    print("="*90)
    L, nv, coords = make_grid(64)
    ev, evc = np.linalg.eigh(L); om = np.sqrt(np.clip(ev,0,None))
    wD = float(om[-1]); xi = np.linspace(-wD, wD, nv)

    # V_phonon (g=1 归一): 去零模+正则化
    inv = 1.0/(om[1:]**2 + 0.2)
    Vp_base = evc[:,1:] @ np.diag(inv) @ evc[:,1:].T
    Vp_base /= np.max(np.abs(Vp_base))   # 归一
    # V_magnetic (J=1): Néel
    sigma = np.array([(-1)**(i+j) for (i,j) in coords.values()])
    Vm_base = -np.outer(sigma, sigma).astype(float); Vm_base /= np.max(np.abs(Vm_base))

    rows = []
    print(f"\nGrid8x8, ω_D={wD:.3f}. 扫 g(声子) vs J(磁):")
    print("  g     J    Δ_mean|Δ|   Δmin   Δmax   负顶点  对称性   主导")
    for g in [1.0]:
        for J in [0.0, 0.3, 0.5, 0.8, 1.0, 1.3, 1.6, 2.0]:
            V = g**2 * Vp_base + J**2 * Vm_base
            D = solve_gap(V, xi)
            nneg = int((D < -1e-6).sum())
            amean = float(np.mean(np.abs(D)))
            if nneg == 0:
                sym = "s-wave"; dom = "声子"
            else:
                sym = "d-wave(混)"; dom = "磁"
            rows.append({"g":g, "J":J, "Delta_abs_mean":round(amean,4),
                         "Delta_min":round(float(D.min()),4), "Delta_max":round(float(D.max()),4),
                         "n_neg":nneg, "sym":sym, "dominant":dom})
            print(f"  {g:.1f}  {J:.2f}   {amean:.4f}   {D.min():.4f}  {D.max():.4f}   {nneg:2d}   {sym:10s}  {dom}")
    out["competition_scan"] = rows

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_compete.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_compete.json")

if __name__ == "__main__":
    main()
