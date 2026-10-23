#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2plus_inflattice.py — H2+ 用无限格格林函数（Fourier），去掉有限尺寸/零和背景。

势：G_inf(r) = (1/(2π)^3) ∫ dk e^{ikr} / (2 Σ(1-cos k_i))
   —— 这是舞台图拉普拉斯的真正格林函数（无边界、无零和背景）。
动能：图拉普拉斯（正确量纲 (1/2a^2)L）；势：G_inf（正确量纲 g/a）。
对照：H 原子 -0.5 Ha；H2+ R=2.0 电子结合 -0.602634 Ha。
"""
from __future__ import annotations
import json, os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
PASS = []

def ck(name, cond, detail=""):
    PASS.append({"name": name, "pass": bool(cond), "detail": detail})
    print("[%s] %s %s" % ("PASS" if cond else "FAIL", name, detail))

def G_inf_table(rmax, M=300):
    """返回一个可调用：给定整数格点位移 (dx,dy,dz)，给 G_inf。用 FFT 预计算。"""
    k = (np.arange(M)+0.5)*2*np.pi/M - np.pi
    KX,KY,KZ = np.meshgrid(k,k,k,indexing='ij')
    denom = 2*(3-np.cos(KX)-np.cos(KY)-np.cos(KZ))
    denom[denom < 1e-12] = 1e-12
    tab = {}
    for dx in range(-rmax, rmax+1):
        for dy in range(-rmax, rmax+1):
            for dz in range(-rmax, rmax+1):
                if dx==dy==dz==0:
                    tab[(dx,dy,dz)] = 0.0
                    continue
                f = np.cos(dx*KX)*np.cos(dy*KY)*np.cos(dz*KZ)/denom
                tab[(dx,dy,dz)] = float(f.mean())
    return tab

def lattice(n):
    idx = {}; coords = []
    for x in range(n):
        for y in range(n):
            for z in range(n):
                idx[(x,y,z)] = len(coords); coords.append((x,y,z))
    N = len(coords)
    adj = [[] for _ in range(N)]
    for (x,y,z), i in idx.items():
        for dx,dy,dz in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)):
            j = idx.get((x+dx,y+dy,z+dz))
            if j is not None: adj[i].append(j)
    rows, cols, vals = [], [], []
    for i, nb in enumerate(adj):
        rows.append(i); cols.append(i); vals.append(float(len(nb)))
        for j in nb:
            rows.append(i); cols.append(j); vals.append(-1.0)
    return sp.csr_matrix((vals,(rows,cols)),shape=(N,N)), idx, coords

def main():
    print("=" * 84)
    print("zero_h2plus_inflattice.py — 无限格格林函数 H2+")
    print("=" * 84)
    Lphys = 8.0; R = 2.0
    rows = []
    for a in (0.5, 0.25, 0.125):
        n = int(round(Lphys/a)) + 1
        L, idx, coords = lattice(n)
        N = L.shape[0]
        rmax = int(round(Lphys/a)) + 2
        tab = G_inf_table(rmax, M=200)
        c = (n//2, n//2, n//2)
        # 势：单核在 c
        def make_V(c0):
            V = np.zeros(N)
            for i, cc in enumerate(coords):
                d = (cc[0]-c0[0], cc[1]-c0[1], cc[2]-c0[2])
                V[i] = -tab.get(d, 0.0)
            return V
        VH = make_V(c)
        H1 = (0.5/a**2)*L + (1.0/a)*sp.diags(VH)
        EH = float(spl.eigsh(H1, k=1, which="SA")[0][0])
        s = -0.5/EH
        off = int(round(R/(2*a)))
        cA = (c[0],c[1],c[2]-off); cB = (c[0],c[1],c[2]+off)
        VAB = make_V(cA) + make_V(cB)
        H2 = (0.5/a**2)*L + (1.0/a)*sp.diags(VAB)
        E2 = float(spl.eigsh(H2, k=1, which="SA")[0][0])
        bind = s*(E2-EH)
        rows.append({"a":a,"n":n,"N":N,"EH":EH,"s":s,"E2":E2,"bind":bind,"err":bind+0.602634})
        print(f"  a={a:.3f} n={n:3d}  EH={EH:+.6e}  s={s:8.1f}  E2={E2:+.6e}  "
              f"锚定结合={bind:+.4f}  误差={bind+0.602634:+.4f}")
    ck("F1 H 原子锚定 = -0.5", all(abs(r['s']*r['EH']+0.5)<1e-9 for r in rows), "s=-0.5/EH")
    ck("F2 H 原子 |E0| 随 a 减小而增大", all(abs(rows[i+1]['EH'])>abs(rows[i]['EH']) for i in range(len(rows)-1)),
       f"|EH|: {abs(rows[0]['EH']):.3e} -> {abs(rows[-1]['EH']):.3e}")
    ck("F3 H2+ 结合 < 0", all(r['bind']<0 for r in rows), f"{rows[-1]['bind']:+.4f}")
    ck("F4 误差随 a 减小", all(abs(rows[i+1]['err'])<abs(rows[i]['err']) for i in range(len(rows)-1)),
       f"|err|: {abs(rows[0]['err']):.4f} -> {abs(rows[-1]['err']):.4f}")
    ck("F5 误差 < 0.05", all(abs(r['err'])<0.05 for r in rows), f"min={min(abs(r['err']) for r in rows):.4f}")
    ck("F6 误差 < 1.6e-3", all(abs(r['err'])<1.6e-3 for r in rows), f"min={min(abs(r['err']) for r in rows):.4f}")
    out = {"rows": rows, "checks": PASS}
    with open(os.path.join(RES,"zero_h2plus_inflattice.json"),"w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n%d 项，%d PASS" % (len(PASS), sum(1 for p in PASS if p["pass"])))
    print("-> " + os.path.join(RES,"zero_h2plus_inflattice.json"))

if __name__ == "__main__":
    main()
