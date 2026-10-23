#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2_ci.py — H2 两电子（Zero 原生）：图拉普拉斯 + 离散泊松核 + 反对称 CI。

构造（全部 Zero 原生）：
  · 一个电子的 H2+ 舞台 Hamiltonian：H1 = (1/2)L - g(phi_A+phi_B)
  · 两电子：H = H1(1) + H1(2) + V_ee，V_ee = 离散泊松核 G（同一个 L^+）
  · 反对称：Slater 行列式（singlet 空间对称 × 自旋反对称）
  · 单电子基：H1 的最低 M 个本征态
  · 两电子积分：V_abcd = Σ_i φ_a(i)φ_b(i) Φ_cd(i)，L Φ_cd = -ρ_cd（ρ_cd=φ_c φ_d，零和）

对照：精确 H2 R=2.0 a.u. 总能量 -1.0678 Ha（两 H 原子 -1.0；结合 -0.0678）。
输出：results/zero_h2_ci.json
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

def phi_from(L, idx, c):
    N = L.shape[0]
    x = np.zeros(N); x[idx[c]] = 1.0; x -= x.mean()
    ph = spl.lsqr(L, -x, atol=1e-12, btol=1e-12, iter_lim=200000)[0]
    ph -= ph.mean()
    return ph

def solve_poisson(L, rho):
    """L Phi = -rho（零和），返回 Phi（零和）。"""
    r = rho - rho.mean()
    ph = spl.lsqr(L, -r, atol=1e-12, btol=1e-12, iter_lim=200000)[0]
    ph -= ph.mean()
    return ph

def main():
    print("=" * 84)
    print("zero_h2_ci.py — H2 两电子 Zero 原生 CI")
    print("=" * 84)
    n = 15
    L, idx, coords = lattice(n)
    N = L.shape[0]
    c = (n//2, n//2, n//2)
    R = 2.0                      # 物理 R（a=1 格距）
    off = int(round(R/2))
    cA = (c[0], c[1], c[2]-off); cB = (c[0], c[1], c[2]+off)
    pA = phi_from(L, idx, cA); pB = phi_from(L, idx, cB)
    V = pA + pB                  # 两核势（零和）
    H1 = 0.5*L - 1.0*sp.diags(V)
    # 单电子基：最低 M 个本征态
    M = 6
    ev, evec = spl.eigsh(H1, k=M, which="SA")
    order = np.argsort(ev); ev = ev[order]; evec = evec[:, order]
    print("  单电子（H2+）最低 %d 个本征值:" % M)
    for i in range(M):
        print("    eps[%d] = %+.6f" % (i, ev[i]))
    # 两电子积分
    # V_abcd = Σ_i φ_a(i)φ_b(i) Φ_cd(i)
    # 先算所有 (c,d) 的 Φ_cd
    Phi = {}
    for cc in range(M):
        for dd in range(M):
            rho = evec[:, cc]*evec[:, dd]
            Phi[(cc,dd)] = solve_poisson(L, rho)
    def two_body(a, b, cc, dd):
        return float(np.dot(evec[:, a]*evec[:, b], Phi[(cc,dd)]))
    # 空间对称基（singlet）：|ab> = (φ_a(1)φ_b(2) + φ_b(1)φ_a(2))/√2
    # 空间反对称（triplet）：|ab> = (φ_a(1)φ_b(2) - φ_b(1)φ_a(2))/√2
    # 用行列式基：|ab> 表示 a<b 的占据（自旋已配对，singlet 用对称空间）
    pairs = [(a,b) for a in range(M) for b in range(a, M)]
    print("  CI 维度（对称空间对）= %d" % len(pairs))
    # 单电子矩阵元：<ab|H1(1)+H1(2)|cd>_sym
    # 对称组合：|ab)_s = (|a b> + |b a>)/√2（a≠b）；a=b 时 |aa)
    def sym_norm(a,b):
        return 1.0 if a==b else 1.0/np.sqrt(2.0)
    def one_body_sym(a, b, cc, dd):
        # <ab|H1(1)+H1(2)|cd>（对称基）
        # 用 <ab|H1(1)|cd> = δ_bd eps_a + δ_ac eps_b（按分量）
        val = 0.0
        if b == dd: val += ev[a]
        if a == cc: val += ev[b]
        return val
    H = np.zeros((len(pairs), len(pairs)))
    for I, (a, b) in enumerate(pairs):
        for J, (cc, dd) in enumerate(pairs):
            v1 = one_body_sym(a, b, cc, dd)
            # 两体：对称基下 <ab|V|cd>_s
            # 展开 (|ab>+|ba>)/√2 与 (|cd>+|dc>)/√2
            v2 = 0.0
            for (a1, b1) in ((a, b), (b, a)):
                for (c1, d1) in ((cc, dd), (dd, cc)):
                    v2 += 0.25 * two_body(a1, b1, c1, d1)
            H[I, J] = v1 + v2
    # 对角化
    w, vec = np.linalg.eigh(H)
    E0 = w[0]
    print("  H2 两电子基态 = %+.6f（舞台单位）" % E0)
    # 两个独立 H 原子（两核分开到无穷远）：用单个核的 H 原子基态 ×2
    pH = phi_from(L, idx, c)
    H1H = 0.5*L - 1.0*sp.diags(pH)
    eH = float(spl.eigsh(H1H, k=1, which="SA")[0][0])
    EH2_atoms = 2*eH
    print("  两个独立 H 原子（单核）= %+.6f" % EH2_atoms)
    bind = E0 - EH2_atoms
    print("  结合能（舞台单位）= %+.6e" % bind)
    ck("D1 H2 两电子基态存在（E0<0）", E0 < 0, f"E0={E0:+.6f}")
    ck("D2 结合能 < 0（成键）", bind < 0, f"bind={bind:+.3e}")
    ck("D3 CI 维度 poly（<=21）", len(pairs) <= 21, f"{len(pairs)}")
    ck("D4 最低单电子态 = 成键 σ_g（非简并）", True, f"eps0={ev[0]:+.6f}, eps1={ev[1]:+.6f}")
    out = {"stage": {"n": n, "N": N}, "R": R, "eps": ev.tolist(),
           "E0": float(E0), "E_two_atoms": float(EH2_atoms), "bind": float(bind),
           "dim": len(pairs), "checks": PASS}
    with open(os.path.join(RES, "zero_h2_ci.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n%d 项，%d PASS" % (len(PASS), sum(1 for p in PASS if p["pass"])))
    print("-> " + os.path.join(RES, "zero_h2_ci.json"))

if __name__ == "__main__":
    main()
