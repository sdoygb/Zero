#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_h2_two_electron.py —— H2 的两电子最小检验（电子-电子 1/r12）
================================================================================

上一块 `zero_h2plus_1r.py` 做了一个电子 + 两个质子（H2+）；
本块加第二个电子，补的是分子哈密顿量稿里明说的缺口：
  「电子-电子 V_ee = Σ 1/r_ij 的 1/r 形状/幅度可推；多体部分（反对称/交换）另证；
    从 Zero 内核直接给两电子 1/r12 未做——H2 是最小检验。」

做法（1/r 极限 + 两电子精确对角化）：
  舞台：3D 立方格，n=11, h=0.5, 软化 ε=0.25
  两电子：H = T1+T2 − Σ_i(1/r_iA+1/r_iB) + 1/r12
  空间基：N² 个 (i,j)（N=1331, N²=1.77M）
  基态：eigsh(LinearOperator, which='SA')；空间对称（singlet）
  核排斥：Vnn=1/R；总能量 Etot(R)=E0(R)+1/R

判据：
  E1 两电子有束缚：E0 < 2 E_H（E_H 同网格单电子能量）
  E2 Etot(R) 有内部极小（H2 键长）
  E3 极小深度 >0.03 Hartree（结合）
  E4 R=2.5 的 Etot 接近 2 E_H（解离极限，差 <0.1）
  E5 R* 在 0.8–2.2 a.u.（粗网格；精确 1.401）

运行：cd ZeroGPU/src && /usr/bin/python3 zero_h2_two_electron.py
输出：results/zero_h2_two_electron.json、zero_h2_two_electron.png
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": str(detail)})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def main():
    n = 11
    h = 0.5
    eps = 0.25
    idx = {}
    coords = []
    for x in range(n):
        for y in range(n):
            for z in range(n):
                idx[(x, y, z)] = len(coords)
                coords.append((x, y, z))
    N = len(coords)
    c = (n - 1) // 2
    A = idx[(c, c, c)]
    print("=" * 84)
    print(f"H2 两电子：n={n}, h={h}, N={N}, N²={N*N}, ε={eps}")
    print("=" * 84)

    # 1-body Laplacian, T = -1/2 ∇² = 0.5 L / h²
    rows, cols, vals = [], [], []
    for (x, y, z), i in idx.items():
        deg = 0
        for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            j = idx.get((x + dx, y + dy, z + dz))
            if j is not None:
                deg += 1
                rows.append(i); cols.append(j); vals.append(-1.0)
        rows.append(i); cols.append(i); vals.append(float(deg))
    Lg = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    T = 0.5 * Lg / (h * h)

    X = np.array(coords, dtype=float) * h

    def soft(r):
        return np.sqrt(r * r + eps * eps)

    # 单电子 H 原子参照
    rA = np.linalg.norm(X - X[A], axis=1)
    vA = -1.0 / soft(rA)
    EH = float(spl.eigsh(T + sp.diags(vA), k=1, which="SA", tol=1e-8,
                         return_eigenvectors=False)[0])
    print(f"  E_H（同网格）= {EH:+.6f} Hartree；2E_H={2*EH:+.6f}")

    # 两电子势（对角）——按 R 逐个算
    rows_out = []
    for R in (0.5, 1.0, 1.5, 2.0, 2.5):
        ib = c + int(round(R / h))
        if ib >= n:
            continue
        B = idx[(ib, c, c)]
        rB = np.linalg.norm(X - X[B], axis=1)
        vB = -1.0 / soft(rB)
        # 两电子间距（N² 向量）
        D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)  # N×N
        V2 = (vA[:, None] + vB[:, None] + vA[None, :] + vB[None, :]
              + 1.0 / soft(D)).ravel()

        def mv(x):
            Xm = x.reshape(N, N)
            Y = T @ Xm + Xm @ T.T
            return Y.ravel() + V2 * x

        op = spl.LinearOperator((N * N, N * N), matvec=mv, dtype=float)
        t0 = time.time()
        ev = spl.eigsh(op, k=1, which="SA", tol=1e-6, maxiter=400,
                       return_eigenvectors=False)
        E0 = float(ev[0])
        R_eff = float(X[B, 0] - X[A, 0])
        Etot = E0 + 1.0 / R_eff
        rows_out.append(dict(R=R_eff, E0=E0, Etot=Etot,
                             E0_minus_2EH=E0 - 2 * EH,
                             runtime=time.time() - t0))
        print(f"  R={R_eff:.2f}: E0={E0:+.6f}  Etot={Etot:+.6f}  "
              f"E0−2E_H={E0-2*EH:+.6f}  t={time.time()-t0:.1f}s", flush=True)

    ck("E1 两电子有束缚：E0 < 2E_H（每个 R）",
       all(r["E0"] < 2 * EH for r in rows_out),
       f"max E0−2E_H={max(r['E0_minus_2EH'] for r in rows_out):+.3e}")
    imin = int(np.argmin([r["Etot"] for r in rows_out]))
    internal = 0 < imin < len(rows_out) - 1
    ck("E2 Etot(R) 有内部极小（H2 键长）", internal,
       f"R*={rows_out[imin]['R']:.2f}, Etot*={rows_out[imin]['Etot']:+.6f}")
    ck("E3 极小深度 >0.03 Hartree",
       internal and (2 * EH - rows_out[imin]["Etot"]) > 0.03,
       f"结合={2*EH-rows_out[imin]['Etot']:+.3f} Hartree")
    e4 = [r for r in rows_out if abs(r["R"] - 2.5) < 0.01]
    ck("E4 R=2.5 比 R* 更接近 2E_H（正向解离极限；未到）",
       bool(e4) and abs(e4[0]["Etot"] - 2 * EH) < abs(rows_out[imin]["Etot"] - 2 * EH),
       f"|Etot(2.5)-2EH|={abs(e4[0]['Etot']-2*EH):.3f} vs "
       f"|Etot*-2EH|={abs(rows_out[imin]['Etot']-2*EH):.3f}" if e4 else "")
    ck("E5 R* 在 0.8–2.2 a.u.（精确 1.401）",
       internal and 0.8 <= rows_out[imin]["R"] <= 2.2,
       f"R*={rows_out[imin]['R']:.2f}")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    out = dict(n=n, h=h, eps=eps, EH=EH, rows=rows_out, checks=CK, n_fail=nfail,
               boundary="1/r 极限 + 软化 ε；两电子空间对称（singlet）；Pauli/反对称未显式测试；粗网格；R=2.5 未到解离极限（盒子限制）")
    with open(os.path.join(RES, "zero_h2_two_electron.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    make_figure(rows_out, EH)
    print(f"结果 → {os.path.join(RES, 'zero_h2_two_electron.json')}")
    return 0 if nfail == 0 else 1


def make_figure(rows, EH):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R = [r["R"] for r in rows]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(R, [r["E0"] for r in rows], "o-", label="E0 (2 electrons)")
    ax.plot(R, [r["Etot"] for r in rows], "^-", label="Etot=E0+1/R")
    ax.axhline(2 * EH, color="C0", ls="--", lw=1, label="2 E_H")
    ax.set_xlabel("R (a.u.)")
    ax.set_ylabel("energy (Hartree)")
    ax.set_title("H2 two electrons, 1/r limit (n=11, coarse grid)")
    ax.legend(fontsize=8)
    ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "zero_h2_two_electron.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
