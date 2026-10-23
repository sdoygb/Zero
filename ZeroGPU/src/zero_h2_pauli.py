#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_h2_pauli.py —— H2 的 Pauli/交换能：单重态 vs 三重态（Z006）
================================================================================

上一块 `zero_h2_two_electron.py` 只取了空间对称基态（singlet），
Pauli/反对称没有显式测试。本块补这一条：

  两电子空间波函数：对称（singlet）／反对称（triplet）
  交换算符 P12： (P12 ψ)(i,j)=ψ(j,i)
  用 H ± λP12 把对称/反对称分开：
    H - λP12 的最低 = 空间对称（singlet）
    H + λP12 的最低 = 空间反对称（triplet；对称态被抬高 λ）
  交换能 J = E_triplet - E_singlet > 0（Pauli/交换）

舞台与势同 H2 两电子：n=11,h=0.5,ε=0.25，1/r 极限：
  H=T1+T2-Σ(1/r_iA+1/r_iB)+1/r12

判据：
  P1 E_triplet > E_singlet（Pauli）
  P2 J=E_triplet-E_singlet > 0
  P3 J(R=2.5) < J(R=2.0)（交换随 R 衰减）
  P4 H-λP12 的最低与 H 的基态一致（singlet 识别）
  P5 H+λP12 的最低本征矢反对称（||ψ+Pψ||/||ψ|| < 0.05）

运行：cd ZeroGPU/src && /usr/bin/python3 zero_h2_pauli.py
输出：results/zero_h2_pauli.json
"""
from __future__ import annotations

import json
import os
import sys

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
    print(f"H2 Pauli：n={n}, h={h}, N={N}, N²={N*N}, ε={eps}")
    print("=" * 84)

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

    rA = np.linalg.norm(X - X[A], axis=1)
    vA = -1.0 / soft(rA)
    lam = 0.5

    rows_out = []
    for R in (2.0, 2.5):
        B = idx[(c + int(round(R / h)), c, c)]
        rB = np.linalg.norm(X - X[B], axis=1)
        vB = -1.0 / soft(rB)
        D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
        V2 = (vA[:, None] + vB[:, None] + vA[None, :] + vB[None, :]
              + 1.0 / soft(D)).ravel()

        def Hx(x):
            Xm = x.reshape(N, N)
            return (T @ Xm + Xm @ T.T).ravel() + V2 * x

        def Px(x):
            return x.reshape(N, N).T.ravel()

        def mv_plus(x):
            return Hx(x) + lam * Px(x)

        def mv_minus(x):
            return Hx(x) - lam * Px(x)

        def solve(mv, k=1, vec=False):
            op = spl.LinearOperator((N * N, N * N), matvec=mv, dtype=float)
            if vec:
                w, v = spl.eigsh(op, k=k, which="SA", tol=1e-6, maxiter=400)
                return float(w[0]), v[:, 0]
            return float(spl.eigsh(op, k=k, which="SA", tol=1e-6, maxiter=400,
                                   return_eigenvectors=False)[0])

        E0 = solve(Hx)
        Ep, vp = solve(mv_plus, vec=True)
        Em = solve(mv_minus)
        Ep_true = Ep + lam
        Em_true = Em + lam
        J = Ep_true - E0
        # 反对称度
        asym = float(np.linalg.norm(vp + Px(vp)) / np.linalg.norm(vp))
        rows_out.append(dict(R=R, E_singlet=E0, E_triplet=Ep_true, J=J,
                             E_minus=Em_true, asym=asym))
        print(f"  R={R}: E_singlet={E0:+.6f}  E_triplet={Ep_true:+.6f}  "
              f"J={J:+.6f}  E_minus={Em_true:+.6f}  asym={asym:.3e}", flush=True)

    s2 = [r for r in rows_out if abs(r["R"] - 2.0) < 0.01][0]
    s25 = [r for r in rows_out if abs(r["R"] - 2.5) < 0.01][0]
    ck("P1 E_triplet > E_singlet（Pauli）", s2["E_triplet"] > s2["E_singlet"],
       f"{s2['E_triplet']:+.6f} > {s2['E_singlet']:+.6f}")
    ck("P2 交换能 J>0", s2["J"] > 0, f"J={s2['J']:+.6f}")
    ck("P3 J(R=2.5) < J(R=2.0)（交换衰减）", s25["J"] < s2["J"],
       f"{s25['J']:+.6f} < {s2['J']:+.6f}")
    ck("P4 H−λP 的最低与 H 基态一致（singlet 识别）",
       abs(s2["E_minus"] - s2["E_singlet"]) < 1e-6,
       f"{s2['E_minus']:+.6f} vs {s2['E_singlet']:+.6f}")
    ck("P5 H+λP 最低本征矢反对称（asym<0.05）",
       s2["asym"] < 0.05, f"asym={s2['asym']:.3e}")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    out = dict(n=n, h=h, eps=eps, lambda_penalty=lam, rows=rows_out,
               checks=CK, n_fail=nfail,
               boundary="1/r 极限 + 软化 ε；空间对称=singlet、反对称=triplet；粗网格")
    with open(os.path.join(RES, "zero_h2_pauli.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {os.path.join(RES, 'zero_h2_pauli.json')}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
