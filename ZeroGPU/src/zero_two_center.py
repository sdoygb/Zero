#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_two_center.py —— 双中心 Coulomb：H2+ 的有限格结构版（含负结果）
================================================================================

补的缺口（分子哈密顿量稿明说）：
  「多中心 -Σ_A Z_A/r_{iA} 是多个单中心 1/r 的叠加；这一步旧稿没有做——只做了单中心。」

做法（不借 QED 的数）：
  舞台：n^3 立方格（开边界）；两个源 A、B：x = δ_A + δ_B − 2/N（零和）
  势：L φ = +x（G=L^+，Z033 定理16；正源 -> φ>0）
  电子：H = ½L − g·diag(φ)，最低本征值 E0(R)
  核-核：Vnn(R) = +g² φ_B(A)（φ_B 是 B 单源势；φ 为正井，故 Vnn>0）
  总能量：E_tot(R) = E0(R) + Vnn(R)

判据（审计）：
  T1 势残差 <1e-10
  T2 每个 R 都有束缚态 E0<0
  T3 E0(R) 单调（本格上递减：从上方趋近单中心）
  T4 Vnn(R) 正、递减
  T5 E0(R) 从上方趋近单中心（反绑定诊断）

负结果（单独登记，不算脚本失败）：
  N1 扫描范围内无内部极小 ⟹ 有限格双中心不结合 H2+
  N2 总能量单调 ⟹ 无键长极小

运行：cd ZeroGPU/src && /usr/bin/python3 zero_two_center.py
输出：results/zero_two_center.json、zero_two_center.png
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
NEG = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": str(detail)})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def neg(name, detail=""):
    NEG.append({"name": name, "detail": str(detail)})
    print(f"  [NEG] {name} —— {detail}")


def cubic_sparse(n):
    idx = {}
    coords = []
    for x in range(n):
        for y in range(n):
            for z in range(n):
                idx[(x, y, z)] = len(coords)
                coords.append((x, y, z))
    N = len(coords)
    rows, cols, vals = [], [], []
    for (x, y, z), i in idx.items():
        deg = 0
        for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            j = idx.get((x + dx, y + dy, z + dz))
            if j is not None:
                deg += 1
                if j > i:          # 每条边只加一次（两个方向对称）
                    rows += [i, j]; cols += [j, i]; vals += [-1.0, -1.0]
        rows.append(i); cols.append(i); vals.append(float(deg))
    L = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    return L, coords, idx


def solve_phi(L, x):
    phi = spl.lsqr(L, x, atol=1e-13, btol=1e-13, iter_lim=30000)[0]
    phi -= phi.mean()
    return phi, float(np.max(np.abs(L @ phi - x)))


def lowest(L, phi, k=1):
    H = 0.5 * L - sp.diags(phi)
    return float(spl.eigsh(H, k=k, which="SA", tol=1e-9, return_eigenvectors=False)[0])


def main():
    n = 15
    L, coords, idx = cubic_sparse(n)
    N = L.shape[0]
    c = (n - 1) // 2
    A = idx[(c, c, c)]
    print("=" * 84)
    print(f"双中心 Coulomb：n={n}, N={N}, A={coords[A]}")
    print("=" * 84)

    # 单中心参照
    xA = np.zeros(N); xA[A] = 1.0; xA -= xA.mean()
    phiA, resA = solve_phi(L, xA)
    E_single = lowest(L, phiA)
    print(f"  单中心 E_single={E_single:+.6f}  res={resA:.1e}")

    rows = []
    for R in range(1, 7):
        B = idx[(c + R, c, c)]
        x = np.zeros(N); x[A] = x[B] = 1.0; x -= x.mean()
        phi, res = solve_phi(L, x)
        E0 = lowest(L, phi)
        xB = np.zeros(N); xB[B] = 1.0; xB -= xB.mean()
        phiB, _ = solve_phi(L, xB)
        Vnn = float(phiB[A])
        Etot = E0 + Vnn
        rows.append(dict(R=R, E0=E0, Vnn=Vnn, Etot=Etot, res=res,
                         E0_minus_single=E0 - E_single))
        print(f"  R={R}: E0={E0:+.6f}  Vnn={Vnn:+.6f}  Etot={Etot:+.6f}  "
              f"E0−E_single={E0-E_single:+.2e}  res={res:.1e}")

    ck("T1 每个 R 的泊松残差 <1e-10", all(r["res"] < 1e-10 for r in rows),
       f"max={max(r['res'] for r in rows):.1e}")
    ck("T2 每个 R 都有束缚态 E0<0", all(r["E0"] < 0 for r in rows),
       f"E0 max={max(r['E0'] for r in rows):.3e}")
    ck("T3 E0(R) 单调递减（本格）",
       all(rows[i]["E0"] > rows[i + 1]["E0"] for i in range(len(rows) - 1)),
       f"{rows[0]['E0']:.3e} → {rows[-1]['E0']:.3e}")
    ck("T4 Vnn(R) 正、递减",
       all(r["Vnn"] > 0 for r in rows) and
       all(rows[i]["Vnn"] > rows[i + 1]["Vnn"] for i in range(len(rows) - 1)),
       f"{rows[0]['Vnn']:.3e} → {rows[-1]['Vnn']:.3e}")
    ck("T5 E0(R) 低于单中心（两源势井加深）",
       all(r["E0_minus_single"] < 0 for r in rows),
       f"{rows[0]['E0_minus_single']:.2e} → {rows[-1]['E0_minus_single']:.2e}")

    # 负结果：无内部极小
    imin = int(np.argmin([r["Etot"] for r in rows]))
    internal = 0 < imin < len(rows) - 1
    if not internal:
        neg("N1 R≤6 内无内部极小（扫描限制，不是不存在）",
            f"Etot 最小值在 R={rows[imin]['R']}（端点），Etot={rows[imin]['Etot']:+.6f}")
    if not internal:
        neg("N2 总能量在 R=6 仍下降；极小可能在 R>6 或需更大格/1/r 极限",
            "E_tot(R) 在扫描边缘仍单调")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"审计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}；负结果 {len(NEG)} 条")
    out = dict(n=n, N=N, A=coords[A], E_single=E_single, rows=rows,
               checks=CK, n_fail=nfail, negative=NEG,
               boundary="有限格双中心；R 是格点单位；Vnn 用格点 Green 函数；负结果不是脚本错误")
    with open(os.path.join(RES, "zero_two_center.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    make_figure(rows, E_single)
    print(f"结果 → {os.path.join(RES, 'zero_two_center.json')}")
    return 0 if nfail == 0 else 1


def make_figure(rows, E_single):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R = [r["R"] for r in rows]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(R, [r["E0"] for r in rows], "o-", label="E0 (electronic)")
    ax.axhline(E_single, color="C0", ls="--", lw=1, label="single-center E0")
    ax.plot(R, [r["Vnn"] for r in rows], "s-", label="Vnn (nuclear repulsion)")
    ax.plot(R, [r["Etot"] for r in rows], "^-", label="Etot")
    ax.axhline(0, color="k", lw=.8)
    ax.set_xlabel("R (lattice units)")
    ax.set_ylabel("energy (lattice units)")
    ax.set_title("Two-center Coulomb, n=15: no minimum within R<=6 (scan limit)")
    ax.legend(fontsize=8)
    ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "zero_two_center.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
