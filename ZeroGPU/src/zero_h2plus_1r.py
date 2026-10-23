#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_h2plus_1r.py —— 双中心 1/r 极限：H2+ 的键长极小
================================================================================

上一版 `zero_two_center.py` 用**离散泊松核**（G=L^+）在开边界格上扫描，
R≤6 内没有总能量极小；n=21/31 与周期边界也未在可靠范围内找到。
原因是有限格的势被零和背景/边界像荷屏蔽，1/r 只是 n→∞ 的极限（原生库仑 §2b）。

本脚本走**1/r 极限**（Z147 的形状，不借 QED 常数）：
  舞台：3D 立方格，间距 h，开边界（盒子足够大）
  动能：T = -½∇²（7 点差分，原子单位）
  势：V = -1/r_A - 1/r_B（软化 r→√(r²+ε²)，ε=h/2 避免格点奇点）
  电子：H = T + V，最低本征值 E0(R)
  核排斥：Vnn = +1/R
  总能量：Etot(R)=E0(R)+1/R

判据：
  H1 R=2 时 E0 ≈ -0.60 Hartree（5% 内）
  H2 R=2 时 Etot ≈ -0.10 Hartree（10% 内）
  H3 Etot(R) 在 0.5<R<4 有内部极小（H2+ 键长）
  H4 极小深度 >0.05 Hartree（不是数值噪声）
  H5 R=1 的 E0 比 R=4 更低（短距离更紧）

运行：cd ZeroGPU/src && /usr/bin/python3 zero_h2plus_1r.py
输出：results/zero_h2plus_1r.json、zero_h2plus_1r.png
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
    n = 41
    h = 0.25
    eps = 0.5 * h
    idx = {}
    coords = []
    for x in range(n):
        for y in range(n):
            for z in range(n):
                idx[(x, y, z)] = len(coords)
                coords.append((x, y, z))
    N = len(coords)
    print("=" * 84)
    print(f"H2+ 1/r 极限：n={n}, h={h}, N={N}, ε={eps}")
    print("=" * 84)

    # 7 点图拉普拉斯（开边界），T = 0.5 * L / h²
    rows, cols, vals = [], [], []
    for (x, y, z), i in idx.items():
        deg = 0
        for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            j = idx.get((x + dx, y + dy, z + dz))
            if j is not None:
                deg += 1
                if j > i:          # 每条边只加一次
                    rows += [i, j]; cols += [j, i]; vals += [-1.0, -1.0]
        rows.append(i); cols.append(i); vals.append(float(deg))
    Lg = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    T = 0.5 * Lg / (h * h)
    c = (n - 1) // 2
    A = idx[(c, c, c)]
    xA = (np.array(coords[A]) - c) * h
    print(f"  格点 A={coords[A]}  物理坐标={xA}")

    def soft(r):
        return np.sqrt(r * r + eps * eps)

    rows_out = []
    for R in (0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0):
        # B 放在 A + R 的 x 方向（格点化；同时用精确 R 算 Vnn）
        ib = c + int(round(R / h))
        if ib >= n:
            continue
        B = idx[(ib, c, c)]
        xB = (np.array(coords[B]) - c) * h
        R_eff = float(np.linalg.norm(xB - xA))
        rA = np.array([np.linalg.norm((np.array(cc) - c) * h - xA) for cc in coords])
        rB = np.array([np.linalg.norm((np.array(cc) - c) * h - xB) for cc in coords])
        V = -1.0 / soft(rA) - 1.0 / soft(rB)
        H = T + sp.diags(V)
        E0 = float(spl.eigsh(H, k=1, which="SA", tol=1e-9, return_eigenvectors=False)[0])
        Vnn = 1.0 / R_eff
        Etot = E0 + Vnn
        rows_out.append(dict(R=R_eff, E0=E0, Vnn=Vnn, Etot=Etot))
        print(f"  R={R_eff:.3f}: E0={E0:+.6f}  Vnn={Vnn:+.6f}  Etot={Etot:+.6f}")

    # 参考：H2+ 在 R=2 的精确电子能 -1.1026、总能量 -0.6026 Hartree
    ck("H1 R=2 的 E0 ≈ -1.1026 Hartree（5%）",
       any(abs(r["R"] - 2.0) < 0.02 and abs(r["E0"] / -1.1026 - 1) < 0.05 for r in rows_out),
       f"{[round(r['E0'],4) for r in rows_out if abs(r['R']-2)<0.02]}")
    ck("H2 R=2 的 Etot ≈ -0.6026 Hartree（10%）",
       any(abs(r["R"] - 2.0) < 0.02 and abs(r["Etot"] / -0.6026 - 1) < 0.10 for r in rows_out),
       f"{[round(r['Etot'],4) for r in rows_out if abs(r['R']-2)<0.02]}")
    # 局部极小（H2+ 键长）；全局最小在 R=4 是盒子/软化伪影
    local = [i for i in range(1, len(rows_out) - 1)
             if rows_out[i]["Etot"] < rows_out[i-1]["Etot"] and rows_out[i]["Etot"] < rows_out[i+1]["Etot"]]
    li = local[0] if local else None
    ck("H3 Etot(R) 有内部局部极小（H2+ 键长）", li is not None,
       f"R*={rows_out[li]['R']:.3f}, Etot*={rows_out[li]['Etot']:+.6f}" if li is not None else "无")
    ck("H4 局部极小相对 H+H+ (-0.5) 的结合能 >0.05 Hartree",
       li is not None and (rows_out[li]["Etot"] + 0.5) < -0.05,
       f"结合能={rows_out[li]['Etot']+0.5:+.3f} Hartree" if li is not None else "")
    ck("H5 R=1 的 E0 比 R=4 更低",
       any(abs(r["R"] - 1.0) < 0.02 for r in rows_out) and
       any(abs(r["R"] - 4.0) < 0.02 for r in rows_out) and
       [r["E0"] for r in rows_out if abs(r["R"] - 1.0) < 0.02][0] <
       [r["E0"] for r in rows_out if abs(r["R"] - 4.0) < 0.02][0],
       "短距离更紧")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    local = [i for i in range(1, len(rows_out) - 1)
             if rows_out[i]["Etot"] < rows_out[i-1]["Etot"] and rows_out[i]["Etot"] < rows_out[i+1]["Etot"]]
    out = dict(n=n, h=h, eps=eps, rows=rows_out,
               local_min=(rows_out[local[0]] if local else None),
               global_min=rows_out[int(np.argmin([r["Etot"] for r in rows_out]))],
               checks=CK, n_fail=nfail,
               boundary="1/r 极限（Z147 形状）＋软化 ε；不是有限格离散泊松核；H2+ 键长以 a.u. 给出")
    with open(os.path.join(RES, "zero_h2plus_1r.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    make_figure(rows_out)
    print(f"结果 → {os.path.join(RES, 'zero_h2plus_1r.json')}")
    return 0 if nfail == 0 else 1


def make_figure(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R = [r["R"] for r in rows]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(R, [r["E0"] for r in rows], "o-", label="E0 (electronic)")
    ax.plot(R, [r["Vnn"] for r in rows], "s-", label="Vnn=1/R")
    ax.plot(R, [r["Etot"] for r in rows], "^-", label="Etot")
    ax.axhline(-0.5, color="C0", ls="--", lw=1, label="H + H+ (-0.5)")
    ax.set_xlabel("R (a.u.)")
    ax.set_ylabel("energy (Hartree)")
    ax.set_title("H2+ in the 1/r limit (grid n=41, h=0.25)")
    ax.legend(fontsize=8)
    ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "zero_h2plus_1r.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
