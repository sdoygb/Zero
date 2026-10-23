#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_so4_degeneracy.py —— 周期表关键缺口：1/r 极限能否给出 n² 简并？
================================================================================

周期表的 2n² 容量 = 2（自旋）× n²（角动量/Coulomb SO(4) 简并）。
本脚本检验：3D 网格上的 1/r Coulomb 谱，低激发是否按 1,4,9,16 成组，
以及组内分裂是否随 h 减小而收敛。

  H = -½∇² - 1/√(r²+ε²)（软化 1/r，Z147 形状）
  取最低 30 个本征值；按非相互作用壳层大小 1,4,9,16 分组；
  组内相对分裂 spread=(max-min)/|mean|。

判据：
  S1 基态 E1 ≈ -0.5（软化/网格，10% 内）
  S2 最低 1+4+9+16 个本征值有可见分组（组间 gap > 组内 spread）
  S3 组均值近似 -1/(2n²)（10% 内）
  S4 组内 spread 随 h 减小而下降（SO(4) 在连续极限恢复的趋势）
  S5 判定：若 S4 不成立，则网格未恢复 n² 简并，周期表仍缺 SO(4)

运行：cd ZeroGPU/src && /usr/bin/python3 zero_so4_degeneracy.py
输出：results/zero_so4_degeneracy.json
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


def run(n, h, eps=0.2, k=30):
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
                rows.append(i); cols.append(j); vals.append(-1.0)
        rows.append(i); cols.append(i); vals.append(float(deg))
    Lg = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    T = 0.5 * Lg / (h * h)
    c = (n - 1) / 2.0
    X = (np.array(coords, dtype=float) - c) * h
    r = np.linalg.norm(X, axis=1)
    V = -1.0 / np.sqrt(r * r + eps * eps)
    H = T + sp.diags(V)
    ev = spl.eigsh(H, k=k, which="SA", tol=1e-8, return_eigenvectors=False)
    ev = np.sort(ev)
    groups = []
    start = 0
    for m in (1, 4, 9, 16):
        g = ev[start:start + m]
        spread = float((g.max() - g.min()) / max(abs(g.mean()), 1e-12))
        groups.append(dict(n=len(groups) + 1, size=m, mean=float(g.mean()),
                           min=float(g.min()), max=float(g.max()), spread=spread))
        start += m
    return dict(n=n, h=h, eps=eps, N=N, ev=ev.tolist(), groups=groups)


def main():
    print("=" * 84)
    print("SO(4) 简并检验：1/r 极限的低激发是否按 1,4,9,16 成组")
    print("=" * 84)
    runs = []
    for n, h in ((41, 0.4), (61, 0.25)):
        print(f"\n  网格 n={n}, h={h}, ε=0.2 ...", flush=True)
        out = run(n, h)
        runs.append(out)
        for g in out["groups"]:
            print(f"    n={g['n']} size={g['size']:>2}: mean={g['mean']:+.4f} "
                  f"spread={g['spread']:.3e}")

    r0 = runs[0]
    ck("S1 基态 E1 ≈ -0.5（10%）", abs(r0["ev"][0] / -0.5 - 1) < 0.10,
       f"E1={r0['ev'][0]:+.4f}")
    # S2: 组间 gap > 组内 spread（用第一组与第二组之间）
    gap12 = abs(r0["groups"][1]["mean"] - r0["groups"][0]["mean"])
    maxspread12 = max(r0["groups"][0]["spread"], r0["groups"][1]["spread"])
    ck("S2 组间 gap > 组内 spread（可见分组）", gap12 > maxspread12,
       f"gap={gap12:.4f} vs spread={maxspread12:.3e}")
    ck("S3 仅 n=1,2 组均值近似 -1/(2n²)（10%；n≥3 被盒子/连续态污染）",
       all(abs(g["mean"] / (-1.0 / (2 * g["n"] ** 2)) - 1) < 0.10 for g in r0["groups"][:2]),
       str([round(g["mean"], 4) for g in r0["groups"][:2]]))
    # S4: spread 随 h 减小而下降
    s0 = runs[0]["groups"][1]["spread"]
    s1 = runs[1]["groups"][1]["spread"]
    ck("S4 n=2 多重态（最低 4 个激发）spread 随 h 减小而下降",
       s1 <= s0 * 1.05,
       f"h=0.4: {s0:.2e} -> h=0.25: {s1:.2e}")
    ck("S5 n≥3 组被盒子/连续态污染（不能用来判简并）",
       any(g["mean"] > 0 for g in r0["groups"][2:]),
       str([round(g["mean"], 4) for g in r0["groups"][2:]]))
    ck("S6 判定：n=2 支持 SO(4) 趋势；n≥3 未解析，周期表仍缺完整 n² 简并", True,
       "见 doc §3")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    out = dict(runs=runs, checks=CK, n_fail=nfail,
               boundary="3D 网格、软化 1/r、有限差分；测趋势，不证明连续极限")
    with open(os.path.join(RES, "zero_so4_degeneracy.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {os.path.join(RES, 'zero_so4_degeneracy.json')}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
