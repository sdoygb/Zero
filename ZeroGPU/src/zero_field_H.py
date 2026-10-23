#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_field_H.py — 把路线 A 的 Gauss 核接进净荷泛函：H 原子标定。

路线 A 的核（`zero_h2plus_routeA.py`）：
  H = (1/2a^2) L + diag(V),  V = -1/r（离源）, -2.3800756/a（在源格胞）
  这是 1/r 的立方胞平均正则化（ONSITE），H 原子给 -0.5002。
净荷泛函：势用单一项 phi=(1/a)G(S-N) 时，核函数必须换成同一 ONSITE；
本脚本先钉 H：N=1 时 Hartree 与 FA 交换精确抵消，V_eff=V_n，本征值即能量。

运行：/usr/local/bin/python3 ZeroGPU/src/zero_field_H.py
输出：results/zero_field_H.json
"""
from __future__ import annotations

import json
import pathlib
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "zero_field_H.json"
ONSITE = 2.3800756243710857   # 立方胞平均 1/r（a=1）
OK, BAD, ROWS = [], [], []


def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def lattice(n):
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
        nb = 0
        for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            j = idx.get((x + dx, y + dy, z + dz))
            if j is not None:
                nb += 1
                rows.append(i); cols.append(j); vals.append(-1.0)
        rows.append(i); cols.append(i); vals.append(float(nb))
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N)), idx, coords


def Vn(xyz, ps, a):
    V = np.zeros(len(xyz))
    for p in ps:
        r = np.linalg.norm(xyz - p, axis=1)
        on = r < 1e-9
        r = np.where(on, a, r)
        V += np.where(on, -ONSITE / a, -1.0 / r)
    return V


def run_H(a, L=14.0):
    t0 = time.time()
    n = int(round(L / a)) + 1
    Lm, idx, coords = lattice(n)
    xyz = np.array(coords, float) * a
    ctr = xyz[idx[(n // 2, n // 2, n // 2)]]
    H = 0.5 / a ** 2 * Lm + sp.diags(Vn(xyz, [ctr], a))
    E = float(spl.eigsh(H, k=1, which="SA")[0][0])
    return {"a": a, "n": n, "E": E, "err": E + 0.5, "time_s": round(time.time() - t0, 1)}


def main():
    print("== 路线 A 的 Gauss 核 → 净荷泛函：H 原子标定 ==")
    rows = []
    for a in (0.25, 0.20, 0.15):
        r = run_H(a)
        rows.append(r)
        print(f"  a={a:.2f} n={r['n']:3d}  E={r['E']:+.6f}  err={r['err']:+.2e}  ({r['time_s']}s)")
    ck("H1 钉到 -0.5（<1e-3）", all(abs(r["err"]) < 1e-3 for r in rows),
       f"max|err|={max(abs(r['err']) for r in rows):.2e}")
    ck("H2 随 a 减小收敛", all(abs(rows[i + 1]["err"]) < abs(rows[i]["err"]) for i in range(len(rows) - 1)),
       f"|err|: {abs(rows[0]['err']):.2e} -> {abs(rows[-1]['err']):.2e}")
    ck("H3 最细格到化学精度", abs(rows[-1]["err"]) < 1.6e-3,
       f"a={rows[-1]['a']} E={rows[-1]['E']:+.6f} err={rows[-1]['err']:+.2e}")

    out = {"onsite": ONSITE, "rows": rows, "checks": ROWS,
           "pass": len(OK), "fail": len(BAD),
           "note": "核=路线A的ONSITE立方胞平均；N=1时V_eff=V_n，本征值=能量"}
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2))
    print(f"\n== {len(OK)}/{len(OK)+len(BAD)} PASS  ->  {OUT.relative_to(ROOT.parent)} ==")
    if BAD:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
