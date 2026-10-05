#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R55 · 量子线：Z* 下投影 + 有向面积倾斜（打破反射对称）

构造（按用户的图像）：
  (1) Z* 的成员 = 长度 L 的闭合 ±1 路径（Σw = 0）
  (2) 下投影    = 取每条闭合路径的**全部前缀**，前缀 σ 是开路径（演化层状态）
  (3) 原生权重  = 挂在 σ 上的闭合路径数（补全数）
  (4) 打破对称  = 用**有向面积** A 加权：每条完整闭合路径权重 ∝ exp(-β A)
                  A = (1/2) Σ_k S_k   （S_k = 部分和；闭曲线下的几何面积）
                  β>0 压低正面积，β<0 压低负面积 —— 这是唯一的漂移项

判据：量子线 S_max(ω) > 2，其中 ω 是前缀类上的诱导测度。

实现：精确整数 DP 算补全的面积生成函数，再对 β 扫描。
"""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
MU = (math.sqrt(5), (5 - math.sqrt(5)) / 2, (5 - math.sqrt(5)) / 2)
OUT = {}


def smax(ws):
    """ws: 降序权重（已归一或未归一皆可）→ S_max = <顶三归一, mu>"""
    v = sorted(ws, reverse=True)[:3]
    s = sum(v)
    if s <= 0:
        return 0.0
    v = [x / s for x in v]
    return sum(a * b for a, b in zip(v, MU))


def area_gf(maxsteps, cutoff):
    """g[n][ds] = Counter{2A : 条数}，n 步、净位移 ds 的路径。
    2A = Σ_k (s + S_k)，S_k 为该段内的部分和（S_0=0）。"""
    g = [dict() for _ in range(maxsteps + 1)]   # g[n][ds] -> dict(2A->count)
    g[0][0] = {0: 1}
    for n in range(1, maxsteps + 1):
        cur = defaultdict(lambda: defaultdict(int))
        for ds, areas in g[n - 1].items():
            for step in (1, -1):
                nds = ds + step
                # 新的一段：从 s 走到 s+step，这段贡献 2*area = (s + (s+step))
                # 但 area 是相对于起点的高度，故贡献 = ds + ds+step = 2ds+step
                add = 2 * ds + step
                for a, c in areas.items():
                    na = a + add
                    if abs(na) <= cutoff:
                        cur[nds][na] += c
        g[n] = {ds: dict(d) for ds, d in cur.items()}
    return g


def run(L, betas, cutoff=None):
    if cutoff is None:
        cutoff = 4 * L * L
    g = area_gf(L, cutoff)
    # 前缀 (k,s)：补全 n=L-k 步，净位移需为 -s
    # 前缀自身贡献 2A_pref = (k+1)*s
    recs = []
    for k in range(0, L + 1):
        n = L - k
        for s in range(-k, k + 1, 2):
            tbl = g[n].get(-s)
            if not tbl:
                continue
            pref = (k + 1) * s
            recs.append((k, s, pref, dict(tbl)))
    out = {"L": L, "n_prefix_classes": len(recs), "sweep": []}
    for beta in betas:
        w = []
        for k, s, pref, tbl in recs:
            tot = 0.0
            for a2, c in tbl.items():
                tot += c * math.exp(-beta * 0.5 * (pref + a2))
            if tot > 1e-300:
                w.append(tot)
        if len(w) < 3:
            continue
        S = smax(w)
        q = sum((x / sum(w)) ** 2 for x in w)
        mx = max(w) / sum(w)
        out["sweep"].append({"beta": beta, "classes": len(w),
                             "S_max": round(S, 6), "q": round(q, 8),
                             "max_share": round(mx, 6),
                             "contextual": bool(S > 2)})
    return out


if __name__ == "__main__":
    BETAS = [0.0, 0.05, -0.05, 0.1, -0.1, 0.2, -0.2, 0.3, -0.3,
             0.5, -0.5, 0.7, -0.7, 1.0, -1.0, 1.5, -1.5, 2.0, -2.0]
    for L in [8, 12, 16]:
        r = run(L, BETAS)
        OUT["L%d" % L] = r
        print("=" * 72)
        print("L = %d   前缀类数 = %d" % (L, r["n_prefix_classes"]))
        print("=" * 72)
        print("   beta   类数   S_max      最大类占比   量子?")
        for row in r["sweep"]:
            print("  %+-6.2f %-5d  %.6f   %.6f    %s"
                  % (row["beta"], row["classes"], row["S_max"],
                     row["max_share"], "是" if row["contextual"] else "否"))
        best = max(r["sweep"], key=lambda x: x["S_max"])
        OUT.setdefault("best", {})["L%d" % L] = best
        print("  -> 最大 S_max = %.6f 在 beta=%+.2f" % (best["S_max"], best["beta"]))
        print()
    with open(os.path.join(HERE, "R55_area_tilt_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("-> R55_area_tilt_results.json")
