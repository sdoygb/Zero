#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R61 · 电阻度规路（甲）：初始段图上的有效电阻与标度指数

原材料（全部来自已证猜想）：
  K9 初始段      → 图的顶点 = 部分和取值（壳层），边 = ±1 步
  K3 量子线      → 测度 = 初始段权重
  K11 局部配平   → 边权 = 该边被"能闭合到零的路径"穿越的计数

构造：
  顶点 v ∈ ℤ（部分和取值）
  边 (v, v+1) 的**电导** = 穿越该边的初始段计数
       w_v = Σ_{k} #{长度 k 的路径, 起点 0, 终点 v, 且可续接成闭路径}
  有效电阻（路径图上的串联）：
       R(0, n) = Σ_{v=0}^{n-1} 1 / w_v

测：
  甲1  R 是否随距离单调增长（非退化）
  甲2  标度指数 ζ = d ln R / d ln r 是否稳定 ⇒ 反推维数
  甲3  与连续极限的对照（1D 自由游走：R ~ ln r；3D 以上：R → 常数）
"""
from __future__ import annotations

import json
import math
import os
from collections import Counter
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}


def edge_conductance(L):
    """w_v = 穿越边 (v, v+1) 的初始段计数（含可续接约束）。
    对起点 0、终点 v 的长度 k 初始段，其续接数 = C(L-k, (L-k-v)/2)。
    边 (v,v+1) 被穿越 ⟺ 路径在某个时刻位于 v 且下一步为 +1。"""
    w = Counter()
    for k in range(L + 1):
        for v in range(-k, k + 1, 2):
            rem = L - k
            if (rem - v) % 2:
                continue
            a = (rem - v) // 2
            if a < 0 or a > rem:
                continue
            W = comb(rem, a)          # 该初始段的续接数
            # 从 v 出发的两步：+1 跨过边 (v,v+1)，-1 跨过边 (v-1,v)
            w[v] += W                 # 跨 (v, v+1)
            w[v - 1] += W             # 跨 (v-1, v)
    return w


def resistance(w, nmax):
    """R(0,n) = Σ_{v=0}^{n-1} 1/w_v（路径图串联）"""
    R = [0.0]
    for n in range(1, nmax + 1):
        v = n - 1
        wv = w.get(v, 0) + w.get(-v - 1, 0)   # 对称性：用 v 与 -v-1 两段之和
        R.append(R[-1] + (1.0 / wv if wv > 0 else float("inf")))
    return R


if __name__ == "__main__":
    print("=" * 74)
    print("甲1/甲2：初始段图上的有效电阻与标度指数")
    print("=" * 74)
    allrows = []
    for L in [8, 12, 16, 20, 24]:
        w = edge_conductance(L)
        nmax = L // 2
        R = resistance(w, nmax)
        print()
        print("L = %d" % L)
        print("  边 (v,v+1) 的电导（前若干）:",
              [w.get(v, 0) for v in range(0, min(nmax, 8))])
        print("  R(0,n) 前若干:", [round(x, 6) for x in R[:min(nmax, 8) + 1]])
        # 标度指数
        rows = []
        for n in range(1, nmax + 1):
            if R[n] > 0 and n > 1:
                zeta = (math.log(R[n]) - math.log(R[n - 1])) / (
                    math.log(n) - math.log(n - 1))
                rows.append((n, R[n], zeta))
        print("  n, R, 局部指数 ζ:")
        for n, r, z in rows[-6:]:
            print("    n=%-3d R=%-12.6f ζ=%.4f" % (n, r, z))
        if rows:
            zs = [z for _, _, z in rows if math.isfinite(z)]
            allrows.append({"L": L, "nmax": nmax, "R_max": R[nmax],
                            "zeta_last": zs[-1] if zs else None,
                            "zeta_mean": sum(zs) / len(zs) if zs else None,
                            "R_diverges": R[nmax] > 5 * R[max(1, nmax // 2)]})
            print("  R_max=%.6f  ζ_last=%.4f  ζ_mean=%.4f"
                  % (R[nmax], zs[-1] if zs else float('nan'),
                     sum(zs) / len(zs) if zs else float('nan')))
    OUT["rows"] = allrows

    print()
    print("=" * 74)
    print("甲3：与已知情形的对照")
    print("=" * 74)
    print("  1D 自由游走（均匀边权）: R ~ n    ⇒ ζ → 1")
    print("  1D 对数发散特例        : R ~ ln n ⇒ ζ → 0")
    print("  3D 以上（暂留）        : R → 常数 ⇒ ζ → 0 但 R 有界")
    print()
    print("  判读规则：")
    print("    ζ 稳定且 ≈1  ⇒ 图给出线性度规，维数 1")
    print("    ζ → 0 且 R 有界 ⇒ 电阻退化（大尺度不可分辨）")
    print("    ζ → 0 且 R 发散 ⇒ 对数型，维数 1 但非标准")

    print()
    print("=" * 74)
    print("ζ 随 L 的稳定性（决定能否反推维数）")
    print("=" * 74)
    print("  L    R_max        ζ_last     ζ_mean")
    for r in allrows:
        print("  %-4d %-12.6f %-10s %s"
              % (r["L"], r["R_max"],
                 ("%.4f" % r["zeta_last"]) if r["zeta_last"] is not None else "—",
                 ("%.4f" % r["zeta_mean"]) if r["zeta_mean"] is not None else "—"))
    zs = [r["zeta_last"] for r in allrows if r["zeta_last"] is not None]
    OUT["zeta_stable"] = bool(len(zs) > 1 and (max(zs) - min(zs)) < 0.05)
    OUT["zeta_span"] = (max(zs) - min(zs)) if zs else None
    print()
    print("  ζ 的跨度 = %s  ⇒ %s"
          % (("%.4f" % OUT["zeta_span"]) if OUT["zeta_span"] is not None else "—",
             "稳定" if OUT["zeta_stable"] else "不稳定"))
    with open(os.path.join(HERE, "R61_resistance_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("\n-> R61_resistance_results.json")
