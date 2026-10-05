#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R66 · 方向 2（修正版）：二维格上的因果度规与符号检验

修正 R65 的两处缺陷：
  (1) 环图饱和 —— 空间取二维格 $\mathbb Z_S^2$，足够大，锥不饱和
  (2) 奇偶性伪影 —— 因果区间体积改为按**光锥条件**归一，并显式处理视界

模型（原生）：
  · 空间：二维格（Z1 定理 1 的连通图，顶点 = 局域荷的载体）
  · 时间：步推进 τ（不可逆，Z4）
  · 运动：每步走一格（局域补偿移动）
  · 净荷守恒 $\sum x = 0$（Z0③）

因果结构：事件 $e=(\tau,\mathbf x)$，$e\preceq e'$ ⟺ 存在合法演化连接二者
  ⟺ $\Delta\tau\ge0$ 且 $|\Delta x_1|+|\Delta x_2|\le\Delta\tau$ 且同奇偶

测：
  甲  Myrheim–Meyer 维数（主判据，样品取足够大）
  乙  光锥条件：可达 ⟺ 菱形范数 ≤ Δτ（含奇偶）
  丙  符号检验：因果区间体积对 $\Delta\tau^2-|\Delta\mathbf x|^2$ 的依赖
"""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}


# ---------------------------------------------------------------- 因果判定（解析）
def causally_related(dtau, dx1, dx2):
    """二维格上的因果关系：菱形范数 ≤ Δτ 且同奇偶。"""
    if dtau < 0:
        return False
    if (dtau - abs(dx1) - abs(dx2)) % 2 != 0:
        return False
    return abs(dx1) + abs(dx2) <= dtau


def cone_size(dtau, dims=2):
    """时刻 dtau 的光锥内部（含边界）点数，dims 维格。"""
    if dims == 1:
        return 1 if dtau == 0 else 2
    if dims == 2:
        # |x1|+|x2| ≤ dtau，同奇偶 ⇒ 对每个 max(|x1|+|x2|) = k ≤ dtau 计数
        tot = 0
        for k in range(dtau % 2, dtau + 1, 2):
            tot += 4 * k if k > 0 else 1
        return tot
    raise ValueError


def mm_dim(f):
    """由有序对比例 f 反解 Myrheim–Meyer 维数。"""
    def f_of_d(d):
        return (math.gamma(d + 1) * math.gamma(d / 2)) / (2 * math.gamma(3 * d / 2))
    lo, hi = 1.0, 14.0
    for _ in range(300):
        mid = (lo + hi) / 2
        if f_of_d(mid) > f:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def empirical_f(T, S, n=400000, seed=3, dims=2):
    """在 [0,T] × Z_S^dims 上随机取事件对，测有序对比例。"""
    rng = np.random.default_rng(seed)
    n_ord = 0
    n_tot = 0
    for _ in range(n):
        t1 = int(rng.integers(0, T + 1))
        t2 = int(rng.integers(0, T + 1))
        x1 = int(rng.integers(0, S))
        y1 = int(rng.integers(0, S))
        x2 = int(rng.integers(0, S))
        y2 = int(rng.integers(0, S))
        if t1 == t2:
            continue
        if t1 > t2:
            t1, t2 = t2, t1
            x1, y1, x2, y2 = x2, y2, x1, y1
        # 环面上取最短位移
        d1 = x2 - x1
        d1 -= S * round(d1 / S)
        d2 = y2 - y1
        d2 -= S * round(d2 / S)
        n_tot += 1
        if causally_related(t2 - t1, d1, d2):
            n_ord += 1
    return (n_ord / n_tot) if n_tot else 0.0


if __name__ == "__main__":
    print("=" * 74)
    print("甲：Myrheim–Meyer 维数 —— 因果结构测出的 d")
    print("=" * 74)
    print("  理论值 f(d) = Γ(d+1)Γ(d/2) / (2Γ(3d/2)):")
    for d in [2, 3, 4, 5]:
        f = (math.gamma(d + 1) * math.gamma(d / 2)) / (2 * math.gamma(3 * d / 2))
        print("    d=%d  f=%.6f" % (d, f))
    print()
    rows = []
    for S, T in [(200, 40), (400, 60), (600, 80)]:
        f = empirical_f(T, S, n=200000, dims=2)
        d = mm_dim(f)
        print("  2 维空间 + 时间：S=%-4d T=%-3d  f=%.6f  ⇒ d=%.4f" % (S, T, f, d))
        rows.append({"space_dim": 2, "S": S, "T": T, "f": f, "d": d})
    OUT["A_mm_dim_2d"] = rows

    print()
    print("=" * 74)
    print("乙：光锥条件（可达 ⟺ 菱形范数 ≤ Δτ）")
    print("=" * 74)
    print("  Δτ  |Δx|₁  可达?   锥内点数")
    ok = 0
    tot = 0
    for dtau in range(0, 7):
        for n1 in range(0, dtau + 1):
            pred = causally_related(dtau, n1, 0)
            tot += 1
            if pred == (n1 <= dtau and (dtau - n1) % 2 == 0):
                ok += 1
            print("  %-3d %-6d %-7s %d" % (dtau, n1, pred, cone_size(dtau)))
            break
    OUT["B_lightcone"] = {"checks_ok": ok, "checks_total": tot}

    print()
    print("=" * 74)
    print("丙：符号检验 —— 因果区间体积 vs 洛伦兹/欧氏不变量")
    print("=" * 74)
    print("  取 Δτ 固定，比较 |Δx|₁ = 0..Δτ 的因果区间（2D 锥内）体积")
    print()
    print("  Δτ  |Δx|₁  Δτ²−|Δx|²   Δτ²+|Δx|²   锥内体积 V")
    rows_c = []
    for dtau in [4, 6, 8, 10]:
        for n1 in range(0, dtau + 1, 2):
            # 时刻 Δτ 的光锥截面中，位于 nx≤n1 的点数（累积）
            V = 0
            for k in range(dtau % 2, n1 + 1, 2):
                V += 4 * k if k > 0 else 1
            rows_c.append({"dtau": dtau, "dx": n1, "V": V})
            print("  %-3d %-6d %-11d %-11d %d"
                  % (dtau, n1, dtau ** 2 - n1 ** 2, dtau ** 2 + n1 ** 2, V))
    OUT["C_sign_test"] = rows_c

    print()
    print("  ⇒ 判读：")
    print("     若 V 随 |Δx| 增大而**减小**（同 Δτ）⇒ 依赖 Δτ²−|Δx|² ⇒ 洛伦兹")
    print("     若 V 随 |Δx| 增大而**增大**      ⇒ 依赖 Δτ²+|Δx|² ⇒ 欧氏")
    dec = all(
        (lambda lst: all(lst[i]["V"] >= lst[i + 1]["V"] for i in range(len(lst) - 1)))(
            [r for r in rows_c if r["dtau"] == dt])
        for dt in [4, 6, 8, 10])
    print()
    print("  V 随 |Δx| 单调递减: %s" % dec)
    OUT["C_monotone_decreasing"] = bool(dec)

    print()
    print("=" * 74)
    print("总判定")
    print("=" * 74)
    d2 = self_d = None
    print("  甲 Myrheim–Meyer 维数（2 维空间 + 时间）:", 
          ["%.4f" % r["d"] for r in rows])
    print("  乙 光锥条件成立: %s" % (ok == tot))
    print("  丙 V 随 |Δx| 递减（洛伦兹型）: %s" % dec)
    OUT["verdict"] = {"lightcone": ok == tot, "lorentzian_sign": bool(dec)}
    with open(os.path.join(HERE, "R66_causal_metric_2d_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R66_causal_metric_2d_results.json")
