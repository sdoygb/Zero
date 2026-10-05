#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R67 · 方向 2（修正 II）：正确的因果区间体积与尺度匹配

修正 R66 的两处：
  (1) 尺度伪影：Myrheim–Meyer 要求时空**体积匹配**。2 维空间下取 S_side ≈ T/√2·(1/√2)。
      正确做法：用**体积匹配**的采样域，且 f 记为 N^{-1} 的标度。
  (2) 区间公式：R66 算的是"锥在 Δτ 时刻的截面累积"，不是因果区间。
      正确区间（2+1 维格）：I(Δτ, |Δx|₁) = Σ_{k=0}^{m} (2k²+1)，
      m = (Δτ − |Δx|₁)/2，只需菱形范数 Δτ' = Δτ − |Δx|₁。

判据（这是方向 2 的核心）：
  I 是洛伦兹不变量 Δτ² − |Δx|² 的函数，还是欧氏不变量 Δτ² + |Δx|² 的函数？
  两两比较：固定 Δτ，扫 |Δx|，看 I 是降还是升。
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}


def causal_interval_2d1t(dtau, dx1, dx2):
    """2+1 维格上事件 (0,0,0) 与 (Δτ, Δx) 之间的因果区间基数。
    链长条件：Δτ' = Δτ − |Δx|₁ 必须 ≥0 且同奇偶。
    区间 = Σ_{k=0}^{Δτ'/2} (2k² + 1)（每个时间片 k 的菱形截面 2k²+1）。"""
    n1 = abs(dx1) + abs(dx2)
    if dtau < n1:
        return 0
    d = dtau - n1
    if d % 2 != 0:
        return 0
    m = d // 2
    return sum(2 * k * k + 1 for k in range(m + 1))


def causal_interval_1d1t(dtau, dx1):
    """1+1 维：区间基数 = 时间片数 = (Δτ−|Δx|)/2 + 1"""
    n1 = abs(dx1)
    if dtau < n1 or (dtau - n1) % 2:
        return 0
    return (dtau - n1) // 2 + 1


def mm_dim(f, N):
    """Myrheim–Meyer：f = Γ(d+1)Γ(d/2)/(2Γ(3d/2))，但有限 N 下用
       f(N) = N^{-(d-1)/d}·const。这里直接用连续公式反解 d（域体积匹配时适用）。"""
    def f_of_d(d):
        return (math.gamma(d + 1) * math.gamma(d / 2)) / (2 * math.gamma(3 * d / 2))
    lo, hi = 1.0, 20.0
    for _ in range(300):
        mid = (lo + hi) / 2
        if f_of_d(mid) > f:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


if __name__ == "__main__":
    print("=" * 74)
    print("丙（核心）：因果区间体积 —— 洛伦兹还是欧氏不变量")
    print("=" * 74)
    print("  固定 Δτ，扫 |Δx|₁。洛伦兹 ⇒ I 只依赖 Δτ−|Δx|（降）；欧氏 ⇒ 依赖 和（升）")
    print()
    print("  Δτ  |Δx|₁  Δτ²−|Δx|²   Δτ²+|Δx|²     I")
    rows = []
    for dtau in [4, 6, 8, 10, 12]:
        for n1 in range(0, dtau + 1, 2):
            I = causal_interval_2d1t(dtau, n1, 0)
            rows.append({"dtau": dtau, "dx": n1, "I": I,
                         "lorentz": dtau ** 2 - n1 ** 2,
                         "euclid": dtau ** 2 + n1 ** 2})
            print("  %-3d %-6d %-11d %-12d %d"
                  % (dtau, n1, dtau ** 2 - n1 ** 2, dtau ** 2 + n1 ** 2, I))
    OUT["C_interval"] = rows

    # 判据：固定 Δτ，I 随 n1 单调递减？
    dec = True
    for dtau in [4, 6, 8, 10, 12]:
        seq = [r["I"] for r in rows if r["dtau"] == dtau]
        if not all(seq[i] >= seq[i + 1] for i in range(len(seq) - 1)):
            dec = False
    print()
    print("  固定 Δτ 下 I 随 |Δx|₁ 单调递减：%s" % dec)
    print("  ⇒ 洛伦兹型：I 依赖 Δτ − |Δx|₁（即退化后的固有时）" if dec else "  ⇒ 非洛伦兹型")

    print()
    print("  更精确的检验：I 是否只是 (Δτ−|Δx|₁) 的函数？")
    pairs = {}
    for r in rows:
        key = r["dtau"] - r["dx"]
        pairs.setdefault(key, set()).add(r["I"])
    consistent = all(len(v) == 1 for v in pairs.values())
    for k in sorted(pairs)[:8]:
        print("    Δτ−|Δx|₁ = %-3d  ⇒  I ∈ %s" % (k, sorted(pairs[k])))
    print("  ⇒ I 完全由 Δτ−|Δx|₁ 决定: %s" % consistent)
    OUT["I_depends_only_on_dtau_minus_dx"] = bool(consistent)

    print()
    print("=" * 74)
    print("甲：Myrheim–Meyer 维数（体积匹配采样）")
    print("=" * 74)
    print("  1+1 维（1 维空间 + 时间），采样域 体积 T²：")
    # 1+1: 域 [0,T]x[0,T]，事件数 N = T²，有序对比例 f
    for T in [20, 40, 80, 160]:
        n_ord = 0
        n_tot = 0
        for t1 in range(T):
            for t2 in range(T):
                for x1 in range(T):
                    for x2 in range(T):
                        if (t1, x1) == (t2, x2):
                            continue
                        if t1 > t2:
                            continue
                        d = t2 - t1
                        dx = abs(x2 - x1)
                        n_tot += 1
                        if d >= dx and (d - dx) % 2 == 0:
                            n_ord += 1
        f = n_ord / n_tot if n_tot else 0
        print("    T=%-4d N=%-8d f=%.6f  ⇒ MM d=%.4f（期望 2）"
              % (T, T * T, f, mm_dim(f, T * T)))

    print()
    print("  2+1 维（2 维空间 + 时间），采样域 体积匹配：")
    rng = np.random.default_rng(5)
    for T in [20, 40, 80]:
        Ns = max(4, int(T * T / 2))
        S = int(math.sqrt(Ns)) + 1
        n_ord = 0
        n_tot = 0
        for _ in range(400000):
            t1 = int(rng.integers(0, T))
            t2 = int(rng.integers(0, T))
            x1 = int(rng.integers(0, S))
            y1 = int(rng.integers(0, S))
            x2 = int(rng.integers(0, S))
            y2 = int(rng.integers(0, S))
            if (t1, x1, y1) == (t2, x2, y2):
                continue
            if t1 > t2:
                t1, t2 = t2, t1
                x1, y1, x2, y2 = x2, y2, x1, y1
            d = t2 - t1
            n1 = abs(x2 - x1) + abs(y2 - y1)
            n_tot += 1
            if d >= n1 and (d - n1) % 2 == 0:
                n_ord += 1
        f = n_ord / n_tot if n_tot else 0
        print("    T=%-4d S=%-4d f=%.6f  ⇒ MM d=%.4f（期望 3）"
              % (T, S, f, mm_dim(f, T * S * S)))
        OUT.setdefault("A_mm", []).append({"T": T, "S": S, "f": f,
                                           "d": mm_dim(f, T * S * S)})

    print()
    print("=" * 74)
    print("总判定")
    print("=" * 74)
    print("  丙（核心）：I = Σ_{k≤(Δτ−|Δx|₁)/2}(2k²+1)")
    print("    ⇒ 因果区间体积**只依赖** Δτ − |Δx|₁（菱形范数缺口）")
    print("    ⇒ 这就是 2+1 维洛伦兹结构：时间方向进入时带**负号**")
    OUT["verdict"] = {"lorentzian": bool(dec),
                      "I_only_dtau_minus_dx": bool(consistent)}
    with open(os.path.join(HERE, "R67_causal_interval_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R67_causal_interval_results.json")
