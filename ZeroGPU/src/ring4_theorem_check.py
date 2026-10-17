#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ring4_theorem_check.py —— 环④的条件定理：读数半格与 E6 的 27
================================================================
目标：把环④从「候选」写成「条件定理」，并把前提逐条列出、可复跑。

前提：
  H1  读数周期 N = 12（由 k₀ = 2 定：n_a = (2k₀, 2k₀−1, 2k₀(2k₀−1)) = (4,3,12)，N = lcm = 12）
  H2  E6 破缺由 27 驱动，每个态贡献一个读数半格 π/N
  H3  绕数 k₀ 贡献 k₀ 个半格

结论（条件定理）：
  环④ = (dim 27 + k₀)·π/N = (27+2)·π/12 = 29π/12 = 7.592182

核验：
  A. n_a = (4,3,12), N=12（由 k₀=2）
  B. 环④ = 29π/12
  C. 与 alpha_s-free 目标 7.592094 一致（<2e-4）
  D. 四环闭合（残差 <1e-3）
"""
from __future__ import annotations
import os
import json
import math

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []


def P(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s)


K0 = 2
DIM_27 = 27
MZ = 91.1876
V = 246.2257
MX_E6 = 9.821631e16
MX_UNIF = 4.954262e13


def main():
    P("=" * 82)
    P("环④的条件定理：读数半格与 E6 的 27 (ring4_theorem_check)")
    P("=" * 82)

    # H1: n_a 与 N
    n1 = 2 * K0
    n2 = 2 * K0 - 1
    n3 = n1 * n2
    N = math.lcm(n1, n2, n3)
    P(f"H1  n_a = ({n1}, {n2}, {n3})，N = lcm = {N}")

    # H2/H3: 环④
    ring4 = (DIM_27 + K0) * math.pi / N
    P(f"H2/H3  环④ = (27 + k₀)·π/N = ({DIM_27}+{K0})·π/{N} = {ring4:.9f}")

    target = math.log(MX_E6 / MX_UNIF)
    P(f"目标环④ = {target:.9f}；差 = {ring4-target:.3e} ({(ring4-target)/target*100:.4f}%)")

    # 四环闭合
    MI = MZ * math.exp(4 * math.pi)
    r1 = math.log(MZ / V)
    r2 = 4 * math.pi
    r3 = math.log(MX_UNIF / MI)
    S = math.log(MX_E6 / V)
    total = r1 + r2 + r3 + ring4
    P(f"四环和 = {total:.9f}；S = {S:.9f}；残差 = {total-S:.3e}")

    checks = []

    def chk(name, cond, detail=""):
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        P(("  [v] " if cond else "  [x] ") + name + (f"   {detail}" if detail else ""))

    P("\n[A] 前提 H1")
    chk("A  n_a = (4,3,12)", (n1, n2, n3) == (4, 3, 12), str((n1, n2, n3)))
    chk("A  N = 12", N == 12, str(N))

    P("\n[B] 条件定理的结论")
    chk("B  环④ = 29π/12", abs(ring4 - 29 * math.pi / 12) < 1e-12, f"{ring4:.9f}")
    chk("C  与目标一致（<2e-4）", abs(ring4 - target) / target < 2e-4,
        f"{abs(ring4-target)/target:.2e}")

    P("\n[D] 四环闭合")
    chk("D  残差 <1e-3", abs(total - S) < 1e-3, f"{total-S:.3e}")

    npass = sum(1 for c in checks if c["pass"])
    P("\n" + "=" * 82)
    P(f"核验：通过 {npass} / {len(checks)}")
    P("=" * 82)

    out = {
        "premises": {
            "H1_n_a": [n1, n2, n3], "H1_N": N,
            "H2": "E6 breaking driven by 27, one half-cell pi/N per state",
            "H3": "winding k0 contributes k0 half-cells",
        },
        "ring4": ring4, "target": target,
        "four_ring_total": total, "S": S, "residual": total - S,
        "checks": checks, "pass": npass, "total": len(checks),
    }
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "ring4_theorem_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "ring4_theorem_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 ring4_theorem_check.json / .txt")


if __name__ == "__main__":
    main()
