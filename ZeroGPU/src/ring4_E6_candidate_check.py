#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ring4_E6_candidate_check.py —— 环④（E6 破缺）的相位候选与命中密度
======================================================================
对象：Z153 §5.6 的四环链中唯一开放的一环 ④ = ln(M_X^{E6}/M_X^{unif})。

目标值（旧库 α_s-free 链的机器精度解）：
  M_X^{E6} = 9.821631e16 GeV,  M_X^{unif} = 4.954262e13 GeV
  ⇒ 环④ = 7.5920942609

候选：
  旧：π(1+√2) = 7.5844756          差 0.1003%（旧库已拒斥）
  新：2π(27+k₀)/24 = 29π/12        差 0.0012%
      其中 27 = dim(E6 基本表示)，k₀ = 2，24 = dim(SU(3)³ 伴随) = 3×8

命中密度：
  在“自然表达式”集 {πp/q (p,q≤30)} ∪ {2π(a+b/n) (a≤7, |b|≤6, n∈{3,4,6,12,24})}
  中，0.02% 内命中的个数。

判据：新候选相对差 < 2e-4；命中个数 = 1；旧候选 > 1e-3。
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


MX_E6 = 9.821631e16
MX_UNIF = 4.954262e13
PI = math.pi
K0 = 2
DIM_27 = 27
DIM_SU3CUBE = 24   # 3 × 8


def main():
    P("=" * 82)
    P("环④（E6 破缺）相位候选与命中密度 (ring4_E6_candidate_check)")
    P("=" * 82)

    target = math.log(MX_E6 / MX_UNIF)
    old = PI * (1 + math.sqrt(2))
    new = 2 * PI * (DIM_27 + K0) / DIM_SU3CUBE

    def rel(v):
        return abs(v - target) / target

    P(f"目标环④ = ln(M_X^E6/M_X^unif) = {target:.10f}")
    P(f"旧候选 π(1+√2)      = {old:.10f}  差 {rel(old)*100:.4f}%")
    P(f"新候选 2π(27+k₀)/24 = {new:.10f}  差 {rel(new)*100:.5f}%")
    P(f"      29π/12        = {29*PI/12:.10f}（同一值）")

    checks = []

    def chk(name, cond, detail=""):
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        P(("  [v] " if cond else "  [x] ") + name + (f"   {detail}" if detail else ""))

    P("\n[A] 候选")
    chk("A  新候选相对差 < 2e-4", rel(new) < 2e-4, f"{rel(new):.2e}")
    chk("A  旧候选相对差 > 1e-3", rel(old) > 1e-3, f"{rel(old):.2e}")
    chk("A  2π(27+k₀)/24 = 29π/12", abs(new - 29 * PI / 12) < 1e-12)

    # 命中密度
    P("\n[B] 命中密度（0.02% 内）")
    hits = []
    for p in range(1, 31):
        for q in range(1, 31):
            v = PI * p / q
            if rel(v) < 2e-4:
                hits.append((f"π·{p}/{q}", v, rel(v)))
    for a in range(0, 8):
        for b in range(-6, 7):
            for n in (3, 4, 6, 12, 24):
                v = 2 * PI * (a + b / n)
                if rel(v) < 2e-4:
                    hits.append((f"2π·({a}+{b}/{n})", v, rel(v)))
    P(f"  命中数 = {len(hits)}（按表达式）")
    uniq = []
    for name, v, r in hits:
        if not any(abs(v - u[1]) < 1e-12 for u in uniq):
            uniq.append((name, v, r))
    P(f"  去重后 = {len(uniq)}（按值；π·29/12 与 2π(1+5/24) 是同一值）")
    for name, v, r in uniq:
        P(f"    {name:18s} = {v:.8f}  差 {r*100:.5f}%")
    chk("B  去重后命中数 = 1", len(uniq) == 1, f"{len(uniq)}")
    chk("B  唯一命中 = 29π/12", len(uniq) == 1 and abs(uniq[0][1] - new) < 1e-12)

    npass = sum(1 for c in checks if c["pass"])
    P("\n" + "=" * 82)
    P(f"核验：通过 {npass} / {len(checks)}")
    P("=" * 82)

    out = {
        "target_ring4": target,
        "old_candidate": {"expr": "pi(1+sqrt2)", "value": old, "rel": rel(old)},
        "new_candidate": {"expr": "2*pi*(27+k0)/24 = 29*pi/12", "value": new, "rel": rel(new),
                          "dim_27": DIM_27, "k0": K0, "dim_SU3cube": DIM_SU3CUBE},
        "hit_density": {"hits": [h[0] for h in hits], "threshold": 2e-4},
        "checks": checks, "pass": npass, "total": len(checks),
    }
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "ring4_E6_candidate_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "ring4_E6_candidate_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 ring4_E6_candidate_check.json / .txt")


if __name__ == "__main__":
    main()
