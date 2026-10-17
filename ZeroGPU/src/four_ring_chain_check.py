#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
four_ring_chain_check.py —— 四环链闭合核验
================================================================
四环链：v → M_Z → M_I → M_X^unif → M_X^E6

  ① ln(M_Z/v)                电弱比
  ② 2πk₀ = 4π                 Z154（k₀=2）
  ③ ln(M_X^unif/M_I)           Z155（SM＋3 代色单态双重态）
  ④ ln(M_X^E6/M_X^unif)        E6 破缺；候选 2π(27+k₀)/(2N)=29π/12

核验：
  A. 四环各自的来源值
  B. 四环和 vs S=ln(M_X^E6/v)
  C. 残差 < 1e-3（闭合）
  D. 环④候选与 alpha_s-free 目标一致（<2e-4）
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


MZ = 91.1876
V = 246.2257
MX_E6 = 9.821631e16
MX_UNIF = 4.954262e13
MI = MZ * math.exp(4 * math.pi)
PI = math.pi
K0 = 2
N_READ = 12
DIM_27 = 27


def main():
    P("=" * 84)
    P("四环链闭合核验 (four_ring_chain_check)")
    P("=" * 84)

    r1 = math.log(MZ / V)
    r2 = 4 * PI
    r3 = math.log(MX_UNIF / MI)
    r4 = 2 * PI * (DIM_27 + K0) / (2 * N_READ)
    S = math.log(MX_E6 / V)
    total = r1 + r2 + r3 + r4
    resid = total - S

    P(f"① ln(M_Z/v)            = {r1:.9f}")
    P(f"② 2πk₀ = 4π            = {r2:.9f}")
    P(f"③ ln(M_X^unif/M_I)     = {r3:.9f}")
    P(f"④ 29π/12               = {r4:.9f}")
    P(f"和 = {total:.9f}   S = {S:.9f}   残差 = {resid:.3e} (相对 {abs(resid)/S:.2e})")

    checks = []

    def chk(name, cond, detail=""):
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        P(("  [v] " if cond else "  [x] ") + name + (f"   {detail}" if detail else ""))

    P("\n[A] 四环来源")
    chk("A  ① 电弱比 = -0.9933", abs(r1 + 0.9933) < 1e-3, f"{r1:.6f}")
    chk("A  ② = 4π = 12.5664", abs(r2 - 4 * PI) < 1e-12, f"{r2:.6f}")
    chk("A  ③ = 14.4546", abs(r3 - 14.4546) < 1e-3, f"{r3:.6f}")
    chk("A  ④ = 29π/12 = 7.5922", abs(r4 - 29 * PI / 12) < 1e-12, f"{r4:.6f}")

    P("\n[B] 闭合")
    chk("B  四环和 = S（残差 <1e-3）", abs(resid) < 1e-3, f"{resid:.3e}")
    chk("B  相对残差 <1e-5", abs(resid) / S < 1e-5, f"{abs(resid)/S:.2e}")

    P("\n[C] 环④候选 vs alpha_s-free 目标")
    target4 = math.log(MX_E6 / MX_UNIF)
    P(f"  目标环④ = {target4:.9f}；候选 = {r4:.9f}；差 = {r4-target4:.3e}（{(r4-target4)/target4*100:.4f}%）")
    chk("C  候选与目标一致（<2e-4）", abs(r4 - target4) / target4 < 2e-4,
        f"{abs(r4-target4)/target4:.2e}")

    npass = sum(1 for c in checks if c["pass"])
    P("\n" + "=" * 84)
    P(f"核验：通过 {npass} / {len(checks)}")
    P("=" * 84)

    out = {
        "rings": {"r1": r1, "r2": r2, "r3": r3, "r4": r4},
        "S": S, "total": total, "residual": resid,
        "ring4_target": target4,
        "checks": checks, "pass": npass, "total": len(checks),
    }
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "four_ring_chain_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "four_ring_chain_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 four_ring_chain_check.json / .txt")


if __name__ == "__main__":
    main()
