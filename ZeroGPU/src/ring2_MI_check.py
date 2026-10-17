#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ring2_MI_check.py —— 环②：M_I = M_Z·e^{2πk₀}，k₀ = A_ℓ² = 2
================================================================
对象：Z153 §5.6 四环链的环②。

链条：
  Z014 定理 28：θ_M = π/2（圆分割唯一解）
  Z014 定理 29：振幅 √2
  弦长恒等式：A_ℓ = |1+e^{iθ_M}| = 2cos(θ_M/2) = √2
  几何不变量：弦 × 割 ≡ 2  ⇒  k₀ ≡ A_ℓ² = 2
  两标度关系：M_I = M_Z·e^{2πk₀} = M_Z·e^{4π}
  环②：ln(M_I/M_Z) = 2πk₀ = 4π = 12.5663706

核验：
  A. 弦长恒等式 |1+e^{iθ}| = 2cos(θ/2)（多角度）
  B. 弦 × 割 ≡ 2（多角度）
  C. θ_M=π/2 时 A_ℓ=√2、k₀=2
  D. M_I = M_Z e^{4π}；环② = 4π
  E. 与登记值 M_I=2.615e7 GeV 相符

判据：A–C 相对差 <1e-12；D 相对差 <1e-12；E 相对差 <1e-2。
"""
from __future__ import annotations
import os
import json
import math
import cmath

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []


def P(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s)


M_Z = 91.1876          # GeV
M_I_REF = 2.615e7      # GeV（登记值）
THETA_M = math.pi / 2  # Z014 定理 28


def rel(a, b):
    return abs(a - b) / abs(b)


def main():
    P("=" * 80)
    P("环②：M_I = M_Z·e^{2πk₀}，k₀ = A_ℓ² = 2 (ring2_MI_check)")
    P("=" * 80)

    checks = []

    def chk(name, a, b, tol):
        r = abs(a) if b == 0 else rel(a, b)
        ok = r < tol
        checks.append({"name": name, "value": a, "ref": b,
                       "rel": r, "pass": bool(ok)})
        P(("  [v] " if ok else "  [x] ") + name +
          f"  值={a:.12g}  参考={b:.12g}  差={r:.2e}")

    # A. 弦长恒等式
    P("\n[A] 弦长恒等式 |1+e^{iθ}| = 2cos(θ/2)")
    A_vals = []
    for th in [math.pi / 3, math.pi / 2, 2 * math.pi / 3, 5 * math.pi / 6]:
        lhs = abs(1 + cmath.exp(1j * th))
        rhs = 2 * math.cos(th / 2)
        A_vals.append(rel(lhs, rhs))
        P(f"  θ={th:.4f}: |1+e^iθ|={lhs:.12f}, 2cos(θ/2)={rhs:.12f}")
    chk("A  弦长恒等式（最大相对差）", max(A_vals), 0.0, 1e-12)

    # B. 弦 × 割 ≡ 2
    P("\n[B] 弦 × 割 ≡ 2")
    B_vals = []
    for th in [math.pi / 3, math.pi / 2, 2 * math.pi / 3, 5 * math.pi / 6]:
        chord = abs(1 + cmath.exp(1j * th))
        secant = 1 / math.cos(th / 2)
        B_vals.append(rel(chord * secant, 2.0))
        P(f"  θ={th:.4f}: 弦×割={chord*secant:.12f}")
    chk("B  弦×割≡2（最大相对差）", max(B_vals), 0.0, 1e-12)

    # C. θ_M=π/2
    P("\n[C] θ_M = π/2 ⇒ A_ℓ=√2, k₀=2")
    A_l = abs(1 + cmath.exp(1j * THETA_M))
    k0 = A_l * A_l
    chk("C  A_ℓ = √2", A_l, math.sqrt(2), 1e-12)
    chk("C  k₀ = A_ℓ² = 2", k0, 2.0, 1e-12)

    # D. M_I
    P("\n[D] M_I = M_Z·e^{2πk₀}")
    ring2 = 2 * math.pi * k0
    M_I = M_Z * math.exp(ring2)
    P(f"  2πk₀ = {ring2:.10f}  (= 4π = {4*math.pi:.10f})")
    P(f"  M_I = {M_I:.6e} GeV")
    chk("D  环② = 4π", ring2, 4 * math.pi, 1e-12)
    chk("D  M_I = 2.615e7 GeV", M_I, M_I_REF, 1e-2)

    npass = sum(1 for c in checks if c["pass"])
    P("\n" + "=" * 80)
    P(f"核验：通过 {npass} / {len(checks)}")
    P("=" * 80)

    out = {
        "theta_M": THETA_M, "A_l": A_l, "k0": k0,
        "M_Z_GeV": M_Z, "ring2": ring2, "M_I_GeV": M_I,
        "checks": checks, "pass": npass, "total": len(checks),
    }
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "ring2_MI_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "ring2_MI_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 ring2_MI_check.json / .txt")


if __name__ == "__main__":
    main()
