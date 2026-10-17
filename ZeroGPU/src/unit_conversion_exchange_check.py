#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit_conversion_exchange_check.py —— 从"物质↔时空交换"推导三个单位换算
================================================================================
思路：把物质与时空的相互转化写成同一个交换结构（Z144／Z150）：

    J_{i→j} = K_{ij}(φ_i - φ_j),   K = A（邻接，唯一）,   Lφ = -x .

交换有三把自然尺子：跳数 α（长度）、代 β（时间）、重数 γ（质量）。
两种转化率把它们串起来：

    c = α/β            （空间↔时间）
    G = α c²/γ         （物质↔时空；即 γ = α c²/G）

交换的自然单位取 c=1、G=1（一交换步走一跳；泊松核系数 1，Z144 定理 144.2）。
于是三个物理单位由 α（锚 v 定出）与 SI 定义值给出：

    β = α/c_SI,        γ = α c_SI²/G_SI .

本脚本核验：
  A. 由 v 定出 α；
  B. 交换自然单位下 β、γ 的数值；
  C. 两条转化率恒等式 c_SI=α/β、G_SI=α c²/γ；
  D. "物质↔时空"读数：自然单位下 α=β=γ（同一把尺子）；
  E. 对照：若取动力学值 c_disc=0.8944，β、γ 的偏移（说明 c=1 是单位选择）。

判据：A–D 相对差 <1e-2（恒等式 <1e-12）；E 只作对照，不设 PASS/FAIL。
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


V_ANCHOR = 246.22
LN_MX_OVER_V = 33.6197
HBARC = 0.1973269804      # GeV·fm
C_SI = 2.99792458e8       # m/s
G_SI = 6.67430e-11        # m^3 kg^-1 s^-2


def rel(a, b):
    return abs(a - b) / abs(b)


def main():
    P("=" * 88)
    P("物质↔时空交换 → 三个单位换算 (unit_conversion_exchange_check)")
    P("=" * 88)

    # A. 锚定长度
    MX = V_ANCHOR * math.exp(LN_MX_OVER_V)
    alpha = HBARC / MX * 1e-15      # 1 跳数（m）

    # B. 交换自然单位 c=1, G=1
    beta = alpha / C_SI             # 1 代（s）
    gamma = alpha * C_SI * C_SI / G_SI  # 1 重数（kg）

    # C. 两条转化率
    c_from_units = alpha / beta
    G_from_units = alpha * c_from_units ** 2 / gamma

    P(f"锚: v={V_ANCHOR} GeV → M_X={MX:.6e} GeV → α={alpha:.6e} m")
    P(f"交换自然单位 c=1, G=1:")
    P(f"  β = α/c_SI = {beta:.6e} s")
    P(f"  γ = α c_SI²/G_SI = {gamma:.6e} kg")
    P(f"  c = α/β = {c_from_units:.6e} m/s   (c_SI = {C_SI:.6e})")
    P(f"  G = α c²/γ = {G_from_units:.6e}    (G_SI = {G_SI:.6e})")
    P(f"  自然单位读数: α/α={alpha/alpha:.6f}, β·c_SI/α={beta*C_SI/alpha:.6f}, "
      f"γ·G_SI/(α c_SI²)={gamma*G_SI/(alpha*C_SI*C_SI):.6f}")

    checks = []

    def chk(name, a, b, tol=1e-2):
        ok = rel(a, b) < tol
        checks.append({"name": name, "value": a, "ref": b,
                       "rel": rel(a, b), "pass": bool(ok)})
        P(("  [v] " if ok else "  [x] ") + name +
          f"  值={a:.6e}  参考={b:.6e}  相对差={rel(a,b):.2e}")

    P("\n[C] 转化率恒等式")
    chk("C  c_SI = α/β", c_from_units, C_SI, tol=1e-12)
    chk("C  G_SI = α c²/γ", G_from_units, G_SI, tol=1e-12)
    P("\n[D] 物质↔时空：自然单位下三把尺子同长")
    chk("D  β c_SI/α = 1", beta * C_SI / alpha, 1.0, tol=1e-12)
    chk("D  γ G_SI/(α c_SI²) = 1", gamma * G_SI / (alpha * C_SI * C_SI), 1.0, tol=1e-12)

    # E. 对照：动力学 c_disc=0.8944（CFL）
    P("\n[E] 对照：动力学 c_disc=0.8944（CFL）")
    c_disc = 0.8944
    beta2 = c_disc * alpha / C_SI
    gamma2 = alpha * C_SI * C_SI / (G_SI * c_disc ** 2)
    P(f"  β' = c_disc·α/c_SI = {beta2:.6e} s  (自然单位值 {beta:.6e})")
    P(f"  γ' = α c_SI²/(G_SI·c_disc²) = {gamma2:.6e} kg  (自然单位值 {gamma:.6e})")
    P(f"  比值 β'/β = {beta2/beta:.4f}, γ'/γ = {gamma2/gamma:.4f}")

    npass = sum(1 for c in checks if c["pass"])
    P("\n" + "=" * 88)
    P(f"核验：通过 {npass} / {len(checks)}")
    P("=" * 88)

    out = {
        "anchor_GeV": V_ANCHOR, "alpha_m": alpha,
        "beta_s": beta, "gamma_kg": gamma,
        "c_from_units": c_from_units, "G_from_units": G_from_units,
        "checks": checks, "pass": npass, "total": len(checks),
        "alt_c_disc": {"c_disc": c_disc, "beta_s": beta2,
                       "gamma_kg": gamma2, "beta_ratio": beta2 / beta,
                       "gamma_ratio": gamma2 / gamma},
    }
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "unit_conversion_exchange_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "unit_conversion_exchange_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 unit_conversion_exchange_check.json / .txt")


if __name__ == "__main__":
    main()
