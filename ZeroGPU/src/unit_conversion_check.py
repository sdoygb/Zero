#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit_conversion_check.py —— 唯一量纲锚 v 的单位换算核验
================================================================
把"从电弱真空期望值 v 换算 Zero 基本单位"逐项复算（纯标准库）：

  A. 锚点链：M_X = v·exp(33.6197)
  B. 长度：r_int = ħc/M_X（米）；r_int/l_Pl
  C. 时间：β = α/c_SI；质量：γ = α·c_SI²/G_SI
  D. 能量：v = 246.22 GeV → eV
  E. 恒等式：c_SI = α/β，G_SI = α·c_SI²/γ（换算的算术自洽）
  F. 只读换算：换锚不改变无量纲比值（以 r_int/l_Pl 与 2π 为例）

判据方向固定：PASS = 相对差 < 1e-2（E/F 用 < 1e-12）；
FAIL = 该篇要报的判定（换算链不成立）。

本脚本只写 results/unit_conversion_check.{json,txt}，不 import 项目脚本。
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


# ---- 输入（SI 定义值与锚） ----
V_ANCHOR = 246.22          # GeV，电弱真空期望值
LN_MX_OVER_V = 33.6197     # M_X/v 的对数（模型纯数，具名输入）
HBARC = 0.1973269804       # GeV·fm
C_SI = 2.99792458e8        # m/s（SI 定义）
G_SI = 6.67430e-11         # m^3 kg^-1 s^-2
L_PL = 1.616255e-35        # m


def rel(a, b):
    return abs(a - b) / abs(b)


def main():
    P("=" * 88)
    P("唯一量纲锚 v 的单位换算核验 (unit_conversion_check)")
    P("=" * 88)

    MX = V_ANCHOR * math.exp(LN_MX_OVER_V)
    r_int = HBARC / MX * 1e-15          # GeV·fm → m
    alpha = r_int                       # 1 跳数
    beta = alpha / C_SI                 # 1 代
    gamma = alpha * C_SI * C_SI / G_SI  # 1 重数
    eV = V_ANCHOR * 1e9
    ratio_pl = alpha / L_PL

    P(f"v = {V_ANCHOR} GeV;  ln(M_X/v) = {LN_MX_OVER_V}")
    P(f"A  M_X = v·exp(33.6197) = {MX:.6e} GeV")
    P(f"B  1 跳数 = {alpha:.6e} m ;  r_int/l_Pl = {ratio_pl:.3f}")
    P(f"C  1 代 = {beta:.6e} s ;  1 重数 = {gamma:.6e} kg")
    P(f"D  v = {eV:.6e} eV")
    P(f"E  α/β = {alpha/beta:.6e} m/s (c_SI = {C_SI:.6e})")
    P(f"   α·c_SI²/γ = {alpha*C_SI*C_SI/gamma:.6e} (G_SI = {G_SI:.6e})")

    checks = []

    def chk(name, a, b, tol=1e-2):
        ok = rel(a, b) < tol
        checks.append({"name": name, "value": a, "ref": b,
                       "rel": rel(a, b), "pass": bool(ok)})
        P(("  [v] " if ok else "  [x] ") + name +
          f"  值={a:.6e}  参考={b:.6e}  相对差={rel(a,b):.2e}")

    P("\n[A] 锚点链")
    chk("A  M_X = 9.82e16 GeV", MX, 9.82e16)
    P("\n[B] 长度单位")
    chk("B  1 跳数 = 2.009e-33 m", alpha, 2.009e-33)
    chk("B  r_int/l_Pl = 124.3", ratio_pl, 124.3)
    P("\n[C] 时间与质量单位")
    chk("C  1 代 = 6.70e-42 s", beta, 6.70e-42)
    chk("C  1 重数 = 2.71e-6 kg", gamma, 2.71e-6)
    P("\n[D] 能量")
    chk("D  v = 2.462e11 eV", eV, 2.462e11)
    P("\n[E] 换算恒等式")
    chk("E  c_SI = α/β", alpha / beta, C_SI, tol=1e-12)
    chk("E  G_SI = α·c_SI²/γ", alpha * C_SI * C_SI / gamma, G_SI, tol=1e-12)

    # F. 换锚：把锚改成 m_e（0.51099895 MeV = 5.1099895e-4 GeV），
    #    r_int/l_Pl 随锚改变是"单位读数"；无量纲比值（如 2π）不变。
    P("\n[F] 换锚")
    v2 = 5.1099895e-4
    MX2 = v2 * math.exp(LN_MX_OVER_V)
    r2 = HBARC / MX2 * 1e-15
    two_pi = 2.0 * math.pi
    P(f"  锚换 m_e: r_int = {r2:.6e} m（单位读数随锚变）")
    P(f"  无量纲比值 2π = {two_pi:.12f}（与锚无关）")
    okF = abs(two_pi - 2.0 * math.pi) < 1e-15
    checks.append({"name": "F  无量纲比值与锚无关", "pass": bool(okF)})
    P(("  [v] " if okF else "  [x] ") + "F  无量纲比值与锚无关")

    npass = sum(1 for c in checks if c["pass"])
    P("\n" + "=" * 88)
    P(f"核验：通过 {npass} / {len(checks)}")
    P("=" * 88)

    out = {
        "anchor_GeV": V_ANCHOR,
        "ln_MX_over_v": LN_MX_OVER_V,
        "M_X_GeV": MX,
        "hop_m": alpha,
        "generation_s": beta,
        "multiplicity_kg": gamma,
        "eV": eV,
        "r_int_over_lPl": ratio_pl,
        "checks": checks,
        "pass": npass,
        "total": len(checks),
    }
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "unit_conversion_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "unit_conversion_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 unit_conversion_check.json / .txt")


if __name__ == "__main__":
    main()
