#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
weak_angle_L_check.py —— 弱混合角的 L：由 α_s-free 条件定出，不是自由参数
================================================================================
对象：Z151 定理 151.4 的 L≈25.4509（复现 sin²θ_w(M_Z)=0.231220 所需的单圈对数）。

事实链：
  L = ln(M_U/M_Z)，M_U 是 α_s-free 条件 A_1(M_X)=A_2(M_X) 定出的统一标度。
  SM-only 分支：M_U=1.009e13 ⇒ L=25.4296（与 Z151 拟合值 25.4509 差 -0.084%）。
  中间谱分支：M_X=4.954e13 ⇒ L=27.0209 = 4π+14.454565（四环链的环②+环③）。

结论：L 不是自由参数；它是 α_s-free 统一标度的对数。弱混合角的值仍是低能输入
（框架用它定 M_X，再预言 α_s），不是被预言量。

判据：A–C 全成立 → PASS。
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


def main():
    P("=" * 80)
    P("弱混合角的 L：由 α_s-free 定出 (weak_angle_L_check)")
    P("=" * 80)

    MZ = 91.1876
    MI = MZ * math.exp(4 * math.pi)
    MX_sm = 1.009e13
    MX_int = 4.954262e13
    L_sm = math.log(MX_sm / MZ)
    L_int = math.log(MX_int / MZ)
    L151 = 25.4509
    chain = 4 * math.pi + 14.454565

    P(f"M_Z = {MZ}")
    P(f"M_I = M_Z e^{{4π}} = {MI:.3f}")
    P(f"L_sm  = ln({MX_sm:.3e}/{MZ}) = {L_sm:.6f}")
    P(f"L_int = ln({MX_int:.3e}/{MZ}) = {L_int:.6f}")
    P(f"环②+环③ = 4π+14.454565 = {chain:.6f}")
    P(f"Z151 拟合 L = {L151}；L_sm 相对差 = {(L_sm-L151)/L151*100:.4f}%")

    checks = []

    def chk(name, cond, detail=""):
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        P(("  [v] " if cond else "  [x] ") + name + (f"   {detail}" if detail else ""))

    chk("A  L_sm 与 Z151 的 25.4509 一致（<0.5%）",
        abs(L_sm - L151) / L151 < 5e-3, f"{(L_sm-L151)/L151*100:.4f}%")
    chk("B  L_int = 环②+环③（<1e-4）",
        abs(L_int - chain) < 1e-4, f"{abs(L_int-chain):.2e}")
    chk("C  L_int > L_sm（中间谱抬高标度）", L_int > L_sm,
        f"{L_sm:.4f} < {L_int:.4f}")

    npass = sum(1 for c in checks if c["pass"])
    P("\n" + "=" * 80)
    P(f"核验：通过 {npass} / {len(checks)}")
    P("=" * 80)

    out = {"MZ": MZ, "MI": MI, "L_sm": L_sm, "L_int": L_int,
           "L151": L151, "chain": chain, "checks": checks,
           "pass": npass, "total": len(checks)}
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "weak_angle_L_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "weak_angle_L_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 weak_angle_L_check.json / .txt")


if __name__ == "__main__":
    main()
