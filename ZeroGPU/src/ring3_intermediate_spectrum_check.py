#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ring3_intermediate_spectrum_check.py —— 环③的中间谱：色单态双重态与 alpha_s-free 标度
=====================================================================================
对象：环③ = ln(M_X^unif/M_I)。SM-only 内容给 M_X=1.031e13、环③=12.88；
      加上「3 代色单态 SU(2) 双重态」（27 的 C 块）后，alpha_s-free 条件给出
      M_X=4.954e13、环③=14.4546、alpha_s(M_Z)=0.1179。

中间谱：
  SM 段（M_Z→M_I）：b_c=-7, b_L=-19/6, b_Y=(5/3)*4.1, b_R=(3/4)(b_Y-b_L/3)
  中间段（M_I→M_X）：b_L += 2, b_Y += 2（3 代色单态双重态），b_c 不变

条件：alpha_R(M_X) = alpha_L(M_X)（alpha_s-free）。
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
K0 = 2
AL_MZ = 29.57
AY_MZ = (5 / 3) * 58.98
AR_MZ = (3 / 4) * (AY_MZ - AL_MZ / 3)
BL, BY = -19 / 6, (5 / 3) * 4.1
BR = (3 / 4) * (BY - BL / 3)
B_SM = (-7.0, BL, BR)
B_INT = (-7.0, BL + 2.0, (3 / 4) * ((BY + 2.0) - (BL + 2.0) / 3))


def main():
    P("=" * 86)
    P("环③的中间谱：色单态双重态与 alpha_s-free 标度 (ring3_intermediate_spectrum_check)")
    P("=" * 86)

    t1 = K0
    numL = AL_MZ - B_SM[1] * t1
    numR = AR_MZ - B_SM[2] * t1
    t2 = (numR - numL) / (B_INT[2] - B_INT[1])
    MI = MZ * math.exp(2 * math.pi * K0)
    MX = MI * math.exp(2 * math.pi * t2)
    ring3 = math.log(MX / MI)
    ivL = AL_MZ - B_SM[1] * t1 - B_INT[1] * t2
    inv_ac = ivL + B_SM[0] * t1 + B_INT[0] * t2
    as_pred = 1 / inv_ac

    P(f"b_SM = ({B_SM[0]:.4f}, {B_SM[1]:.6f}, {B_SM[2]:.6f})")
    P(f"b_int= ({B_INT[0]:.4f}, {B_INT[1]:.6f}, {B_INT[2]:.6f})")
    P(f"t1=k0={t1},  t2={t2:.6f}")
    P(f"MI={MI:.5e}  MX={MX:.5e}  ring3={ring3:.6f}  alpha_s={as_pred:.5f}")

    checks = []

    def chk(name, got, want, tol):
        r = abs(got - want) / abs(want)
        ok = r < tol
        checks.append({"name": name, "value": got, "ref": want, "rel": r, "pass": bool(ok)})
        P(("  [v] " if ok else "  [x] ") + name +
          f"  值={got:.8g}  参考={want:.8g}  相对差={r:.2e}")

    P("\n[A] 中间谱复现")
    chk("A  M_X = 4.954262e13 GeV", MX, 4.954262e13, 1e-3)
    chk("A  环③ = 14.4546", ring3, 14.4546, 1e-3)
    chk("A  alpha_s(M_Z) = 0.1179", as_pred, 0.1179, 1e-2)
    chk("A  t2 = 2.3005", t2, 2.3005, 1e-3)

    P("\n[B] 对照：SM-only 内容（(Y,L) 基，b1=41/10）")
    def MX_sm_only():
        # SM-only alpha_1=alpha_2 crossing (GUT normalized)
        A1 = (3 / 5) * AY_MZ
        A2 = AL_MZ
        L = 2 * math.pi * (A1 - A2) / (4.1 - BL)
        return MZ * math.exp(L)
    MXsm = MX_sm_only()
    P(f"  SM-only M_X = {MXsm:.4e} GeV（中间谱高 {MX/MXsm:.2f} 倍）")
    chk("B  SM-only M_X = 1.009e13", MXsm, 1.009e13, 1e-2)

    npass = sum(1 for c in checks if c["pass"])
    P("\n" + "=" * 86)
    P(f"核验：通过 {npass} / {len(checks)}")
    P("=" * 86)

    out = {
        "b_SM": list(B_SM), "b_int": list(B_INT), "t1": t1, "t2": t2,
        "M_I_GeV": MI, "M_X_GeV": MX, "ring3": ring3, "alpha_s": as_pred,
        "M_X_SM_only_GeV": MXsm, "checks": checks, "pass": npass, "total": len(checks),
    }
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "ring3_intermediate_spectrum_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "ring3_intermediate_spectrum_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 ring3_intermediate_spectrum_check.json / .txt")


if __name__ == "__main__":
    main()
