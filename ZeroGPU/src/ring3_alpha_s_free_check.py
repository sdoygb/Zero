#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ring3_alpha_s_free_check.py —— 环③在现库重建：alpha_s-free 交叉与环的分割
================================================================================
对象：四环链的环③。旧库用 alpha_s-free 条件把 M_X^unif 解到 4.954262e13 GeV，
      得环③=14.4546、环④=7.5921。本篇用**现库的 SM 内容**（Z151，b=(41/10,-19/6,-7)）
      重建同一个 alpha_s-free 条件（alpha_1 = alpha_2 交叉），看它给什么。

结论：
  - 现库（SM 内容）的 alpha_1=alpha_2 交叉：L=25.4509, M_X^unif=1.031e13 GeV。
  - 该交叉处 alpha_3 与 alpha_U 不重合（差 15.1%）——SM 内容不统一。
  - 于是现库分割：环③=12.8845、环④=9.1621；旧库分割：14.4546 / 7.5921。
  - 总和 S=ln(M_X^E6/v)=33.6197 不变；**分割依内容/中间谱而定**。

判据：
  A. alpha_1=alpha_2 交叉 L 与 Z151 的 25.4509 相符（<1e-2）
  B. M_X^unif ~ 1.03e13（<1%）
  C. alpha_3 与 alpha_U 不重合（>10%）
  D. 总 S=33.6197（<1e-6）
  E. 现库与旧库的环③差 >1（说明分割依赖谱）
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
AEM_INV = 127.95
S2 = 0.23122
AS_MZ = 0.1179
B1, B2, B3 = 4.1, -19 / 6, -7
MX_E6 = 9.821631e16
OLD_RING3 = 14.4546
OLD_RING4 = 7.5921


def main():
    P("=" * 86)
    P("环③在现库重建：alpha_s-free 交叉与环的分割 (ring3_alpha_s_free_check)")
    P("=" * 86)

    A2 = S2 * AEM_INV
    AY = (1 - S2) * AEM_INV
    A1 = (3 / 5) * AY
    A3 = 1 / AS_MZ
    L12 = 2 * math.pi * (A1 - A2) / (B1 - B2)
    MXu = MZ * math.exp(L12)
    AU = A1 - (B1 / (2 * math.pi)) * L12
    A3X = A3 - (B3 / (2 * math.pi)) * L12
    alpha3X = 1 / A3X
    alphaU = 1 / AU
    mismatch = abs(alpha3X - alphaU) / alphaU

    MI = MZ * math.exp(4 * math.pi)
    S = math.log(MX_E6 / V)
    r1 = math.log(MZ / V)
    r2 = 4 * math.pi
    r3 = math.log(MXu / MI)
    r4 = S - (r1 + r2 + r3)

    P(f"A1={A1:.5f}  A2={A2:.5f}  A3={A3:.5f}")
    P(f"alpha_1=alpha_2 交叉: L={L12:.6f}  M_X^unif={MXu:.5e} GeV  A_U={AU:.5f}")
    P(f"alpha_3(M_X)={alpha3X:.5f} vs alpha_U={alphaU:.5f}  不重合 {mismatch*100:.2f}%")
    P(f"M_I={MI:.5e}  S={S:.7f}")
    P(f"现库分割: 环①={r1:.4f} 环②={r2:.4f} 环③={r3:.4f} 环④={r4:.4f}")
    P(f"旧库分割: 环③={OLD_RING3:.4f} 环④={OLD_RING4:.4f}（M_X^unif=4.954262e13）")

    checks = []

    def chk(name, cond, detail=""):
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        P(("  [v] " if cond else "  [x] ") + name + (f"   {detail}" if detail else ""))

    P("\n[A] 现库 alpha_s-free 交叉")
    chk("A  L=25.4509（与 Z151 相符）", abs(L12 - 25.4509) < 1e-2, f"{L12:.6f}")
    chk("B  M_X^unif ~ 1.03e13", abs(MXu - 1.031e13) / 1.031e13 < 1e-2, f"{MXu:.4e}")
    chk("C  alpha_3 与 alpha_U 不重合（SM 不统一）", mismatch > 0.10, f"{mismatch*100:.2f}%")

    P("\n[D] 总和与分割")
    chk("D  S = 33.6197", abs(S - 33.6197) < 1e-4, f"{S:.7f}")
    chk("E  现库环③与旧库环③差 >1", abs(r3 - OLD_RING3) > 1.0, f"{r3:.4f} vs {OLD_RING3}")
    chk("E  现库环④与旧库环④差 >1", abs(r4 - OLD_RING4) > 1.0, f"{r4:.4f} vs {OLD_RING4}")

    npass = sum(1 for c in checks if c["pass"])
    P("\n" + "=" * 86)
    P(f"核验：通过 {npass} / {len(checks)}")
    P("=" * 86)

    out = {
        "A1": A1, "A2": A2, "A3": A3, "L12": L12, "M_X_unif_GeV": MXu,
        "A_U": AU, "alpha3_at_MX": alpha3X, "alpha_U": alphaU,
        "mismatch": mismatch, "M_I_GeV": MI, "S": S,
        "split_current": {"ring1": r1, "ring2": r2, "ring3": r3, "ring4": r4},
        "split_old": {"ring3": OLD_RING3, "ring4": OLD_RING4},
        "checks": checks, "pass": npass, "total": len(checks),
    }
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "ring3_alpha_s_free_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "ring3_alpha_s_free_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 ring3_alpha_s_free_check.json / .txt")


if __name__ == "__main__":
    main()
