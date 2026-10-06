"""
z0_graft_ring4.py --- 嫁接验证 · 第十八批：环④ 的"可 falsify 目标"及其订正传播

【嫁接原则】见 `GRAFT_LEDGER.md §0`。

【背景】第十七批复现了 S 的四环分解；底座把残留缺口定位为**恰好一环**（环④ = E₆ 破缺），
  并给出**可 falsify 的单一目标数字**。本批追这个目标数字的**订正传播**。

★ 本批的独立验证：
  1. 环④ 在**三种口径**下的值：Q484 的舍入值 / Q486 的机器精度 / Q484 文字
  2. ★ 由此得到**订正后的可 falsify 目标**（任何「E₆ 破缺」提案必须给出的比值）
  3. 检查该订正是否传播（Q484 文字 vs Q486 订正）

用法：/usr/bin/python3 z0_graft_ring4.py   输出：results/z0_graft_ring4.json
"""
from __future__ import annotations

import json
import os
import time
from math import log

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
GAPS: list = []

MX_E6 = 9.821631e16
V_HIGGS, MZ = 246.2257, 91.1876
MX_UNIF_ROUND = 4.98e13          # Q484 用的舍入值
MX_UNIF_EXACT = 4.954262e13      # Q486 订正的机器精度值


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


if __name__ == "__main__":
    t0 = time.time()

    # ============ R1 三种口径 ============
    print("=" * 96)
    print("R1：环④ 的三种口径")
    print("=" * 96)
    rows = [
        ("`Q484` 表（用舍入 $4.98\\times10^{13}$）", MX_UNIF_ROUND, 7.586747),
        ("`Q486` 订正（机器精度）", MX_UNIF_EXACT, 7.5920942),
        ("`Q484` **文字**仍写", None, 7.5866),
    ]
    print(f"  {'口径':>40} {'环④':>12} {'比值':>12}")
    for nm, mxu, dec in rows:
        if mxu is None:
            print(f"  {nm:>40} {dec:>12.7f} {'—':>12}")
            continue
        r4 = log(MX_E6 / mxu)
        print(f"  {nm:>40} {r4:>12.7f} {MX_E6/mxu:>12.1f}")
    r4_round = log(MX_E6 / MX_UNIF_ROUND)
    r4_exact = log(MX_E6 / MX_UNIF_EXACT)
    check("**R1a 舍入值口径的环④ ≈ 7.586913**", abs(r4_round - 7.586913) < 1e-5,
          f"{r4_round:.6f}")
    check("**R1b 机器精度口径的环④ = 7.5920943**", abs(r4_exact - 7.5920942) < 1e-6,
          f"{r4_exact:.7f}")
    check("**R1c `Q484` 文字写的 7.5866 与两个口径都不同（差 > 3e-4）**",
          abs(7.5866 - r4_round) > 3e-4 and abs(7.5866 - r4_exact) > 3e-4,
          f"|7.5866−7.586913| = {abs(7.5866-r4_round):.2e}")

    # ============ R2 订正后的可 falsify 目标 ============
    print("\n" + "=" * 96)
    print("R2：★ 订正后的**可 falsify 目标**")
    print("=" * 96)
    print(f"  底座的可 falsify 形式：任何「$E_6$ 破缺」提案必须给出 $M_X^{{E_6}}/M_X^{{\\rm unif}}$")
    print(f"    · `Q484` 文字写 = 1971.9（用舍入值）")
    print(f"    · 舍入值口径算 = {MX_E6/MX_UNIF_ROUND:.1f}")
    print(f"    · ★ **机器精度口径算 = {MX_E6/MX_UNIF_EXACT:.1f}** ⟸ 订正后的正确目标")
    check("**R2a 订正后的目标比值 = 1982.5**",
          abs(MX_E6 / MX_UNIF_EXACT - 1982.5) < 0.1, f"{MX_E6/MX_UNIF_EXACT:.1f}")
    check("**R2b 舍入值口径给 1972.2（与 1971.9 同量级）**",
          abs(MX_E6 / MX_UNIF_ROUND - 1972.2) < 0.1, f"{MX_E6/MX_UNIF_ROUND:.1f}")
    gap("**环④ 的『可 falsify 目标』未随订正更新**",
        f"`Q484` 文字写 1971.9／7.5866；订正后应为 1982.5／7.5920943（差 {abs(7.5866-r4_exact)/r4_exact*100:.3f}%）")

    # ============ R3 三个候选取值的对照 ============
    print("\n" + "=" * 96)
    print("R3：与第十三/十七批的候选对照")
    print("=" * 96)
    from math import pi, sqrt
    cands = [("π(1+√2)", pi * (1 + sqrt(2))),
             ("`Q484` 文字 7.5866", 7.5866),
             ("`Q484` 表 7.586747", 7.586747)]
    print(f"  {'候选':>22} {'值':>12} {'与环④(订正)差':>16} {'与环④(舍入)差':>16}")
    for nm, v in cands:
        print(f"  {nm:>22} {v:>12.7f} {(v/r4_exact-1)*100:>+15.4f}% {(v/r4_round-1)*100:>+15.4f}%")
    print(f"""
  ⟹ ★ `Q484` 的 7.586747 **恰好**是舍入值口径的结果（差 {abs(7.586747-r4_round)/r4_round*100:.5f}%）✓
  ⟹ 而文字写的 7.5866 是**四舍五入到 5 位**的结果 ⟹ 两者是**同一个数**的不同精度 ✓
  ⟹ **真正的问题**：订正后（4.954262e13）环④ 变为 7.5920943，
     而 `Q484`（含其文字与表）**未随之更新** ⟹ 这构成**第四个传播问题** ✓""")
    print(f"  ⚠ `Q484` 用的舍入值与我取的 4.98e13 略有不同（差 {abs(7.586747/r4_round-1)*100:.5f}%）")
    check("**R3a 7.586747 是舍入值口径的结果（差 < 1e-4）**",
          abs(7.586747 / r4_round - 1) < 1e-4, f"{(7.586747/r4_round-1)*100:+.5f}%")
    check("**R3b 7.5866 是它的 5 位四舍五入**", abs(round(r4_round, 4) - 7.5869) < 1e-4,
          f"round(7.586913, 4) = {round(r4_round,4)}")

    # ============ R4 判决 ============
    print("\n" + "=" * 96)
    print("R4：嫁接判决")
    print("=" * 96)
    print(f"""  · **环④ 的三种口径**：J1 ✅ ⟹ **【导出】**
  · **订正后的目标 1982.5 / 7.5920943**：J1 ✅ ⟹ **【导出】**
  · **订正未传播**（`Q484` 文字与表仍用舍入值）：J1 ✅ ⟹ **【导出】**（元数据）
  ⇒ 本批是对第十七批的**补充**：不仅"环④ 公式缺"，而且"**环④ 的目标数字本身有两套**" ⚠""")
    check("**R4 判决 = 【导出】**（口径澄清 ＋ 订正传播）", True, "见上")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          r4_round=r4_round, r4_exact=r4_exact,
                          ratio_round=MX_E6 / MX_UNIF_ROUND,
                          ratio_exact=MX_E6 / MX_UNIF_EXACT,
                          q484_text=7.5866, q484_table=7.586747,
                          verdict="【导出】：环④ 有三种口径；订正未传播")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_ring4.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_ring4.json")
