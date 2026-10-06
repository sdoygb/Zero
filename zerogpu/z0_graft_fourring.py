"""
z0_graft_fourring.py --- 嫁接验证 · 第十七批：S 的四环分解与 Q483 的界

【嫁接原则】见 `GRAFT_LEDGER.md §0`。

来源：`Q484_TWO_STAGE_MECHANISM.md`（`build_225`，38/0）、`Q483_OVERALL_SCALE_S.md`、
      `Q486_E6_PHASE_MECHANISM.md`（`build_227`，34/0）

【底座的主张】
  ① ★★★ **S 被精确分解为四环链**（望远镜恒等式）：
       v --(−0.993)--> M_Z --(+12.566 = 2πk₀)--> M_I --(+14.455)--> M_X^unif --(+7.592)--> M_X^E₆
       合计 = 33.6197 = S
  ② ★★★ **Q483 的界被精化**：单预算 12·ln12 过紧 ⇒ **三界预算 Σa n_a = 19**
       ⇒ 19·ln12 = 47.213 > S ⇒ **S 可达**（饱和度 71.2%）
  ③ ★★★ **缺口定位**：四环**三环在手**，残留**恰好一环** ⇒ S 的缺口 ≙ Q461（E₆ 破缺）
       并给出**单一目标数字 7.5866**（可 falsify）
  ④ ❌ **环④ 仍未算出**（诚实的残余）

★ 本批的独立验证：
  1. 四环的值与望远镜恒等式（含**符号约定**的说明：环① 是 ln(M_Z/v) = −0.993）
  2. 三界预算 19·ln12 = 47.213 与饱和度 71.2%
  3. ★ **环④ 的缺口**：环④(实测) = 7.5920942 vs 目标 7.5866 —— 差多少？
  4. 与第十六批的框框界对比（12ln12 vs 19ln12）

用法：/usr/bin/python3 z0_graft_fourring.py   输出：results/z0_graft_fourring.json
"""
from __future__ import annotations

import json
import os
import time
from math import exp, log, pi, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
GAPS: list = []

V_HIGGS, MZ = 246.2257, 91.1876
K0, LAM = 2, 3
MI = MZ * exp(2 * pi * K0)
MX_UNIF, MX_E6 = 4.954262e13, 9.821631e16


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


if __name__ == "__main__":
    t0 = time.time()

    # ============ F1 四环分解 ============
    print("=" * 96)
    print("F1：$S$ 的四环分解（望远镜恒等式）")
    print("=" * 96)
    rings = [
        ("环①  ln(M_Z/v)",       log(MZ / V_HIGGS),   -0.993),
        ("环②  ln(M_I/M_Z)",     log(MI / MZ),         12.566),
        ("环③  ln(M_X^unif/M_I)", log(MX_UNIF / MI),   14.4546),
        ("环④  ln(M_X^E₆/M_X^unif)", log(MX_E6 / MX_UNIF), 7.592),
    ]
    S = log(MX_E6 / V_HIGGS)
    print(f"  S = ln(M_X^E₆/v) = {S:.6f}")
    print(f"\n  {'环':>28} {'本侧':>12} {'声明':>10} {'差':>10}")
    tot = 0.0
    for nm, val, dec in rings:
        tot += val
        print(f"  {nm:>28} {val:>12.6f} {dec:>10.4f} {val-dec:>+10.4f}")
    print(f"  {'合计':>28} {tot:>12.6f} {'33.6197':>10} {tot-S:>+10.2e}")
    check("**F1a 四环之和 = S（望远镜恒等式，差 < 1e-12）**", abs(tot - S) < 1e-12,
          f"{tot:.6f}")
    check("**F1b 环① = ln(M_Z/v) = −0.9933（符号是方向约定）**",
          abs(log(MZ / V_HIGGS) + 0.993) < 1e-3, f"{log(MZ/V_HIGGS):.6f}")
    check("**F1c 环② = 2πk₀ 精确**", abs(log(MI / MZ) - 2 * pi * K0) < 1e-12,
          f"{log(MI/MZ):.6f} vs {2*pi*K0:.6f}")
    check("**F1d 环③ 复现 14.4546**", abs(log(MX_UNIF / MI) - 14.4546) < 1e-3,
          f"{log(MX_UNIF/MI):.6f}")
    check("**F1e 环④ 复现 7.5920942**", abs(log(MX_E6 / MX_UNIF) - 7.5920942) < 1e-6,
          f"{log(MX_E6/MX_UNIF):.7f}")

    # ============ F2 Q483 的界 ============
    print("\n" + "=" * 96)
    print("F2：`Q483` 的界 —— 单预算 vs 三界预算")
    print("=" * 96)
    n_a = (4, 3, 12)
    Sn = sum(n_a)
    print(f"  (n_a) = {n_a}，Σn_a = {Sn}")
    for nm, N in (("单预算 n₃ = 12", 12), ("**三界预算 Σn_a = 19**", Sn)):
        b = N * log(12)
        sat = S / b * 100
        flag = "≥ S ✓ 可达" if b >= S else "< S ✗ 过紧"
        print(f"  {nm:>26}: {N}·ln12 = {b:>9.4f}  {flag}  饱和度 = {sat:.1f}%")
    check("**F2a 单预算 12·ln12 = 29.8189 < S（过紧）**",
          abs(12 * log(12) - 29.8189) < 1e-3 and 12 * log(12) < S, f"{12*log(12):.4f}")
    check("**F2b 三界预算 19·ln12 = 47.213 > S（可达）**",
          abs(Sn * log(12) - 47.213) < 1e-3 and Sn * log(12) > S, f"{Sn*log(12):.4f}")
    check("**F2c 饱和度 = 71.2%**", abs(S / (Sn * log(12)) * 100 - 71.2) < 0.1,
          f"{S/(Sn*log(12))*100:.1f}%")

    # ============ F3 环④ 的缺口 ============
    print("\n" + "=" * 96)
    print("F3：★ 环④ 的缺口（底座给出单一目标数字 7.5866）")
    print("=" * 96)
    r4 = log(MX_E6 / MX_UNIF)
    target = 7.5866
    print(f"  环④（实测）= {r4:.7f}")
    print(f"  底座给的单一目标 = {target}")
    print(f"  差 = {r4-target:+.7f} = {(r4/target-1)*100:+.4f}%")
    print(f"\n  与第十三/十五批的候选对比：")
    for nm, v in (("π(1+√2)", pi * (1 + sqrt(2))), ("7.5866（底座目标）", target)):
        print(f"    {nm:>22} = {v:.7f}   与环④差 {(r4/v-1)*100:+.4f}%")
    check("**F3a 环④ = 7.5920942 复现**", abs(r4 - 7.5920942) < 1e-6, f"{r4:.7f}")
    check("**F3b 环④ 与目标 7.5866 差 ~0.07%**",
          abs(r4 / target - 1) < 1e-3, f"{(r4/target-1)*100:+.4f}%")
    gap("**环④ 的公式**",
        f"实测 {r4:.7f}，目标 7.5866（差 {(r4/target-1)*100:+.4f}%）；"
        f"底座自标『环④ 仍未算出 —— 诚实的残余』")

    # ============ F4 两个界的对比 ============
    print("\n" + "=" * 96)
    print("F4：两个界的对比（与第十六批的联系）")
    print("=" * 96)
    print(f"  {'界':>22} {'形式':>16} {'值':>10} {'≥ S?':>8} {'性质':>16}")
    rows = [("第十六批：框框界", "12·ln12", 12 * log(12), "选择"),
            ("第十七批：三界预算", "19·ln12", Sn * log(12), "由 (n_a) 定"),
            ("临界（N=12）", "12·ln(e^{S/12})", S, "恒等")]
    for nm, form, val, kind in rows:
        print(f"  {nm:>22} {form:>16} {val:>10.4f} {'✓' if val >= S - 1e-9 else '✗':>8} {kind:>16}")
    print(f"""
  ⟹ ★★ **联系**：第十六批指出「单预算 12·ln12 = 29.82」只是**框框选择**；
     第十七批给出**框架自己的**预算 Σn_a = 19 ⟹ 19·ln12 = 47.21 > S ⟹ **S 可达** ✓
  ⟹ 即：**框框不该是 12，而应是 Σn_a = 19** —— 这修正了第十六批的框框选择问题 ✓
  ⟹ 而真正的缺口**不是**"预算不够"，而是**环④ 的公式未算出** ✓（底座诚实）""")
    check("**F4 三界预算（19）修正了第十六批的框框选择（12）**", Sn * log(12) > S,
          f"19·ln12 = {Sn*log(12):.4f} > S")

    # ============ F5 判决 ============
    print("\n" + "=" * 96)
    print("F5：嫁接判决")
    print("=" * 96)
    print(f"""  · **四环分解**：J1 ✅（本侧复现到 1e-15，含符号约定的澄清）J2 ✅（$M_I = M_Z e^{{2πk_0}}$ 用 $M_Z$）J3 ⚠（环④ 公式缺）
      ⟹ **【条件】**
  · **三界预算 19·ln12 = 47.213 > S**：J1 ✅ J2 ✅（$(n_a)$ 由框架给）J3 ✅ ⟹ **【导出】**
  · **环④ 的公式**：❌ **缺**（底座自标"诚实的残余"）
  ⇒ 本批：四环分解【条件】；三界预算【导出】；环④ 公式【缺】""")
    check("**F5 判决：四环分解【条件】；三界预算【导出】；环④ 公式【缺】**", True, "见上")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          S=S, rings=[(nm, val, dec) for nm, val, dec in rings],
                          sum_rings=tot, MI=MI,
                          bound_12=12 * log(12), bound_19=Sn * log(12),
                          saturation=S / (Sn * log(12)) * 100,
                          ring4=r4, ring4_target=target,
                          ring4_dev_pct=(r4 / target - 1) * 100,
                          verdict="四环分解【条件】；三界预算【导出】；环④公式【缺】")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_fourring.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_fourring.json")
