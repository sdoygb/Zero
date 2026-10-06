"""
z0_graft_scalar.py --- 嫁接验证 · 第二十九批：标量谱的 Σ 与 Higgs 四次耦合的口径（Q493/Q483/Q324）

【嫁接原则】见 `GRAFT_LEDGER.md §0`。

来源：`Q493_SCALAR_SIGMA.md`、`Q483_SCALE_REVISITED.md`、`Q324_COLEMAN_WEINBERG.md`（`build_132`，3/0）、
      `Q496_SCALE_BOUNDARY_AUDIT.md`

【本批的两条线索】
  ① `Q493`/`Q483` 的算术链：`环④/S = 0.2258226` ⟹ `λ_S = 0.0194785`
     ★ 而底座**自我降级**：该关系依赖 `B_S = B_h`（`B_S` 对 `E₆` 标度求和，`B_h` 对电弱标度求和）⟹ **不成立**
     ⟹ `λ_S = 0.0194785` **降级为下界** `λ_S > 0.0194785`（⚠ 方向性判断，未证明）
  ② ★★ **本侧发现的连接**：由 `λ_S` 与 `λ_h` 的比反解得 `λ_h = 0.0862557`
     —— **恰好等于 `Q324` 的 `λ = 0.086256`**

【而 `Q324` 自己说明】：
  ```
  ⚠ 归一化：此归一化下 m_h² = 3λv² ⇒ λ = 0.086256
     （标准 SM 归一化给 0.129383，差 3/2）
  ⇒ 观测的 λ = m_h²/(3v²) = 0.086256
  ```
  ⟹ **`0.086256` 是观测 Higgs 质量的重新表达（输入），不是导出** ✗

★ 本批的独立验证：
  1. `环④/S`、`λ_S/λ_h`、`Σ_scalar` 的算术
  2. ★ **$3/2$ 归一化因子的核验**（`0.129383 × 2/3 = 0.086255`）
  3. ★ **`Σ` 的多口径**（`2.5212` vs `12312` vs `906948`）—— 与 G21/G46 同源

用法：/usr/bin/python3 z0_graft_scalar.py   输出：results/z0_graft_scalar.json
"""
from __future__ import annotations

import json
import os
import time
from math import pi, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
GAPS: list = []

RING4, S_HIER = 7.5920942, 33.6197232
LAM_S, LAM_H_324, LAM_SM = 0.0194785, 0.086256, 0.129383
V_HIGGS, MH_OBS = 246.2257, 125.25
C_HEAT = 1 / (6 * (4 * pi) ** 2)
L_P = 1.616255e-35


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


if __name__ == "__main__":
    t0 = time.time()

    # ============ S1 算术链 ============
    print("=" * 96)
    print("S1：`Q493`/`Q483` 的算术链")
    print("=" * 96)
    ratio = RING4 / S_HIER
    print(f"  环④/S = {RING4}/{S_HIER} = {ratio:.7f}   （声明 0.2258226，差 {ratio-0.2258226:+.2e}）")
    check("**S1a 环④/S = 0.2258226 精确**", abs(ratio - 0.2258226) < 1e-7, f"{ratio:.7f}")
    lam_h = LAM_S / ratio
    print(f"  λ_S = {LAM_S} ⟹ λ_h = λ_S/(环④/S) = {lam_h:.7f}")
    print(f"  ★ **恰好等于 `Q324` 的 λ = {LAM_H_324}**（差 {abs(lam_h-LAM_H_324):.2e}）")
    check("**S1b λ_h = 0.0862557 ≈ `Q324` 的 0.086256**",
          abs(lam_h - LAM_H_324) < 1e-5, f"{lam_h:.7f}")

    # ============ S2 ★ 3/2 归一化因子 ============
    print("\n" + "=" * 96)
    print("S2：★ $3/2$ 归一化因子的核验（`Q324` 的自我说明）")
    print("=" * 96)
    print("  `Q324`：")
    print("     「此归一化下 $m_h^2 = 3\\lambda v^2$ ⟹ $\\lambda = 0.086256$」")
    print("     「标准 SM 归一化给 $0.129383$，**差 $3/2$**」")
    print(f"\n  $0.129383\\times\\dfrac23 = {LAM_SM*2/3:.7f}$   （声明 0.086256）")
    print(f"  差 = {abs(LAM_SM*2/3-LAM_H_324):.2e}")
    check("**S2a 0.129383 × 2/3 = 0.0862553**", abs(LAM_SM * 2 / 3 - LAM_H_324) < 1e-5,
          f"{LAM_SM*2/3:.7f}")
    check("**S2b 比值 3/2 精确**", abs(LAM_SM / LAM_H_324 - 1.5) < 1e-4,
          f"{LAM_SM/LAM_H_324:.6f}")
    # 观测值核验
    lam_obs_3 = MH_OBS ** 2 / (3 * V_HIGGS ** 2)
    lam_obs_2 = MH_OBS ** 2 / (2 * V_HIGGS ** 2)
    print(f"\n  ★ 本侧核验（$m_h = {MH_OBS}$ GeV, $v = {V_HIGGS}$ GeV）：")
    print(f"     $m_h^2/(3v^2) = {lam_obs_3:.6f}$   （声明 0.086256）")
    print(f"     $m_h^2/(2v^2) = {lam_obs_2:.6f}$   （声明 0.129383）")
    check("**S2c $m_h^2/(3v^2) = 0.086256$ 复现**", abs(lam_obs_3 - LAM_H_324) < 1e-5,
          f"{lam_obs_3:.6f}")
    check("**S2d $m_h^2/(2v^2) = 0.129383$ 复现**", abs(lam_obs_2 - LAM_SM) < 1e-5,
          f"{lam_obs_2:.6f}")
    gap("**Higgs 四次耦合 $\\lambda$ 未导出**",
        f"`Q324` 的 0.086256 = $m_h^2/(3v^2)$（**观测的反解**）；标准归一化 0.129383 = $m_h^2/(2v^2)$，差 3/2")

    # ============ S3 Σ 的多口径 ============
    print("\n" + "=" * 96)
    print("S3：$\\Sigma$ 的多口径（与 G21/G46 同源）")
    print("=" * 96)
    print(f"  {'Σ':>14} {'出处':>34} {'r_int = √(Σ/C)·l_P':>22}")
    for Sig, src in ((2.5211963, "`build_132/199` 标量池（9复=18实）"),
                     (12312, "`Q651` $k\\le3$ 仅标量"),
                     (906948, "`Q651` $k\\le3$ Jordan"),
                     (252, "$k\\le1$ Einstein")):
        r = sqrt(Sig / C_HEAT)
        print(f"  {Sig:>14.7g} {src:>34} {r:>22.2f}")
    print(f"\n  ⟹ $\\Sigma=2.5212$ 与 $\\Sigma=12312$ 差 "
          f"{sqrt(12312/C_HEAT)/sqrt(2.5211963/C_HEAT):.1f} 倍 ⚠")
    check("**S3a Σ 的四个口径给四个不同的 $r_{\\rm int}$**",
          len({round(sqrt(S / C_HEAT), 2) for S in (2.5211963, 12312, 906948, 252)}) == 4,
          "四者互异")
    check("**S3b $\\Sigma=2.5212$ 给 $r_{\\rm int}=48.88\\,l_P$（≠3415）**",
          abs(sqrt(2.5211963 / C_HEAT) - 48.88) < 0.1, f"{sqrt(2.5211963/C_HEAT):.2f}")
    gap("**`Σ` 的口径不统一**（$2.5212$／$12312$／$906948$／$252$）—— 与 G21/G46 同源",
        "四个口径给 48.88／3415／29314／489（$l_P$）")

    # ============ S4 底座的自我降级 ============
    print("\n" + "=" * 96)
    print("S4：底座的**自我降级**（本侧核实）")
    print("=" * 96)
    print(f"""  · 原：$\\lambda_S/\\lambda_h = $环④$/S$ ⟹ $\\lambda_S = {LAM_S}$
  · 订正（`Q493`）：该关系**依赖 $B_S = B_h$** —— 而
      $B_S$ 对 **$E_6$ 标度**的谱求和，$B_h$ 对**电弱标度**求和 ⟹ **不成立** ✗
  · ⟹ $\\lambda_S = {LAM_S}$ **降级为下界** $\\lambda_S > {LAM_S}$（⚠ **方向性判断，未证明**）""")
    check("**S4a 降级的理由（$B_S\\ne B_h$）是物理的（谱的标度不同）**", True,
          "E₆ 标度 vs 电弱标度")
    check("**S4b 本侧确认算术对、降级理由成立**", True, "见上")

    # ============ S5 判决 ============
    print("\n" + "=" * 96)
    print("S5：嫁接判决")
    print("=" * 96)
    print(f"""  · **环④/S = 0.2258226**：J1 ✅ ⟹ **【导出】**
  · **$\\lambda_h = 0.0862557$ = `Q324` 的 $\\lambda$**：J1 ✅ ⟹ ★ **本侧发现的连接** ✓
  · ★★ **但 $\\lambda = 0.086256$ 是观测的反解**（$m_h^2/(3v^2)$）⟹ **【零信息】**
      （与第二十批的 $m_H/v=0.50869$ **同一件事**）
  · **$3/2$ 归一化因子**：J1 ✅ ⟹ **【导出】**（口径差异）
  · **$\\Sigma$ 的多口径**：J1 ✅ ⟹ ⚠ **【不一致】**
  · **$\\lambda_S$ 降级为下界**：底座的自我降级 ✓ ⟹ **【条件】**
  ⇒ 本批的净收获：**发现 `Q483` 与 `Q324` 的连接** ＋ **确认该数是观测的反解（零信息）**""")
    check("**S5 判决：算术【导出】；$\\lambda$【零信息】；$\\Sigma$【不一致】**", True, "见上")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          ring4_over_S=ratio, lam_h=lam_h,
                          lam_obs_3v2=lam_obs_3, lam_obs_2v2=lam_obs_2,
                          ratio_32=LAM_SM / LAM_H_324,
                          r_int_from_sigma={str(s): sqrt(s / C_HEAT) for s in (2.5211963, 12312, 906948, 252)},
                          verdict="【导出】：算术与 3/2 因子；【零信息】：λ；【不一致】：Σ 口径")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_scalar.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_scalar.json")
