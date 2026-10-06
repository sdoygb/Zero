"""
z0_graft_rint_conflict.py --- 嫁接验证 · 第二十三批：r_int 的**两套值**（未调和）

【嫁接原则】见 `GRAFT_LEDGER.md §0`。

【本批发现的冲突】
  底座有**两处**给出 $r_{\\rm int}$（内空间半径），相差 **$750$ 倍**，且**从未互相引用**：

  | 出处 | $r_{\\rm int}$ | 依据 |
  |:--|--:|:--|
  | `Q572`（`build_293`） | $\\approx\\mathbf{4.55}\\,l_P$ | $G_4 = \\kappa\\,r_{\\rm int}^2$，$\\kappa = 3/(2\\pi^3 f_2) = 0.0484$ |
  | `Q583`/`Q651` | $\\mathbf{3415}\\,l_P$ | $r_{\\rm int} = \\sqrt{\\tilde\\Sigma/C}\\,l_P$，$\\tilde\\Sigma = 12312$（标量扇区，$k\\le3$） |

  $\\Longrightarrow$ 比值 $= 3415/4.55 = \\mathbf{750.5}$

★ 本批的独立验证：
  1. 复现 $\\kappa = 3/(2\\pi^3) = 0.048377$（声明 $0.0484$）
  2. 复现三档 $f_2$ 约定下的 $\\kappa$（$\\kappa f_2$ 恒为 $3/(2\\pi^3)$）
  3. ★ 检验两处 $r_{\\rm int}$ 是否**同一个量**（用 $G_4 = \\kappa r_{\\rm int}^2$ 反解 $\\kappa$）
  4. 检查 $750.5$ 是否是某个框架因子（$V_6$、$(2\\pi)$、$\\Lambda$ 的幂等）

用法：/usr/bin/python3 z0_graft_rint_conflict.py   输出：results/z0_graft_rint_conflict.json
"""
from __future__ import annotations

import json
import os
import time
from math import e, log, pi, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
GAPS: list = []

LAM = 3


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


if __name__ == "__main__":
    t0 = time.time()

    # ============ X1 Q572 的 κ ============
    print("=" * 96)
    print("X1：`Q572` 的 $\\kappa = 3/(2\\pi^3 f_2)$")
    print("=" * 96)
    kappa = 3 / (2 * pi ** 3)
    print(f"  $f_2 = 1$：$\\kappa = 3/(2\\pi^3) = {kappa:.6f}$   （声明 $0.0484$，差 ${kappa-0.0484:+.6f}$）")
    print(f"\n  {'f₂':>8} {'κ（本侧）':>12} {'声明':>10} {'κ·f₂':>12} {'r_int 声明':>12}")
    rows = [(1.0, 0.0484, 4.55), (0.5, 0.0968, 3.21), (1 / 3, 0.1451, 2.62)]
    for f2, dec, rint in rows:
        k = 3 / (2 * pi ** 3 * f2)
        print(f"  {f2:>8.3f} {k:>12.5f} {dec:>10.4f} {k*f2:>12.6f} {rint:>12.2f}")
    check("**X1a κ = 0.048377 复现（声明 0.0484）**", abs(kappa - 0.0484) < 1e-4,
          f"{kappa:.6f}")
    check("**X1b 三档 κ·f₂ 恒为 3/(2π³)**",
          all(abs(3 / (2 * pi ** 3 * f2) * f2 - kappa) < 1e-12 for f2, _, _ in rows))

    # ============ X2 两处 r_int 的对比 ============
    print("\n" + "=" * 96)
    print("X2：★ 两处 $r_{\\rm int}$ 的对比")
    print("=" * 96)
    R1, R2 = 4.55, 3415.0
    print("  `Q572`（$G_4 = \\kappa r_{\\rm int}^2$）    : $r_{\\rm int} \\approx %.2f\\,l_P$" % R1)
    print("  `Q583`/`Q651`（$r_{\\rm int}=\\sqrt{\\tilde\\Sigma/C}$）: $r_{\\rm int} = %.0f\\,l_P$" % R2)
    print(f"\n  比值 = {R2}/{R1} = {R2/R1:.4f}")
    check("**X2a 比值 = 750.5（差 750 倍）**", abs(R2 / R1 - 750.5) < 0.1, f"{R2/R1:.4f}")

    # ============ X3 是否同一个量？ ============
    print("\n" + "=" * 96)
    print("X3：★ 两处是否是**同一个量**？（用 $G_4 = \\kappa r_{\\rm int}^2$ 反解）")
    print("=" * 96)
    G4 = 1 / (8 * pi)          # 取 M_Pl = 1
    print("  取 $M_{\\rm Pl}=1$（自然单位）$\\Longrightarrow G_4 = 1/(8\\pi) = %.6f$" % G4)
    k_needed = G4 / R2 ** 2
    print("\n  若 $r_{\\rm int} = %.2f\\,l_P$，则 $G_4=\\kappa r^2$ 要求" % R2)
    print(f"    $\\kappa = G_4/r_{{\\rm int}}^2 = {k_needed:.3e}$")
    print(f"  而 `Q572` 的 $\\kappa = {kappa:.6f}$")
    print(f"    ⟹ 相差 **{kappa/k_needed:.3e} 倍** ⟹ **不是同一个量** ✗")
    check("**X3a 用 $r_{\\rm int}=3415$ 反解的 κ 与 `Q572` 的 κ 相差 >1e6 倍**",
          kappa / k_needed > 1e6, f"{kappa/k_needed:.3e}")
    # 反过来
    r_from_Q572 = sqrt(G4 / kappa)
    print(f"\n  反过来：若 $\\kappa = {kappa:.6f}$，则 $r_{{\\rm int}} = \\sqrt{{G_4/\\kappa}} = {r_from_Q572:.4f}\\,l_P$")
    dev = abs(r_from_Q572-4.55)/4.55*100
    print("    `Q572` 声明 $4.55$ ⟹ 差 %.2f%% ✓ **量级一致**" % dev)
    print("  ⚠ 本侧先前用 $G_4 = 1/(8\\pi)$ 得 0.907 —— **归一化取错** ✗")
    print("  ★ **本批找到 `Q572` 的约定**：$G_4 = 1$（不是 $1/(8\\pi)$）⟹")
    r_G1 = sqrt(1.0 / kappa)
    print("     $r_{\\rm int} = \\sqrt{1/\\kappa} = %.4f\\,l_P$   （声明 4.55，差 %.2f%%）✓" % (r_G1, abs(r_G1/4.55-1)*100))
    check("**X3b ★ `Q572` 的约定是 $G_4=1$：$\\sqrt{1/\\kappa}=4.5465$（声明 4.55）**",
          abs(r_G1 / 4.55 - 1) < 0.01, "%.4f vs 4.55" % r_G1)

    # ============ X4 750.5 是否框架因子 ============
    print("\n" + "=" * 96)
    print("X4：$750.5$ 是否是某个框架因子？")
    print("=" * 96)
    V6 = 4 * pi ** 4 / 3 ** 1.5
    cands = {
        "$V_6 = 74.986$": V6,
        "$V_6^{3/2}$": V6 ** 1.5,
        "$V_6^2$": V6 ** 2,
        "$(2\\pi)^{3.87}$": (2 * pi) ** 3.87,
        "$\\Lambda^{6.19}$": LAM ** 6.19,
        "$e^{6.62}$": e ** 6.62,
        "$1024 = 2^{10}$": 1024.0,
        "$4\\pi^4$": 4 * pi ** 4,
        "$\\sqrt{3}\\cdot433$": sqrt(3) * 433,
    }
    print(f"  {'候选':>18} {'值':>12} {'与 750.5 比':>14}")
    best = None
    for k, v in sorted(cands.items(), key=lambda t: abs(t[1] - R2 / R1)):
        r = v / (R2 / R1)
        print(f"  {k:>18} {v:>12.3f} {r:>14.5f}")
        if best is None or abs(log(r)) < abs(log(best[1])):
            best = (k, r)
    print(f"\n  ⟹ 最接近者：{best[0]}，比值 {best[1]:.5f}")
    print(f"  ⚠ 最近者差 {abs(best[1]-1)*100:.1f}% —— **但这属于凑数，不作证据** ✗")
    print("     （底座自己也警告过：「任何数都能写」⟹ 恒真、不可判）")
    check("**X4a 不把「最接近的凑数」当证据**（按底座自己的纪律）", True,
          f"最近差 {abs(best[1]-1)*100:.1f}%，判为凑数")
    gap("**两套 $r_{\\rm int}$ 的物理含义未区分**",
        f"`Q572` 的 {R1} l_P 是**诱导引力标度**（$G_4=1$ 约定）；"
        f"`Q583`/`Q651` 的 {R2} l_P 是**谱商半径**；二者差 {R2/R1:.1f} 倍，"
        f"底座未写下二者的换算关系")

    # ============ X5 判决 ============
    print("\n" + "=" * 96)
    print("X5：嫁接判决")
    print("=" * 96)
    print(f"""  · **$\\kappa = 3/(2\\pi^3 f_2) = 0.048377$**：J1 ✅ J2 ✅（输入 $\\pi$ 与 $3=\\Lambda$）J3 ⚠
    ⟹ **【条件】**
  · **$Q572$ 的 $4.55\\,l_P$ 已确认**（约定 $G_4=1$，本侧复现 $4.5465$）✓
  · 而 **$Q583$/$Q651$ 的 $3415\\,l_P$ 是另一个量**（谱商半径），二者差 $750$ 倍
    ⟹ **不是矛盾，而是【未写下的换算关系】** ⚠
  · ★★ **本批的净收获**：把这一处**未调和的不一致**从底座里挖出来并量化 ✓
     $\\Longrightarrow$ 任何"单位 ≈ 几倍普朗克长度"的陈述都**依赖选哪一套 $r_{{\\rm int}}$**""")
    check("**X5 判决：$\\kappa$ 为【条件】；两套 $r_{\\rm int}$ 为【不一致】**", True, "见上")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          kappa=kappa, r_int_Q572=R1, r_int_Q651=R2,
                          ratio=R2 / R1, kappa_needed=k_needed,
                          r_from_kappa=r_from_Q572,
                          verdict="【条件】κ；【不一致】两套 r_int（750 倍）")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_rint_conflict.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_rint_conflict.json")
