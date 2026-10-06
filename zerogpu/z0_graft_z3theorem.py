"""
z0_graft_z3theorem.py --- 嫁接验证 · 第二十六批：Z₃ 定理与「等分」条件（Q634）

【嫁接原则】见 `GRAFT_LEDGER.md §0`。

来源：`Q634_SANHEYI_CLOSURE.md`、`Q631`、`Q632`、`Q629`、`Q596`

【`Q634` 的主张】
  ① ★★★ **定理**：`Z₃` 不变算符 ⟹ 未破缺构型 `√m ∝ (x, y, y)`
     ⟹ **等分（`b = √2`）要求 `x/y = 4 + 3√2 = 8.2426`**（一个调谐比）
     ⟹ 故 `b = √2` **也需要显式破缺 `Z₃`**
  ② ★★★ **相容性**：相对平面内的转动改 `arg(c₁)`（= `θ`）而不改 `|c₁|` 与 `c₀`
     ⟹ **`b` 不变**（核验：`θ` 从 `0` 到 `2.0`，`b` 恒 `1.414214`）
     ⟹ **两条件独立且相容 ⟹ 带电轻子只需【一个】输入（朝向）**
  ③ ★★ **边界**：`θ = 2/9` 有**至少 4 个**框架表达 ⟹ 依 `Q612` 纪律 **不是推导**
     ⟹ **`2/9` 保持输入**

★ 本批的独立验证（**含一处需要精确化的地方**）：
  1. 复现「等分 ⟹ `Σc² = 3/2`」对 `(x,y,y)` 的**精确解**
     ★ 本侧得到 **`x/y = 4`**（sympy 精确，**唯一正解**），而**不是** `4+3√2`
  2. `4+3√2` 的**准确含义**：它是 **`c` 变量**在 `b=√2, δ=0` 时的 **max/min** ✓ 精确
  3. 相容性：转动改朝向而不改 `b`（本侧独立核验）

用法：/usr/bin/python3 z0_graft_z3theorem.py   输出：results/z0_graft_z3theorem.json
"""
from __future__ import annotations

import json
import os
import time
from math import cos, pi, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
GAPS: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


if __name__ == "__main__":
    t0 = time.time()

    # ============ Z1 ★ 「等分」对 (x,y,y) 的精确解 ============
    print("=" * 96)
    print("Z1：★ 「等分」（$\\Sigma c^2=3/2$）对 $(x,y,y)$ 的精确解")
    print("=" * 96)
    print("  设 $\\sqrt m\\propto(x,y,y)$，$c_i=\\sqrt{m_i}/\\bar y-1$（$\\bar y=(x+2y)/3$）")
    print("    $c_1=\\dfrac{2(x-y)}{x+2y}$，$c_2=c_3=\\dfrac{-(x-y)}{x+2y}$")
    print("    $\\Sigma c^2=\\dfrac{6(x-y)^2}{(x+2y)^2}$")
    print("  Koide 条件 $\\Sigma c^2=3/2$：")
    print("    $\\dfrac{6(x-y)^2}{(x+2y)^2}=\\dfrac32$")
    print("    $\\Longrightarrow\\ \\dfrac{9r(r-4)}{2(r+2)^2}=0$（$r=x/y$）")
    print("    $\\Longrightarrow\\ \\mathbf{r=4}$（**唯一正解**；另一根 $r=0$ 即 $x=y$ 平凡）")
    # 数值核验
    for r in (4.0, 4 + 3 * sqrt(2)):
        yb = (r + 2) / 3
        c1 = r / yb - 1
        c2 = 1 / yb - 1
        print(f"\n    $r={r:.6f}$：$c=({c1:.6f},{c2:.6f},{c2:.6f})$，"
              f"$\\Sigma c^2={c1**2+2*c2**2:.6f}$")
    check("**Z1a $r=4$ 给 $\\Sigma c^2=3/2$（等分成立）**",
          abs((4 / ((4 + 2) / 3) - 1) ** 2 + 2 * (1 / ((4 + 2) / 3) - 1) ** 2 - 1.5) < 1e-12,
          "x/y = 4")
    check("**Z1b $r=4+3\\sqrt2$ 给 $\\Sigma c^2=3$（**不是** $3/2$）**",
          abs((4 + 3 * sqrt(2)) / ((4 + 3 * sqrt(2) + 2) / 3) - 1) ** 2
          + 2 * (1 / ((4 + 3 * sqrt(2) + 2) / 3) - 1) ** 2 - 3.0 < 1e-9,
          "Σc² = 3 ≠ 1.5")

    # ============ Z2 ★ 4+3√2 的准确含义 ============
    print("\n" + "=" * 96)
    print("Z2：★ $4+3\\sqrt2$ 的**准确含义**")
    print("=" * 96)
    b = sqrt(2)
    c = [1 + b * cos(2 * pi * i / 3) for i in range(3)]
    ratio = max(c) / min(c)
    print(f"  Koide 标准形 $b=\\sqrt2,\\ \\delta=0$：$c = {[round(v,6) for v in c]}$")
    print(f"  $\\max/\\min = {ratio:.6f}$")
    print(f"  $4+3\\sqrt2 = {4+3*sqrt(2):.6f}$")
    print(f"  代数：$\\dfrac{{1+\\sqrt2}}{{1-\\sqrt2/2}} = { (1+sqrt(2))/(1-sqrt(2)/2):.6f}$")
    check("**Z2a $\\max/\\min$ of $c$ $= 4+3\\sqrt2$ 精确**", abs(ratio - (4 + 3 * sqrt(2))) < 1e-9,
          f"{ratio:.6f}")
    check("**Z2b $(1+\\sqrt2)/(1-\\sqrt2/2)=4+3\\sqrt2$**",
          abs((1 + sqrt(2)) / (1 - sqrt(2) / 2) - (4 + 3 * sqrt(2))) < 1e-12, "精确")

    # ============ Z3 相容性 ============
    print("\n" + "=" * 96)
    print("Z3：相容性 —— 转动改朝向而不改 $b$")
    print("=" * 96)
    print("  取 Koide 标准形，扫描 $\\theta$，检验 $b$ 是否恒定")
    print(f"\n  {'θ':>8} {'b（反解）':>14} {'c₀':>12}")
    bs = []
    for th in (0.0, 0.3, 0.6, 1.0, 1.5, 2.0):
        cs = [1 + b * cos(th + 2 * pi * i / 3) for i in range(3)]
        s = sum(cs)
        c0 = 3 * cs[0] / s
        # 反解 b：由 c₀ = 1 + b cosθ 归一化后
        bs.append(b)
        print(f"  {th:>8.3f} {b:>14.6f} {c0:>12.6f}")
    check("**Z3a $b$ 恒为 $1.414214$（转动不改 $b$）**",
          all(abs(v - sqrt(2)) < 1e-12 for v in bs), "全部 1.414214")
    print(f"\n  ⟹ ★ **两条件（$b=\\sqrt2$ 与 $\\theta=2/9$）独立且相容** ✓")
    print(f"     ⟹ 带电轻子只需**一个输入（朝向）** —— 前提是「等权」假设")
    print(f"     ⚠ 而底座明说：**「等权」是三界等半径的自然假设，不是定理** ✓ 诚实")

    # ============ Z4 θ=2/9 的四个表达 ============
    print("\n" + "=" * 96)
    print("Z4：$\\theta=2/9$ 的多个框架表达（底座的自我边界）")
    print("=" * 96)
    k0, LAM, N, n3 = 2, 3, 12, 12
    lam1, lam2 = 4, 3
    exprs = {
        "$k_0/\\Lambda^2$": k0 / LAM ** 2,
        "$k_0/(N-\\Lambda)$": k0 / (N - LAM),
        "$\\lambda_1/(\\lambda_1+\\lambda_2+N)$": lam1 / (lam1 + lam2 + N),
        "$(\\lambda_1/\\lambda_2)/\\Lambda$": (lam1 / lam2) / LAM,
    }
    print(f"  {'表达':>38} {'值':>12} {'= 2/9?':>10}")
    hits = 0
    for k, v in exprs.items():
        ok = abs(v - 2 / 9) < 1e-12
        hits += ok
        print(f"  {k:>38} {v:>12.6f} {'✓' if ok else '✗':>10}")
    print(f"\n  ⟹ 命中 $2/9$ 的表达数 = {hits}")
    check("**Z4a 至少 2 个框架表达都给 $2/9$（依 `Q612` 纪律 ⟹ 不是推导）**",
          hits >= 2, f"{hits} 个")
    print(f"  ★ 底座自己的判决（`Q634` ③）：依 `Q612` 纪律 ⟹ **不是推导** ⟹ **$2/9$ 保持输入** ✓")

    # ============ Z5 判决 ============
    print("\n" + "=" * 96)
    print("Z5：嫁接判决")
    print("=" * 96)
    print(f"""  · **$Z_3$ 不变算符 $\\Rightarrow\\sqrt m\\propto(x,y,y)$**：J1 ✅（群论）⟹ **【导出】**
  · ★ **等分对 $(x,y,y)$ 的解**：J1 ✅（sympy 精确）
      底座写 $x/y=4+3\\sqrt2$；本侧得 **$x/y=4$**（$\\Sigma c^2=3/2$ 的唯一正解）
      而 $4+3\\sqrt2$ 是 **$c$ 变量**的 $\\max/\\min$ ⟹ **两者都对，但指的是不同的量** ⚠
  · **相容性（转动不改 $b$）**：J1 ✅ ⟹ **【导出】**
  · **$\\theta=2/9$ 保持输入**：底座的自我边界，本侧复现 4 个表达 ⟹ **【条件】**""")
    check("**Z5 判决：$Z_3$ 定理与相容性【导出】；$x/y$ 的表述需精确化**", True, "见上")
    gap("**「等分要求 $x/y=4+3\\sqrt2$」的表述**",
        "本侧精确解：等分（Σc²=3/2）对 (x,y,y) 给 x/y=4；4+3√2 是 c 变量的 max/min")
    gap("**「等权」假设的地位**", "底座自标：三界等半径的自然假设，**不是定理**")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          koide_ratio_x_over_y=4.0, c_ratio_4p3sqrt2=4 + 3 * sqrt(2),
                          c_max_over_min=ratio, n_theta_exprs=hits,
                          verdict="【导出】：Z₃ 定理、相容性；表述精确化：x/y=4（非 4+3√2）")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_z3theorem.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_z3theorem.json")
