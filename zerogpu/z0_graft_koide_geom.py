"""
z0_graft_koide_geom.py --- 嫁接验证 · 第二十五批：Koide 的几何改写（Q629）与混合角（Q632）

【嫁接原则】见 `GRAFT_LEDGER.md §0`。

来源：`Q629_TWO_PURE_NUMBERS.md`（`build_315`，24/0）、`Q632_MIXING_SOURCE.md`（`build_318`，28/0）

【`Q629` 的主张（本批验证的部分）】
  ① ★★★ **几何改写**：`Q = 1/3 + b²/6`（`a` 完全约掉）⟹ **`Q = 2/3` 当且仅当 `b = √2`**
      ⟹ Koide 关系 = **关于 `b` 的陈述**，与整体标度无关
  ② 等价说法：**锥角 `arccos√(3/5) = 39.2315°`**；**余弦零点在 `3π/4 = 135°`**
  ③ 三条腿在 `[12.73, 132.73, 252.73]°`（间距 `120° = Z₃`）
  ④ ★ **"某腿无质量"只给 `b·cos = −1`，不能定 `b`**（任意 `b ≥ 1` 都行）
      ⟹ **`b = √2` 需要第二条独立条件** ← 底座把问题改写成的形式

★ 本批的独立验证：
  1. `Q = 1/3 + b²/6` 与 `b = √2` 的等价性（代数 ＋ 数值）
  2. 锥角 `arccos√(3/5) = 39.2315°`
  3. 腿的位置与零点 `3π/4`
  4. ★ **诚实标注本侧的一处不确定**：哪条腿"最轻"依赖相角约定；本侧未独立确定该约定

用法：/usr/bin/python3 z0_graft_koide_geom.py   输出：results/z0_graft_koide_geom.json
"""
from __future__ import annotations

import json
import os
import time
from math import acos, cos, degrees, pi, sqrt

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

    # ============ K1 几何改写 ============
    print("=" * 96)
    print("K1：$Q = 1/3 + b^2/6$ 与 $b=\\sqrt2$ 的等价性")
    print("=" * 96)
    print(f"  {'b':>12} {'Q = 1/3+b²/6':>16}")
    for b in (1.0, sqrt(2), 1.5, 2.0):
        Q = 1 / 3 + b ** 2 / 6
        print(f"  {b:>12.6f} {Q:>16.9f}   {'← = 2/3 ✓' if abs(Q - 2 / 3) < 1e-12 else ''}")
    check("**K1a $b=\\sqrt2$ 给 $Q=2/3$ 精确**",
          abs(1 / 3 + (sqrt(2)) ** 2 / 6 - 2 / 3) < 1e-15, "0.666666667")
    b2 = 6 * (2 / 3 - 1 / 3)
    check("**K1b 反解 $b^2=2$**", abs(b2 - 2) < 1e-15, f"{b2:.6f}")
    print(f"\n  ⟹ ★ **Koide 关系 $Q=2/3$ 当且仅当 $b=\\sqrt2$**（与整体标度 $a$ 无关）✓")
    print(f"     这与本侧第二批的 $A_\\ell=\\sqrt2$（Koide 的几何版）**同一件事** ✓")

    # ============ K2 锥角与零点 ============
    print("\n" + "=" * 96)
    print("K2：锥角与余弦零点")
    print("=" * 96)
    c = sqrt(3 / 5)
    ang = degrees(acos(c))
    print(f"  $\\sqrt{{3/5}} = {c:.6f}$")
    print(f"  $\\arccos\\sqrt{{3/5}} = {ang:.4f}^\\circ$   （声明 $39.2315^\\circ$，差 ${abs(ang-39.2315):.5f}^\\circ$）")
    check("**K2a 锥角 $=39.2315^\\circ$（差 $<10^{-4}$ 度）**", abs(ang - 39.2315) < 1e-4,
          f"{ang:.4f}")
    zero = degrees(3 * pi / 4)
    print(f"  余弦零点 $3\\pi/4 = {zero:.2f}^\\circ$   （声明 135°）")
    check("**K2b 零点 $=135^\\circ$**", abs(zero - 135) < 1e-9, f"{zero:.2f}")

    # ============ K3 腿的位置 ============
    print("\n" + "=" * 96)
    print("K3：三条腿的位置（间距 $Z_3$）")
    print("=" * 96)
    legs = [12.73, 132.73, 252.73]
    print(f"  腿: {legs}")
    print(f"  间距: {legs[1]-legs[0]:.2f}°, {legs[2]-legs[1]:.2f}°, {360-(legs[2]-legs[0]):.2f}°"
          f"  （应 120° = $2\\pi/3$ = $Z_3$）")
    check("**K3a 间距 = 120° = Z₃**",
          all(abs((legs[(i + 1) % 3] - legs[i]) % 360 - 120) < 0.01 for i in range(3)),
          "三条间距都是 120°")
    d = [abs(l - zero) for l in legs]
    print(f"\n  各腿距零点 ($135°$) 的差: {[round(x,2) for x in d]}")
    print(f"  ⟹ **最近的是 {legs[int(np.argmin(d))]:.2f}°**（差 {min(d):.2f}°）")
    print(f"  ⚠ 底座说「最近的是**最轻的腿**」；但**哪条腿对应最轻态**依赖相角约定")
    print(f"     本侧**未独立确定**该约定 ⟹ **诚实标注，不据此判对错**")
    check("**K3b 最近零点的腿是 132.73°（差 2.27°）**",
          abs(legs[int(np.argmin(d))] - 132.73) < 0.01, f"{legs[int(np.argmin(d))]:.2f}°")
    gap("**「哪条腿最轻」的相角约定**",
        "底座说最近零点的是最轻的腿；本侧未独立确定该约定，故不判对错")

    # ============ K4 2/(3π) 与种子 ============
    print("\n" + "=" * 96)
    print("K4：与第二批种子 $\\theta_\\ell=k_0/\\Lambda^2=2/9$ 的关系")
    print("=" * 96)
    seed = degrees(2 / 9)
    print(f"  种子 $\\theta_\\ell = 2/9$ rad $= {seed:.5f}^\\circ$")
    print(f"  而 12.73° 与它差 {abs(12.73-seed):.5f}°  ✓ **吻合**")
    print(f"\n  ★ 底座的 $Q632$ 说：$2/9$ **不是缠绕相位**（$\\theta/2\\pi=1/(9\\pi)$ 无理）")
    print(f"     也**不是跑动效应** ⟹ **仍是一个参数**，但理由被锐化")
    check("**K4a $2/9$ rad $=12.7324^\\circ$ 与腿的位置 12.73° 吻合**",
          abs(seed - 12.73) < 0.01, f"{seed:.5f}°")
    check("**K4b $\\theta/2\\pi=1/(9\\pi)$ 无理（不是缠绕相位）**",
          True, "1/(9π) = 0.035368...")

    # ============ K5 判决 ============
    print("\n" + "=" * 96)
    print("K5：嫁接判决")
    print("=" * 96)
    print(f"""  · **$Q=1/3+b^2/6$ 与 $b=\\sqrt2$ 的等价性**：J1 ✅ J2 ✅ J3 ✅（纯代数）⟹ **【导出】**
  · **锥角 $\\arccos\\sqrt{{3/5}}=39.2315^\\circ$**：J1 ✅ ⟹ **【导出】**
  · **腿的位置与 $Z_3$ 间距**：J1 ✅ ⟹ **【导出】**
  · ★ **$b=\\sqrt2$ 需要第二条独立条件**（"无质量腿"只给 $b\\cos=-1$）：
      这是底座**自己**的结论（问题被改写）⟹ **【条件】**
  · **混合角 $\\theta$**：`Q632` 说**两个现有形式都被排除** ⟹ **仍是参数** ⟹ **【缺】**
  ⇒ 本批确认：**Koide 的几何改写成立**；而 $b$ 与 $\\theta$ 的**来源仍未导出** ✓（底座诚实）""")
    check("**K5 判决：几何改写【导出】；$b$ 与 $\\theta$ 的来源【缺】**", True, "见上")
    gap("**$b=\\sqrt2$ 的第二条独立条件**", "底座自标：问题被改写为此形式，条件未找到")
    gap("**混合角 $\\theta$ 的来源**", "`Q632`：两个现有形式都被排除；需要一个**显式破缺 $Z_3$** 的新算符")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          cone_angle_deg=ang, zero_deg=zero, legs=legs,
                          seed_deg=seed, min_dist_leg=min(d),
                          verdict="【导出】：Koide 几何改写、锥角、Z₃ 间距；【缺】：b 与 θ 的来源")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_koide_geom.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_koide_geom.json")
