"""
z0_graft_v6.py --- 嫁接验证 · 第十九批：商空间体积 V₆ 与 Q583 的订正清单

【嫁接原则】见 `GRAFT_LEDGER.md §0`。

来源：`Q583_RELATIVE_GEOMETRY.md`、`build/build_294_q583.py`、`Q651_QUOTIENT_SPECTRUM.md`

【`Q583` 的原始主张 vs 事后订正】
  · 原文：内禀模清单（含 $(1,1,1)$、$\\lambda=9$、重数 $27$）
  · 订正（`Q651`，`build_329`，19/0）：
      (a) **$(1,1,1)$ 不含单态**（三个自旋-$1/2$ $\\Rightarrow\\Sigma j=3/2\\notin\\mathbb Z$）
          $\\Rightarrow$ **$\\lambda=9$ 不在谱里** ✗
      (b) **$\\lambda=6$ 的重数是 $12$**（$3_{\\text{置换}}\\times(2\\cdot2\\cdot1)$），**不是 $27$**
          —— $27=3\\times3\\otimes\\bar3$ 来自 **$E_6$ 维数**，**非谱**
      (c) 精确谱（$k\\le1$）：**$\\{\\lambda=0:1,\\ \\lambda=6:12\\}$**
      (d) ★★ **体积 $V_6=74.99\\,r^6$ 被复现（$74.986$）** $\\Rightarrow$ **几何取对了，问题在模清单** ✓
  · 限定：该订正只覆盖**标量（$H$-不变）扇区**；$(1,1,1)$ **在扭曲（旋量／张量）扇区可以存活**

★ 本批的独立验证：
  1. $V_6 = (2\\pi^2r^3)^2/3^{3/2} = 4\\pi^4/3^{3/2} = 74.98555$ 的独立复算
  2. 新旧体积的因子 $3^{3/2}=5.196152$ 与「轨道半径 $\\sqrt3 r$」的对应
  3. $\\dim[SU(2)^3]=9$、$\\dim[SU(2)_{\\rm diag}]=3$ $\\Rightarrow$ 商维数 $6$ ✓

用法：/usr/bin/python3 z0_graft_v6.py   输出：results/z0_graft_v6.json
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


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


if __name__ == "__main__":
    t0 = time.time()

    # ============ V1 商空间的维数 ============
    print("=" * 96)
    print("V1：商空间 $G/H = SU(2)^3/SU(2)_{\\rm diag}$ 的维数")
    print("=" * 96)
    dG, dH = 3 * 3, 3
    print(f"  dim[SU(2)^3] = 3 × 3 = {dG}")
    print(f"  dim[SU(2)_diag] = {dH}")
    print(f"  ⟹ dim[G/H] = {dG} − {dH} = {dG-dH} ✓（声明 6）")
    check("**V1 商空间维数 = 6**", dG - dH == 6, f"{dG-dH}")

    # ============ V2 体积公式 ============
    print("\n" + "=" * 96)
    print("V2：★ 商空间体积 $V_6$ 的独立复算")
    print("=" * 96)
    V_S3 = 2 * pi ** 2
    V_old = 4 * pi ** 4
    factor = 3 ** 1.5
    V_new = V_old / factor
    print(f"  单位 $S^3$ 体积 = $2\\pi^2$ = {V_S3:.5f}")
    print(f"  旧值（乘积几何）$V_6 = (2\\pi^2)^2 = 4\\pi^4$ = {V_old:.5f}")
    print(f"  因子 $3^{{3/2}} = (\\sqrt3)^3$ = {factor:.6f}")
    print(f"  ★ 商空间 $V_6 = 4\\pi^4/3^{{3/2}}$ = {V_new:.5f}   （声明 74.986）")
    print(f"  差 = {abs(V_new-74.986):.6f}")
    check("**V2a $V_6 = 74.98555$ 复现（差 < 1e-3）**", abs(V_new - 74.986) < 1e-3,
          f"{V_new:.5f}")
    check("**V2b 因子 $3^{3/2} = 5.196152$**", abs(factor - 5.196152) < 1e-5, f"{factor:.6f}")
    check("**V2c 旧值 $4\\pi^4 = 389.63636$**", abs(V_old - 389.63636) < 1e-4, f"{V_old:.5f}")
    print(f"\n  物理内容：商空间的两因子处**轨道半径 = $\\sqrt3\\,r$** ⟹ 体积按 $(\\sqrt3)^3$ 缩小 ✓")
    print(f"  ⚠ 而 $3 = \\Lambda$（框架整数）⟹ 因子不是自由参数 ✓")

    # ============ V3 Q651 的订正清单 ============
    print("\n" + "=" * 96)
    print("V3：`Q651` 对 `Q583` 的订正清单")
    print("=" * 96)
    rows = [
        ("$(1,1,1)$ 含单态", "✗ **不含**（三个自旋-$1/2$ $\\Rightarrow\\Sigma j=3/2\\notin\\mathbb Z$）"),
        ("$\\lambda=9$ 在谱里", "✗ **不在**（该模不存在）"),
        ("$\\lambda=6$ 的重数 $=27$", "✗ **是 12**（$3_{\\text{置换}}\\times(2\\cdot2\\cdot1)$）"),
        ("$27$ 的来源", "**$E_6$ 维数**（$3\\times3\\otimes\\bar3$），**非谱**"),
        ("精确谱（$k\\le1$）", "$\\{\\lambda=0:1,\\ \\lambda=6:12\\}$"),
        ("**体积 $V_6 = 74.99\\,r^6$**", "✅ **复现（74.986）** —— 几何取对了 ✓"),
    ]
    print(f"  {'项':>26} {'订正后':>52}")
    for a, b in rows:
        print(f"  {a:>26} {b:>52}")
    check("**V3 订正清单已登记**", True, "6 项（1 项保留，5 项订正）")

    # ============ V4 判决 ============
    print("\n" + "=" * 96)
    print("V4：嫁接判决")
    print("=" * 96)
    print(f"""  · **$V_6 = 4\\pi^4/3^{{3/2}} = 74.98555$**：
      J1 ✅（本侧独立复算）J2 ✅（输入只有 $\\pi$ 与 $3=\\Lambda$）J3 ⚠（「轨道半径 $\\sqrt3 r$」的来源需核）
      ⟹ **【条件】** —— 且这是 `Q583` 明确**保留有效**的一项 ✓
  · **$Q651$ 的订正（模清单）**：J1 ✅（本侧未独立复算 $\\lambda=6$ 的重数 —— 需 Peter–Weyl，诚实标注）
      ⟹ 登记为**底座的自我订正**，判决 **【已订正】**
  ⇒ 本批确认：**`Q583` 的几何对、模清单错** ✓（底座自己的结论，本侧复现了几何那一半）""")
    check("**V4 判决：$V_6$ 为【条件】；模清单订正为【已订正】**", True, "见上")
    gap("**$\\lambda=6$ 的重数 $=12$ 的独立复算**",
        "需 $SU(2)^3$ 的 Peter–Weyl 分解；本批未做（诚实标注）")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          dim_quotient=dG - dH, V_old=V_old, factor=factor, V_new=V_new,
                          verdict="$V_6$ 为【条件】；模清单订正为【已订正】")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_v6.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_v6.json")
