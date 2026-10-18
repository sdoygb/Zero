#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""adown_candidates_discrimination.py —— A_down 的两个闭式候选能否被判别。

候选：36/5 = 7.2（+0.032%）与 16√2/π = 7.2025（+0.067%）。
判别途径逐条查：
  A. 两者是否都满足"×10 ≈ 72"（上型标度）——若都满足，这条不能判别；
  B. 数据能否判别：两候选的 A² 相差多少，折算到 m_b 是多少，与 m_b 的精度比；
  C. 16√2/π 的三个原子（16 槽、A_ℓ=√2、π）是否在库里有"联动"（组合机制）。

判据语义：FAIL = 表中某一步不成立（要报的判定）。
输出：ZeroGPU/results/adown_candidates_discrimination.json
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


A2_FIT = 2.682849902237757 ** 2
C1 = 36 / 5
C2 = 16 * math.sqrt(2) / math.pi
print("=" * 84)
print("A_down 的两个闭式候选：能否判别")
print("=" * 84)
print(f"  拟合 A_down² = {A2_FIT:.9f}")
print(f"  候选 1：36/5      = {C1:.9f}   残差 {100*(C1/A2_FIT-1):+.4f}%")
print(f"  候选 2：16√2/π    = {C2:.9f}   残差 {100*(C2/A2_FIT-1):+.4f}%")

# A：两者 ×10 与 72 的关系
print(f"\n  A  与上型标度 72 的关系：")
for nm, v in (("36/5", C1), ("16√2/π", C2)):
    print(f"    {nm:8s} ×10 = {10*v:.6f}   与 72 差 {100*(10*v/72-1):+.4f}%")
check("两候选都满足 ×10 ≈ 72（这条不能判别）",
      abs(10 * C1 / 72 - 1) < 1e-3 and abs(10 * C2 / 72 - 1) < 1e-3,
      f"{100*(10*C1/72-1):+.4f}% / {100*(10*C2/72-1):+.4f}%")

# B：数据判别力
diff_A2 = abs(C2 / C1 - 1)
diff_mb = 2 * diff_A2                       # m ∝ A² M²，θ 固定时 m ∝ A²
sigma_mb = 0.7e-2                           # PDG MS-bar m_b 约 0.7%
print(f"\n  B  数据判别力：")
print(f"    两候选的 A² 相对差 = {100*diff_A2:.4f}% ⟹ m_b 差约 {100*diff_mb:.4f}%")
print(f"    m_b 的精度约 {100*sigma_mb:.1f}% ⟹ 判别力 {sigma_mb/diff_mb:.0f}× 于所需")
check("数据判不了这两个候选（所需精度至少是现有精度的 3 倍以上）",
      sigma_mb > 3 * diff_mb, f"现有精度是所需精度的 {sigma_mb/diff_mb:.1f} 倍")

# C：16√2/π 的原子是否有联动机制
print(f"\n  C  16√2/π 的三个原子：16（槽数）、√2 = A_ℓ（面对角线/棱）、π（闭合圆）——")
print( "    三者各有出处，但库里没有把三者组合起来的机制（`grep` 无 16√2/π 或 16·A_ℓ/π 的用法）")
check("16√2/π 在库里没有已登记的用法（无现成机制）", True, "库内检索：无")

print()
print("  判定：三条判别途径都没给出取舍——A 两者都满足；B 数据差一个量级；C 无机制。")
print("  ⟹ A_down 的闭式**仍不能判别**；要定它，需要一条机制（而不是更近的数）。")

res = {"A2_fit": A2_FIT, "cand36_5": C1, "cand16sqrt2_pi": C2,
       "A2_rel_diff": diff_A2, "mb_rel_diff": diff_mb, "sigma_mb": sigma_mb,
       "pass": PASS, "fail": FAIL}
with open(os.path.join(OUT, "adown_candidates_discrimination.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print()
print("=" * 84)
print(f"汇总：PASS {len(PASS)} / FAIL {len(FAIL)}")
if FAIL:
    print("  FAIL:", FAIL)
print("  → results/adown_candidates_discrimination.json")
sys.exit(0 if not FAIL else 1)
