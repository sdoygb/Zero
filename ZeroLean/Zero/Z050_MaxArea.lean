/-
  Z050（|a|_max(h) = −h² + (L+1)h − L/2）的机器检查版本。

  审计修订（2026-10-14）：原版只有两个**数值代入**（`-1^2 + (4+1)*1 - 4/2 = 2` 之类），
  公式本身没被陈述。现在改成对**任意 L** 的真命题：
    · 该二次型在 `h = (L+1)/2` 处取最大（配方上界）；
    · 在顶层 `h = H = L/2` 处它恰等于 `H² = L²/4`（与 Z050 §注 的
      "|a|_max(H) = L²/4" 一致）。
-/

import Mathlib

namespace Zero

/-- **Z050 · 最大面积公式的上界**：`−h² + (L+1)h − L/2 ≤ ((L+1)/2)² − L/2`，
    等号在 `h = (L+1)/2`（配方：`= 顶点值 − (h − (L+1)/2)²`）。 -/
theorem z050_max_area_upper_bound (L h : ℝ) :
    -h ^ 2 + (L + 1) * h - L / 2 ≤ ((L + 1) / 2) ^ 2 - L / 2 := by
  nlinarith [sq_nonneg (h - (L + 1) / 2)]

/-- **Z050 · 顶层值 = H²**：在 `h = H = L/2` 处，`−h² + (L+1)h − L/2 = (L/2)² = L²/4`
    （即原文 §注 的 `|a|_max(H) = L²/4 = H²`，对任意 L 成立）。 -/
theorem z050_max_area_at_top (L : ℝ) :
    -(L / 2) ^ 2 + (L + 1) * (L / 2) - L / 2 = (L / 2) ^ 2 := by ring

#print axioms z050_max_area_upper_bound
#print axioms z050_max_area_at_top

end Zero
