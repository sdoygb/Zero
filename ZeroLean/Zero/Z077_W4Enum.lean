/-
  Z077（W₄ 的划分给 q = 5/9）。

  审计修订（2026-10-14）：原版写 `(2/3)² + (1/3)² = 5/9`（纯算术，没有"划分"这个对象），
  另一条 `z077_nflip_q` 更是 `True := by trivial`。现在改成**真的用 W₄ 的层划分算出来**：
  层计数 (b₁,b₂) = (4,2) 来自 `words4` 与层高函数 `layer`（Z079_Ledger），q = Σ(b_h/6)² 。
-/

import Mathlib
import Zero.Z079_Ledger

namespace Zero

/-- **Z077 · hmax 划分的 q = 5/9**：用 `words4`（W₄ 的 6 个词）与层高 `layer` 真的数出
    层计数 (b₁,b₂) = (4,2)，再算 `q = (b₁/6)² + (b₂/6)² = 5/9`。
    与 Z079 的 `z079_kappa1` 一致（同一划分、同一数值，这里从划分本身算出）。 -/
theorem z077_layer_partition_q :
    (((words4.filter (fun w => layer w = 1)).length : ℚ) / 6) ^ 2
      + (((words4.filter (fun w => layer w = 2)).length : ℚ) / 6) ^ 2 = 5 / 9 := by
  have h1 : (words4.filter (fun w => layer w = 1)).length = 4 := by decide
  have h2 : (words4.filter (fun w => layer w = 2)).length = 2 := by decide
  rw [h1, h2]
  norm_num

#print axioms z077_layer_partition_q

end Zero
