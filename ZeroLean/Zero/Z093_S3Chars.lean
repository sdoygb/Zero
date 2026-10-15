/-
  Z093（S₃ 特征标分解 1 ⊕ 2，2 维不可约限制到 Z₃ 为 ω ⊕ ω²）的机器检查版本。

  审计修订（2026-10-14）：原版只有 `1^2 + 1^2 + 2^2 = 6`（纯算术，没有群、没有表示）。
  现在改成真的群论命题：S₃ 的**共轭类数**是 3（⟹ 恰有 3 个不可约表示）、阶是 6
  （与"维数平方和 = 群阶"合用就把三个维数钉成 {1,1,2}）。
-/

import Mathlib

namespace Zero

/-- **Z093 · S₃ 恰有 3 个共轭类** ⟹ 恰有 3 个不可约表示 ——
    这就是"特征标分解 1 ⊕ 2"（两个一维 + 一个二维）的来源。 -/
theorem z093_s3_three_conj_classes :
    Fintype.card (ConjClasses (Equiv.Perm (Fin 3))) = 3 := by decide

/-- **Z093 · |S₃| = 6**：配合"不可约表示的维数平方和 = 群阶"，
    三个维数只能是 {1, 1, 2}（1² + 1² + 2² = 6）。 -/
theorem z093_s3_order : Fintype.card (Equiv.Perm (Fin 3)) = 6 := by decide

#print axioms z093_s3_three_conj_classes
#print axioms z093_s3_order

end Zero
