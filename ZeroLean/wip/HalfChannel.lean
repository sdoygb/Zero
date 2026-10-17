/-
  D0（半槽规则）的算术核心。

  通道计数：$n$ 个满通道 + 一个半通道（权重 $1/2$），总权重 $n+1/2$；
  耦合只计满通道，故矩阵元 = 满通道占比 $n/(n+1/2)$。

  实例：$n=784$（质量单元的两体通道）与 $n=16$（$M_R$ 的单态通道）。
-/

import Mathlib

namespace Zero

/-- **半槽占比**：满通道占比 $n/(n+\tfrac12)=\dfrac{2n}{2n+1}$。 -/
theorem half_channel_fraction (n : ℕ) (hn : 0 < n) :
    ((n : ℝ) / (n + 1 / 2)) = (2 * n) / (2 * n + 1) := by
  have h : (n : ℝ) + 1 / 2 ≠ 0 := by positivity
  field_simp

/-- **质量单元实例**：$784/(784+\tfrac12)=1568/1569$。 -/
theorem half_channel_784 : ((784 : ℝ) / (784 + 1 / 2)) = 1568 / 1569 := by
  norm_num

/-- **Majorana 实例**：$16/(16+\tfrac12)=32/33$。 -/
theorem half_channel_16 : ((16 : ℝ) / (16 + 1 / 2)) = 32 / 33 := by
  norm_num

#print axioms half_channel_fraction
#print axioms half_channel_784
#print axioms half_channel_16

end Zero
