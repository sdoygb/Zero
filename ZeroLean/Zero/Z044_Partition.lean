/-
  Z044（能谱之桥的判定）的机器检查版本 —— 配分函数与"模能均值 = Shannon 熵"。

  Z044 的 F 类两条：
    · `Σ_h e^{−K_h} = Σ_h ω_h = 1`（配分函数归一，K_h = −log ω_h 的直接后果）；
    · `⟨K⟩_ω = H(w)`（模能均值恰是层权 Shannon 熵）。

  两条都写成**一般权重**（不只是 L=4 的 (2/3,1/3)），因为它们的"内容"正是
  "这些等式是定义展开，不依赖具体数值" —— 这比逐个代数值更有信息量。
-/

import Mathlib

namespace Zero

/-- **Z044 · 配分函数归一**：`Σ_h e^{−K_h} = Σ_h ω_h`，而后者为 1。
    即 `K_h = −log ω_h` 时配分函数自动归一（不需要额外假设）。 -/
theorem z044_partition_normalized (ω : Fin 2 → ℝ) (hpos : ∀ i, 0 < ω i) :
    ∑ i, Real.exp (-(-Real.log (ω i))) = ∑ i, ω i := by
  exact Finset.sum_congr rfl fun i _ => by
    rw [neg_neg, Real.exp_log (hpos i)]

/-- **Z044 · 配分函数 = 1**（在权重归一的前提下）。 -/
theorem z044_partition_eq_one (ω : Fin 2 → ℝ) (hpos : ∀ i, 0 < ω i)
    (hsum : ∑ i, ω i = 1) :
    ∑ i, Real.exp (-(-Real.log (ω i))) = 1 := by
  rw [z044_partition_normalized ω hpos, hsum]

/-- **Z044 · 模能均值 = Shannon 熵**：`Σ_h ω_h K_h = −Σ_h ω_h log ω_h`（`K_h = −log ω_h`）。
    —— 这条说明"模能均值等于熵"是定义展开的后果，不是数值巧合。 -/
theorem z044_entropy_mean (ω : Fin 2 → ℝ) :
    ∑ i, ω i * (-Real.log (ω i)) = -∑ i, ω i * Real.log (ω i) := by
  rw [← Finset.sum_neg_distrib]
  exact Finset.sum_congr rfl fun i _ => by ring

#print axioms z044_partition_normalized
#print axioms z044_partition_eq_one
#print axioms z044_entropy_mean

end Zero
