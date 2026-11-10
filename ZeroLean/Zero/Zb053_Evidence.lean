/-
  Zb053 的机器检查版本：证据的计数判据。

  正文 §1–§2 的数学核：
    · 期望命中 `E[N] = M p`（有限和，只需线性性）；
    · 门槛 `ε < 1/(2M) ⟹ M ε < 1/2`（Markov 界的具体算术形式）；
    · 两个具体门槛：`M=10^4, ε=10^{-6} ⟹ Mε=10^{-2}`。

  未覆盖：把「候选值在目标附近无信息地分布」编码成概率模型（那是建模选择，
  正文已记为记账约定），以及 Markov 不等式本身的测度论版本（Mathlib 已有，
  但本文件只取有限/算术核）。
-/

import Mathlib

namespace Zero

/-- **Zb053 · 有限期望的线性性**：`M` 个候选、每个命中概率 `p`，
    则期望命中数 `Σ_{i<M} p = M p`。 -/
theorem zb053_expectation (M : ℕ) (p : ℝ) :
    (∑ _ : Fin M, p) = (M : ℝ) * p := by
  simp

/-- **Zb053 · 门槛的算术形式**：`M > 0` 且 `ε < 1/(2M)` ⟹ `M ε < 1/2`。 -/
theorem zb053_markov_threshold (M : ℕ) (ε : ℝ) (hM : 0 < M)
    (hε : ε < 1 / (2 * (M : ℝ))) :
    (M : ℝ) * ε < 1 / 2 := by
  have hM' : (0 : ℝ) < (M : ℝ) := by exact_mod_cast hM
  have hmul : (M : ℝ) * ε < (M : ℝ) * (1 / (2 * (M : ℝ))) :=
    mul_lt_mul_of_pos_left hε hM'
  have hhalf : (M : ℝ) * (1 / (2 * (M : ℝ))) = 1 / 2 := by
    field_simp
  linarith

/-- **Zb053 · 具体门槛一**：`M = 10^4`、`ε = 10^{-6}` 给 `M ε = 10^{-2}`。 -/
theorem zb053_example_M1e4 : (10000 : ℝ) * (1 / 1000000) = 1 / 100 := by
  norm_num

/-- **Zb053 · 具体门槛二**：`M = 10^{10}`、`ε = 10^{-12}` 给 `M ε = 10^{-2}`。 -/
theorem zb053_example_M1e10 : ((10 : ℝ) ^ 10) * (1 / (10 : ℝ) ^ 12) = 1 / 100 := by
  norm_num

/-- **Zb053 · 精度不构成证据**：同一个 `ε = 10^{-5}` 在 `M = 10^2` 与 `M = 10^5`
    下给出 `Mε = 10^{-3}` 与 `Mε = 1`——前者是强证据，后者与偶然无法区分。 -/
theorem zb053_precision_relative :
    (100 : ℝ) * (1 / 100000) = 1 / 1000 ∧
      (100000 : ℝ) * (1 / 100000) = 1 := by
  norm_num

#print axioms zb053_expectation
#print axioms zb053_markov_threshold
#print axioms zb053_example_M1e4
#print axioms zb053_example_M1e10
#print axioms zb053_precision_relative

end Zero
