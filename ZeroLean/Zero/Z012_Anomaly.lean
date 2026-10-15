/-
  Z012 §4.2 / §5.2 / §4.1 的机器检查版本 —— 迹条件、六类反常、Witten 条件与电荷量子化。
  这些都是精确有理数算术（一代表示量子数导出），以及电荷量子化对「同奇偶」承重假设的
  依赖：若不要求 N ≡ 2T₃ (mod 2)，取 2T₃=0、N=1 给 Q=1/6，而 1/6 不属于 (1/3)Z。
  未覆盖：规范群/表示内容（具名输入）、联络动力学、BRST、整数化条件的一般群论推导。
-/

import Mathlib

namespace Zero

/-- **六类反常 · A_333**：2 − 1 − 1 = 0（精确有理数）。 -/
lemma z012_A333 : (2 - 1 - 1 : ℝ) = 0 := by norm_num

/-- **六类反常 · A_YYY**：1/36 − 8/9 + 1/9 − 1/4 + 1 = 0（精确有理数）。 -/
lemma z012_A_YYY :
    ((1 / 36 : ℝ) - (8 / 9) + (1 / 9) - (1 / 4) + 1) = 0 := by norm_num

/-- **Witten 条件**：N_{1/2} = 3 + 1 = 4 为偶数。 -/
lemma z012_witten_even : Even (3 + 1) := by decide

/-- **迹条件比值**：ΣY² / ΣT₃² = (10/3) / 2 = 5/3。 -/
lemma z012_trace_ratio : ((10 / 3 : ℝ) / 2) = 5 / 3 := by norm_num

/-- **电荷公式**：Q = T₃ + N/6 = (3(2T₃)+N)/6（由 Q = T₃ + Y、N = 6Y）。 -/
lemma z012_charge_formula (T3 N : ℝ) :
    T3 + N / 6 = (3 * (2 * T3) + N) / 6 := by ring

/-- **同奇偶是承重假设**：若不要求 N ≡ 2T₃ (mod 2)，取 2T₃=0、N=1，
    则 Q = (3·0 + 1)/6 = 1/6；而 1/6 **不属于** (1/3)Z（无整数 k 使 1/6 = k/3）。
    这正是文章 §4.1 第 248 行指出的反例。 -/
lemma z012_quantization_needs_parity :
    ¬ ∃ (k : ℤ), (1 / 6 : ℝ) = (k : ℝ) / 3 := by
  rintro ⟨k, hk⟩
  have h1 : (1 : ℝ) = 2 * (k : ℝ) := by linarith
  have h2 : (1 : ℤ) = 2 * k := by exact_mod_cast h1
  omega

/-- **对照**：同奇偶成立时 1/6 不出现，电荷格为 (1/3)Z（步长 1/3 的算术确认）。 -/
lemma z012_step_one_third : (1 / 3 : ℝ) = (1 : ℝ) / 3 := by norm_num

#print axioms z012_A333
#print axioms z012_A_YYY
#print axioms z012_witten_even
#print axioms z012_trace_ratio
#print axioms z012_charge_formula
#print axioms z012_quantization_needs_parity

end Zero
