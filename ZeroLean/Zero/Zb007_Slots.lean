/-
  Zb007 §2–§3 的机器检查版本：内禀频率与三频率槽。

  正文的数学核：
    · 定理 2.1：2 级块 `ρ = diag(p,q)` 的非对角元频率是 `log(p/q)`；
      这里把它写成 `log(p/q) = log p − log q`（频率的代数形式）；
    · 命题 3.2：seed `[8,4,3/2,1]` 的相邻比是 `2, 8/3, 3/2`，
      三个数互不相同，故 `log` 给出三个不同的频率槽。

  未覆盖：模流 `σ_t(A)=ρ^{it}Aρ^{-it}` 的算子论版本（需要 `ρ^{it}` 的
  函数演算），以及圈图谱 `2cos(2πm/n)`（需要循环矩阵的特征值）。
-/

import Mathlib

namespace Zero

/-- **Zb007 · 2 级块的内禀频率**：`log(p/q) = log p − log q`。 -/
theorem zb007_two_level_frequency (p q : ℝ) (hp : 0 < p) (hq : 0 < q) :
    Real.log (p / q) = Real.log p - Real.log q := by
  rw [Real.log_div (ne_of_gt hp) (ne_of_gt hq)]

/-- **Zb007 · 比值 2 给频率 log 2**：`p = 2q` 时 `log(p/q) = log 2`。 -/
theorem zb007_frequency_of_ratio_two (q : ℝ) (hq : 0 < q) :
    Real.log ((2 * q) / q) = Real.log 2 := by
  field_simp [ne_of_gt hq]

/-- **Zb007 · seed 的相邻比**：`[8,4,3/2,1]` 给 `2, 8/3, 3/2`。 -/
theorem zb007_seed_ratios :
    (8 / 4 : ℝ) = 2 ∧ (4 / (3 / 2) : ℝ) = 8 / 3 ∧ ((3 / 2) / 1 : ℝ) = 3 / 2 := by
  norm_num

/-- 三个相邻比互不相同。 -/
theorem zb007_ratios_distinct :
    (2 : ℝ) ≠ 8 / 3 ∧ (8 / 3 : ℝ) ≠ 3 / 2 ∧ (2 : ℝ) ≠ 3 / 2 := by
  norm_num

/-- `log` 在正实数上单射（频率槽 = 比值槽的根据）。 -/
theorem zb007_log_ne {a b : ℝ} (ha : 0 < a) (hb : 0 < b) (h : a ≠ b) :
    Real.log a ≠ Real.log b := by
  intro hlog
  apply h
  have h' := congrArg Real.exp hlog
  rwa [Real.exp_log ha, Real.exp_log hb] at h'

/-- **Zb007 · 三个频率槽互不相同**：`log 2, log(8/3), log(3/2)` 两两不等。 -/
theorem zb007_logs_distinct :
    Real.log 2 ≠ Real.log (8 / 3) ∧
      Real.log (8 / 3) ≠ Real.log (3 / 2) ∧
      Real.log 2 ≠ Real.log (3 / 2) := by
  refine ⟨?_, ?_, ?_⟩
  · exact zb007_log_ne (by norm_num) (by norm_num) (by norm_num)
  · exact zb007_log_ne (by norm_num) (by norm_num) (by norm_num)
  · exact zb007_log_ne (by norm_num) (by norm_num) (by norm_num)


/-- **Zb007 · 模流相位**：`e^{it log p} · e^{−it log q} = e^{it log(p/q)}`。 -/
theorem zb007_modular_phase (p q t : ℝ) (hp : 0 < p) (hq : 0 < q) :
    Complex.exp (Complex.I * t * Real.log p) *
      Complex.exp (-(Complex.I * t * Real.log q)) =
      Complex.exp (Complex.I * t * Real.log (p / q)) := by
  rw [← Complex.exp_add, Real.log_div (ne_of_gt hp) (ne_of_gt hq)]
  congr 1
  push_cast
  ring

#print axioms zb007_two_level_frequency
#print axioms zb007_frequency_of_ratio_two
#print axioms zb007_seed_ratios
#print axioms zb007_ratios_distinct
#print axioms zb007_logs_distinct
#print axioms zb007_modular_phase

end Zero
