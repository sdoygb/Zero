/-
  Z086（Σs_i = 3⟨s⟩ ⟹ Σw_i = 0 ⟹ 偏离空间仅 2 连续自由度）的机器检查版本。

  审计修订（2026-10-14）：原版第二条是 `(3 : ℕ) - 1 = 2`（纯算术，没有"超平面"这个对象）。
  现在改成真的线性代数：零和向量构成 `Fin 3 → ℝ` 的一个子空间，其**维数**是 2
  （秩–零化度定理：`finrank (range) + finrank (ker) = finrank M`，而 range = ⊤ 的
  finrank 是 1，M 的是 3）。
-/

import Mathlib

namespace Zero

open Finset

/-- **Z086 · 均值偏离和为零**：若 s₁+s₂+s₃ = 3μ，则 (s₁−μ)+(s₂−μ)+(s₃−μ) = 0。 -/
theorem z086_mean_deviation_sum_zero (s1 s2 s3 μ : ℝ)
    (h : s1 + s2 + s3 = 3 * μ) :
    (s1 - μ) + (s2 - μ) + (s3 - μ) = 0 := by linarith

/-- 求和线性映射 `w ↦ Σ w i`。 -/
def sumLM3 : (Fin 3 → ℝ) →ₗ[ℝ] ℝ where
  toFun w := ∑ i, w i
  map_add' x y := by simp [Finset.sum_add_distrib]
  map_smul' c x := by simp [Finset.mul_sum]

/-- **Z086 · 偏离空间是 2 维超平面**：零和子空间 `{w : Fin 3 → ℝ // Σ w = 0}`
    作为 ℝ-向量空间的维数恰为 `2`（三个分量减一个线性约束）。
    这是"三代偏离只有两个连续自由度"的机器检查版本。 -/
theorem z086_hyperplane_finrank :
    Module.finrank ℝ (LinearMap.ker sumLM3) = 2 := by
  have h := LinearMap.finrank_range_add_finrank_ker sumLM3
  have htop : LinearMap.range sumLM3 = ⊤ := by
    rw [LinearMap.range_eq_top]
    intro y
    refine ⟨fun _ => y / 3, ?_⟩
    simp only [sumLM3, LinearMap.coe_mk, AddHom.coe_mk, Fin.sum_univ_three]
    ring
  rw [htop] at h
  have htopfin : Module.finrank ℝ (⊤ : Submodule ℝ ℝ) = 1 := by simp
  have hM : Module.finrank ℝ (Fin 3 → ℝ) = 3 := by simp
  rw [htopfin, hM] at h
  omega

#print axioms z086_mean_deviation_sum_zero
#print axioms z086_hyperplane_finrank

end Zero
