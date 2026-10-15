/-
  Z131（原生演化律）的机器检查版本 —— 模流 σ_k(A) = ρ^k A ρ^{-k} 的三条承重性质。

  背景：Z093 标下的缺口是"Zero 无 Z005 时间演化律"。Z131 的判定是：**原生演化律就是模流**，
  它在中心态上平凡（Z005 定理 3），在自发非中心态上非平凡、以内禀频率 log(p/q) 振荡；
  K₄ 的三频率槽 [0.4055, 0.6931, 0.9808] 正是稳态记录相邻权重比的 log。

  这里把「能量标度／时间演化」这条最承重的线里**可形式化的部分**钉死：
    · 模流对**对角元**（= 与 ρ 交换的中心部分）逐点不动 —— 即"中心态 ⇒ 平凡"；
    · 对**非对角元**，模流按与步数无关的固定比值 ρ i / ρ j 逐次缩放 —— 这个比值就是内禀频率
      （它的 log = log(ρ i) − log(ρ j)，即 Z131 的 log(p/q)）；
    · 模流满足**复合律** σ_l ∘ σ_k = σ_{k+l}（"演化律"该有的群作用性质）；
    · K₄ 稳态记录 [8,4,1.5,1] 的相邻权重比恰为 2、8/3、3/2 —— 三槽是其 log。

  注意（诚实边界）："取哪张图、seed 首项为何是 8、为何用这个模流"都是**具名输入**，
  本文件不证也证不了；能被 kernel 检查的是"给定这些输入后，推论对不对"。
-/

import Mathlib

namespace Zero

/-- **模流（离散步）**：对角密度 `ρ`（有限指标集 `Fin n` → ℝ）上的模流把矩阵元 `A i j`
    乘以 `(ρ i / ρ j) ^ k`。连续参数 t 的版本是 `ρ^{it} A ρ^{-it}`，这里取整数步 `k : ℤ`
    以避免复指数幂的额外记号；频率结论与步长无关（见 `z131_offdiagonal_ratio`）。 -/
noncomputable def modularFlow {n : ℕ} (ρ : Fin n → ℝ) (A : Matrix (Fin n) (Fin n) ℝ) (k : ℤ) :
    Matrix (Fin n) (Fin n) ℝ :=
  fun i j => (ρ i / ρ j) ^ k * A i j

/-- **Z131 · 中心态 ⇒ 模流平凡**：对角元（与 ρ 交换的那部分）在模流下逐点不动。
    这是 Z005 定理 3 的有限维离散版本：中心态上 σ = id。 -/
theorem z131_flow_fixes_diagonal {n : ℕ} (ρ : Fin n → ℝ) (A : Matrix (Fin n) (Fin n) ℝ)
    (k : ℤ) (i : Fin n) (h : ρ i ≠ 0) :
    modularFlow ρ A k i i = A i i := by
  simp [modularFlow, div_self h]

/-- **Z131 · 非中心态 ⇒ 内禀频率**：非对角元的模流／原值之比是 `(ρ i / ρ j) ^ k` ——
    一个**与 k 无关的固定比值** `ρ i / ρ j` 的 k 次幂。这就是"内禀频率"的精确含义：
    振荡因子由权重比决定，不需要任何外部时间单位。 -/
theorem z131_offdiagonal_ratio {n : ℕ} (ρ : Fin n → ℝ) (A : Matrix (Fin n) (Fin n) ℝ)
    (k : ℤ) {i j : Fin n} (hA : A i j ≠ 0) :
    modularFlow ρ A k i j / A i j = (ρ i / ρ j) ^ k := by
  simp only [modularFlow]
  rw [mul_div_assoc, div_self hA, mul_one]

/-- **Z131 · 内禀频率的对数形式**：`log(ρ i / ρ j) = log(ρ i) − log(ρ j)`。
    正权重下，"频率 = log 权重比"（即 Z131 的 log(p/q)）。 -/
theorem z131_frequency_is_log_ratio {n : ℕ} (ρ : Fin n → ℝ) (i j : Fin n)
    (hi : 0 < ρ i) (hj : 0 < ρ j) :
    Real.log (ρ i / ρ j) = Real.log (ρ i) - Real.log (ρ j) :=
  Real.log_div hi.ne' hj.ne'

/-- **Z131 · 模流的复合律**：σ_l ∘ σ_k = σ_{k+l} —— 模流是 ℤ 在可观测量上的作用，
    这正是"演化律"该满足的性质（时间可加）。 -/
theorem z131_flow_composition {n : ℕ} (ρ : Fin n → ℝ) (hρ : ∀ i, ρ i ≠ 0)
    (A : Matrix (Fin n) (Fin n) ℝ) (k l : ℤ) :
    modularFlow ρ (modularFlow ρ A k) l = modularFlow ρ A (k + l) := by
  ext i j
  simp only [modularFlow]
  rw [zpow_add₀ (div_ne_zero (hρ i) (hρ j)) k l]
  ring

/-- **Z131 · K₄ 稳态记录的相邻权重比**：记录 = seed/Σseed，seed = `[8, 4, 1.5, 1]`，
    相邻比恰为 `2`、`8/3`、`3/2`（与归一化常数无关 —— 这是"槽位是比值"的关键）。 -/
theorem z131_k4_slot_ratios :
    (8 : ℝ) / 4 = 2 ∧ (4 : ℝ) / 1.5 = 8 / 3 ∧ (1.5 : ℝ) / 1 = 3 / 2 := by
  norm_num

/-- **Z131 · 三频率槽 = 相邻权重比的 log**：`[log(3/2), log 2, log(8/3)]`
    = `[0.4055, 0.6931, 0.9808]`（Z131 §摘要的数值对证的形式版本）。 -/
theorem z131_k4_slots_are_logs :
    Real.log ((8 : ℝ) / 4) = Real.log 2 ∧
      Real.log ((4 : ℝ) / 1.5) = Real.log ((8 : ℝ) / 3) ∧
      Real.log ((1.5 : ℝ) / 1) = Real.log ((3 : ℝ) / 2) := by
  refine ⟨?_, ?_, ?_⟩ <;> norm_num

#print axioms z131_flow_fixes_diagonal
#print axioms z131_offdiagonal_ratio
#print axioms z131_frequency_is_log_ratio
#print axioms z131_flow_composition
#print axioms z131_k4_slot_ratios
#print axioms z131_k4_slots_are_logs

end Zero
