/-
  Z114（配对对称性与介导竞争）的**线性代数核心**：配对势的对称性 ⟹ 能隙的符号结构。

  Z114 §二 给出两种配对势的核：
    · 图上磁介导（二部图 Néel 序，σ_k = ±1）：`V^{mag}_{kk'} = −J²·σ_k·σ_{k'}`  ← `z114NeelKernel`
    · 常数势（BCS 情形）：`V_{kk'} = −V₀`（全 1 核）                            ← `z114ConstKernel`

  Z114 §二 的判据（原句：「对称性判定：s-wave=能隙全同号；d-wave=符号交替（节点过零）」）
  在这里被逐条翻成关于核的纯线性代数命题：

    · `z114_neel_eigen`  ：`K·σ = −J²·|ι|·σ`，即 Néel 核把 σ 映成自己的倍数（本征向量 ∝ σ）；
    · `z114_neel_image`  ：`K·x = (−J²·Σ_j σ_j x_j)·σ`，像恒落在 span{σ} 里（秩 ≤ 1），
      于是任何本征向量都 ∝ σ，σ 的符号结构就是能隙的符号结构；
    · `z114_gap_alternates`：`Δ = c·σ`（c ≠ 0）在 Néel 边（σ 反号的两点）上符号交替 —— 这是
      Z114 §三 那条「磁介导 ⟹ d-wave，节点在跨 Néel 边」的代数骨架（符号交替 = d-wave）；
    · `z114_const_image` / `z114_const_gap_constant_of_ne_zero`：常数势的像恒落在 span{1} 里
      （`J·Δ = (ΣΔ)·1`），故非零本征值（μ ≠ 0）的本征向量必为常向量 ⟹ 全同号 —— 这是
      Z114 §三 表里「常数势（BCS）⟹ 全同号 ⟹ 纯 s-wave」那一行的代数内容（s-wave）。

  ⚠️ 任务清单里最后一条 `z114_const_gap_constant`（假设 `μ ≠ |ι|`）**按字面为假**，故未收录。
  成因：`μ = 0` 是 J 的本征值，其本征空间是 `{Δ | ΣΔ = 0}`，远大于常向量；`ΣΔ = 0` 并不蕴含
  `Δ = 0`（只有「Δ 是常向量」时才成立）。机器核验的反例见
  `z114_const_gap_constant_counterexample`（`ι = Fin 2`、`Δ = (1, −1)`、`μ = 0`：`J·Δ = (ΣΔ)·1 = 0 = 0·Δ`，
  且 `0 ≠ |ι| = 2`，但 `Δ 0 ≠ Δ 1`）。把该假设换成 `μ ≠ 0` 后陈述为真
  （`z114_const_gap_constant_of_ne_zero`）——s-wave 判据的正确代数内容。

  机器检查：`lake env lean Zero/Z114_NeelKernel.lean`（本文件只依赖 Mathlib；无 sorry/axiom 引入）。
-/

import Mathlib

namespace Zero

open scoped Matrix

/-! ## 一、Néel 核（磁介导势）与 d-wave 符号结构 -/

/-- Néel 核：`K_ij = −J² · σ_i · σ_j`（`σ` 在每个点取 ±1）
    —— 即 Z114 §二 的 `V^{mag}_{kk'} = −J²·σ_k·σ_{k'}`。 -/
noncomputable def z114NeelKernel {ι : Type*} [Fintype ι] (σ : ι → ℝ) (J : ℝ) :
    Matrix ι ι ℝ := fun i j => -(J ^ 2) * (σ i * σ j)

/-- **Néel 核的本征向量 ∝ σ**：`K · σ = −J²·|ι|·σ`。
    逐分量算：`(K·σ)_i = Σ_j (−J²σ_iσ_j)σ_j = −J²σ_i·Σ_j σ_j² = −J²·|ι|·σ_i`（用 `σ_j² = 1`）。 -/
theorem z114_neel_eigen {ι : Type*} [Fintype ι] (σ : ι → ℝ) (J : ℝ)
    (hσ : ∀ i, σ i * σ i = 1) :
    z114NeelKernel σ J *ᵥ σ = (-(J ^ 2) * (Fintype.card ι : ℝ)) • σ := by
  ext i
  simp only [Matrix.mulVec, dotProduct, z114NeelKernel, Pi.smul_apply, smul_eq_mul]
  have hcongr : ∑ j, (-(J ^ 2) * (σ i * σ j)) * σ j = ∑ _j : ι, (-(J ^ 2) * σ i) :=
    Finset.sum_congr rfl (fun j _ => by
      calc -(J ^ 2) * (σ i * σ j) * σ j
          = -(J ^ 2) * σ i * (σ j * σ j) := by ring
        _ = -(J ^ 2) * σ i * 1 := by rw [hσ j]
        _ = -(J ^ 2) * σ i := by ring)
  rw [hcongr, Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
  ring

/-- **像恒落在 span{σ} 里**（rank-one 结构）：`K · x = (−J²·Σ_j σ_j x_j)·σ`。
    这条把「K 的秩 ≤ 1 且值域 = span{σ}」写成了等式：无论 x 是什么，K·x 都是 σ 的倍数，
    故 K 的每个本征向量都 ∝ σ（Néel 核不产生 σ 之外的符号结构）。 -/
theorem z114_neel_image {ι : Type*} [Fintype ι] (σ : ι → ℝ) (J : ℝ) (x : ι → ℝ) :
    z114NeelKernel σ J *ᵥ x = (-(J ^ 2) * ∑ j, σ j * x j) • σ := by
  ext i
  simp only [Matrix.mulVec, dotProduct, z114NeelKernel, Pi.smul_apply, smul_eq_mul]
  have hcongr : ∑ j, (-(J ^ 2) * (σ i * σ j)) * x j = ∑ j, σ i * (-(J ^ 2) * (σ j * x j)) :=
    Finset.sum_congr rfl (fun j _ => by ring)
  rw [hcongr, ← Finset.mul_sum, ← Finset.mul_sum]
  ring

/-- **d-wave 符号交替**：若 `Δ = c·σ`（`c ≠ 0`）且两点 `σ i = 1`、`σ j = −1`，
    则 `Δ i` 与 `Δ j` 异号 —— 即 Z114 §二 判据里的「d-wave = 符号交替（节点过零）」：
    Néel 序反号的两点之间必过零（节点落在跨 Néel 边上）。 -/
theorem z114_gap_alternates {ι : Type*} (σ : ι → ℝ) (c : ℝ) (hc : c ≠ 0) (i j : ι)
    (hi : σ i = 1) (hj : σ j = -1) :
    (c • σ) i * (c • σ) j < 0 := by
  have hi' : (c • σ) i = c := by simp [hi]
  have hj' : (c • σ) j = -c := by simp [hj]
  rw [hi', hj']
  nlinarith [mul_self_pos.mpr hc]

/-! ## 二、常数势（BCS / s-wave 情形）的核 -/

/-- 全 1 核（常数势）：`J_ij = 1` —— Z114 §二 的常数势 `V_{kk'} = −V₀`（取 V₀ = −1 的归一化）。 -/
noncomputable def z114ConstKernel {ι : Type*} [Fintype ι] : Matrix ι ι ℝ := fun _ _ => 1

/-- **常数势的像恒落在 span{1} 里**：`J · Δ = (Σ_j Δ_j) · 1`。
    这是 s-wave 判据的无条件代数内容：常数势不可能生成非平凡的非 s-wave 结构 ——
    任何本征向量只要在 J 的像里（即 μ ≠ 0），就只能是常向量。 -/
theorem z114_const_image {ι : Type*} [Fintype ι] (Δ : ι → ℝ) :
    z114ConstKernel (ι := ι) *ᵥ Δ = (∑ j, Δ j) • (fun _ : ι => (1 : ℝ)) := by
  ext i
  simp only [z114ConstKernel, Matrix.mulVec, dotProduct, one_mul, Pi.smul_apply, smul_eq_mul,
    mul_one]

/-- **s-wave：常数势 ⟹ 能隙的分量全相等**（故全同号）。
    由 `z114_const_image`：`(ΣΔ)·1 = μ·Δ` 逐分量给出 `Σ_j Δ_j = μ·Δ_i`；当 `μ ≠ 0` 时
    左端与 i 无关，故 `μ·Δ_i = μ·Δ_j`，约去 μ ≠ 0 得 `Δ_i = Δ_j`。
    对照 Z114 §三 表：「常数势（BCS）⟹ 能隙全同号（+15.9，负顶点 0，节点 0）⟹ 纯 s-wave」。 -/
theorem z114_const_gap_constant_of_ne_zero {ι : Type*} [Fintype ι]
    (Δ : ι → ℝ) (μ : ℝ)
    (h : z114ConstKernel (ι := ι) *ᵥ Δ = μ • Δ) (hμ : μ ≠ 0) :
    ∀ i j : ι, Δ i = Δ j := by
  have hS : ∀ i, ∑ j, Δ j = μ * Δ i := by
    intro i
    have := congrFun h i
    simpa [z114ConstKernel, Matrix.mulVec, dotProduct, Pi.smul_apply, smul_eq_mul] using this
  intro i j
  have hmul : μ * Δ i = μ * Δ j := by rw [← hS i, ← hS j]
  exact mul_left_cancel₀ hμ hmul

/-- **任务清单那条 `z114_const_gap_constant` 的反例**（故该陈述按字面为假、未收录；
    此处只记录反例，不改弱任何被要求的陈述）。
    取 `ι = Fin 2`、`Δ = (1, −1)`、`μ = 0`：`J·Δ = (ΣΔ)·1 = 0 = 0·Δ` 且 `μ = 0 ≠ |ι| = 2`，
    但 `Δ 0 ≠ Δ 1`。即：`μ = 0` 的本征空间是 `{Δ | ΣΔ = 0}`，远比常向量集大。 -/
theorem z114_const_gap_constant_counterexample :
    ∃ (Δ : Fin 2 → ℝ) (μ : ℝ),
      z114ConstKernel (ι := Fin 2) *ᵥ Δ = μ • Δ ∧ μ ≠ (Fintype.card (Fin 2) : ℝ) ∧
        ¬ (∀ i j : Fin 2, Δ i = Δ j) := by
  refine ⟨![1, -1], 0, ?_, ?_, ?_⟩
  · ext i
    fin_cases i <;>
      simp [z114ConstKernel, Matrix.mulVec, dotProduct]
  · norm_num
  · intro hc
    have := hc 0 1
    norm_num at this

#print axioms z114_neel_eigen
#print axioms z114_neel_image
#print axioms z114_gap_alternates
#print axioms z114_const_image
#print axioms z114_const_gap_constant_of_ne_zero
#print axioms z114_const_gap_constant_counterexample

end Zero
