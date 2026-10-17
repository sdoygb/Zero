/-
  Z117（Zero 量子双缝：图上的量子行走）的机器检查版本。

  Z117 的推导骨架只有三步，这里逐步机器核验：

    ① **反厄米性**：`H` 厄米 ⟹ `A := −i t H` 反厄米（`Aᴴ = −A`）。
       对应的定理是 `z117_skew_adjoint`。图上量子行走的步进算子写成 `U = e^{−iHt}`，
       而 `−iHt` 反厄米正是 `U` 酉的**唯一来源**（不是额外假设）。
    ② **酉性**：反厄米 `A` 的指数 `e^A` 满足 `Uᴴ U = 1` 且 `U Uᴴ = 1`。
       对应的定理是 `z117_evolution_unitary`；算术内核是
       `Uᴴ U = e^{Aᴴ} e^{A} = e^{Aᴴ + A} = e^{0} = 1`（`Aᴴ = −A` 使 `Aᴴ` 与 `A` 交换）。
    ③ **概率守恒**：`U` 酉 ⟹ 模方和（概率和）不变，Born 规则的算术内容。
       对应的定理是 `z117_normSq_preserved` 与 `z117_probability_sum_one`。

  实现说明：`exp` 用 `NormedSpace.exp`。矩阵上的 `exp_add_of_commute` 用
  `Matrix.exp_add_of_commute`（该引理把矩阵范数的选择藏在证明内部，只要求标量 `ℂ`
  是 `NormedRing`），矩阵共轭转置与指数交换用 `Matrix.exp_conjTranspose`。
  这样无需 `open scoped Matrix.Norms.L2Operator` 引入任何矩阵范数约定。
-/

import Mathlib

namespace Zero

open Matrix

/-- 概率/模方（自给自足，不用希尔伯特空间范数） -/
noncomputable def z116NormSq {n : ℕ} (ψ : Fin n → ℂ) : ℝ := ∑ i, Complex.normSq (ψ i)

/-- **Z117 ① 的算术内核**：反厄米 `A`（`Aᴴ = −A`）的指数是酉的。

    证明就是 Z117 §演化的那三行：
    `Uᴴ U = e^{Aᴴ} e^{A} = e^{Aᴴ + A} = e^{0} = 1`（`Aᴴ + A = 0` 用到反厄米），
    右乘版同理。两条都证，因为"酉"要求双侧逆。 -/
theorem z117_exp_conjTranspose_mul_self {n : ℕ} (A : Matrix (Fin n) (Fin n) ℂ)
    (hA : A.conjTranspose = -A) :
    (NormedSpace.exp A).conjTranspose * NormedSpace.exp A = 1 ∧
      NormedSpace.exp A * (NormedSpace.exp A).conjTranspose = 1 := by
  -- `Aᴴ = −A` 使得 `Aᴴ` 与 `A` 交换（左/右各一次），这是用指数加法性的前提
  have h₁ : Commute A.conjTranspose A := by rw [hA]; exact (Commute.refl A).neg_left
  have h₂ : Commute A A.conjTranspose := by rw [hA]; exact (Commute.refl A).neg_right
  constructor
  · -- `Uᴴ U = e^{Aᴴ} · e^{A} = e^{Aᴴ + A} = e^{0} = 1`
    -- `Matrix.exp_conjTranspose` 把目标里的 `(e^A)ᴴ` 换成 `e^{Aᴴ}`，再接 `Aᴴ = −A`
    have h := Matrix.exp_add_of_commute A.conjTranspose A h₁
    rw [hA, neg_add_cancel, NormedSpace.exp_zero] at h
    rw [← Matrix.exp_conjTranspose A, hA]
    exact h.symm
  · -- `U Uᴴ = e^{A} · e^{Aᴴ} = e^{A + Aᴴ} = e^{0} = 1`
    have h := Matrix.exp_add_of_commute A A.conjTranspose h₂
    rw [hA, add_neg_cancel, NormedSpace.exp_zero] at h
    rw [← Matrix.exp_conjTranspose A, hA]
    exact h.symm

/-- **辅助引理（模方和的点积写法）**：模方和是"共轭点积"的实部，
    即 `∑ i, normSq (ψ i) = (ψ̄ ⬝ᵥ ψ).re`。这一步只是把 Born 规则的
    `|ψ_i|²` 翻译成矩阵语言，便于用 `Uᴴ U = 1` 传递。 -/
theorem z117_normSq_eq_re_dotProduct {n : ℕ} (ψ : Fin n → ℂ) :
    (star ψ ⬝ᵥ ψ).re = z116NormSq ψ := by
  have h : star ψ ⬝ᵥ ψ = ((z116NormSq ψ : ℝ) : ℂ) := by
    simp only [z116NormSq, dotProduct, Pi.star_apply, Complex.star_def,
      Complex.normSq_eq_conj_mul_self, Complex.ofReal_sum]
  rw [h, Complex.ofReal_re]

/-- **Z117 ② → ③ 的传递步骤**：`Uᴴ U = 1`（左逆即酉）⟹ 模方和不变。

    证明：`(Uψ)‾ ⬝ᵥ (Uψ) = ψ̄ ⬝ᵥ Uᴴ(Uψ) = ψ̄ ⬝ᵥ (UᴴU)ψ = ψ̄ ⬝ᵥ ψ`，
    两边取实部并用上面的点积写法即得。 -/
theorem z117_normSq_mulVec_of_conjTranspose_mul_self {n : ℕ} (M : Matrix (Fin n) (Fin n) ℂ)
    (hM : M.conjTranspose * M = 1) (ψ : Fin n → ℂ) :
    z116NormSq (M *ᵥ ψ) = z116NormSq ψ := by
  have key : star (M *ᵥ ψ) ⬝ᵥ (M *ᵥ ψ) = star ψ ⬝ᵥ ψ := by
    rw [Matrix.star_mulVec, ← Matrix.dotProduct_mulVec, Matrix.mulVec_mulVec, hM,
      Matrix.one_mulVec]
  have hkey := congrArg Complex.re key
  rwa [z117_normSq_eq_re_dotProduct, z117_normSq_eq_re_dotProduct] at hkey

/-- **Z117 · 反厄米性**：`H` 厄米 ⟹ `−i t H` 反厄米（`Aᴴ = −A`），这是酉性的来源 -/
theorem z117_skew_adjoint {n : ℕ} (H : Matrix (Fin n) (Fin n) ℂ) (hH : H.IsHermitian) (t : ℝ) :
    ((-(t : ℂ) * Complex.I) • H).conjTranspose = -((-(t : ℂ) * Complex.I) • H) := by
  -- 标量层面：`star (−t·i) = −(−t·i)`（实系数 t 取共轭不变，i 变号两次）
  have hstar : star (-(t : ℂ) * Complex.I) = -(-(t : ℂ) * Complex.I) := by
    simp only [Complex.star_def, map_mul, map_neg, Complex.conj_I, Complex.conj_ofReal]
    ring
  rw [Matrix.conjTranspose_smul, hH, hstar, ← neg_smul]

/-- **Z117 · 厄米 H 给酉演化**：`U = exp(−i t H)` 满足 `Uᴴ U = 1` 且 `U Uᴴ = 1` -/
theorem z117_evolution_unitary {n : ℕ} (H : Matrix (Fin n) (Fin n) ℂ) (hH : H.IsHermitian) (t : ℝ) :
    let U := NormedSpace.exp ((-(t : ℂ) * Complex.I) • H)
    U.conjTranspose * U = 1 ∧ U * U.conjTranspose = 1 :=
  z117_exp_conjTranspose_mul_self _ (z117_skew_adjoint H hH t)

/-- **Z117 · 模方守恒**：酉演化不改变模方和 -/
theorem z117_normSq_preserved {n : ℕ} (H : Matrix (Fin n) (Fin n) ℂ) (hH : H.IsHermitian) (t : ℝ)
    (ψ : Fin n → ℂ) :
    z116NormSq ((NormedSpace.exp ((-(t : ℂ) * Complex.I) • H)) *ᵥ ψ) = z116NormSq ψ :=
  z117_normSq_mulVec_of_conjTranspose_mul_self _
    (z117_exp_conjTranspose_mul_self _ (z117_skew_adjoint H hH t)).1 ψ

/-- **Z117 · 概率和归一**（Born 规则的算术内容）：初态模方和为 1 ⟹ 末态概率和仍为 1

    注：`*ᵥ` 在 Mathlib 里是 `scoped infixr:73`（`Matrix.mulVec`），优先级低于函数应用，
    因此"取第 i 个分量"必须显式加括号写成 `((… ) *ᵥ ψ) i`；不加括号时 `… *ᵥ ψ i`
    会被解析成 `… *ᵥ (ψ i)` 而类型不符。这里只补了这对括号，命题内容不变。 -/
theorem z117_probability_sum_one {n : ℕ} (H : Matrix (Fin n) (Fin n) ℂ) (hH : H.IsHermitian) (t : ℝ)
    (ψ : Fin n → ℂ) (hψ : z116NormSq ψ = 1) :
    ∑ i, Complex.normSq (((NormedSpace.exp ((-(t : ℂ) * Complex.I) • H)) *ᵥ ψ) i) = 1 := by
  have h := z117_normSq_preserved H hH t ψ
  rw [hψ] at h
  simpa only [z116NormSq] using h

#print axioms z117_exp_conjTranspose_mul_self
#print axioms z117_normSq_eq_re_dotProduct
#print axioms z117_normSq_mulVec_of_conjTranspose_mul_self
#print axioms z117_skew_adjoint
#print axioms z117_evolution_unitary
#print axioms z117_normSq_preserved
#print axioms z117_probability_sum_one

end Zero