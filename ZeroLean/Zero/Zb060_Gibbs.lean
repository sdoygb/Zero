/-
  Zb060 §1–§2 的对角核：Gibbs 配分函数与第一定律。

  正文用矩阵 `ρ_β = e^{−βH}/Z`；这里取 `H` 对角（`E : Fin n → ℝ`）的交换情形，
  把同一组代数恒等式写成有限和：
    · `p_i = e^{−βE_i}/Z`、`Z = Σ_i e^{−βE_i}`；
    · 归一：`Σ_i p_i = 1`；
    · 第一定律：`−Σ_i p_i log p_i = β Σ_i p_i E_i + log Z`。

  未覆盖：矩阵版 Gibbs 变分（Klein 不等式/最大熵）、pinching 第二定律、
  `T = 1/β` 与模温度 `β = 1` 的区分——那些需要矩阵泛函演算与 majorization。
-/

import Mathlib

namespace Zero

open Finset

/-- 配分函数 `Z(β) = Σ_i e^{−β E_i}`。 -/
noncomputable def zb060Z {n : ℕ} (β : ℝ) (E : Fin n → ℝ) : ℝ :=
  ∑ i, Real.exp (-β * E i)

/-- Gibbs 权重 `p_i = e^{−β E_i} / Z`。 -/
noncomputable def zb060p {n : ℕ} (β : ℝ) (E : Fin n → ℝ) (i : Fin n) : ℝ :=
  Real.exp (-β * E i) / zb060Z β E

/-- **Zb060 · 配分函数归一**：`Z ≠ 0` 时 `Σ_i p_i = 1`。 -/
theorem zb060_partition_normalized {n : ℕ} (β : ℝ) (E : Fin n → ℝ)
    (hZ : zb060Z β E ≠ 0) :
    ∑ i, zb060p β E i = 1 := by
  have hnum : (∑ i, Real.exp (-(β * E i))) = zb060Z β E := by
    rw [zb060Z]
    exact Finset.sum_congr rfl (fun i _ => by congr 1; ring)
  have hsum : (∑ i, zb060p β E i) = (∑ i, Real.exp (-(β * E i))) / zb060Z β E := by
    simp [zb060p, ← Finset.sum_div]
  rw [hsum, hnum, div_self hZ]

/-- **Zb060 · 第一定律（对角/Gibbs 情形）**：
    `−Σ_i p_i log p_i = β Σ_i p_i E_i + log Z`。 -/
theorem zb060_first_law {n : ℕ} (β : ℝ) (E : Fin n → ℝ) (hZ : 0 < zb060Z β E) :
    -∑ i, zb060p β E i * Real.log (zb060p β E i) =
      β * ∑ i, zb060p β E i * E i + Real.log (zb060Z β E) := by
  have hZne : zb060Z β E ≠ 0 := ne_of_gt hZ
  have hlog : ∀ i, Real.log (zb060p β E i) = -β * E i - Real.log (zb060Z β E) := by
    intro i
    rw [zb060p, Real.log_div (Real.exp_ne_zero _) hZne, Real.log_exp]
  have hsum : ∑ i, zb060p β E i = 1 := zb060_partition_normalized β E hZne
  have hterm : ∀ i, zb060p β E i * Real.log (zb060p β E i)
      = -β * (zb060p β E i * E i) - Real.log (zb060Z β E) * zb060p β E i := by
    intro i
    rw [hlog i]
    ring
  calc
    -∑ i, zb060p β E i * Real.log (zb060p β E i)
        = -∑ i, (-β * (zb060p β E i * E i)
            - Real.log (zb060Z β E) * zb060p β E i) := by
          rw [Finset.sum_congr rfl (fun i _ => hterm i)]
    _ = β * ∑ i, zb060p β E i * E i + Real.log (zb060Z β E) := by
          rw [Finset.sum_sub_distrib, neg_sub]
          have h1 : (∑ i, -β * (zb060p β E i * E i))
              = -β * ∑ i, zb060p β E i * E i := by
            rw [← Finset.mul_sum]
          have h2 : (∑ i, Real.log (zb060Z β E) * zb060p β E i)
              = Real.log (zb060Z β E) := by
            rw [← Finset.mul_sum, hsum, mul_one]
          rw [h1, h2]
          ring

#print axioms zb060_partition_normalized

/-- 有限概率向量上的 KL 散度（对角情形）。 -/
noncomputable def zb060KL {n : ℕ} (p q : Fin n → ℝ) : ℝ :=
  ∑ i, p i * Real.log (p i / q i)

/-- **Zb060 · KL 非负（对角 Gibbs 不等式）**：
    `p,q` 为正概率向量时 `KL(p‖q) ≥ 0`。 -/
theorem zb060_kl_nonneg {n : ℕ} (p q : Fin n → ℝ)
    (hp : ∀ i, 0 < p i) (hq : ∀ i, 0 < q i)
    (hsp : ∑ i, p i = 1) (hsq : ∑ i, q i = 1) :
    0 ≤ zb060KL p q := by
  have hlog : ∀ i, Real.log (p i / q i) = -Real.log (q i / p i) := by
    intro i
    rw [Real.log_div (ne_of_gt (hp i)) (ne_of_gt (hq i)),
      Real.log_div (ne_of_gt (hq i)) (ne_of_gt (hp i))]
    ring
  have hbound : (∑ i, p i * Real.log (q i / p i)) ≤ 0 := by
    have h1 : (∑ i, p i * Real.log (q i / p i))
        ≤ ∑ i, p i * (q i / p i - 1) := by
      exact Finset.sum_le_sum fun i _ =>
        mul_le_mul_of_nonneg_left
          (Real.log_le_sub_one_of_pos (div_pos (hq i) (hp i)))
          (le_of_lt (hp i))
    have h2 : (∑ i, p i * (q i / p i - 1)) = ∑ i, (q i - p i) := by
      apply Finset.sum_congr rfl
      intro i _
      field_simp [ne_of_gt (hp i)]
    have h3 : (∑ i, (q i - p i)) = 0 := by
      rw [Finset.sum_sub_distrib, hsp, hsq]
      ring
    linarith
  have hkl : zb060KL p q = -∑ i, p i * Real.log (q i / p i) := by
    simp only [zb060KL, hlog, Finset.sum_neg_distrib, mul_neg]
  rw [hkl]
  linarith




/-- **Zb060 · KL 的等号情形**：正概率向量上 `KL(p‖q)=0 ⟹ p=q`。
    这是 Gibbs 变分「唯一性」的对角核心。 -/
theorem zb060_kl_eq_zero {n : ℕ} (p q : Fin n → ℝ)
    (hp : ∀ i, 0 < p i) (hq : ∀ i, 0 < q i)
    (hsp : ∑ i, p i = 1) (hsq : ∑ i, q i = 1)
    (hKL : zb060KL p q = 0) : p = q := by
  have hlog : ∀ i, Real.log (p i / q i) = -Real.log (q i / p i) := by
    intro i
    rw [Real.log_div (ne_of_gt (hp i)) (ne_of_gt (hq i)),
      Real.log_div (ne_of_gt (hq i)) (ne_of_gt (hp i))]
    ring
  have hsum_f : (∑ i, p i * Real.log (q i / p i)) = 0 := by
    have hkl : zb060KL p q = -∑ i, p i * Real.log (q i / p i) := by
      simp only [zb060KL, hlog, Finset.sum_neg_distrib, mul_neg]
    linarith
  have hsum_g : (∑ i, p i * (q i / p i - 1)) = 0 := by
    have h2 : (∑ i, p i * (q i / p i - 1)) = ∑ i, (q i - p i) := by
      apply Finset.sum_congr rfl
      intro i _
      field_simp [ne_of_gt (hp i)]
    have h3 : (∑ i, (q i - p i)) = 0 := by
      rw [Finset.sum_sub_distrib, hsp, hsq]
      ring
    rw [h2, h3]
  have hle : ∀ i, p i * Real.log (q i / p i) ≤ p i * (q i / p i - 1) := by
    intro i
    exact mul_le_mul_of_nonneg_left
      (Real.log_le_sub_one_of_pos (div_pos (hq i) (hp i)))
      (le_of_lt (hp i))
  have hsums : (∑ i, p i * Real.log (q i / p i))
      = ∑ i, p i * (q i / p i - 1) := by
    rw [hsum_f, hsum_g]
  have hle' : ∀ i ∈ (Finset.univ : Finset (Fin n)),
      p i * Real.log (q i / p i) ≤ p i * (q i / p i - 1) := by
    intro i _
    exact hle i
  have hterm : ∀ i, p i * Real.log (q i / p i) = p i * (q i / p i - 1) := by
    have h := (Finset.sum_eq_sum_iff_of_le hle').mp hsums
    exact fun i => h i (Finset.mem_univ i)
  funext i
  have hlog_eq : Real.log (q i / p i) = q i / p i - 1 := by
    exact mul_left_cancel₀ (ne_of_gt (hp i)) (hterm i)
  have hx : q i / p i = 1 := by
    by_contra hne
    have hlt := Real.log_lt_sub_one_of_pos (div_pos (hq i) (hp i)) hne
    linarith
  field_simp [ne_of_gt (hp i)] at hx
  linarith


/-- **Zb060 · 2×2 对角 Gibbs 变分**（`Fin 2` 的显式情形）：
    自由能差 = `(1/β) KL(p‖p_β) ≥ 0`。 -/
theorem zb060_gibbs_variational_two (β : ℝ) (hβ : 0 < β) (E : Fin 2 → ℝ)
    (p : Fin 2 → ℝ) (hp : ∀ i, 0 < p i) (hsp : ∑ i, p i = 1) :
    (∑ i, p i * E i) - (1 / β) * (-∑ i, p i * Real.log (p i))
      ≥ (∑ i, zb060p β E i * E i)
          - (1 / β) * (-∑ i, zb060p β E i * Real.log (zb060p β E i)) := by
  have hZ : 0 < zb060Z β E := by
    rw [zb060Z, Fin.sum_univ_two]
    positivity
  have hZne : zb060Z β E ≠ 0 := ne_of_gt hZ
  have hpβ : ∀ i, 0 < zb060p β E i := by
    intro i
    rw [zb060p]
    exact div_pos (Real.exp_pos _) hZ
  have hspβ : ∑ i, zb060p β E i = 1 := zb060_partition_normalized β E hZne
  have hKL := zb060_kl_nonneg p (zb060p β E) hp hpβ hsp hspβ
  have hp01 : p 0 + p 1 = 1 := by simpa [Fin.sum_univ_two] using hsp
  have hlogβ : ∀ i, Real.log (zb060p β E i) = -β * E i - Real.log (zb060Z β E) := by
    intro i
    rw [zb060p, Real.log_div (Real.exp_ne_zero _) hZne, Real.log_exp]
  have hFp : (∑ i, p i * E i) - (1 / β) * (-∑ i, p i * Real.log (p i))
      = (∑ i, p i * E i) + (1 / β) * ∑ i, p i * Real.log (p i) := by
    ring
  have hFβ : (∑ i, zb060p β E i * E i)
        - (1 / β) * (-∑ i, zb060p β E i * Real.log (zb060p β E i))
      = -Real.log (zb060Z β E) / β := by
    have hlaw := zb060_first_law β E hZ
    rw [hlaw]
    field_simp [ne_of_gt hβ]
    ring
  have hKL_eq : zb060KL p (zb060p β E)
      = (∑ i, p i * Real.log (p i)) + β * (∑ i, p i * E i)
        + Real.log (zb060Z β E) := by
    simp only [zb060KL, Fin.sum_univ_two]
    rw [Real.log_div (ne_of_gt (hp 0)) (ne_of_gt (hpβ 0)),
      Real.log_div (ne_of_gt (hp 1)) (ne_of_gt (hpβ 1)),
      hlogβ 0, hlogβ 1]
    ring_nf
    linear_combination (Real.log (zb060Z β E)) * hp01
  have hdiff : (∑ i, p i * E i) + (1 / β) * ∑ i, p i * Real.log (p i)
      - (-Real.log (zb060Z β E) / β) = (1 / β) * zb060KL p (zb060p β E) := by
    rw [hKL_eq]
    field_simp [ne_of_gt hβ]
    ring
  rw [hFp, hFβ]
  have hnonneg : 0 ≤ (1 / β) * zb060KL p (zb060p β E) :=
    mul_nonneg (le_of_lt (one_div_pos.mpr hβ)) hKL
  linarith

#print axioms zb060_gibbs_variational_two

#print axioms zb060_kl_eq_zero

#print axioms zb060_kl_nonneg

#print axioms zb060_first_law

end Zero
