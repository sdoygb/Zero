/-
  Zb014 §1 的矩阵核：Hermitian 矩阵 + `i b I`（`b ≠ 0`）的核是平凡的。

  文章里的 `M = c²L − f² I + 2iγf I` 是 `A + i b I` 的形状，其中
  `A = c²L − f² I` 是 Hermitian、`b = 2γf`。这里只证一般引理：
    `Aᴴ = A`、`b ≠ 0` ⟹ `(A + (b·I) • 1) *ᵥ v = 0 → v = 0`。

  证明：对 `M v = 0` 取 `⟨v, M v⟩ = 0`；实部 `⟨v, A v⟩` 是实数（Hermitian），
  虚部给 `b · ‖v‖² = 0`；`b ≠ 0` 与 `v ≠ 0` 矛盾。
-/

import Mathlib

namespace Zero

open Matrix

/-- **Zb014 · Hermitian + `ibI` 的核平凡**：
    `Aᴴ = A`、`b ≠ 0` 时，`A + (b·I) • 1` 的核只有 `0`。 -/
theorem zb014_hermitian_add_i_smul_one_kernel
    {n : Type*} [Fintype n] [DecidableEq n]
    (A : Matrix n n ℂ) (hA : A.IsHermitian) (b : ℝ) (hb : b ≠ 0) :
    ∀ v : n → ℂ, (A + (b * Complex.I) • (1 : Matrix n n ℂ)) *ᵥ v = 0 → v = 0 := by
  intro v hv
  by_cases hv0 : v = 0
  · exact hv0
  · exfalso
    have hinner : dotProduct (star v)
        (Matrix.mulVec (A + (b * Complex.I) • (1 : Matrix n n ℂ)) v) = 0 := by
      simpa using congrArg (fun x => dotProduct (star v) x) hv
    rw [Matrix.add_mulVec] at hinner
    rw [Matrix.smul_mulVec] at hinner
    rw [Matrix.one_mulVec] at hinner
    rw [dotProduct_add] at hinner
    rw [dotProduct_smul] at hinner
    set z : ℂ := dotProduct (star v) (Matrix.mulVec A v) with hzdef
    set w : ℂ := dotProduct (star v) v with hwdef
    have hz : z = star z := by
      have h := hA.star_dotProduct_mulVec_comm v v
      simpa [z] using h.symm
    have hw : star w = w := by
      simp [w, dotProduct, mul_comm]
    have hz_im : z.im = 0 := by
      have h := congrArg Complex.im hz
      simp at h
      linarith
    have hw_im : w.im = 0 := by
      have h := congrArg Complex.im hw
      simp at h
      linarith
    have hinner' : z + (b * Complex.I) * w = 0 := by
      simpa [z, w, smul_eq_mul] using hinner
    have him : b * w.re = 0 := by
      have h := congrArg Complex.im hinner'
      simp only [Complex.add_im, Complex.mul_im, Complex.mul_re,
        Complex.ofReal_re, Complex.ofReal_im, Complex.I_re, Complex.I_im,
        hz_im, hw_im, mul_zero, zero_mul, add_zero, zero_add] at h
      simpa using h
    have hwre : w.re = 0 := by
      rcases mul_eq_zero.mp him with hb0 | hw0
      · exact absurd hb0 hb
      · exact hw0
    have hw0 : w = 0 := by
      apply Complex.ext
      · simp [hwre]
      · simp [hw_im]
    have hw_re0 : w.re = 0 := by rw [hw0]; simp
    have hw_sum : w.re = ∑ i, Complex.normSq (v i) := by
      rw [show w = ∑ i, star (v i) * v i by rfl]
      rw [Complex.re_sum]
      apply Finset.sum_congr rfl
      intro i _
      rw [Complex.star_def, Complex.mul_re, Complex.conj_re, Complex.conj_im, Complex.normSq_apply]
      ring
    have hsum0 : ∑ i, Complex.normSq (v i) = 0 := by
      rw [← hw_sum]
      exact hw_re0
    have hterm : ∀ i, Complex.normSq (v i) = 0 := by
      intro i
      exact (Finset.sum_eq_zero_iff_of_nonneg (fun j _ => Complex.normSq_nonneg (v j))).mp
        hsum0 i (Finset.mem_univ i)
    have hv_eq : v = 0 := by
      funext i
      exact Complex.normSq_eq_zero.mp (hterm i)
    exact hv0 hv_eq

#print axioms zb014_hermitian_add_i_smul_one_kernel

end Zero
