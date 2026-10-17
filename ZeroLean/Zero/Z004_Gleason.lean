/-
  Z004 里那条「Gleason 定理的假设在该模型的格上不成立」的**最小形式**。

  Gleason 定理要的是"投影格上的测度都由态诱导"。二维（ℂ²）时这条前提坏掉：Born 规则
  给出的 {0,1}-赋值必然**集中在一根射线上** —— 也就是说，两个给概率 1 的方向必须是同一
  根射线，于是这种赋值无法在投影格上"连续地"铺开。下面这条定理就是那个刚性步骤。

  完整的有限反例（二维投影格上存在**非**态诱导的有限可加测度）依赖超滤子/选择，
  留作后续；此处只做可判定的刚性子步骤，并在 Z 篇里如实标注范围。
-/

import Mathlib

namespace Zero

/-- **Z004 · Born 赋值的刚性**：单位态 `ψ` 对单位向量 `u`、`v` 的 Born 概率都是 1 时，
    `u` 与 `v` 必是彼此的相位倍数（Cauchy–Schwarz 取等情形）。
    推论：{0,1}-赋值不可能区分两条不同射线 —— Gleason 前提在二维失效的那一步。 -/
theorem z004_born_one_forces_parallel {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℂ E]
    (ψ u v : E) (hψ : ‖ψ‖ = 1) (hu : ‖u‖ = 1) (hv : ‖v‖ = 1)
    (hu1 : ‖inner ℂ ψ u‖ = 1) (hv1 : ‖inner ℂ ψ v‖ = 1) :
    ∃ r : ℂ, r ≠ 0 ∧ v = r • u := by
  have hψ0 : ψ ≠ 0 := by intro h; rw [h, norm_zero] at hψ; norm_num at hψ
  have hu0 : u ≠ 0 := by intro h; rw [h, norm_zero] at hu; norm_num at hu
  have hv0 : v ≠ 0 := by intro h; rw [h, norm_zero] at hv; norm_num at hv
  -- u 与 ψ 同向
  have hpu : ∃ r : ℂ, r ≠ 0 ∧ u = r • ψ := by
    refine (norm_inner_eq_norm_iff (𝕜 := ℂ) hψ0 hu0).mp ?_
    rw [hψ, hu, mul_one]
    exact hu1
  -- v 与 ψ 同向
  have hpv : ∃ s : ℂ, s ≠ 0 ∧ v = s • ψ := by
    refine (norm_inner_eq_norm_iff (𝕜 := ℂ) hψ0 hv0).mp ?_
    rw [hψ, hv, mul_one]
    exact hv1
  obtain ⟨r, hr0, hur⟩ := hpu
  obtain ⟨s, hs0, hsv⟩ := hpv
  refine ⟨s * r⁻¹, mul_ne_zero hs0 (inv_ne_zero hr0), ?_⟩
  have hψu : ψ = r⁻¹ • u := by
    rw [hur, smul_smul, inv_mul_cancel₀ hr0, one_smul]
  rw [hsv, hψu, smul_smul]

#print axioms z004_born_one_forces_parallel

/-!
## 命题 10.1 的主体：频率 6 与频率 2 的 Fourier 正交性

`f(θ) = (1 + cos(6θ))/2` 是一个显式反例：它满足 (i) `f ≥ 0`、(ii) `f(θ)+f(θ+π/2)=1`，
但它的频率 6 分量非零；而任何态诱导的 Born 概率都有形式 `a + b cos(2θ) + c sin(2θ)`，
频率 6 分量为零。故 `f` 不是态诱导的。
-/

open Real intervalIntegral MeasureTheory

theorem z004_gleason_int_cos_nat_mul (k : ℕ) (hk : k ≠ 0) :
    ∫ θ in (0 : ℝ)..(2 * π), Real.cos ((k : ℝ) * θ) = 0 := by
  have hk' : (k : ℝ) ≠ 0 := by exact_mod_cast hk
  have h := intervalIntegral.integral_comp_mul_left (f := Real.cos) (a := (0 : ℝ))
    (b := (2 * π)) (c := (k : ℝ)) hk'
  rw [h, integral_cos]
  have h2 : (k : ℝ) * (2 * π) = ((2 * k : ℕ) : ℝ) * π := by push_cast; ring
  rw [h2, Real.sin_nat_mul_pi]
  simp

theorem z004_gleason_int_sin_nat_mul (k : ℕ) (hk : k ≠ 0) :
    ∫ θ in (0 : ℝ)..(2 * π), Real.sin ((k : ℝ) * θ) = 0 := by
  have hk' : (k : ℝ) ≠ 0 := by exact_mod_cast hk
  have h := intervalIntegral.integral_comp_mul_left (f := Real.sin) (a := (0 : ℝ))
    (b := (2 * π)) (c := (k : ℝ)) hk'
  rw [h, integral_sin]
  have h2 : (k : ℝ) * (2 * π) = ((2 * k : ℕ) : ℝ) * π := by push_cast; ring
  rw [h2]
  have hc : Real.cos (((2 * k : ℕ) : ℝ) * π) = 1 := by
    rw [Real.cos_nat_mul_pi]; norm_num
  rw [hc]; simp

theorem z004_gleason_int_cos6_cos2 :
    ∫ θ in (0 : ℝ)..(2 * π), Real.cos (6 * θ) * Real.cos (2 * θ) = 0 := by
  have h : ∀ θ : ℝ, Real.cos (6 * θ) * Real.cos (2 * θ)
      = (Real.cos (4 * θ) + Real.cos (8 * θ)) / 2 := by
    intro θ
    have := Real.two_mul_cos_mul_cos (6 * θ) (2 * θ)
    have e1 : 6 * θ - 2 * θ = 4 * θ := by ring
    have e2 : 6 * θ + 2 * θ = 8 * θ := by ring
    rw [e1, e2] at this
    linarith
  simp_rw [h]
  rw [intervalIntegral.integral_div]
  have c4 : Continuous (fun θ : ℝ => Real.cos (4 * θ)) := by fun_prop
  have c8 : Continuous (fun θ : ℝ => Real.cos (8 * θ)) := by fun_prop
  have i4 : IntervalIntegrable (fun θ : ℝ => Real.cos (4 * θ)) volume 0 (2 * π) :=
    c4.intervalIntegrable _ _
  have i8 : IntervalIntegrable (fun θ : ℝ => Real.cos (8 * θ)) volume 0 (2 * π) :=
    c8.intervalIntegrable _ _
  have h4 : ∫ θ in (0 : ℝ)..(2 * π), Real.cos (4 * θ) = 0 := by
    simpa using z004_gleason_int_cos_nat_mul 4 (by norm_num)
  have h8 : ∫ θ in (0 : ℝ)..(2 * π), Real.cos (8 * θ) = 0 := by
    simpa using z004_gleason_int_cos_nat_mul 8 (by norm_num)
  rw [intervalIntegral.integral_add i4 i8, h4, h8]
  norm_num

theorem z004_gleason_int_cos6_sin2 :
    ∫ θ in (0 : ℝ)..(2 * π), Real.cos (6 * θ) * Real.sin (2 * θ) = 0 := by
  have h : ∀ θ : ℝ, Real.cos (6 * θ) * Real.sin (2 * θ)
      = (Real.sin (8 * θ) - Real.sin (4 * θ)) / 2 := by
    intro θ
    have := Real.two_mul_sin_mul_cos (2 * θ) (6 * θ)
    have e1 : 2 * θ - 6 * θ = -(4 * θ) := by ring
    have e2 : 2 * θ + 6 * θ = 8 * θ := by ring
    rw [e1, e2, Real.sin_neg] at this
    linarith
  simp_rw [h]
  rw [intervalIntegral.integral_div]
  have c8 : Continuous (fun θ : ℝ => Real.sin (8 * θ)) := by fun_prop
  have c4 : Continuous (fun θ : ℝ => Real.sin (4 * θ)) := by fun_prop
  have i8 : IntervalIntegrable (fun θ : ℝ => Real.sin (8 * θ)) volume 0 (2 * π) :=
    c8.intervalIntegrable _ _
  have i4 : IntervalIntegrable (fun θ : ℝ => Real.sin (4 * θ)) volume 0 (2 * π) :=
    c4.intervalIntegrable _ _
  have h8 : ∫ θ in (0 : ℝ)..(2 * π), Real.sin (8 * θ) = 0 := by
    simpa using z004_gleason_int_sin_nat_mul 8 (by norm_num)
  have h4 : ∫ θ in (0 : ℝ)..(2 * π), Real.sin (4 * θ) = 0 := by
    simpa using z004_gleason_int_sin_nat_mul 4 (by norm_num)
  rw [intervalIntegral.integral_sub i8 i4, h8, h4]
  norm_num

theorem z004_gleason_int_cos6_sq :
    ∫ θ in (0 : ℝ)..(2 * π), Real.cos (6 * θ) ^ 2 = π := by
  have h : ∀ θ : ℝ, Real.cos (6 * θ) ^ 2 = (1 / 2 : ℝ) + Real.cos (12 * θ) / 2 := by
    intro θ
    rw [Real.cos_sq]
    have e : 2 * (6 * θ) = 12 * θ := by ring
    rw [e]
  simp_rw [h]
  have c1 : Continuous (fun _ : ℝ => (1 / 2 : ℝ)) := continuous_const
  have c2 : Continuous (fun θ : ℝ => Real.cos (12 * θ) / 2) := by fun_prop
  have i1 : IntervalIntegrable (fun _ : ℝ => (1 / 2 : ℝ)) volume 0 (2 * π) :=
    c1.intervalIntegrable _ _
  have i2 : IntervalIntegrable (fun θ : ℝ => Real.cos (12 * θ) / 2) volume 0 (2 * π) :=
    c2.intervalIntegrable _ _
  rw [intervalIntegral.integral_add i1 i2]
  have hc : ∫ θ in (0 : ℝ)..(2 * π), Real.cos (12 * θ) / 2 = 0 := by
    rw [intervalIntegral.integral_div]
    have h12 : ∫ θ in (0 : ℝ)..(2 * π), Real.cos (12 * θ) = 0 := by
      simpa using z004_gleason_int_cos_nat_mul 12 (by norm_num)
    rw [h12]; norm_num
  rw [hc, intervalIntegral.integral_const]
  ring_nf

theorem z004_gleason_int_state_form (a b c : ℝ) :
    ∫ θ in (0 : ℝ)..(2 * π),
      (a + b * Real.cos (2 * θ) + c * Real.sin (2 * θ)) * Real.cos (6 * θ) = 0 := by
  have h : ∀ θ : ℝ, (a + b * Real.cos (2 * θ) + c * Real.sin (2 * θ)) * Real.cos (6 * θ)
      = a * Real.cos (6 * θ) + b * (Real.cos (2 * θ) * Real.cos (6 * θ))
        + c * (Real.sin (2 * θ) * Real.cos (6 * θ)) := by
    intro θ; ring
  simp_rw [h]
  have ca : Continuous (fun θ : ℝ => a * Real.cos (6 * θ)) := by fun_prop
  have cb : Continuous (fun θ : ℝ => b * (Real.cos (2 * θ) * Real.cos (6 * θ))) := by fun_prop
  have cc : Continuous (fun θ : ℝ => c * (Real.sin (2 * θ) * Real.cos (6 * θ))) := by fun_prop
  have ia : IntervalIntegrable (fun θ : ℝ => a * Real.cos (6 * θ)) volume 0 (2 * π) :=
    ca.intervalIntegrable _ _
  have ib : IntervalIntegrable (fun θ : ℝ => b * (Real.cos (2 * θ) * Real.cos (6 * θ)))
      volume 0 (2 * π) := cb.intervalIntegrable _ _
  have ic : IntervalIntegrable (fun θ : ℝ => c * (Real.sin (2 * θ) * Real.cos (6 * θ)))
      volume 0 (2 * π) := cc.intervalIntegrable _ _
  rw [intervalIntegral.integral_add (by exact ia.add ib) ic]
  rw [intervalIntegral.integral_add ia ib]
  rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul]
  have h6 : ∫ θ in (0 : ℝ)..(2 * π), Real.cos (6 * θ) = 0 := by
    simpa using z004_gleason_int_cos_nat_mul 6 (by norm_num)
  rw [h6]
  have hb : ∫ x in (0 : ℝ)..(2 * π), Real.cos (2 * x) * Real.cos (6 * x) = 0 := by
    simpa [mul_comm] using z004_gleason_int_cos6_cos2
  have hc2 : ∫ x in (0 : ℝ)..(2 * π), Real.sin (2 * x) * Real.cos (6 * x) = 0 := by
    simpa [mul_comm] using z004_gleason_int_cos6_sin2
  rw [hb, hc2]
  ring

noncomputable def z004_gleason_f (θ : ℝ) : ℝ := (1 + Real.cos (6 * θ)) / 2

theorem z004_gleason_f_nonneg (θ : ℝ) : 0 ≤ z004_gleason_f θ := by
  unfold z004_gleason_f
  have h : -1 ≤ Real.cos (6 * θ) := Real.neg_one_le_cos _
  linarith

theorem z004_gleason_f_add (θ : ℝ) :
    z004_gleason_f θ + z004_gleason_f (θ + π / 2) = 1 := by
  unfold z004_gleason_f
  have h : Real.cos (6 * (θ + π / 2)) = - Real.cos (6 * θ) := by
    have e : 6 * (θ + π / 2) = 6 * θ + 3 * π := by ring
    rw [e, Real.cos_add]
    have h3 : Real.cos (3 * π) = -1 := by
      have e3 : 3 * π = π + 2 * π := by ring
      rw [e3, Real.cos_add_two_pi, Real.cos_pi]
    have h3s : Real.sin (3 * π) = 0 := by
      have e3 : 3 * π = π + 2 * π := by ring
      rw [e3, Real.sin_add_two_pi, Real.sin_pi]
    rw [h3, h3s]
    ring
  rw [h]; ring

theorem z004_gleason_f_int_cos6 :
    ∫ θ in (0 : ℝ)..(2 * π), z004_gleason_f θ * Real.cos (6 * θ) = π / 2 := by
  have h : ∀ θ : ℝ, z004_gleason_f θ * Real.cos (6 * θ)
      = (Real.cos (6 * θ) + Real.cos (6 * θ) ^ 2) / 2 := by
    intro θ; unfold z004_gleason_f; ring
  simp_rw [h]
  rw [intervalIntegral.integral_div]
  have c1 : Continuous (fun θ : ℝ => Real.cos (6 * θ)) := by fun_prop
  have c2 : Continuous (fun θ : ℝ => Real.cos (6 * θ) ^ 2) := by fun_prop
  have i1 : IntervalIntegrable (fun θ : ℝ => Real.cos (6 * θ)) volume 0 (2 * π) :=
    c1.intervalIntegrable _ _
  have i2 : IntervalIntegrable (fun θ : ℝ => Real.cos (6 * θ) ^ 2) volume 0 (2 * π) :=
    c2.intervalIntegrable _ _
  rw [intervalIntegral.integral_add i1 i2]
  have h6 : ∫ θ in (0 : ℝ)..(2 * π), Real.cos (6 * θ) = 0 := by
    simpa using z004_gleason_int_cos_nat_mul 6 (by norm_num)
  rw [h6, z004_gleason_int_cos6_sq]
  ring

/-- **命题 10.1 的主体（机器版）**：显式函数 `f(θ)=(1+cos6θ)/2` 不是态诱导的——
    态诱导的 Born 概率只有 `a+b cos2θ+c sin2θ` 的形式（频率 6 分量为零），
    而 `∫ f·cos6θ = π/2 ≠ 0`。 -/
theorem z004_gleason_not_state_induced :
    ¬ ∃ a b c : ℝ, ∀ θ : ℝ,
      z004_gleason_f θ = a + b * Real.cos (2 * θ) + c * Real.sin (2 * θ) := by
  rintro ⟨a, b, c, h⟩
  have h0 := z004_gleason_int_state_form a b c
  have h1 := z004_gleason_f_int_cos6
  have heq : (fun θ : ℝ => z004_gleason_f θ * Real.cos (6 * θ))
      = (fun θ : ℝ => (a + b * Real.cos (2 * θ) + c * Real.sin (2 * θ)) * Real.cos (6 * θ)) := by
    funext θ; rw [h θ]
  rw [heq, h0] at h1
  linarith [Real.pi_pos]

#print axioms z004_gleason_not_state_induced

end Zero
