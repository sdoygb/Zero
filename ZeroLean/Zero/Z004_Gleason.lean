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

end Zero
