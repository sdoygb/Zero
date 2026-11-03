/-
  Zb005 §1–§3 的机器检查版本：格点色散与逆映射。

  对应正文：
    · 定理 1.2：`(Lφ)_j = 4 sin²(q/2) φ_j`，其纯数学核是
      `2 − 2 cos q = 4 sin²(q/2)`；
    · 定理 3.2：`0 < f < 2c` 时 `ω(q) = f` 在 `(0, π)` 上有唯一解
      `q = 2 arcsin(f/(2c))`；下面的版本取 `c = 1`（无量纲化后一般 c 只是缩放）。

  证明策略（与正文 §1 的“倍角 + 平方关系”同一条）：
    cos q = 1 − 2 sin²(q/2)，代入 2 − 2cos q 即得恒等式；
    逆映射用 sin 在 [0, π/2] 上的严格单调性（由 `arcsin_sin` / `sin_arcsin` 取逆）。

  说明：本文件不引入任何 Zero 侧公理；依赖清单应只有
  `propext` / `Quot.sound` / `Classical.choice`（见文件末 `#print axioms`）。
-/

import Mathlib

namespace Zero

/-- **Zb005 · 色散恒等式**：`2 − 2cos q = 4 sin²(q/2)`。 -/
theorem zb005_dispersion_identity (q : ℝ) :
    2 - 2 * Real.cos q = 4 * Real.sin (q / 2) ^ 2 := by
  have h2 : Real.cos q = 1 - 2 * Real.sin (q / 2) ^ 2 := by
    have h := Real.cos_two_mul' (q / 2)
    rw [show 2 * (q / 2) = q by ring, Real.cos_sq'] at h
    linarith
  linarith

/-- **Zb005 · 逆映射（正向）**：`q ∈ [0, π]`、`0 ≤ f`、`2 − 2cos q = f²`
    ⟹ `q = 2·arcsin(f/2)`。 -/
theorem zb005_inverse_of_dispersion (q f : ℝ) (hq0 : 0 ≤ q) (hqpi : q ≤ Real.pi)
    (hf0 : 0 ≤ f) (h : 2 - 2 * Real.cos q = f ^ 2) :
    q = 2 * Real.arcsin (f / 2) := by
  have hq2_nonneg : 0 ≤ q / 2 := by linarith
  have hq2_le : q / 2 ≤ Real.pi / 2 := by linarith
  have hsin_nonneg : 0 ≤ Real.sin (q / 2) :=
    Real.sin_nonneg_of_nonneg_of_le_pi hq2_nonneg (by linarith)
  have hdisp := zb005_dispersion_identity q
  have hsq : (2 * Real.sin (q / 2)) ^ 2 = f ^ 2 := by nlinarith [hdisp, h]
  have hsin_eq : Real.sin (q / 2) = f / 2 := by
    have hprod : (2 * Real.sin (q / 2) - f) * (2 * Real.sin (q / 2) + f) = 0 := by
      nlinarith [hsq]
    rcases mul_eq_zero.mp hprod with hcase | hcase
    · linarith
    · linarith
  have harc : Real.arcsin (Real.sin (q / 2)) = q / 2 := Real.arcsin_sin (by linarith) hq2_le
  rw [hsin_eq] at harc
  linarith

/-- **Zb005 · 逆映射（逆向）**：`0 ≤ f ≤ 2`、`q = 2·arcsin(f/2)`
    ⟹ `2 − 2cos q = f²`。 -/
theorem zb005_dispersion_of_inverse (q f : ℝ)
    (hf0 : 0 ≤ f) (hf2 : f ≤ 2)
    (h : q = 2 * Real.arcsin (f / 2)) :
    2 - 2 * Real.cos q = f ^ 2 := by
  have hq2 : q / 2 = Real.arcsin (f / 2) := by linarith
  have hsin : Real.sin (q / 2) = f / 2 := by
    rw [hq2, Real.sin_arcsin (by linarith) (by linarith)]
  have hdisp := zb005_dispersion_identity q
  nlinarith [hdisp, hsin]

/-- **Zb005 · 长波极限的代数核**：`q ≠ 0` 时
    `(2 sin(q/2))/q = sin(q/2)/(q/2)`，故 `q → 0` 的极限就是 `sin x / x → 1`。
    这里只证代数改写；极限本身是 Mathlib 的 `Real.tendsto_sin_div` 一类引理。 -/
theorem zb005_longwave_algebra (q : ℝ) (hq : q ≠ 0) :
    2 * Real.sin (q / 2) / q = Real.sin (q / 2) / (q / 2) := by
  field_simp [hq]


/-- **Zb005 · 共扼闭包（反射分支）**：`q ∈ [−π,0]` 的逆色散。
    由原定理作用到 `−q`（`cos(−q)=cos q`）再取负。 -/
theorem zb005_inverse_of_dispersion_neg (q f : ℝ)
    (hq0 : -Real.pi ≤ q) (hqpi : q ≤ 0) (hf0 : 0 ≤ f)
    (h : 2 - 2 * Real.cos q = f ^ 2) :
    q = -2 * Real.arcsin (f / 2) := by
  have hq0' : 0 ≤ -q := by linarith
  have hqpi' : -q ≤ Real.pi := by linarith
  have hcos : 2 - 2 * Real.cos (-q) = f ^ 2 := by
    simpa [Real.cos_neg] using h
  have hneg := zb005_inverse_of_dispersion (-q) f hq0' hqpi' hf0 hcos
  linarith

/-- **Zb005 · 共扼闭包（全区间）**：`q ∈ [−π,π]` 的逆色散（正负两支）。 -/
theorem zb005_inverse_of_dispersion_full (q f : ℝ)
    (hq0 : -Real.pi ≤ q) (hqpi : q ≤ Real.pi) (hf0 : 0 ≤ f)
    (h : 2 - 2 * Real.cos q = f ^ 2) :
    q = 2 * Real.arcsin (f / 2) ∨ q = -2 * Real.arcsin (f / 2) := by
  rcases le_total 0 q with hq | hq
  · left
    exact zb005_inverse_of_dispersion q f hq hqpi hf0 h
  · right
    exact zb005_inverse_of_dispersion_neg q f hq0 hq hf0 h

#print axioms zb005_dispersion_identity
#print axioms zb005_inverse_of_dispersion
#print axioms zb005_dispersion_of_inverse
#print axioms zb005_longwave_algebra
#print axioms zb005_inverse_of_dispersion_neg
#print axioms zb005_inverse_of_dispersion_full

end Zero
