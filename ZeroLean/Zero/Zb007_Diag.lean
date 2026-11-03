/-
  Zb007 §2 的一般 `n` 对角情形：对角 `ρ = diag(λ)` 的模流矩阵元。

  与 `Zb007_ModularFlow.lean` 的 2×2 块相比，这里把矩阵元推广到任意有限维：
    `D(t) = diag(e^{it log λ_i})`、`E_{ab}` 是 (a,b) 处的矩阵单位；
    则 `(D(t) E_{ab} D(−t))_{ab} = e^{it log(λ_a/λ_b)}`。
  这正是「非对角元的模流频率 = log(λ_a/λ_b)」的逐元版本。
-/

import Zero.Zb007_ModularFlow

namespace Zero

open Matrix

/-- 对角 `ρ^{it}`：`diag(e^{it log λ_i})`。 -/
noncomputable def zb007Diag {n : Type*} [DecidableEq n] (lam : n → ℝ) (t : ℝ) :
    Matrix n n ℂ :=
  Matrix.diagonal (fun i => Complex.exp (Complex.I * t * Real.log (lam i)))

/-- `(a,b)` 处的矩阵单位。 -/
def zb007Unit {n : Type*} [DecidableEq n] (a b : n) : Matrix n n ℂ :=
  Matrix.of (fun i j => if i = a ∧ j = b then (1 : ℂ) else 0)

/-- **Zb007 · 一般 `n` 对角的模流矩阵元**：
    `(D(t) E_{ab} D(−t))_{ab} = e^{it log(λ_a/λ_b)}`。 -/
theorem zb007_diag_modular_flow_entry
    {n : Type*} [Fintype n] [DecidableEq n]
    (lam : n → ℝ) (hlam : ∀ i, 0 < lam i) (t : ℝ) (a b : n) :
    (zb007Diag lam t * zb007Unit a b * zb007Diag lam (-t)) a b
      = Complex.exp (Complex.I * t * Real.log (lam a / lam b)) := by
  simp [zb007Diag, zb007Unit, Matrix.diagonal, Matrix.mul_apply, Finset.sum_ite_eq',
    zb007_phase_mul (lam a) (lam b) t (hlam a) (hlam b)]

#print axioms zb007_diag_modular_flow_entry

end Zero
