/-
  Zb007 §2 的矩阵形式：2 级块的模流频率。

  取 `ρ = diag(p,q)`（p,q>0），`U(t) = diag(e^{it log p}, e^{it log q})`（即 `ρ^{it}` 的对角写法），
  `E = [[0,1],[0,0]]`。则
      `U(t) E U(t)^{-1} = e^{it log(p/q)} E`。
  这把「非对角元的模流频率 = log(p/q)」写成 2×2 矩阵恒等式；不需要一般矩阵的 `cpow`。
-/

import Mathlib

namespace Zero

open Matrix

/-- `ρ^{it}` 的 2 级对角写法。 -/
noncomputable def zb007U (p q t : ℝ) : Matrix (Fin 2) (Fin 2) ℂ :=
  !![Complex.exp (Complex.I * t * Real.log p), 0;
     0, Complex.exp (Complex.I * t * Real.log q)]

/-- `U(t)` 的逆（对角矩阵的逆仍对角）。 -/
noncomputable def zb007Uinv (p q t : ℝ) : Matrix (Fin 2) (Fin 2) ℂ :=
  !![Complex.exp (-(Complex.I * t * Real.log p)), 0;
     0, Complex.exp (-(Complex.I * t * Real.log q))]

/-- 非对角矩阵元 `E = |0><1|`。 -/
def zb007E : Matrix (Fin 2) (Fin 2) ℂ :=
  !![0, 1; 0, 0]

/-- **Zb007 · 2 级块的相位因子**：`e^{it log p}·e^{−it log q}=e^{it log(p/q)}`。 -/
theorem zb007_phase_mul (p q t : ℝ) (hp : 0 < p) (hq : 0 < q) :
    Complex.exp (Complex.I * t * Real.log p) *
      Complex.exp (-(Complex.I * t * Real.log q)) =
      Complex.exp (Complex.I * t * Real.log (p / q)) := by
  rw [← Complex.exp_add, Real.log_div (ne_of_gt hp) (ne_of_gt hq)]
  congr 1
  push_cast
  ring

/-- **Zb007 · U 与 U⁻¹**：`U(t) U(t)⁻¹ = 1`。 -/
theorem zb007U_mul_inv (p q t : ℝ) :
    zb007U p q t * zb007Uinv p q t = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [zb007U, zb007Uinv, Matrix.mul_apply, Fin.sum_univ_two, Complex.exp_neg,
      Complex.exp_ne_zero]

/-- **Zb007 · 2 级块的模流**：`U(t) E U(t)⁻¹ = e^{it log(p/q)} E`。 -/
theorem zb007_modular_flow_E (p q t : ℝ) (hp : 0 < p) (hq : 0 < q) :
    zb007U p q t * zb007E * zb007Uinv p q t =
      Complex.exp (Complex.I * t * Real.log (p / q)) • zb007E := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [zb007U, zb007Uinv, zb007E, Matrix.mul_apply, Fin.sum_univ_two,
      zb007_phase_mul p q t hp hq]

#print axioms zb007_phase_mul
#print axioms zb007U_mul_inv
#print axioms zb007_modular_flow_E

end Zero
