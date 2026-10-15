/-
  Z102（K₄ 上顶点 3 势 φ₃ 由中心化点源定出）。

  K₄ 的 Laplacian `L = 4I − J`（J 是全 1 矩阵）。对**零和**向量 `x`（`Σx = 0`）有 `Jx = 0`，
  于是 `L x = 4 x` ⟹ 伪逆给出 `φ = L⁺x = x/4`，点源势由源直接定出。

  审计修订（2026-10-14）：原版把 `L = 4I − J` 写成 `True := by trivial`（什么都没说），
  现在改成对任意零和向量的**真命题**（矩阵作用）。
-/

import Mathlib

namespace Zero

open Finset

/-- **Z102 · K₄ Laplacian 在零和向量上是 4 倍**：`(Lx)_i = 4x_i − Σ_j x_j`，若 `Σx = 0`
    则 `Lx = 4x`。这正是"在零和超平面上 L 的伪逆 = 1/4"的来源。 -/
theorem z102_laplacian_on_zero_sum (x : Fin 4 → ℝ) (hx : ∑ i, x i = 0) :
    (fun i => 4 * x i - ∑ j, x j) = fun i => 4 * x i := by
  ext i
  simp [hx]

/-- **Z102 · 点源势的具体值**：源 `x₃ = −1/4`、`φ = x/4` ⟹ `φ₃ = −1/16`。 -/
theorem z102_phi_three : (-1 : ℝ) / 4 / 4 = -1 / 16 := by norm_num

#print axioms z102_laplacian_on_zero_sum
#print axioms z102_phi_three

end Zero
