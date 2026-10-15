/-
  Z024（离散连续性方程 ⟹ 总荷守恒）的机器检查版本。

  审计修订（2026-10-14）：原版**从未编译通过**（求和交换的证明留下未解决目标）。
  这里重写：先把 `Σ_i Σ_j` 交换成 `Σ_j Σ_i`，再用"每列和为 0"把每一项归零。
-/

import Mathlib

namespace Zero

open Finset

/-- **Z024 · 总荷守恒**：若关联矩阵 `B` 的**每列和为 0**（离散连续性方程），
    则总荷 `Σ_i Σ_j B(i,j) f(j) = 0` —— 即荷的净产生量为零。
    证明只用求和交换与列和为零，不需要作用量或 Noether 定理。 -/
theorem z024_charge_conservation {n m : ℕ} (B : Fin n → Fin m → ℤ) (f : Fin m → ℤ)
    (hcol : ∀ j : Fin m, ∑ i : Fin n, B i j = 0) :
    ∑ i : Fin n, ∑ j : Fin m, B i j * f j = 0 := by
  rw [Finset.sum_comm]
  refine Finset.sum_eq_zero fun j _ => ?_
  rw [← Finset.sum_mul, hcol j, zero_mul]

#print axioms z024_charge_conservation

end Zero
