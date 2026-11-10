/-
  Zb010 §2 的机器检查版本：支撑锥（最近邻传输）。

  正文定理 2.2：`n` 步最近邻传输后，支撑严格落在 `[-n, n]` 内。
  一维模型：每一步 `+1` 或 `−1`（`Bool` 编码）；`n` 步后的位移是
  `Σ step_i`，每条 `|step_i| = 1`，故 `|Σ step_i| ≤ n`。

  未覆盖：图上的支撑锥（需要把 `Walk` 编码进来）、`O(a²)` Taylor 余项
  （需要 Mathlib 的渐近展开），以及锥边概率的测度论版本（这里只给计数核）。
-/

import Mathlib

namespace Zero

/-- 一步的位移：`true ↦ +1`、`false ↦ −1`。 -/
def zb010Step (b : Bool) : ℤ := if b then 1 else -1

/-- `n` 步路径的净位移。 -/
def zb010WalkSum {n : ℕ} (w : Fin n → Bool) : ℤ :=
  ∑ i : Fin n, zb010Step (w i)

/-- 每一步的绝对值都是 1。 -/
theorem zb010_step_abs (b : Bool) : |zb010Step b| = 1 := by
  cases b <;> rfl

/-- **Zb010 · 支撑锥**：`n` 步最近邻路径的净位移必在 `[-n, n]` 内。 -/
theorem zb010_support {n : ℕ} (w : Fin n → Bool) :
    -(n : ℤ) ≤ zb010WalkSum w ∧ zb010WalkSum w ≤ (n : ℤ) := by
  have h : |zb010WalkSum w| ≤ (n : ℤ) := by
    calc
      |zb010WalkSum w| ≤ ∑ i : Fin n, |zb010Step (w i)| :=
        Finset.abs_sum_le_sum_abs (fun i => zb010Step (w i)) Finset.univ
      _ = ∑ _ : Fin n, (1 : ℤ) := by
        exact Finset.sum_congr rfl (fun i _ => zb010_step_abs (w i))
      _ = (n : ℤ) := by simp
  exact abs_le.mp h

/-- **Zb010 · 路径计数**：`n` 步的 `±1` 路径共 `2^n` 条（锥边概率 `2^{-n}` 的计数核）。 -/
theorem zb010_paths_count (n : ℕ) : Fintype.card (Fin n → Bool) = 2 ^ n := by
  simp

#print axioms zb010_step_abs
#print axioms zb010_support
#print axioms zb010_paths_count

end Zero
