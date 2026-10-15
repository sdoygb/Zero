/-
  Z066（对 Z024 §3.3 的核查）的机器检查版本：C₄ 上组态集 `|R| = 19` 的**两种自洽口径**，
  以及 Z024 原句「`|xᵢ| ≤ 2` 且 `Σx = 0` ⟹ `|R| = 19`」的**否证**（实为 85）。

  数值互证（本轮同口径复算，`ZeroGPU` 侧枚举）：
    · `{Σx = 0, |xᵢ| ≤ 1}` 的态数 = **19**；
    · `{Σx = 0, max x − min x ≤ 2}` 的态数 = **19**；
    · `{Σx = 0, |xᵢ| ≤ 2}` 的态数 = **85**（故 Z024 那句按字面错）；
    · 能量 `E(x) = Σ_edges (x_u − x_v)²` 在 19 态上的分布 = `{0:1, 4:4, 6:8, 8:4, 16:2}`（**非常数**）。

  口径说明：`E` 必须取"逐边差平方和"（等价于 `xᵀ(2I − A)x`，C₄ 是 2-正则）。
  若误取 `xᵀ(4I − A)x` 会得到 `{0:1, 8:4, 10:8, 16:4, 24:2}` —— 与篇内分布不符，故口径唯一。
-/

import Mathlib

namespace Zero

/-- 用 `Fin 3` 编码 `{−1, 0, 1}`。 -/
def z066Int1 {n : ℕ} (v : Fin n → Fin 3) : Fin n → ℤ := fun i => (v i : ℤ) - 1

/-- 用 `Fin 5` 编码 `{−2, …, 2}`。 -/
def z066Int2 {n : ℕ} (v : Fin n → Fin 5) : Fin n → ℤ := fun i => (v i : ℤ) - 2

/-- **口径一**：`|xᵢ| ≤ 1` 且 `Σx = 0` 的 C₄ 组态集。 -/
def z066R1 : Finset (Fin 4 → Fin 3) :=
  Finset.univ.filter (fun v => ∑ i, z066Int1 v i = 0)

/-- **口径二**：格位窗宽 `max − min ≤ 2` 且 `Σx = 0`（用"任意两点差 ≤ 2"表达，避免 `sup/inf`）。 -/
def z066R2 : Finset (Fin 4 → Fin 5) :=
  Finset.univ.filter (fun v => (∑ i, z066Int2 v i = 0) ∧
    ∀ i j, z066Int2 v i - z066Int2 v j ≤ 2)

/-- **Z024 §3.3 的朴素读法**：`|xᵢ| ≤ 2` 且 `Σx = 0`（不设窗宽）的态集。 -/
def z066R3 : Finset (Fin 4 → Fin 5) :=
  Finset.univ.filter (fun v => ∑ i, z066Int2 v i = 0)

/-- **Z066 · 口径一给 19**：`|{Σx = 0, |xᵢ| ≤ 1}| = 19`。 -/
theorem z066_card_R1 : z066R1.card = 19 := by decide

/-- **Z066 · 口径二给 19**：`|{Σx = 0, 窗宽 ≤ 2}| = 19`。 -/
theorem z066_card_R2 : z066R2.card = 19 := by decide

/-- **Z066 · 否证 Z024 的朴素读法**：`|{Σx = 0, |xᵢ| ≤ 2}| = 85 ≠ 19`。 -/
theorem z066_naive_reading_gives_85 : z066R3.card = 85 ∧ z066R3.card ≠ 19 := by decide

/-- C₄ 上的能量：逐边差平方和 `E(x) = Σ_edges (x_u − x_v)²`。 -/
def z066Energy (x : Fin 4 → ℤ) : ℤ :=
  (x 0 - x 1) ^ 2 + (x 1 - x 2) ^ 2 + (x 2 - x 3) ^ 2 + (x 3 - x 0) ^ 2

/-- **Z066 · 能量的值集**：在 19 态上 `E` 取值恰为 `{0, 4, 6, 8, 16}`。 -/
theorem z066_energy_values :
    (z066R1.image (fun v => z066Energy (z066Int1 v))) = ({0, 4, 6, 8, 16} : Finset ℤ) := by
  decide

/-- **Z066 · 能量非常数**（这就是"Z024 的最终结论仍成立、但理由须改"的那一步：
    能量在组态集上不是常值，故不能用"能量守恒 ⟹ 唯一态"来论证）。 -/
theorem z066_energy_not_constant :
    ∃ v ∈ z066R1, ∃ w ∈ z066R1, z066Energy (z066Int1 v) ≠ z066Energy (z066Int1 w) := by
  decide

/-- **Z066 · 能量分布（重数）**：`{0:1, 4:4, 6:8, 8:4, 16:2}`。 -/
theorem z066_energy_multiplicities :
    (z066R1.filter (fun v => z066Energy (z066Int1 v) = 0)).card = 1 ∧
    (z066R1.filter (fun v => z066Energy (z066Int1 v) = 4)).card = 4 ∧
    (z066R1.filter (fun v => z066Energy (z066Int1 v) = 6)).card = 8 ∧
    (z066R1.filter (fun v => z066Energy (z066Int1 v) = 8)).card = 4 ∧
    (z066R1.filter (fun v => z066Energy (z066Int1 v) = 16)).card = 2 := by
  decide

/-- **Z066 · 两种口径自洽**：两种读法的态数与能量值集都一致（"两个解释都好"）。 -/
theorem z066_two_readings_agree :
    z066R1.card = z066R2.card ∧
      (z066R1.image (fun v => z066Energy (z066Int1 v)))
        = (z066R2.image (fun v => z066Energy (z066Int2 v))) := by
  decide

#print axioms z066_card_R1
#print axioms z066_card_R2
#print axioms z066_naive_reading_gives_85
#print axioms z066_energy_values
#print axioms z066_energy_not_constant
#print axioms z066_energy_multiplicities
#print axioms z066_two_readings_agree

end Zero