/-
  Zb048 §1–§3 的机器检查版本：K₄ 的边、三个显式完美匹配、自同构群阶。

  形式化口径：
    · 点集 `Fin 4`；边是 `Finset (Fin 4)` 的二元子集；
    · 三个完美匹配 `M1, M2, M3` 显式给出，并逐条核验它们是完美匹配；
    · `|Aut(K₄)| = |S₄| = 24`、`|S₃| = 6`。

  未覆盖（下一批）：
    · 「恰有 3 个完美匹配」的**枚举版**：对 16 个边子集做 `decide` 会触发内核递归深度；
      改用 `native_decide` 会引入未登记的编译信任公理（R2 拒绝）。
      正确做法是写「含 0 的边只有 3 种，其余两点被迫配对」的结构性证明，或把枚举
      编码成更小的显式表；本轮先给显式见证，不冒充「恰有 3 个」。
    · `S₄ → S₃` 同态与核 `V₄` 的显式构造。
-/

import Mathlib

namespace Zero

/-- K₄ 的一条边：二元子集。 -/
def zb048Edge (a b : Fin 4) : Finset (Fin 4) := {a, b}

/-- 完美匹配：两条互不相交、覆盖全部四点的边。 -/
def zb048IsPerfectMatching (M : Finset (Finset (Fin 4))) : Prop :=
  M.card = 2 ∧
    (∀ e ∈ M, e.card = 2) ∧
    (∀ e1 ∈ M, ∀ e2 ∈ M, e1 ≠ e2 → Disjoint e1 e2) ∧
    (∀ v : Fin 4, ∃ e ∈ M, v ∈ e)

instance zb048Decidable (M : Finset (Finset (Fin 4))) : Decidable (zb048IsPerfectMatching M) := by
  unfold zb048IsPerfectMatching
  infer_instance

/-- **Zb048 · 边数**：K₄ 有 `C(4,2) = 6` 条边。 -/
theorem zb048_edge_count : Nat.choose 4 2 = 6 := by
  decide

/-- 三个完美匹配的显式列表。 -/
def zb048M1 : Finset (Finset (Fin 4)) := {zb048Edge 0 1, zb048Edge 2 3}
def zb048M2 : Finset (Finset (Fin 4)) := {zb048Edge 0 2, zb048Edge 1 3}
def zb048M3 : Finset (Finset (Fin 4)) := {zb048Edge 0 3, zb048Edge 1 2}

/-- `M1` 是完美匹配。 -/
theorem zb048_M1_perfect : zb048IsPerfectMatching zb048M1 := by
  decide

/-- `M2` 是完美匹配。 -/
theorem zb048_M2_perfect : zb048IsPerfectMatching zb048M2 := by
  decide

/-- `M3` 是完美匹配。 -/
theorem zb048_M3_perfect : zb048IsPerfectMatching zb048M3 := by
  decide

/-- **Zb048 · 自同构群阶**：`|Aut(K₄)| = |S₄| = 24`。 -/
theorem zb048_aut_card : Fintype.card (Equiv.Perm (Fin 4)) = 24 := by
  simp [Fintype.card_perm, Nat.factorial]

/-- **Zb048 · 三个匹配的置换群阶**：`|S₃| = 6`（商群 `S₄/V₄` 的阶）。 -/
theorem zb048_s3_card : Fintype.card (Equiv.Perm (Fin 3)) = 6 := by
  simp [Fintype.card_perm, Nat.factorial]

#print axioms zb048_edge_count
#print axioms zb048_M1_perfect
#print axioms zb048_M2_perfect
#print axioms zb048_M3_perfect
#print axioms zb048_aut_card
#print axioms zb048_s3_card

end Zero
