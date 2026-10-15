/-
  Z080（M⁺ 自旋实现恒等式）的机器检查版本。

  Z080 §3.2 的命题：对平衡词 `w ∈ W_L`（`±1` 序列、和为 0），设
  `M⁺(w)` = 前缀和的最大值（含空前缀，即"正内程"），则

      m_j := #{w ∈ W_L : M⁺(w) = j} = C(L, L/2−j) − C(L, L/2−j−1)

  这是经典的**反射原理／ballot 定理**。Mathlib 里**没有**现成的 ballot 定理
  （只有 Dyck 词与 Catalan 数，`Combinatorics/Enumerative/DyckWord.lean`），
  所以本篇分两层：

    · **一般引理**（与 L 无关）：翻转前 k 个符号把词的和改变 `2 × 前缀和`——
      这正是反射双射的代数核心（反射论证 = "首次达到 j+1 的位置翻转之"）；
    · **有限验证**：把 Z080 §3.1 的表格从 L=4 扩到 L=6/8，逐项与公式对齐。

  诚实边界：**一般 L 的完整 ballot 定理仍未证**（需要构造"首次达到"位置上的双射，
  是本队列里最重的一项）。这里的一般引理是它的核心部件，但不等于它；
  L=10 因内核归约规模（1024 个词）未纳入，只做到 L=8。

  记号约定：`true = +1`、`false = −1`。
-/

import Mathlib

namespace Zero

/-! ### 一、计数层：长度 L 的词（`Fin L → Bool`，有 Fintype 可枚举） -/

/-- 符号：`true ↦ +1`、`false ↦ −1`。 -/
def bSgn : Bool → ℤ := fun b => if b then 1 else -1

/-- 词的和（`Σ w`）。 -/
def wSum {L : ℕ} (w : Fin L → Bool) : ℤ := Finset.sum Finset.univ (fun i => bSgn (w i))

/-- 前 `k` 项的部分和。 -/
def wPrefix {L : ℕ} (w : Fin L → Bool) (k : ℕ) : ℤ :=
  Finset.sum (Finset.univ.filter (fun i : Fin L => (i : ℕ) < k)) (fun i => bSgn (w i))

/-- `M⁺`：前缀和的最大值（含 `k = 0`，故恒 ≥ 0）。 -/
def wMplus {L : ℕ} (w : Fin L → Bool) : ℤ :=
  (Finset.range (L + 1)).sup'
    ⟨0, Finset.mem_range.mpr (Nat.lt_succ_of_le (Nat.zero_le _))⟩ (wPrefix w)

/-- `W_L`：长度 `L`、和为 0 的平衡词。 -/
def balancedOf (L : ℕ) : Finset (Fin L → Bool) :=
  Finset.univ.filter (fun w => wSum w = 0)

/-- Z080 的计数公式（`j > L/2` 取 0；`j = L/2` 时第二项按 0 记，即 `C(L,0) = 1`）。 -/
def ballotFormula (L j : ℕ) : ℕ :=
  if L / 2 < j then 0
  else if j = L / 2 then Nat.choose L 0
  else Nat.choose L (L / 2 - j) - Nat.choose L (L / 2 - j - 1)

/-- **Z080 · L=6 的完整 M⁺ 分布**（20 个平衡词）：`m_j` 与公式逐项相符。 -/
theorem z080_ballot_L6 :
    (List.range 4).map (fun j => ((balancedOf 6).filter (fun w => wMplus w = j)).card)
      = (List.range 4).map (ballotFormula 6) := by decide

/-- **Z080 · L=8 的完整 M⁺ 分布**（70 个平衡词）。 -/
theorem z080_ballot_L8 :
    (List.range 5).map (fun j => ((balancedOf 8).filter (fun w => wMplus w = j)).card)
      = (List.range 5).map (ballotFormula 8) := by decide

/-- **Z080 · 重数和 = 中央二项式**（L=6）：M⁺ 的值把 `W_L` 恰好分完。 -/
theorem z080_L6_mass :
    Finset.sum (Finset.range 4) (fun j => ((balancedOf 6).filter (fun w => wMplus w = j)).card)
      = Nat.choose 6 3 := by decide

/-! ### 二、一般引理：翻转前缀（反射双射的代数核心） -/

/-- 符号（`List` 版，与 `bSgn` 同一约定）。 -/
def lSgn : Bool → ℤ := fun b => if b then 1 else -1

/-- 词的和（`List` 版）。 -/
def lSum (w : List Bool) : ℤ := (w.map lSgn).sum

/-- 前 `k` 项的部分和（`List` 版）。 -/
def lPrefix (w : List Bool) (k : ℕ) : ℤ := ((w.take k).map lSgn).sum

/-- **翻转前 `k` 项**：反射论证里"把首次达到 `j+1` 之前的那段翻过来"。 -/
def lFlip (w : List Bool) (k : ℕ) : List Bool := (w.take k).map not ++ w.drop k

/-- 翻转不改变长度。 -/
theorem z080_flip_length (w : List Bool) (k : ℕ) : (lFlip w k).length = w.length := by
  simp [lFlip, List.length_take, List.length_drop]
  omega

/-- **逐个翻转符号使和变号**：`Σ (−w_i) = −Σ w_i`。 -/
theorem z080_map_not_sum (l : List Bool) :
    ((l.map not).map lSgn).sum = -((l.map lSgn).sum) := by
  induction l with
  | nil => rfl
  | cons b t ih =>
    simp only [List.map_cons, List.sum_cons, ih]
    cases b <;> simp [lSgn] <;> ring

/-- **Z080 · 反射的核心恒等式**：`lSum (lFlip w k) = lSum w − 2 · lPrefix w k`。
    即"翻转前 k 项"把和改变 `2 × 前缀和` —— 这正是反射原理里
    "首次达到 `j+1` 的位置翻转后，和变成 `−2(j+1)`" 的那一步。 -/
theorem z080_flip_prefix_sum (w : List Bool) (k : ℕ) :
    lSum (lFlip w k) = lSum w - 2 * lPrefix w k := by
  have hsplit : lSum w = lPrefix w k + ((w.drop k).map lSgn).sum := by
    simp only [lSum, lPrefix, ← List.sum_append, ← List.map_append, List.take_append_drop]
  have hnot := z080_map_not_sum (w.take k)
  have hmain : lSum (lFlip w k) = -lPrefix w k + ((w.drop k).map lSgn).sum := by
    simp only [lFlip, lSum, List.map_append, List.sum_append, lPrefix, hnot]
  rw [hmain, hsplit]
  ring

#print axioms z080_ballot_L6
#print axioms z080_ballot_L8
#print axioms z080_L6_mass
#print axioms z080_flip_length
#print axioms z080_map_not_sum
#print axioms z080_flip_prefix_sum

end Zero
