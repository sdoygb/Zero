/-
  Z021 §1／§3 的**维数字典裁决**（ZeroGPU/doc/dict_verdict.md）的机器检查版本。

  内容三件（全部是算术，不含物理）：
    · **恒等式**：`C(m,2) − (m−1) = C(m−1,2)`（边上链商 = 独立成对身份数）；
    · **唯一性**：`m ≥ 3`、`D ≥ 2` 时 `C(D,2) = C(m,2) − (m−1) ⟺ D = m−1`，
      故错字典 `D = m` 全表失败（商严格小于 `C(m,2)`）；
    · **四维窗与唯一峰**：`F_D(q) = C(D,2)·q^D`，峰在 `D` 的窗是 `((D−2)/D, (D−1)/(D+1))`；
      `D = 4` 给 `(1/2, 3/5)`，且 `q = 5/9` 落在窗内、**不**落在三维窗内；
      `q = 5/9` 时全局唯一峰为 `D = 4`。

  数值互证：`ZeroGPU/src/dict_verdict.py`（13 条断言）。

  工程注记：`Nat.choose_succ_succ n 1` 的右端展开成 `n.choose (Nat.succ 1)`，
  omega 会把 `n.choose 2` 与 `n.choose (Nat.succ 1)` 当成两个**不同**的原子 —— 所以这里
  先把 Pascal 递推写成显式的 `have`（用字面 `2` 陈述），再交给 omega。
-/

import Mathlib

namespace Zero

/-! ## 一、边上链商的恒等式与唯一性 -/

/-- **Z021 · 字典恒等式**：`C(m,2) − (m−1) = C(m−1,2)`。
    来源：`K_m` 的 `C(m,2)` 条边减去顶点相位的 `m−1` 维规范自由度，剩 `C(m−1,2)` 个独立成对身份。
    证明：Pascal 递推 `C(m,2) = C(m−1,1) + C(m−1,2)`，而 `C(m−1,1) = m−1`。 -/
theorem z021_dict_identity (m : ℕ) : m.choose 2 - (m - 1) = (m - 1).choose 2 := by
  rcases Nat.eq_zero_or_pos m with h | h
  · subst h; decide
  · obtain ⟨n, rfl⟩ : ∃ n, m = n + 1 := ⟨m - 1, by omega⟩
    have hsub : n + 1 - 1 = n := by omega
    have hpascal : (n + 1).choose 2 = n.choose 1 + n.choose 2 := Nat.choose_succ_succ n 1
    have hone : n.choose 1 = n := Nat.choose_one_right n
    rw [hsub, hpascal, hone]
    omega

/-- **Z021 · `C(·,2)` 在 `≥ 2` 上严格单调**（唯一性裁决的算术引理；
    Mathlib 的 `Nat.choose_le_choose` 只给单调不减，这里要严格）。 -/
theorem z021_choose_two_lt {a b : ℕ} (ha : 2 ≤ a) (hab : a < b) :
    a.choose 2 < b.choose 2 := by
  have step : ∀ c : ℕ, 2 ≤ c → c.choose 2 < (c + 1).choose 2 := by
    intro c hc
    have hpascal : (c + 1).choose 2 = c.choose 1 + c.choose 2 := Nat.choose_succ_succ c 1
    have hone : c.choose 1 = c := Nat.choose_one_right c
    rw [hpascal, hone]
    omega
  obtain ⟨k, hk⟩ := Nat.exists_eq_add_of_le (Nat.succ_le_of_lt hab)
  have hmono : ∀ k : ℕ, a.choose 2 < (a + 1 + k).choose 2 := by
    intro k
    induction k with
    | zero => simpa using step a ha
    | succ k ih => exact lt_trans ih (step (a + 1 + k) (by omega))
  rw [hk]
  exact hmono k

/-- **Z021 · 字典唯一性**：`m ≥ 3`、`D ≥ 2` 时 `C(D,2) = C(m,2) − (m−1) ⟺ D = m−1`。
    这就是 `dict_verdict.md` 的"代数唯一解"——不是偏好，是方程的解唯一。 -/
theorem z021_dict_unique {m D : ℕ} (hm : 3 ≤ m) (hD : 2 ≤ D) :
    D.choose 2 = m.choose 2 - (m - 1) ↔ D = m - 1 := by
  rw [z021_dict_identity]
  constructor
  · intro h
    rcases lt_trichotomy D (m - 1) with hlt | heq | hgt
    · exact absurd h (ne_of_lt (z021_choose_two_lt hD hlt))
    · exact heq
    · exact absurd h.symm (ne_of_lt (z021_choose_two_lt (by omega) hgt))
  · intro h; rw [h]

/-- **Z021 · 错字典 `D = m` 全表失败**（`m ≥ 2`）：商严格小于 `C(m,2)`，故 `D = m` 不是解。 -/
theorem z021_dict_m_fails (m : ℕ) (hm : 2 ≤ m) : m.choose 2 - (m - 1) < m.choose 2 :=
  Nat.sub_lt (Nat.choose_pos (by omega)) (by omega)

-- 审计修订（2026-10-14 夜）：删去 `z021_dict_consequences`（纯数字合取，且已被
-- `z021_dict_identity` + `z021_dict_unique` 一般定理覆盖）。

/-! ## 二、四维窗与唯一峰：`F_D(q) = C(D,2)·q^D` -/

/-- 峰函数 `F_D(q) = C(D,2)·q^D`（`q` 取有理数，与 `dict_verdict.py` 同口径）。 -/
def z021Peak (D : ℕ) (q : ℚ) : ℚ := (D.choose 2 : ℚ) * q ^ D

/-- `2 · C(n,2) = n(n−1)`（清分母用的算术恒等式）。 -/
theorem z021_two_mul_choose_two (n : ℕ) : 2 * n.choose 2 = n * (n - 1) := by
  rw [Nat.choose_two_right]
  exact Nat.mul_div_cancel' (Nat.even_mul_pred_self n).two_dvd

/-- 上式在 `ℚ` 里的形式（`n ≥ 1`）。 -/
theorem z021_two_mul_choose_two_rat {n : ℕ} (hn : 1 ≤ n) :
    (2 : ℚ) * (n.choose 2 : ℚ) = (n : ℚ) * ((n : ℚ) - 1) := by
  have h : ((2 * n.choose 2 : ℕ) : ℚ) = ((n * (n - 1) : ℕ) : ℚ) := by
    rw [z021_two_mul_choose_two]
  push_cast at h
  rw [Nat.cast_sub hn, Nat.cast_one] at h
  linarith

/-- **Z021 · 四维窗**：`q ∈ (1/2, 3/5)` ⟹ `D = 4` 是峰（`F₃ < F₄` 且 `F₅ < F₄`）。
    窗的两端就是峰条件的两个比值：`F₄/F₃ = 2q > 1`、`F₅/F₄ = (5/3)q < 1`。 -/
theorem z021_window_four (q : ℚ) (h1 : 1/2 < q) (h2 : q < 3/5) :
    z021Peak 3 q < z021Peak 4 q ∧ z021Peak 5 q < z021Peak 4 q := by
  have hq : 0 < q := by linarith
  constructor
  · rw [z021Peak, z021Peak]
    norm_num [Nat.choose_two_right]
    nlinarith [pow_pos hq 3]
  · rw [z021Peak, z021Peak]
    norm_num [Nat.choose_two_right]
    nlinarith [pow_pos hq 4]

/-- **Z021 · 三维窗**：`q ∈ (1/3, 1/2)` ⟹ `D = 3` 是峰。 -/
theorem z021_window_three (q : ℚ) (h1 : 1/3 < q) (h2 : q < 1/2) :
    z021Peak 2 q < z021Peak 3 q ∧ z021Peak 4 q < z021Peak 3 q := by
  have hq : 0 < q := by linarith
  constructor
  · rw [z021Peak, z021Peak]
    norm_num [Nat.choose_two_right]
    nlinarith [pow_pos hq 2]
  · rw [z021Peak, z021Peak]
    norm_num [Nat.choose_two_right]
    nlinarith [pow_pos hq 3]

/-- **Z021 · 四维窗的充要刻画**：`q > 0` 时 `D = 4` 是峰 ⟺ `q ∈ (1/2, 3/5)`。
    这比"5/9 落在窗内"强：窗就是峰条件的**解集**（5/9 ∈ (1/2,3/5) 只是它的一个点）。 -/
theorem z021_peak_four_iff (q : ℚ) (hq : 0 < q) :
    (z021Peak 3 q < z021Peak 4 q ∧ z021Peak 5 q < z021Peak 4 q) ↔ (1/2 < q ∧ q < 3/5) := by
  have h3 : z021Peak 3 q = 3 * q ^ 3 := by rw [z021Peak]; norm_num [Nat.choose_two_right]
  have h4 : z021Peak 4 q = 6 * q ^ 4 := by rw [z021Peak]; norm_num [Nat.choose_two_right]
  have h5 : z021Peak 5 q = 10 * q ^ 5 := by rw [z021Peak]; norm_num [Nat.choose_two_right]
  constructor
  · rintro ⟨h34, h54⟩
    rw [h3, h4] at h34
    rw [h5, h4] at h54
    have hq3 : (0 : ℚ) < q ^ 3 := pow_pos hq 3
    have hq4 : (0 : ℚ) < q ^ 4 := pow_pos hq 4
    constructor <;> nlinarith
  · rintro ⟨h1, h2⟩
    exact z021_window_four q h1 h2

/-- **Z021 · 三维窗的充要刻画**：`q > 0` 时 `D = 3` 是峰 ⟺ `q ∈ (1/3, 1/2)`。
    推论：`q = 5/9` **不**在其中 —— 旧说法"峰是 D = 3"由此被排除。 -/
theorem z021_peak_three_iff (q : ℚ) (hq : 0 < q) :
    (z021Peak 2 q < z021Peak 3 q ∧ z021Peak 4 q < z021Peak 3 q) ↔ (1/3 < q ∧ q < 1/2) := by
  have h2 : z021Peak 2 q = 1 * q ^ 2 := by rw [z021Peak]; norm_num [Nat.choose_two_right]
  have h3 : z021Peak 3 q = 3 * q ^ 3 := by rw [z021Peak]; norm_num [Nat.choose_two_right]
  have h4 : z021Peak 4 q = 6 * q ^ 4 := by rw [z021Peak]; norm_num [Nat.choose_two_right]
  constructor
  · rintro ⟨h23, h43⟩
    rw [h2, h3] at h23
    rw [h4, h3] at h43
    have hq2 : (0 : ℚ) < q ^ 2 := pow_pos hq 2
    have hq3 : (0 : ℚ) < q ^ 3 := pow_pos hq 3
    constructor <;> nlinarith
  · rintro ⟨h1, h2⟩
    exact z021_window_three q h1 h2

/-- **Z021 · `q = 5/9` 时 `F` 在 `D ≥ 4` 上严格递减**：核心是 `5(D+1) < 9(D−1)`（即 `D > 3.5`）。 -/
theorem z021_five_ninths_decay {D : ℕ} (hD : 4 ≤ D) :
    z021Peak (D + 1) (5/9) < z021Peak D (5/9) := by
  have hA := z021_two_mul_choose_two_rat (n := D + 1) (by omega)
  have hB := z021_two_mul_choose_two_rat (n := D) (by omega)
  push_cast at hA hB
  have hDq : (4 : ℚ) ≤ (D : ℚ) := by exact_mod_cast hD
  -- 清分母后的核心不等式：`5(D+1) < 9(D−1)`，由 `D ≥ 4` 给出
  have hprod : (0 : ℚ) < (D : ℚ) * (2 * (D : ℚ) - 7) :=
    mul_pos (by linarith) (by linarith)
  have hkey : 5 * ((D : ℚ) + 1) * (D : ℚ) < 9 * (D : ℚ) * ((D : ℚ) - 1) := by
    nlinarith [hprod]
  have h5A : 10 * ((D + 1).choose 2 : ℚ) = 5 * ((D : ℚ) + 1) * (D : ℚ) := by linarith
  have h9B : 18 * (D.choose 2 : ℚ) = 9 * (D : ℚ) * ((D : ℚ) - 1) := by linarith
  have hker : ((D + 1).choose 2 : ℚ) * (5/9) < (D.choose 2 : ℚ) := by
    nlinarith [h5A, h9B, hkey]
  have hqpos : (0 : ℚ) < (5/9) ^ D := by positivity
  calc z021Peak (D + 1) (5/9)
      = ((D + 1).choose 2 : ℚ) * ((5/9) ^ D * (5/9)) := by rw [z021Peak, pow_succ]
    _ = (((D + 1).choose 2 : ℚ) * (5/9)) * (5/9) ^ D := by ring
    _ < (D.choose 2 : ℚ) * (5/9) ^ D := mul_lt_mul_of_pos_right hker hqpos
    _ = z021Peak D (5/9) := by rw [z021Peak]

/-- **Z021 · 唯一峰**：`q = 5/9` 时 `F_D` 的全局唯一峰在 `D = 4`（`D ≥ 2`）。 -/
theorem z021_five_ninths_peak (D : ℕ) (hD : 2 ≤ D) (hne : D ≠ 4) :
    z021Peak D (5/9) < z021Peak 4 (5/9) := by
  rcases lt_trichotomy D 4 with hlt | heq | hgt
  · interval_cases D <;> norm_num [z021Peak, Nat.choose_two_right]
  · exact absurd heq hne
  · obtain ⟨k, hk⟩ := Nat.exists_eq_add_of_le (show 5 ≤ D by omega)
    have hmono : ∀ k : ℕ, z021Peak (5 + k) (5/9) < z021Peak 4 (5/9) := by
      intro k
      induction k with
      | zero => norm_num [z021Peak, Nat.choose_two_right]
      | succ k ih => exact lt_trans (z021_five_ninths_decay (D := 5 + k) (by omega)) ih
    rw [hk]
    exact hmono k

#print axioms z021_dict_identity
#print axioms z021_choose_two_lt
#print axioms z021_dict_unique
#print axioms z021_dict_m_fails
#print axioms z021_two_mul_choose_two
#print axioms z021_two_mul_choose_two_rat
#print axioms z021_window_four
#print axioms z021_peak_four_iff
#print axioms z021_peak_three_iff
#print axioms z021_window_three
#print axioms z021_five_ninths_decay
#print axioms z021_five_ninths_peak

end Zero