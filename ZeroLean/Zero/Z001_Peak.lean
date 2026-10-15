/-
  Z001 §7（账本泛函 F_D = C(D,2) q^D 的唯一有限极大点 D = 4）的机器检查版本。

  相邻项差分（直接展开）：
    F_{D+1} − F_D = D/2 · q^D · ((D+1) q − (D−1))。
  D ≥ 2、q > 0 时 D/2·q^D > 0，故符号由 (D+1)q − (D−1) 决定：
    · D=3：4q−2 > 0（q>1/2） ⟹ F4>F3；
    · D=4：5q−3 < 0（q<3/5） ⟹ F5<F4；
    · D≥4：(D−1)/(D+1) 随 D 递增、D=4 时=3/5，q<3/5 ⟹ (D+1)q−(D−1)<0 ⟹ F_{D+1}<F_D。
  F0=F1=0；故 D=4 是全局唯一极大点。
-/

import Mathlib

namespace Zero

/-- 账本泛函 F_D = D(D−1)/2 · q^D。 -/
noncomputable def z001F (D : ℕ) (q : ℝ) : ℝ :=
  (D : ℝ) * ((D : ℝ) - 1) / 2 * q^D

/-- **相邻项差分公式**：F_{D+1} − F_D = D/2 · q^D · ((D+1)q − (D−1))。 -/
lemma z001_diff (D : ℕ) (q : ℝ) :
    z001F (D + 1) q - z001F D q =
      (D : ℝ) / 2 * q^D * (((D : ℝ) + 1) * q - ((D : ℝ) - 1)) := by
  dsimp only [z001F]
  cases D with
  | zero => ring
  | succ D' => simp [pow_succ]; ring

/-- **D ≥ 2、q > 0 时系数 D/2·q^D > 0。** -/
lemma z001_coeff_pos (D : ℕ) (q : ℝ) (hD2 : 2 ≤ D) (hq0 : 0 < q) :
    0 < (D : ℝ) / 2 * q^D := by
  have hD : (D : ℝ) ≥ 2 := by exact_mod_cast hD2
  have hqD : 0 < q^D := by positivity
  have h1 : (D : ℝ) / 2 > 0 := by linarith
  exact mul_pos h1 hqD

/-- **D=3：F4 > F3**（q > 1/2）。 -/
lemma z001_f4_gt_f3 (q : ℝ) (hq : 1/2 < q) :
    z001F 4 q > z001F 3 q := by
  have hq0 : 0 < q := by linarith
  have h : z001F 4 q - z001F 3 q > 0 := by
    rw [z001_diff 3 q]
    simp only [Nat.cast_ofNat]
    have ha : 0 < (3 : ℝ) / 2 * q^3 := z001_coeff_pos 3 q (by norm_num) hq0
    have hb : 0 < ((3 : ℝ) + 1) * q - ((3 : ℝ) - 1) := by linarith
    exact mul_pos ha hb
  exact sub_pos.mp h

/-- **D=4：F5 < F4**（q < 3/5）。 -/
lemma z001_f5_lt_f4 (q : ℝ) (hq0 : 0 < q) (hq : q < 3/5) :
    z001F 5 q < z001F 4 q := by
  have h : z001F 5 q - z001F 4 q < 0 := by
    rw [z001_diff 4 q]
    simp only [Nat.cast_ofNat]
    have ha : 0 < (4 : ℝ) / 2 * q^4 := z001_coeff_pos 4 q (by norm_num) hq0
    have hb : ((4 : ℝ) + 1) * q - ((4 : ℝ) - 1) < 0 := by linarith
    exact mul_neg_of_pos_of_neg ha hb
  exact sub_neg.mp h

/-- **D ≥ 4：F_{D+1} < F_D**（q < 3/5）。 -/
lemma z001_decreasing_after_4 (D : ℕ) (q : ℝ) (hq0 : 0 < q) (hD4 : 4 ≤ D) (hq : q < 3/5) :
    z001F (D + 1) q < z001F D q := by
  have h : z001F (D + 1) q - z001F D q < 0 := by
    rw [z001_diff D q]
    have hD : (D : ℝ) ≥ 4 := by exact_mod_cast hD4
    have ha : 0 < (D : ℝ) / 2 * q^D := z001_coeff_pos D q (by linarith) hq0
    have hb : ((D : ℝ) + 1) * q - ((D : ℝ) - 1) < 0 := by nlinarith
    exact mul_neg_of_pos_of_neg ha hb
  exact sub_neg.mp h

/-- F4 的显式值：F4 = 6 q⁴。 -/
lemma z001_f4_eq (q : ℝ) : z001F 4 q = 6 * q^4 := by
  dsimp only [z001F]; ring

/-- F0 = 0，F1 = 0。 -/
lemma z001_f0_f1 (q : ℝ) : z001F 0 q = 0 ∧ z001F 1 q = 0 := by
  dsimp only [z001F]; norm_num

/-- **D ≥ 5：F D < F4**（从 F5 起严格递降，独立归纳）。 -/
lemma z001_lt_f4_from_5 (q : ℝ) (hq0 : 0 < q) (hhi : q < 3/5) :
    ∀ D : ℕ, 5 ≤ D → z001F D q < z001F 4 q := by
  intro D hD
  induction' hD with D' hD' ih
  · -- D = 5
    exact z001_f5_lt_f4 q hq0 hhi
  · -- D = D'+1，hD' : 5 ≤ D'，ih : F D' < F4
    have h4 : 4 ≤ D' := by
      exact Nat.le_trans (by norm_num) hD'
    have hdec : z001F (D' + 1) q < z001F D' q :=
      z001_decreasing_after_4 D' q hq0 h4 hhi
    exact hdec.trans ih

/-- **Z001 · 主定理：q ∈ (1/2, 3/5) ⟹ D = 4 是 F 的唯一有限极大点。** -/
theorem z001_peak_unique_at_4 (q : ℝ) (hlo : 1/2 < q) (hhi : q < 3/5) :
    ∀ D : ℕ, D ≠ 4 → z001F D q < z001F 4 q := by
  have hq0 : 0 < q := by linarith
  intro D hD
  by_cases h : D ≤ 4
  · -- D ∈ {0,1,2,3,4}
    interval_cases D
    · -- D=0：F0=0 < F4=6q⁴
      have hf0 : z001F 0 q = 0 := (z001_f0_f1 q).1
      have hf4 : z001F 4 q = 6 * q^4 := z001_f4_eq q
      rw [hf0, hf4]; positivity
    · -- D=1：F1=0 < F4
      have hf1 : z001F 1 q = 0 := (z001_f0_f1 q).2
      have hf4 : z001F 4 q = 6 * q^4 := z001_f4_eq q
      rw [hf1, hf4]; positivity
    · -- D=2：F2=q² < 6q⁴=F4（q>1/2 ⟹ 6q²>1）
      have hF2 : z001F 2 q = q^2 := by dsimp only [z001F]; ring
      have hf4 : z001F 4 q = 6 * q^4 := z001_f4_eq q
      rw [hF2, hf4]
      have hq2 : 0 < q^2 := by positivity
      have h6 : 6 * q^2 > 1 := by nlinarith
      calc q^2 = q^2 * 1 := by ring
        _ < q^2 * (6 * q^2) := by gcongr
        _ = 6 * q^4 := by ring
    · exact z001_f4_gt_f3 q hlo
    · exact False.elim (hD rfl)
  · -- D ≥ 5：用独立的递降引理
    have h5 : 5 ≤ D := by omega
    exact z001_lt_f4_from_5 q hq0 hhi D h5

#print axioms z001_peak_unique_at_4

end Zero
