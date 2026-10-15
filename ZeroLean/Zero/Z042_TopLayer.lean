/-
  Z042（顶层恰两个词）的机器检查版本：**一般偶长度 L**。

  层高与 Z001/Z047/Z079 一致：`h(w) = max_{k≤L} |前缀和_k|`，壳 `wordShell L h`
  = 平衡（Σ = 0）且层高恰为 h 的词集。平衡词的层高最大为 `L/2`，而**顶层
  （h = L/2）恰有两个词**：「先全上后全下」及其取负。数值侧 B_H = 2 见
  ZeroGPU/src/z047_layer_parity.py（Z047 的 DP 计数）。

  证明骨架（全部 kernel，无枚举、无 sorry）：

  ① 步长恒为 ±1，故 `|前缀和_k| ≤ min(k, L−k)`（后半段的和 = −前缀，因平衡）；
  ② 若某步前缀达 `±L/2`：由 ① 得 `L/2 ≤ k` 且 `L/2 ≤ L−k`，配合 `L` 偶即 `k = L/2`；
     再用「k 项各 ≤ 1 且和 = k ⟹ 每项 = 1」得前一半全 +1、后一半全 −1；
  ③ 层高 `= L/2` 意味着某步达到 `±L/2`（`max'` 取到值），于是 w 只能有两种。
  ④ 反过来这两个词确实在壳里（`Σ = 0` 靠计数，层高达 `L/2` 靠 `prefix_{L/2} = L/2`），
     故壳的基数恰为 2。
-/

import Mathlib
import Zero.Z047_Involution

namespace Zero

open Finset

/-! ### 0. 目标词 -/

/-- "先全上后全下"的词（长度 L 偶） -/
def z042UpThenDown (L : ℕ) : Fin L → Bool := fun i => decide (i.val < L / 2)

/-! ### 1. 步长的上下界与等号 -/

/-- 步长 ≤ 1。 -/
lemma z042_step_le_one (b : Bool) : wordStep b ≤ 1 := by
  cases b <;> norm_num [wordStep]

/-- 步长 ≥ −1。 -/
lemma z042_neg_one_le_step (b : Bool) : -1 ≤ wordStep b := by
  cases b <;> norm_num [wordStep]

/-- `wordStep b = 1 ↔ b = true`。 -/
lemma z042_step_eq_one_iff {b : Bool} : wordStep b = 1 ↔ b = true := by
  cases b <;> simp [wordStep]

/-- `wordStep b = −1 ↔ b = false`。 -/
lemma z042_step_eq_neg_one_iff {b : Bool} : wordStep b = -1 ↔ b = false := by
  cases b <;> simp [wordStep]

/-! ### 2. 前缀和的整体界 -/

/-- 前缀和上界：`prefix_k ≤ k`（每步 ≤ 1）。 -/
lemma z042_prefix_le_self {L : ℕ} (w : Fin L → Bool) (k : ℕ) : wordPrefix w k ≤ (k : ℤ) := by
  unfold wordPrefix
  calc ∑ i ∈ range k, (if h : i < L then wordStep (w ⟨i, h⟩) else 0)
      ≤ ∑ _i ∈ range k, (1 : ℤ) := by
        apply Finset.sum_le_sum
        intro i _
        by_cases h : i < L
        · rw [dite_eq_left h]; exact z042_step_le_one _
        · rw [dite_eq_right h]; norm_num
    _ = (k : ℤ) := by simp

/-- 前缀和下界：`−k ≤ prefix_k`。 -/
lemma z042_neg_le_prefix_self {L : ℕ} (w : Fin L → Bool) (k : ℕ) :
    -((k : ℕ) : ℤ) ≤ wordPrefix w k := by
  unfold wordPrefix
  calc -((k : ℕ) : ℤ) = ∑ _i ∈ range k, (-1 : ℤ) := by simp
    _ ≤ ∑ i ∈ range k, (if h : i < L then wordStep (w ⟨i, h⟩) else 0) := by
        apply Finset.sum_le_sum
        intro i _
        by_cases h : i < L
        · rw [dite_eq_left h]; exact z042_neg_one_le_step _
        · rw [dite_eq_right h]; norm_num

/-- 前缀和到 `L` 就是总和（补零项不存在）。 -/
private lemma z042_prefix_self {L : ℕ} (w : Fin L → Bool) : wordPrefix w L = wordSum w := by
  unfold wordPrefix wordSum
  rw [Finset.sum_fin_eq_sum_range (fun i : Fin L => wordStep (w i))]

/-- 前缀和可拆：`prefix_{k+m} = prefix_k + 后 m 步的和`（用 `k + m = L`）。 -/
lemma z042_prefix_split {L : ℕ} (w : Fin L → Bool) {k m : ℕ} (_hkm : k + m = L) :
    wordPrefix w (k + m) =
      wordPrefix w k +
        ∑ i ∈ range m, (if h : k + i < L then wordStep (w ⟨k + i, h⟩) else 0) := by
  unfold wordPrefix
  rw [Finset.sum_range_add (fun i : ℕ => if h : i < L then wordStep (w ⟨i, h⟩) else 0) k m]

/-- 后 `m` 步的和落在 `[−m, m]`（`k + m = L`）。 -/
lemma z042_tail_bounds {L : ℕ} (w : Fin L → Bool) {k m : ℕ} (_hkm : k + m = L) :
    -((m : ℕ) : ℤ) ≤
        (∑ i ∈ range m, (if h : k + i < L then wordStep (w ⟨k + i, h⟩) else 0)) ∧
      (∑ i ∈ range m, (if h : k + i < L then wordStep (w ⟨k + i, h⟩) else 0)) ≤ ((m : ℕ) : ℤ) := by
  have hlt : ∀ i ∈ range m, k + i < L := by
    intro i hi
    have := Finset.mem_range.mp hi
    omega
  constructor
  · calc -((m : ℕ) : ℤ) = ∑ _i ∈ range m, (-1 : ℤ) := by simp
      _ ≤ ∑ i ∈ range m, (if h : k + i < L then wordStep (w ⟨k + i, h⟩) else 0) := by
          apply Finset.sum_le_sum
          intro i hi
          rw [dite_eq_left (hlt i hi)]
          exact z042_neg_one_le_step _
  · calc (∑ i ∈ range m, (if h : k + i < L then wordStep (w ⟨k + i, h⟩) else 0))
        ≤ ∑ _i ∈ range m, (1 : ℤ) := by
          apply Finset.sum_le_sum
          intro i hi
          rw [dite_eq_left (hlt i hi)]
          exact z042_step_le_one _
      _ = ((m : ℕ) : ℤ) := by simp

/-- 由整数上下界得 `natAbs` 的界。 -/
lemma z042_natAbs_le {x : ℤ} {a : ℕ} (h1 : -((a : ℕ) : ℤ) ≤ x) (h2 : x ≤ ((a : ℕ) : ℤ)) :
    x.natAbs ≤ a := by
  have h : ((x.natAbs : ℕ) : ℤ) ≤ (a : ℤ) := by
    rcases Int.natAbs_eq x with hx | hx
    · rw [← hx]; exact h2
    · rw [hx] at h1; omega
  exact_mod_cast h

/-- 平衡词的前缀和界：`prefix_k ≤ L−k` 且 `−(L−k) ≤ prefix_k`
    （后半段的和 = −前缀）。 -/
lemma z042_prefix_bounds {L : ℕ} (w : Fin L → Bool) (hsum : wordSum w = 0) {k : ℕ} (hk : k ≤ L) :
    wordPrefix w k ≤ ((L - k : ℕ) : ℤ) ∧ -((L - k : ℕ) : ℤ) ≤ wordPrefix w k := by
  have hkm : k + (L - k) = L := Nat.add_sub_of_le hk
  have hsplit := z042_prefix_split w hkm
  have htb := z042_tail_bounds w hkm
  have hpre : wordPrefix w (k + (L - k)) = 0 := by
    rw [hkm, z042_prefix_self, hsum]
  exact ⟨by linarith, by linarith⟩

/-- 平衡词：`|prefix_k| ≤ min k (L−k)`。 -/
lemma z042_prefix_natAbs_le_min {L : ℕ} (w : Fin L → Bool) (hsum : wordSum w = 0) {k : ℕ}
    (hk : k ≤ L) : (wordPrefix w k).natAbs ≤ min k (L - k) := by
  have hb := z042_prefix_bounds w hsum hk
  have h1 : (wordPrefix w k).natAbs ≤ k :=
    z042_natAbs_le (z042_neg_le_prefix_self w k) (z042_prefix_le_self w k)
  have h2 : (wordPrefix w k).natAbs ≤ L - k := z042_natAbs_le hb.2 hb.1
  exact le_min h1 h2

/-- 平衡词的层高不超过 `L/2`（L 偶）。 -/
lemma z042_layer_le_half {L : ℕ} (h2 : L / 2 + L / 2 = L) {w : Fin L → Bool}
    (hsum : wordSum w = 0) : wordLayer w ≤ L / 2 := by
  unfold wordLayer
  apply Finset.max'_le
  intro y hy
  rcases Finset.mem_image.mp hy with ⟨k, hk, rfl⟩
  have hkL : k ≤ L := by
    have := Finset.mem_range.mp hk
    omega
  have h := z042_prefix_natAbs_le_min w hsum hkL
  have hmin : min k (L - k) ≤ L / 2 := by omega
  omega

/-- 层高 `= a` 意味着某步 `|prefix| = a`（`max'` 取到值）。 -/
lemma z042_exists_reach {L : ℕ} {w : Fin L → Bool} {a : ℕ} (h : wordLayer w = a) :
    ∃ k, k ≤ L ∧ (wordPrefix w k).natAbs = a := by
  have hmem : wordLayer w ∈ (range (L + 1)).image (fun k : ℕ => (wordPrefix w k).natAbs) := by
    unfold wordLayer
    exact Finset.max'_mem _ _
  rw [h] at hmem
  rcases Finset.mem_image.mp hmem with ⟨k, hk, hkabs⟩
  exact ⟨k, by have := Finset.mem_range.mp hk; omega, hkabs⟩

/-! ### 3. 达到 `±L/2` 的结构 -/

/-- **核心**：偶长度平衡词若在某步达到 `+L/2`，则必是「先全上后全下」。 -/
lemma z042_reach_pos {L : ℕ} (h2 : L / 2 + L / 2 = L) {w : Fin L → Bool}
    (hsum : wordSum w = 0) {k : ℕ} (hk : k ≤ L) (hpre : wordPrefix w k = ((L / 2 : ℕ) : ℤ)) :
    w = z042UpThenDown L := by
  -- ① 定出 k = L/2（用 |prefix| ≤ min(k, L−k)）
  have hk_le : L / 2 ≤ k := by
    have h := z042_prefix_le_self w k
    rw [hpre] at h
    exact_mod_cast h
  have hLk_le : L / 2 ≤ L - k := by
    have h := (z042_prefix_bounds w hsum hk).1
    rw [hpre] at h
    exact_mod_cast h
  have hk_eq : k = L / 2 := by omega
  have hpre_half : wordPrefix w (L / 2) = ((L / 2 : ℕ) : ℤ) := hk_eq ▸ hpre
  -- ② 前 L/2 步全 +1：和 = L/2 且每项 ≤ 1
  have hhead_sum :
      (∑ i ∈ range (L / 2), (if h : i < L then wordStep (w ⟨i, h⟩) else 0)) =
        ((L / 2 : ℕ) : ℤ) := hpre_half
  have hhead_le : ∀ i ∈ range (L / 2),
      (if h : i < L then wordStep (w ⟨i, h⟩) else 0) ≤ (1 : ℤ) := by
    intro i _
    by_cases h : i < L
    · rw [dite_eq_left h]; exact z042_step_le_one _
    · rw [dite_eq_right h]; norm_num
  have hhead_eq : ∀ i ∈ range (L / 2),
      (if h : i < L then wordStep (w ⟨i, h⟩) else 0) = (1 : ℤ) := by
    have h1 : (∑ i ∈ range (L / 2), (if h : i < L then wordStep (w ⟨i, h⟩) else 0)) =
        ∑ _i ∈ range (L / 2), (1 : ℤ) := by
      rw [hhead_sum]; simp
    exact (Finset.sum_eq_sum_iff_of_le hhead_le).mp h1
  have hhead_true : ∀ i : Fin L, i.val < L / 2 → w i = true := by
    intro i hi
    have hmem : i.val ∈ range (L / 2) := Finset.mem_range.mpr hi
    have h := hhead_eq i.val hmem
    rw [dite_eq_left i.isLt] at h
    exact z042_step_eq_one_iff.mp h
  -- ③ 后 L/2 步全 −1：后半段的和 = −L/2 且每项 ≥ −1
  have hsplit_at := z042_prefix_split w (k := L / 2) (m := L / 2) h2
  have h0 : wordPrefix w (L / 2 + L / 2) = 0 := by
    rw [h2, z042_prefix_self, hsum]
  have htail_sum :
      (∑ i ∈ range (L / 2), (if h : L / 2 + i < L then wordStep (w ⟨L / 2 + i, h⟩) else 0)) =
        -((L / 2 : ℕ) : ℤ) := by
    linarith [hsplit_at, hpre_half, h0]
  have htail_le : ∀ i ∈ range (L / 2),
      (-1 : ℤ) ≤ (if h : L / 2 + i < L then wordStep (w ⟨L / 2 + i, h⟩) else 0) := by
    intro i hi
    rw [dite_eq_left (by have := Finset.mem_range.mp hi; omega)]
    exact z042_neg_one_le_step _
  have htail_eq : ∀ i ∈ range (L / 2),
      (if h : L / 2 + i < L then wordStep (w ⟨L / 2 + i, h⟩) else 0) = (-1 : ℤ) := by
    have h1 : (∑ _i ∈ range (L / 2), (-1 : ℤ)) =
        (∑ i ∈ range (L / 2), (if h : L / 2 + i < L then wordStep (w ⟨L / 2 + i, h⟩) else 0)) := by
      rw [htail_sum]; simp
    intro i hi
    exact (((Finset.sum_eq_sum_iff_of_le htail_le).mp h1) i hi).symm
  have htail_false : ∀ j : Fin L, L / 2 ≤ j.val → w j = false := by
    intro j hj
    have hmem : j.val - L / 2 ∈ range (L / 2) := Finset.mem_range.mpr (by omega)
    have h := htail_eq (j.val - L / 2) hmem
    have hlt : L / 2 + (j.val - L / 2) < L := by omega
    rw [dite_eq_left hlt] at h
    have hfin : (⟨L / 2 + (j.val - L / 2), hlt⟩ : Fin L) = j :=
      Fin.ext (Nat.add_sub_of_le hj)
    have hstep : wordStep (w j) = -1 := by
      rw [← hfin]
      exact h
    exact z042_step_eq_neg_one_iff.mp hstep
  -- ④ 拼起来
  funext j
  rcases Nat.lt_or_ge j.val (L / 2) with hj | hj
  · rw [hhead_true j hj]
    simp only [z042UpThenDown]
    rw [decide_eq_true_iff.mpr hj]
  · rw [htail_false j hj]
    simp only [z042UpThenDown]
    rw [decide_eq_false_iff_not.mpr (Nat.not_lt.mpr hj)]

/-- 达到 `−L/2` 的情形：取负归到 `z042_reach_pos`。 -/
lemma z042_reach_neg {L : ℕ} (h2 : L / 2 + L / 2 = L) {w : Fin L → Bool}
    (hsum : wordSum w = 0) {k : ℕ} (hk : k ≤ L) (hpre : wordPrefix w k = -((L / 2 : ℕ) : ℤ)) :
    w = wordNeg (z042UpThenDown L) := by
  have hsum' : wordSum (wordNeg w) = 0 := by rw [z047_sum_neg, hsum, neg_zero]
  have hpre' : wordPrefix (wordNeg w) k = ((L / 2 : ℕ) : ℤ) := by
    rw [z047_prefix_neg, hpre, neg_neg]
  have h := z042_reach_pos h2 hsum' hk hpre'
  calc w = wordNeg (wordNeg w) := by funext i; simp [wordNeg]
    _ = wordNeg (z042UpThenDown L) := by rw [h]

/-! ### 4. 主定理一：顶层恰两个词 -/

/-- **Z042 · 顶层词恰两个**：`L` 偶、`w` 平衡且层高 `= L/2` ⟹ `w` 是先全上后全下或其取负 -/
theorem z042_top_layer_two {L : ℕ} (hL : 0 < L) (hLeven : Even L) (w : Fin L → Bool)
    (hw : w ∈ wordShell L (L / 2)) :
    w = z042UpThenDown L ∨ w = wordNeg (z042UpThenDown L) := by
  have h2 : L / 2 + L / 2 = L := by
    obtain ⟨m, rfl⟩ := hLeven
    omega
  rw [wordShell, Finset.mem_filter] at hw
  obtain ⟨hsum, hlayer⟩ := hw.2
  obtain ⟨k, hkL, hkabs⟩ := z042_exists_reach hlayer
  rcases Int.natAbs_eq_iff.mp hkabs with hpos | hneg
  · exact Or.inl (z042_reach_pos h2 hsum hkL hpos)
  · exact Or.inr (z042_reach_neg h2 hsum hkL hneg)

/-! ### 5. 两个顶层词确实在壳里（第二条定理用） -/

/-- 「先全上后全下」的词和为 0（`L = m + m`）。 -/
lemma z042_sum_utd (m : ℕ) : wordSum (z042UpThenDown (m + m)) = 0 := by
  have ha : (m + m) / 2 = m := by omega
  have hpad :
      (∑ i ∈ range (m + m),
          (if h : i < m + m then wordStep (z042UpThenDown (m + m) ⟨i, h⟩) else 0)) =
        ∑ i ∈ range (m + m), (if i < m then (1 : ℤ) else -1) := by
    apply Finset.sum_congr rfl
    intro i hi
    rw [dite_eq_left (Finset.mem_range.mp hi)]
    simp only [z042UpThenDown, ha]
    by_cases h : i < m <;> simp [wordStep, h]
  have hsplit :
      (∑ i ∈ range (m + m), (if i < m then (1 : ℤ) else -1)) = 0 := by
    rw [Finset.sum_range_add]
    have h1 : (∑ i ∈ range m, (if i < m then (1 : ℤ) else -1)) = ((m : ℕ) : ℤ) := by
      rw [show (∑ i ∈ range m, (if i < m then (1 : ℤ) else -1)) =
            ∑ _i ∈ range m, (1 : ℤ) from
        Finset.sum_congr rfl (fun i hi => by simp [Finset.mem_range.mp hi])]
      simp
    have h2 : (∑ i ∈ range m, (if m + i < m then (1 : ℤ) else -1)) = -((m : ℕ) : ℤ) := by
      rw [show (∑ i ∈ range m, (if m + i < m then (1 : ℤ) else -1)) =
            ∑ _i ∈ range m, (-1 : ℤ) from
        Finset.sum_congr rfl (fun i _ => by simp [show ¬ m + i < m by omega])]
      simp
    rw [h1, h2]; ring
  unfold wordSum
  rw [Finset.sum_fin_eq_sum_range]
  rw [hpad, hsplit]

/-- 「先全上后全下」在第 `L/2` 步的前缀和恰为 `L/2`（`L = m + m`）。 -/
lemma z042_prefix_utd_top (m : ℕ) : wordPrefix (z042UpThenDown (m + m)) m = ((m : ℕ) : ℤ) := by
  have ha : (m + m) / 2 = m := by omega
  unfold wordPrefix
  rw [show (∑ i ∈ range m,
        (if h : i < m + m then wordStep (z042UpThenDown (m + m) ⟨i, h⟩) else 0)) =
      ∑ _i ∈ range m, (1 : ℤ) from
    Finset.sum_congr rfl (fun i hi => by
      rw [dite_eq_left (by have := Finset.mem_range.mp hi; omega)]
      simp only [z042UpThenDown, ha]
      simp [wordStep, Finset.mem_range.mp hi])]
  simp

/-- 「先全上后全下」的层高恰为 `L/2`（`L = m + m`）。 -/
lemma z042_layer_utd (m : ℕ) : wordLayer (z042UpThenDown (m + m)) = (m + m) / 2 := by
  apply le_antisymm
  · have h := z042_layer_le_half (L := m + m) (by omega) (z042_sum_utd m)
    omega
  · rw [show (m + m) / 2 = m by omega]
    unfold wordLayer
    apply Finset.le_max'
    refine Finset.mem_image.mpr ⟨m, Finset.mem_range.mpr (by omega), ?_⟩
    rw [z042_prefix_utd_top m]
    simp

/-- 「先全上后全下」在顶层壳里（`(m+m)/2 = m`）。 -/
lemma z042_utd_mem (m : ℕ) : z042UpThenDown (m + m) ∈ wordShell (m + m) m := by
  rw [wordShell, Finset.mem_filter]
  refine ⟨Finset.mem_univ _, z042_sum_utd m, ?_⟩
  rw [z042_layer_utd m]
  omega

/-- 其取负也在顶层壳里。 -/
lemma z042_neg_utd_mem (m : ℕ) :
    wordNeg (z042UpThenDown (m + m)) ∈ wordShell (m + m) m := by
  rw [wordShell, Finset.mem_filter]
  refine ⟨Finset.mem_univ _, ?_, ?_⟩
  · rw [z047_sum_neg, z042_sum_utd m, neg_zero]
  · rw [z047_layer_neg, z042_layer_utd m]
    omega

/-- 两个顶层词互异（需要 `m > 0`，即长度非零）。 -/
lemma z042_utd_ne_neg {m : ℕ} (hm : 0 < m) :
    z042UpThenDown (m + m) ≠ wordNeg (z042UpThenDown (m + m)) := by
  intro h
  have h0 := congrFun h ⟨0, by omega⟩
  have hpos : 0 < (m + m) / 2 := by omega
  simp [z042UpThenDown, wordNeg, hpos] at h0

/-! ### 6. 主定理二：顶层词数 = 2 -/

/-- **Z042 · 顶层词数 = 2** -/
theorem z042_top_layer_card {L : ℕ} (hL : 0 < L) (hLeven : Even L) :
    (wordShell L (L / 2)).card = 2 := by
  obtain ⟨m, rfl⟩ := hLeven
  have hm : 0 < m := by omega
  have ha : (m + m) / 2 = m := by omega
  rw [ha]
  have hmem_pair :
      ({z042UpThenDown (m + m), wordNeg (z042UpThenDown (m + m))} :
          Finset (Fin (m + m) → Bool)) ⊆ wordShell (m + m) m := by
    intro x hx
    simp only [Finset.mem_insert, Finset.mem_singleton] at hx
    rcases hx with rfl | rfl
    · exact z042_utd_mem m
    · exact z042_neg_utd_mem m
  have hsub :
      wordShell (m + m) m ⊆
        ({z042UpThenDown (m + m), wordNeg (z042UpThenDown (m + m))} :
          Finset (Fin (m + m) → Bool)) := by
    intro w hw
    rcases z042_top_layer_two (L := m + m) (by omega) ⟨m, rfl⟩ w (by rwa [ha]) with h | h
    · rw [h]; simp
    · rw [h]; simp
  apply le_antisymm
  · calc (wordShell (m + m) m).card
        ≤ ({z042UpThenDown (m + m), wordNeg (z042UpThenDown (m + m))} :
            Finset (Fin (m + m) → Bool)).card := Finset.card_le_card hsub
      _ ≤ 2 := by
          have h := Finset.card_insert_le (z042UpThenDown (m + m))
            ({wordNeg (z042UpThenDown (m + m))} : Finset (Fin (m + m) → Bool))
          simp only [Finset.card_singleton] at h
          omega
  · have h := Finset.card_le_card hmem_pair
    rwa [Finset.card_pair (z042_utd_ne_neg hm)] at h

#print axioms z042_step_le_one
#print axioms z042_prefix_split
#print axioms z042_tail_bounds
#print axioms z042_prefix_bounds
#print axioms z042_layer_le_half
#print axioms z042_reach_pos
#print axioms z042_sum_utd
#print axioms z042_layer_utd
#print axioms z042_top_layer_two
#print axioms z042_top_layer_card

end Zero