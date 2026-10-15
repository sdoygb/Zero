/-
  Z048（面积能否定出 ε_h？—— 不能，但面积确实再分了层）的机器检查版本。

  本篇机器核验 Z048 §1.2 的判据 (1.4)：

      零面积词存在  ⟺  L ≡ 0 (mod 4) 。

  · 必要性（奇偶论证）：每步 ±1 ⟹ 前缀和 `s_k ≡ k (mod 2)`（对 `k ≤ L`），
    故 `a(w) ≡ Σ_{k≤L} k = L(L+1)/2 (mod 2)`；`L ≡ 2 (mod 4)` 时该数为**奇**
    ⟹ 一切词（含平衡词）的面积都是奇数，无零面积词
    （数值：L=6,10,14,18 零面积词数 = 0，见 ZeroGPU/src/zero_sim.py）。
  · 充分性（显式词族）：`L = 4k`、`w = (−)^k (+)^{2k} (−)^k`，其前缀和关于中点 `2k`
    反对称（`s_{4k−j} = −s_j`）⟹ 平衡且面积 = 0
    （数值：L=4 给 `−++−`、L=8 给 `−−++++−−`、L=12 给 `−−−++++++−−−`）。

  ⚠ 任务书原第 1 条（对**一切** `k` 的 `z048_prefix_odd_iff`）为**假**：
    `wordPrefix` 在 `k ≥ L` 后恒定（越界步补 0），于是 `k > L` 且 `k ≢ L (mod 2)` 时奇偶不符。
    机器核验反例：`L = 1`、`w = (true)`、`k = 2` 给 `s_2 = 1`（奇）而 `2` 是偶
    —— 见 `z048_prefix_odd_iff_refuted`（全称陈述的否证）。
    真命题是加 `k ≤ L` 的 `z048_prefix_even_iff_of_le`：**新名**，未复用原名，
    也未在原名下写任何被削弱的陈述。
-/

import Mathlib
import Zero.Z047_Involution

namespace Zero

open Finset

/-- 面积：`a(w) = Σ_{k ≤ L} 前缀和_k` -/
def z048Area {L : ℕ} (w : Fin L → Bool) : ℤ := ∑ k ∈ Finset.range (L + 1), wordPrefix w k

/-- **Z048 · 前缀和的递推**：`s_{k+1} = s_k + 步_k`；`k ≥ L` 时步为 0，故 `s` 越界后恒定。 -/
lemma z048_prefix_succ {L : ℕ} (w : Fin L → Bool) (k : ℕ) :
    wordPrefix w (k + 1) = wordPrefix w k + (if h : k < L then wordStep (w ⟨k, h⟩) else 0) := by
  unfold wordPrefix
  rw [Finset.sum_range_succ]

/-- 词步的两个具体取值（技术引理）。 -/
private lemma z048_wordStep_false : wordStep false = -1 := rfl

private lemma z048_wordStep_true : wordStep true = 1 := rfl

/-- **Z048 · 前缀和的奇偶（`k ≤ L`）**：`s_k − k` 恒为偶，即 `s_k ≡ k (mod 2)`。 -/
lemma z048_prefix_sub_even {L : ℕ} (w : Fin L → Bool) (k : ℕ) (hk : k ≤ L) :
    Even (wordPrefix w k - (k : ℤ)) := by
  revert hk
  induction k with
  | zero => intro _; exact ⟨0, by simp [wordPrefix]⟩
  | succ k ih =>
      intro hk
      have hk' : k ≤ L := by omega
      have hklt : k < L := by omega
      have key : wordPrefix w k + wordStep (w ⟨k, hklt⟩) - ((k + 1 : ℕ) : ℤ)
          = (wordPrefix w k - (k : ℤ)) + (wordStep (w ⟨k, hklt⟩) - 1) := by
        push_cast
        ring
      have hstep : Even (wordStep (w ⟨k, hklt⟩) - 1) := by
        cases hb : w ⟨k, hklt⟩
        · rw [z048_wordStep_false]
          exact ⟨-1, by ring⟩
        · rw [z048_wordStep_true]
          exact ⟨0, by ring⟩
      rw [z048_prefix_succ w k, dite_eq_left hklt, key]
      exact (ih hk').add hstep

/-- 显式词族：`(−)^k (+)^{2k} (−)^k` —— 前 `k` 个 `false`、中间 `2k` 个 `true`、后 `k` 个 `false`。 -/
def z048ZeroAreaWord (k : ℕ) : Fin (4 * k) → Bool :=
  fun i => if i.val < k then false else if i.val < 3 * k then true else false

/-- 词族的取值规则（技术引理）。 -/
private lemma z048_zeroAreaWord_apply {k : ℕ} (i : Fin (4 * k)) :
    z048ZeroAreaWord k i
      = if i.val < k then false else if i.val < 3 * k then true else false := rfl

/-- **Z048 · 词族的三段前缀和**：对 `j ≤ 4k`，
    `s_j = 2·(#中段前 j 位) − j = 2(min(j,3k) − min(j,k)) − j`
    （即 `−j`（j ≤ k）、`j − 2k`（k ≤ j ≤ 3k）、`4k − j`（j ≥ 3k）三段）。 -/
lemma z048_family_prefix (k : ℕ) (j : ℕ) (hj : j ≤ 4 * k) :
    wordPrefix (z048ZeroAreaWord k) j
      = 2 * (((min j (3 * k) : ℕ) : ℤ) - ((min j k : ℕ) : ℤ)) - (j : ℤ) := by
  revert hj
  induction j with
  | zero => intro _; simp [wordPrefix]
  | succ j ih =>
      intro hj
      have hjle : j ≤ 4 * k := by omega
      have hjlt : j < 4 * k := by omega
      rw [z048_prefix_succ (z048ZeroAreaWord k) j, dite_eq_left hjlt, ih hjle,
        z048_zeroAreaWord_apply (k := k) ⟨j, hjlt⟩]
      by_cases h1 : j < k
      · have hm1 : min j (3 * k) = j := Nat.min_eq_left (by omega)
        have hm2 : min j k = j := Nat.min_eq_left (le_of_lt h1)
        have hm3 : min (j + 1) (3 * k) = j + 1 := Nat.min_eq_left (by omega)
        have hm4 : min (j + 1) k = j + 1 := Nat.min_eq_left (by omega)
        rw [ite_eq_left h1, hm1, hm2, hm3, hm4, z048_wordStep_false]
        push_cast
        ring
      · by_cases h2 : j < 3 * k
        · have hm1 : min j (3 * k) = j := Nat.min_eq_left (le_of_lt h2)
          have hm2 : min j k = k := Nat.min_eq_right (by omega)
          have hm3 : min (j + 1) (3 * k) = j + 1 := Nat.min_eq_left (by omega)
          have hm4 : min (j + 1) k = k := Nat.min_eq_right (by omega)
          rw [ite_eq_right h1, ite_eq_left h2, hm1, hm2, hm3, hm4, z048_wordStep_true]
          push_cast
          ring
        · have hm1 : min j (3 * k) = 3 * k := Nat.min_eq_right (by omega)
          have hm2 : min j k = k := Nat.min_eq_right (by omega)
          have hm3 : min (j + 1) (3 * k) = 3 * k := Nat.min_eq_right (by omega)
          have hm4 : min (j + 1) k = k := Nat.min_eq_right (by omega)
          rw [ite_eq_right h1, ite_eq_right h2, hm1, hm2, hm3, hm4, z048_wordStep_false]
          push_cast
          ring

/-- **Z048 · 词族前缀和的反对称性**：`s_{4k−j} = −s_j`（`j ≤ 4k`）—— 面积为零的承重步骤。 -/
lemma z048_family_antisymm (k : ℕ) (j : ℕ) (hj : j ≤ 4 * k) :
    wordPrefix (z048ZeroAreaWord k) (4 * k - j) = -wordPrefix (z048ZeroAreaWord k) j := by
  rw [z048_family_prefix k (4 * k - j) (by omega), z048_family_prefix k j hj]
  by_cases h1 : j ≤ k
  · rw [Nat.min_eq_left (by omega : j ≤ 3 * k), Nat.min_eq_left h1,
      Nat.min_eq_right (by omega : 3 * k ≤ 4 * k - j),
      Nat.min_eq_right (by omega : k ≤ 4 * k - j),
      Nat.cast_sub (by omega : j ≤ 4 * k)]
    push_cast
    ring
  · have h1' : k < j := by omega
    by_cases h2 : j ≤ 3 * k
    · rw [Nat.min_eq_left h2, Nat.min_eq_right (le_of_lt h1'),
        Nat.min_eq_left (by omega : 4 * k - j ≤ 3 * k),
        Nat.min_eq_right (by omega : k ≤ 4 * k - j),
        Nat.cast_sub (by omega : j ≤ 4 * k)]
      push_cast
      ring
    · have h2' : 3 * k < j := by omega
      rw [Nat.min_eq_right (le_of_lt h2'), Nat.min_eq_right (by omega : k ≤ j),
        Nat.min_eq_left (by omega : 4 * k - j ≤ 3 * k),
        Nat.min_eq_left (by omega : 4 * k - j ≤ k),
        Nat.cast_sub (by omega : j ≤ 4 * k)]
      push_cast
      ring

/-- **Z048 · `L(L+1)/2` 的奇偶**：`L ≡ 2 (mod 4) ⟹ L(L+1)/2` 为奇（写 `L = 4m+2` 后
    该式为 `(4m+3)(2m+1)`，两个奇数之积）。 -/
lemma z048_half_prod_odd (L : ℕ) (hL : L % 4 = 2) : Odd ((L + 1) * L / 2) := by
  obtain ⟨m, rfl⟩ : ∃ m, L = 4 * m + 2 := ⟨L / 4, by omega⟩
  have h2 : (4 * m + 2 + 1) * (4 * m + 2) = ((4 * m + 3) * (2 * m + 1)) * 2 := by ring
  rw [h2]
  have hdiv : ((4 * m + 3) * (2 * m + 1)) * 2 / 2 = (4 * m + 3) * (2 * m + 1) := by omega
  rw [hdiv]
  exact Odd.mul ⟨2 * m + 1, by ring⟩ ⟨m, by ring⟩

/-- **Z048 · 反例（机器核验）**：任务书原第 1 条的逐字陈述
    （`∀ L w k, Even (wordPrefix w k) ↔ Even k`）为**假**。
    见证：`L = 1`、`w = (true)`、`k = 2` —— 越界补 0 使 `s_2 = s_1 = 1`（奇），而 `2` 是偶。 -/
theorem z048_prefix_odd_iff_refuted :
    ¬ (∀ (L : ℕ) (w : Fin L → Bool) (k : ℕ), Even (wordPrefix w k) ↔ Even k) := by
  intro h
  have h2 : Even (wordPrefix (fun _ : Fin 1 => true) 2) ↔ Even 2 := h 1 (fun _ => true) 2
  have h1 : wordPrefix (fun _ : Fin 1 => true) 2 = 1 := by decide
  rw [h1] at h2
  have hEven : Even (1 : ℤ) := h2.mpr ⟨1, by norm_num⟩
  have hmod : (1 : ℤ) % 2 = 0 := Int.even_iff.mp hEven
  norm_num at hmod

/-- **Z048 · 前缀和的奇偶（修正版）**：`k ≤ L` 时 `s_k ≡ k (mod 2)`。
    这是原第 1 条的真命题形式（**新名**；原名下的全称陈述已被 `z048_prefix_odd_iff_refuted` 否证）。 -/
theorem z048_prefix_even_iff_of_le {L : ℕ} (w : Fin L → Bool) (k : ℕ) (hk : k ≤ L) :
    (Even (wordPrefix w k) ↔ Even k) :=
  (Int.even_sub.mp (z048_prefix_sub_even w k hk)).trans (Int.even_coe_nat k)

/-- **Z048 · 必要性**：`L ≡ 2 (mod 4)` 时一切词（含平衡词）的面积都是**奇数**，故无零面积词 -/
theorem z048_area_odd_of_mod_four_two {L : ℕ} (hL : L % 4 = 2) (w : Fin L → Bool) :
    Odd (z048Area w) := by
  have hA : Even (z048Area w - ∑ k ∈ Finset.range (L + 1), (k : ℤ)) := by
    rw [z048Area, ← Finset.sum_sub_distrib]
    exact Finset.even_sum _ (fun k hk =>
      z048_prefix_sub_even w k (by have := Finset.mem_range.mp hk; omega))
  have hSU : (∑ k ∈ Finset.range (L + 1), (k : ℤ)) = (((L + 1) * L / 2 : ℕ) : ℤ) := by
    rw [← Nat.cast_sum, Nat.cast_inj]
    rw [Finset.sum_range_id, show (L + 1) - 1 = L from by omega]
  have hB : Odd (∑ k ∈ Finset.range (L + 1), (k : ℤ)) := by
    rw [hSU]
    exact (Int.odd_coe_nat _).mpr (z048_half_prod_odd L hL)
  have hsplit : z048Area w
      = (z048Area w - ∑ k ∈ Finset.range (L + 1), (k : ℤ))
        + ∑ k ∈ Finset.range (L + 1), (k : ℤ) := by ring
  rw [hsplit]
  exact hA.add_odd hB

/-- **Z048 · 充分性的显式词族**：`L = 4k` 时 `(−)^k (+)^{2k} (−)^k` 平衡且面积为零
    （原陈述带的 `hk : 1 ≤ k` 在证明中未被用到，且 `k = 0` 时命题亦真，
    已按审计删除该假设 = 加强命题）。 -/
theorem z048_zero_area_exists (k : ℕ) :
    ∃ w : Fin (4 * k) → Bool, wordSum w = 0 ∧ z048Area w = 0 := by
  -- 全长前缀和 = 词和（与 `Z042_TopLayer.lean` 的 `z042_prefix_self` 逐字同陈述；
  -- 按审计两处只保留 Z042 的那份**声明**，此处改为局部 `have`，不再新立重复引理）
  have hlen : wordPrefix (z048ZeroAreaWord k) (4 * k) = wordSum (z048ZeroAreaWord k) := by
    have h : wordSum (z048ZeroAreaWord k) = ∑ i ∈ Finset.range (4 * k),
        (fun j : ℕ => if h : j < 4 * k then wordStep (z048ZeroAreaWord k ⟨j, h⟩) else 0) i := by
      rw [wordSum]
      rw [← Fin.sum_univ_eq_sum_range
        (fun j : ℕ => if h : j < 4 * k then wordStep (z048ZeroAreaWord k ⟨j, h⟩) else 0) (4 * k)]
      exact Finset.sum_congr rfl (fun i _ => by simp)
    rw [h]
    rfl
  refine ⟨z048ZeroAreaWord k, ?_, ?_⟩
  · rw [← hlen, z048_family_prefix k (4 * k) (le_refl _),
      Nat.min_eq_right (by omega : 3 * k ≤ 4 * k), Nat.min_eq_right (by omega : k ≤ 4 * k)]
    push_cast
    ring
  · have hzero : (∑ j ∈ Finset.range (4 * k + 1), wordPrefix (z048ZeroAreaWord k) j) = 0 := by
      have hreflect : (∑ j ∈ Finset.range (4 * k + 1), wordPrefix (z048ZeroAreaWord k) j)
          = ∑ j ∈ Finset.range (4 * k + 1),
              wordPrefix (z048ZeroAreaWord k) (4 * k - j) := by
        conv_lhs => rw [← Finset.sum_range_reflect
          (fun j => wordPrefix (z048ZeroAreaWord k) j) (4 * k + 1)]
        simp only [Nat.add_sub_cancel]
      have hself : (∑ j ∈ Finset.range (4 * k + 1), wordPrefix (z048ZeroAreaWord k) j)
          = -(∑ j ∈ Finset.range (4 * k + 1), wordPrefix (z048ZeroAreaWord k) j) := by
        have h2 : (∑ j ∈ Finset.range (4 * k + 1),
              wordPrefix (z048ZeroAreaWord k) (4 * k - j))
            = -(∑ j ∈ Finset.range (4 * k + 1), wordPrefix (z048ZeroAreaWord k) j) := by
          rw [Finset.sum_congr rfl (fun j hj => z048_family_antisymm k j (by
            have := Finset.mem_range.mp hj
            omega)), Finset.sum_neg_distrib]
        exact hreflect.trans h2
      omega
    show (∑ j ∈ Finset.range (4 * k + 1), wordPrefix (z048ZeroAreaWord k) j) = 0
    exact hzero

#print axioms z048_prefix_succ
#print axioms z048_prefix_sub_even
#print axioms z048_family_prefix
#print axioms z048_family_antisymm
#print axioms z048_half_prod_odd
#print axioms z048_prefix_odd_iff_refuted
#print axioms z048_prefix_even_iff_of_le
#print axioms z048_area_odd_of_mod_four_two
#print axioms z048_zero_area_exists

end Zero