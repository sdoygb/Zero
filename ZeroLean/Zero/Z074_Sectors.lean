/-
  Z074（扇区：环上正段数 = 1 的类的计数与最低层高）的机器检查版本。

  口径（Zero 库 Z074；数值由 ZeroGPU/src/zero_sim.py 复算 L ≤ 18）：
    · 扇区 `z074Sectors w` = 把词首尾相接后 `+1` 段的个数
      = `#{i : w i = +1 ∧ w (i−1) = −1}`（环指标：`i − 1` 取模 `L`）；
    · `|扇区₁|`（平衡词里恰有 1 个环上正段的词数）`= L`（`L` 偶，L = 4..18 逐位相符）；
    · 该类词恰是 `(+)^{L/2}(−)^{L/2}` 的 `L` 个旋转（互异），故等于 `L`。

  本文件机器核验：
    1. `z074_sector_one_card`：`|扇区₁| = L`（一般偶数 `L ≥ 2`，kernel 证明，非枚举）；
    2. `z074_L4_sectors_indistinguishable`：`L = 4` 时 `r = 1` 与 `r = 2` 两类的最低层高都是 1。

  **未写入**：第 2 条的一般 `L` 版本「扇区 1 类的层高最小值 `= ⌈L/4⌉`」（原定名 `z074_sector_one_hmin`）。
  盒内未能证完，理由与卡点见文件末尾注释（不写弱化版、不改名）。

  证明思路（第 1 条）：
    扇区数 = 环上「假后真」的转向数。关键结构事实：恰有 1 个扇区的词，其真位集合是
    **从唯一转向位起沿后继方向的连续弧**——因为转向位唯一 ⟹ 除它之外的每个真位的前驱仍为真
    （`z074_blk_of_sectors_one` 里的 `hpush`），于是「第一个假位之后不再有真位」，
    真位集合 = `{succ^[s] i0 : s < a}`，再与真位数 `L/2` 比较得 `a = L/2`。
    于是词由转向位唯一确定，`k ↦ z074Blk L k` 是 `Fin L` 到扇区 1 类上的双射（`Finset.card_bij`）。
-/

import Mathlib
import Zero.Z047_Involution

namespace Zero

open Finset

/-! ## 一、定义 -/

/-- `Fin L` 上取模的界（`0 < L` 由 `i : Fin L` 保证）。 -/
theorem z074_mod_lt {L : ℕ} (i : Fin L) (a : ℕ) : a % L < L :=
  Nat.mod_lt a (Fin.pos i)

/-- 环上正段数 -/
def z074Sectors {L : ℕ} (w : Fin L → Bool) : ℕ :=
  ((Finset.univ : Finset (Fin L)).filter
    (fun i => w i = true ∧ w ⟨(i.val + L - 1) % L, z074_mod_lt i _⟩ = false)).card

/-- 环上后继：`i ↦ i + 1 (mod L)`。 -/
def z074Succ {L : ℕ} (i : Fin L) : Fin L :=
  ⟨(i.val + 1) % L, z074_mod_lt i _⟩

/-- 环上前驱：`i ↦ i − 1 (mod L)`（即 `z074Sectors` 定义里的那个位置）。 -/
def z074Pred {L : ℕ} (i : Fin L) : Fin L :=
  ⟨(i.val + L - 1) % L, z074_mod_lt i _⟩

/-- 起段集合：`+1` 位且其前驱是 `−1` 的位置集合。 -/
def z074StartSet {L : ℕ} (w : Fin L → Bool) : Finset (Fin L) :=
  (Finset.univ : Finset (Fin L)).filter (fun i => w i = true ∧ w (z074Pred i) = false)

/-- 真位集合：`+1` 的位置集合。 -/
def z074TrueSet {L : ℕ} (w : Fin L → Bool) : Finset (Fin L) :=
  (Finset.univ : Finset (Fin L)).filter (fun i => w i = true)

/-- 从 `k` 沿后继走 `s` 步（`s < c`）扫出的环上弧。 -/
def z074Arc {L : ℕ} (k : Fin L) (c : ℕ) : Finset (Fin L) :=
  (Finset.range c).image (fun s => (z074Succ^[s]) k)

/-- 弧词 `z074Blk L k`：弧 `z074Arc k (L/2)` 上取 `+1`，弧外取 `−1`。
    这是 `(+)^{L/2}(−)^{L/2}` 的旋转的等价写法（弧从 `k` 起沿后继方向）。 -/
def z074Blk (L : ℕ) (k : Fin L) : Fin L → Bool :=
  fun i => decide (i ∈ z074Arc k (L / 2))

/-- `z074Sectors` 就是起段集合的基数（定义等价）。 -/
theorem z074Sectors_eq_startSet_card {L : ℕ} (w : Fin L → Bool) :
    z074Sectors w = (z074StartSet w).card := rfl

/-- 弧词的取值：弧上为 `true`。 -/
theorem z074Blk_eq_true {L : ℕ} (k i : Fin L) :
    z074Blk L k i = true ↔ i ∈ z074Arc k (L / 2) := by
  simp [z074Blk]

/-- 弧词的取值：弧外为 `false`。 -/
theorem z074Blk_eq_false {L : ℕ} (k i : Fin L) :
    z074Blk L k i = false ↔ i ∉ z074Arc k (L / 2) := by
  constructor
  · intro h hmem
    rw [(z074Blk_eq_true k i).mpr hmem] at h
    exact Bool.noConfusion h
  · intro h
    cases hh : z074Blk L k i with
    | false => rfl
    | true => exact absurd ((z074Blk_eq_true k i).mp hh) h

/-- 布尔的两种取值形式（避免 `cases` 破坏目标形状）。 -/
theorem z074_bool_eq_false_of_ne_true {b : Bool} (h : b ≠ true) : b = false := by
  cases b with
  | false => rfl
  | true => exact absurd rfl h

/-- 布尔的两种取值形式（避免 `cases` 破坏目标形状）。 -/
theorem z074_bool_eq_true_of_ne_false {b : Bool} (h : b ≠ false) : b = true := by
  cases b with
  | false => exact absurd rfl h
  | true => rfl

/-! ## 二、`Fin L` 上的模算术：后继迭代 -/

/-- 模加一：`(a % L + 1) % L = (a + 1) % L`（`L = 0` 时两边都是 `a + 1`；`L = 1` 时都是 0）。 -/
theorem z074_mod_add_one (L a : ℕ) : (a % L + 1) % L = (a + 1) % L := by
  rcases Nat.eq_zero_or_pos L with hL | hL
  · subst hL; simp
  · rcases Nat.lt_or_ge 1 L with h2 | h2
    · rw [Nat.add_mod a 1 L, Nat.mod_eq_of_lt h2]
    · have h1 : L = 1 := by omega
      subst h1; simp

/-- 后继迭代的值：`(succ^[s] i).val = (i.val + s) % L`。 -/
theorem z074Succ_iter_val {L : ℕ} (i : Fin L) (s : ℕ) :
    ((z074Succ^[s]) i).val = (i.val + s) % L := by
  induction s with
  | zero => simp [Nat.mod_eq_of_lt i.isLt]
  | succ s ih =>
    rw [Function.iterate_succ_apply']
    simp only [z074Succ, Fin.val_mk, ih, z074_mod_add_one]
    exact congrArg (fun n => n % L) (by omega)

/-- 迭代步与模：`s < L`、`t < L` 时 `(v + s) % L = (v + t) % L ⟹ s = t`。 -/
theorem z074_mod_eq_of_lt {L v s t : ℕ} (hs : s < L) (ht : t < L)
    (h : (v + s) % L = (v + t) % L) : s = t := by
  have hm : (v + s) ≡ (v + t) [MOD L] := h
  exact Nat.ModEq.eq_of_lt_of_lt (Nat.ModEq.add_left_cancel' v hm) hs ht

/-- 迭代不重复：`s < L`、`t < L` 时 `succ^[s] i = succ^[t] i ⟹ s = t`。 -/
theorem z074_iter_val_eq {L : ℕ} (i : Fin L) {s t : ℕ} (hs : s < L) (ht : t < L)
    (h : ((z074Succ^[s]) i).val = ((z074Succ^[t]) i).val) : s = t := by
  rw [z074Succ_iter_val, z074Succ_iter_val] at h
  exact z074_mod_eq_of_lt hs ht h

/-- 后继迭代在 `range L` 上单射。 -/
theorem z074_iter_injOn {L : ℕ} (k : Fin L) :
    Set.InjOn (fun s : ℕ => (z074Succ^[s]) k) ↑(Finset.range L) := by
  intro s hs t ht h
  exact z074_iter_val_eq k (Finset.mem_range.mp hs) (Finset.mem_range.mp ht)
    (congrArg Fin.val h)

/-- 弧的基数：`c ≤ L` 时 `#(z074Arc k c) = c`。 -/
theorem z074_arc_card {L : ℕ} (k : Fin L) {c : ℕ} (hc : c ≤ L) : (z074Arc k c).card = c := by
  have hsub : ((Finset.range c : Finset ℕ) : Set ℕ) ⊆ ((Finset.range L : Finset ℕ) : Set ℕ) := by
    intro s hs
    exact Finset.mem_range.mpr (lt_of_lt_of_le (Finset.mem_range.mp hs) hc)
  rw [z074Arc, Finset.card_image_of_injOn (Set.InjOn.mono hsub (z074_iter_injOn k)),
    Finset.card_range]

/-- 弧含起点（`c ≥ 1`）。 -/
theorem z074_mem_arc_self {L : ℕ} (k : Fin L) {c : ℕ} (hc : 1 ≤ c) : k ∈ z074Arc k c := by
  rw [z074Arc, Finset.mem_image]
  exact ⟨0, Finset.mem_range.mpr hc, by simp⟩

/-- 弧在后继方向封闭（`s + 1 < c`）。 -/
theorem z074_mem_arc_succ {L : ℕ} (k : Fin L) {c s : ℕ} (hs : s + 1 < c) :
    (z074Succ^[s]) k ∈ z074Arc k c := by
  rw [z074Arc, Finset.mem_image]
  exact ⟨s, Finset.mem_range.mpr (by omega), rfl⟩

/-- 迭代遍历整个 `Fin L`：每个 `i` 都是 `succ^[s] k`（`s < L`）。 -/
theorem z074_iter_surj {L : ℕ} (k : Fin L) :
    ∀ i : Fin L, ∃ s : ℕ, s < L ∧ (z074Succ^[s]) k = i := by
  intro i
  have hcard : (z074Arc k L).card = L := z074_arc_card k le_rfl
  have huniv : z074Arc k L = Finset.univ := by
    apply Finset.eq_univ_of_card
    simpa using hcard
  have hmem : i ∈ z074Arc k L := by rw [huniv]; exact Finset.mem_univ i
  obtain ⟨s, hs, hsk⟩ := Finset.mem_image.mp (by rwa [z074Arc] at hmem)
  exact ⟨s, Finset.mem_range.mp hs, hsk⟩

/-! ## 三、前驱与后继互逆 -/

/-- `pred (succ i) = i`。 -/
theorem z074_pred_succ {L : ℕ} (i : Fin L) : z074Pred (z074Succ i) = i := by
  apply Fin.ext
  simp only [z074Pred, z074Succ, Fin.val_mk]
  by_cases h : i.val + 1 < L
  · rw [Nat.mod_eq_of_lt h]
    have h2 : i.val + 1 + L - 1 = i.val + L := by omega
    rw [h2, Nat.add_mod, Nat.mod_self, Nat.add_zero, Nat.mod_mod, Nat.mod_eq_of_lt i.isLt]
  · have h1 : i.val + 1 = L := by omega
    rw [h1, Nat.mod_self]
    rw [show 0 + L - 1 = L - 1 from by omega, Nat.mod_eq_of_lt (by omega : L - 1 < L)]
    omega

/-- 前驱与迭代：`pred (succ^[s] k) = succ^[s−1] k`（`s ≥ 1`）。 -/
theorem z074_pred_succ_iter {L : ℕ} (k : Fin L) {s : ℕ} (hs : 0 < s) :
    z074Pred ((z074Succ^[s]) k) = (z074Succ^[s - 1]) k := by
  obtain ⟨s', rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : s ≠ 0)
  rw [Function.iterate_succ_apply', z074_pred_succ, Nat.succ_sub_one]

/-- 末位：`succ^[L−1] k = pred k`。 -/
theorem z074_iter_last {L : ℕ} (k : Fin L) : (z074Succ^[L - 1]) k = z074Pred k := by
  apply Fin.ext
  rw [z074Succ_iter_val]
  show (k.val + (L - 1)) % L = (k.val + L - 1) % L
  congr 1
  omega

/-! ## 四、和 = 2·（真位数）− L -/

/-- 词和与真位数：`Σ(w) = 2·#{i : w i = +1} − L`。 -/
theorem z074_sum_eq {L : ℕ} (w : Fin L → Bool) :
    wordSum w = 2 * ((z074TrueSet w).card : ℤ) - L := by
  have hfalseFilter : ((Finset.univ : Finset (Fin L)).filter (fun i => w i = false))
      = (Finset.univ : Finset (Fin L)).filter (fun i => ¬(w i = true)) := by
    apply Finset.filter_congr
    intro i _
    cases w i <;> simp
  have hsplit := Finset.sum_filter_add_sum_filter_not
    (s := (Finset.univ : Finset (Fin L))) (p := fun i => w i = true)
    (f := fun i => wordStep (w i))
  have hA : (∑ i ∈ (Finset.univ : Finset (Fin L)).filter (fun i => w i = true),
        wordStep (w i)) = ((z074TrueSet w).card : ℤ) := by
    have h1 : (∑ i ∈ (Finset.univ : Finset (Fin L)).filter (fun i => w i = true),
          wordStep (w i))
        = ∑ _i ∈ (Finset.univ : Finset (Fin L)).filter (fun i => w i = true), (1 : ℤ) :=
      Finset.sum_congr rfl (fun i hi => by
        have := (Finset.mem_filter.mp hi).2
        simp [wordStep, this])
    rw [h1, z074TrueSet]
    simp
  have hB : (∑ i ∈ (Finset.univ : Finset (Fin L)).filter (fun i => ¬(w i = true)),
        wordStep (w i))
      = -(((Finset.univ : Finset (Fin L)).filter (fun i => w i = false)).card : ℤ) := by
    have h1 : (∑ i ∈ (Finset.univ : Finset (Fin L)).filter (fun i => ¬(w i = true)),
          wordStep (w i))
        = ∑ _i ∈ (Finset.univ : Finset (Fin L)).filter (fun i => w i = false), (-1 : ℤ) := by
      rw [← hfalseFilter]
      exact Finset.sum_congr rfl (fun i hi => by
        have := (Finset.mem_filter.mp hi).2
        simp [wordStep, this])
    rw [h1]
    simp
  have hcard : ((Finset.univ : Finset (Fin L)).filter (fun i => w i = false)).card
      = L - (z074TrueSet w).card := by
    have h := Finset.card_filter_add_card_filter_not
      (s := (Finset.univ : Finset (Fin L))) (p := fun i => w i = true)
    rw [← hfalseFilter] at h
    have h' : (z074TrueSet w).card
        + ((Finset.univ : Finset (Fin L)).filter (fun i => w i = false)).card = L := by
      simpa [z074TrueSet] using h
    omega
  have hTle : (z074TrueSet w).card ≤ L := by
    have := Finset.card_le_univ (z074TrueSet w)
    simpa using this
  unfold wordSum
  rw [← hsplit, hA, hB, hcard, Nat.cast_sub hTle]
  ring

/-- 平衡词的真位数：`Σ(w) = 0 ⟹ #{i : w i = +1} = L/2`。 -/
theorem z074_card_true_of_sum_zero {L : ℕ} (w : Fin L → Bool) (hsum : wordSum w = 0) :
    (z074TrueSet w).card = L / 2 := by
  have h := z074_sum_eq w
  rw [hsum] at h
  have h2 : 2 * ((z074TrueSet w).card : ℤ) = (L : ℤ) := by omega
  have h3 : 2 * (z074TrueSet w).card = L := by omega
  omega

/-- 偶数 `L`：`2·(L/2) = L`。 -/
theorem z074_two_mul_div_two {L : ℕ} (heven : Even L) : 2 * (L / 2) = L := by
  obtain ⟨m, hm⟩ := heven
  have h2 : 2 ∣ L := ⟨m, by omega⟩
  calc 2 * (L / 2) = L / 2 * 2 := by rw [Nat.mul_comm]
    _ = L := Nat.div_mul_cancel h2

/-! ## 五、弧词的扇区数 = 1（基元侧） -/

/-- **弧词的起段集合是单点 `{k}`**：弧上的位置里，只有起点 `k` 的前驱落在弧外。 -/
theorem z074_startSet_blk {L : ℕ} (hL : 2 ≤ L) (k : Fin L) :
    z074StartSet (z074Blk L k) = {k} := by
  have hc : 1 ≤ L / 2 := by omega
  ext i
  simp only [z074StartSet, Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_singleton]
  constructor
  · rintro ⟨htrue, hfalse⟩
    obtain ⟨s, hs, rfl⟩ := Finset.mem_image.mp ((z074Blk_eq_true k i).mp htrue)
    rw [Finset.mem_range] at hs
    rcases Nat.eq_zero_or_pos s with rfl | hspos
    · rfl
    · exfalso
      have hnotmem : (z074Succ^[s - 1]) k ∉ z074Arc k (L / 2) := by
        intro hmem
        rw [← z074_pred_succ_iter k hspos] at hmem
        exact (z074Blk_eq_false k _).mp hfalse hmem
      exact hnotmem (z074_mem_arc_succ k (c := L / 2) (s := s - 1) (by omega))
  · intro hik
    rw [hik]
    refine ⟨(z074Blk_eq_true k k).mpr (z074_mem_arc_self k hc), ?_⟩
    refine (z074Blk_eq_false k (z074Pred k)).mpr fun hmem => ?_
    obtain ⟨s, hs, heq⟩ := Finset.mem_image.mp hmem
    rw [Finset.mem_range] at hs
    have h1 : s = L - 1 := by
      refine z074_iter_val_eq k (s := s) (t := L - 1) (by omega) (by omega) ?_
      rw [heq, z074_iter_last]
    omega

/-- **弧词的和为 0**（`L` 偶）：真位数 = `L/2`，故 `Σ = 2·(L/2) − L = 0`。 -/
theorem z074_blk_sum_zero {L : ℕ} (heven : Even L) (k : Fin L) :
    wordSum (z074Blk L k) = 0 := by
  have hcard : (z074TrueSet (z074Blk L k)).card = L / 2 := by
    have hset : z074TrueSet (z074Blk L k) = z074Arc k (L / 2) := by
      ext i
      simp only [z074TrueSet, Finset.mem_filter, Finset.mem_univ, true_and, z074Blk_eq_true]
    rw [hset, z074_arc_card k (by omega : L / 2 ≤ L)]
  rw [z074_sum_eq, hcard]
  have hmul : (2 : ℤ) * ((L / 2 : ℕ) : ℤ) = (L : ℤ) := by
    have h : ((2 * (L / 2) : ℕ) : ℤ) = (L : ℤ) := by rw [z074_two_mul_div_two heven]
    simpa using h
  omega

/-! ## 六、结构定理：扇区 1 的词 = 弧词 -/

/-- **结构定理**：平衡且恰有 1 个环上正段的词，必是某个弧词 `z074Blk L k`。

    证明要点：起段集合是单点 `{i0}` ⟹ 除 `i0` 外每个真位的前驱仍为真（`hpush`）；
    设 `a` 是对 `i0` 沿后继走时第一个假位，则 `a` 之后不再有真位
    （否则取最小的「`a` 之后的真位」`s0`，由 `hpush` 其前驱也为真，与最小性矛盾），
    于是真位集合恰是 `{succ^[s] i0 : s < a}`，再与真位数 `L/2` 比较得 `a = L/2`。 -/
theorem z074_blk_of_sectors_one {L : ℕ} (hL : 2 ≤ L) {w : Fin L → Bool}
    (hsum : wordSum w = 0) (hsect : z074Sectors w = 1) :
    ∃ k : Fin L, w = z074Blk L k := by
  have hAcard : (z074StartSet w).card = 1 := hsect
  obtain ⟨i0, hi0⟩ := Finset.card_eq_one.mp hAcard
  have hi0mem : i0 ∈ z074StartSet w := by rw [hi0]; exact Finset.mem_singleton_self i0
  have hi0true : w i0 = true := (Finset.mem_filter.mp hi0mem).2.1
  have hi0false : w (z074Pred i0) = false := (Finset.mem_filter.mp hi0mem).2.2
  have huniq : ∀ j : Fin L, w j = true → w (z074Pred j) = false → j = i0 := by
    intro j h1 h2
    have hmem : j ∈ z074StartSet w := Finset.mem_filter.mpr ⟨Finset.mem_univ j, h1, h2⟩
    rw [hi0] at hmem
    exact Finset.mem_singleton.mp hmem
  have hpush : ∀ j : Fin L, w j = true → j ≠ i0 → w (z074Pred j) = true := by
    intro j hj hne
    by_contra hcon
    exact hne (huniq j hj (z074_bool_eq_false_of_ne_true hcon))
  have hTcard : (z074TrueSet w).card = L / 2 := z074_card_true_of_sum_zero w hsum
  -- 第一个假位 `a`：`succ^[s] i0` 首次离开真位集合处
  let S : Finset ℕ := (Finset.Icc 1 (L - 1)).filter (fun s => w ((z074Succ^[s]) i0) = false)
  have hSne : S.Nonempty := by
    refine ⟨L - 1, ?_⟩
    rw [Finset.mem_filter, Finset.mem_Icc]
    exact ⟨⟨by omega, le_rfl⟩, by rw [z074_iter_last]; exact hi0false⟩
  let a : ℕ := S.min' hSne
  have hamem : a ∈ S := Finset.min'_mem S hSne
  have ha1 : 1 ≤ a := (Finset.mem_Icc.mp (Finset.mem_filter.mp hamem).1).1
  have haL : a ≤ L - 1 := (Finset.mem_Icc.mp (Finset.mem_filter.mp hamem).1).2
  have hafalse : w ((z074Succ^[a]) i0) = false := (Finset.mem_filter.mp hamem).2
  have hmin : ∀ s ∈ S, a ≤ s := fun s hs => Finset.min'_le S s hs
  -- `a` 之前全是真位
  have hbefore : ∀ s, s < a → w ((z074Succ^[s]) i0) = true := by
    intro s hs
    rcases Nat.eq_zero_or_pos s with rfl | hspos
    · simpa using hi0true
    · have hnotin : s ∉ S := fun hmem => absurd (hmin s hmem) (by omega)
      have hnf : w ((z074Succ^[s]) i0) ≠ false := fun hf =>
        hnotin (by
          rw [Finset.mem_filter, Finset.mem_Icc]
          exact ⟨⟨hspos, by omega⟩, hf⟩)
      exact z074_bool_eq_true_of_ne_false hnf
  -- `a` 之后全是假位
  have hfalse_all : ∀ s, a ≤ s → s < L → w ((z074Succ^[s]) i0) = false := by
    intro s hs1 hs2
    by_contra hcon
    have hcon' : w ((z074Succ^[s]) i0) ≠ false := hcon
    have hex : ∃ t, a ≤ t ∧ t ≤ L - 1 ∧ w ((z074Succ^[t]) i0) = true :=
      ⟨s, hs1, by omega, z074_bool_eq_true_of_ne_false hcon'⟩
    let Q : Finset ℕ := (Finset.Icc a (L - 1)).filter (fun t => w ((z074Succ^[t]) i0) = true)
    obtain ⟨t0, ht0⟩ := hex
    have hQne : Q.Nonempty := ⟨t0, by
      rw [Finset.mem_filter, Finset.mem_Icc]
      exact ⟨⟨ht0.1, ht0.2.1⟩, ht0.2.2⟩⟩
    let s0 : ℕ := Q.min' hQne
    have hs0mem : s0 ∈ Q := Finset.min'_mem Q hQne
    have hs0min : ∀ t ∈ Q, s0 ≤ t := fun t ht => Finset.min'_le Q t ht
    have hs0a : a ≤ s0 := (Finset.mem_Icc.mp (Finset.mem_filter.mp hs0mem).1).1
    have hs0L : s0 ≤ L - 1 := (Finset.mem_Icc.mp (Finset.mem_filter.mp hs0mem).1).2
    have hs0true : w ((z074Succ^[s0]) i0) = true := (Finset.mem_filter.mp hs0mem).2
    have hs0ne : s0 ≠ a := by
      intro h
      rw [h] at hs0true
      rw [hs0true] at hafalse
      exact Bool.false_ne_true hafalse.symm
    have hpred_le : a ≤ s0 - 1 := by omega
    have hpred_notin : s0 - 1 ∉ Q := fun hmem => absurd (hs0min (s0 - 1) hmem) (by omega)
    have hpred_false : w ((z074Succ^[s0 - 1]) i0) = false := by
      have hmemIcc : s0 - 1 ∈ Finset.Icc a (L - 1) := Finset.mem_Icc.mpr ⟨hpred_le, by omega⟩
      by_contra hc
      exact hpred_notin (by
        rw [Finset.mem_filter]
        exact ⟨hmemIcc, z074_bool_eq_true_of_ne_false hc⟩)
    have hne_i0 : (z074Succ^[s0]) i0 ≠ i0 := by
      intro h
      have h0 : ((z074Succ^[0]) i0).val = ((z074Succ^[s0]) i0).val := by
        rw [h]
        simp
      have := z074_iter_val_eq i0 (s := 0) (t := s0) (by omega) (by omega) h0
      omega
    have hstep := hpush ((z074Succ^[s0]) i0) hs0true hne_i0
    rw [z074_pred_succ_iter i0 (s := s0) (by omega)] at hstep
    rw [hpred_false] at hstep
    exact Bool.false_ne_true hstep
  -- 真位集合 = 从 `i0` 起走 `a` 步的弧
  have hT_eq : z074TrueSet w = (Finset.range a).image (fun s => (z074Succ^[s]) i0) := by
    ext i
    simp only [z074TrueSet, Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_image,
      Finset.mem_range]
    constructor
    · intro hi
      obtain ⟨s, hsL, hsk⟩ := z074_iter_surj i0 i
      have hsa : s < a := by
        by_contra hcon
        have hfa := hfalse_all s (by omega) hsL
        rw [hsk] at hfa
        exact absurd hi (by rw [hfa]; exact Bool.false_ne_true)
      exact ⟨s, hsa, hsk⟩
    · rintro ⟨s, hs, hsk⟩
      rw [← hsk]
      exact hbefore s hs
  have ha_eq : a = L / 2 := by
    have hcard_img : ((Finset.range a).image (fun s => (z074Succ^[s]) i0)).card = a := by
      rw [Finset.card_image_of_injOn]
      · exact Finset.card_range a
      · refine Set.InjOn.mono ?_ (z074_iter_injOn i0)
        intro s hs
        exact Finset.mem_range.mpr (lt_of_lt_of_le (Finset.mem_range.mp hs) (by omega))
    have h1 : (z074TrueSet w).card = a := by rw [hT_eq, hcard_img]
    omega
  -- 收尾：词由起段位唯一确定
  refine ⟨i0, ?_⟩
  have hT_arc : z074TrueSet w = z074Arc i0 (L / 2) := by rw [hT_eq, ha_eq, z074Arc]
  have hiff : ∀ i : Fin L, (w i = true) ↔ i ∈ z074Arc i0 (L / 2) := by
    intro i
    rw [← hT_arc]
    exact ⟨fun h => Finset.mem_filter.mpr ⟨Finset.mem_univ i, h⟩,
      fun h => (Finset.mem_filter.mp h).2⟩
  funext i
  by_cases h : i ∈ z074Arc i0 (L / 2)
  · rw [(hiff i).mpr h, (z074Blk_eq_true i0 i).mpr h]
  · have hw : w i = false := z074_bool_eq_false_of_ne_true fun hh => h ((hiff i).mp hh)
    rw [hw, (z074Blk_eq_false i0 i).mpr h]

/-! ## 七、主定理 -/

/-- **Z074 · 扇区 1 的计数 = L**（平衡词、L ≥ 2 偶） -/
theorem z074_sector_one_card {L : ℕ} (hL : 2 ≤ L) (heven : Even L) :
    ((Finset.univ : Finset (Fin L → Bool)).filter
      (fun w => wordSum w = 0 ∧ z074Sectors w = 1)).card = L := by
  have hbij : (Finset.univ : Finset (Fin L)).card
      = ((Finset.univ : Finset (Fin L → Bool)).filter
          (fun w => wordSum w = 0 ∧ z074Sectors w = 1)).card := by
    refine Finset.card_bij (fun k (_ : k ∈ (Finset.univ : Finset (Fin L))) => z074Blk L k) ?_ ?_ ?_
    · intro k _
      rw [Finset.mem_filter]
      exact ⟨Finset.mem_univ _, z074_blk_sum_zero heven k,
        by rw [z074Sectors_eq_startSet_card, z074_startSet_blk hL k, Finset.card_singleton]⟩
    · intro k _ k' _ h
      have h1 : z074StartSet (z074Blk L k) = z074StartSet (z074Blk L k') := by rw [h]
      rw [z074_startSet_blk hL k, z074_startSet_blk hL k'] at h1
      simpa using h1
    · intro w hw
      rw [Finset.mem_filter] at hw
      obtain ⟨k, hk⟩ := z074_blk_of_sectors_one hL hw.2.1 hw.2.2
      exact ⟨k, Finset.mem_univ k, hk.symm⟩
  simpa using hbij.symm

/-- **Z074 · L=4 两扇区不可分辨**：r=1 与 r=2 两类的最低层高都是 1 -/
theorem z074_L4_sectors_indistinguishable :
    (∀ w : Fin 4 → Bool, wordSum w = 0 → z074Sectors w = 1 → 1 ≤ wordLayer w) ∧
    (∀ w : Fin 4 → Bool, wordSum w = 0 → z074Sectors w = 2 → 1 ≤ wordLayer w) ∧
    (∃ w : Fin 4 → Bool, wordSum w = 0 ∧ z074Sectors w = 1 ∧ wordLayer w = 1) ∧
    (∃ w : Fin 4 → Bool, wordSum w = 0 ∧ z074Sectors w = 2 ∧ wordLayer w = 1) := by
  decide

/-! ## 八、未写入的条目（诚实边界）

  **`z074_sector_one_hmin`（扇区 1 类的层高最小值 `= ⌈L/4⌉`，一般偶数 `L ≥ 4`）本轮未给出。**

  可行的一半与卡点（复现用）：`z074_sector_one_card` 已给出该类的完全刻画——
  扇区 1 类 = `{z074Blk L k : k : Fin L}`（互异）。于是 `∃w`（达到 `(L+3)/4`）与下界
  `∀w` 都化归为同一个计算：`wordLayer (z074Blk L k) = max(d, L/2 − d)`，
  其中 `d` 是基准块词 `(+)^{L/2}(−)^{L/2}` 的前缀和在旋转起点处的取值（范围 `0..L/2`）；
  对 `d` 取最小即 `⌈(L/2)/2⌉ = (L+3)/4`。

  卡点是把这一步形式化：`wordLayer` 经 `wordPrefix`（`range k` 上的 `dite` 和）定义，
  需要（i）旋转下前缀和集合的平移公式 `P^{(rot k)}_t = P_{t+k} − P_k`（含环绕段），
  与（ii）基准块词前缀和的显式值 `P_t = min(t, L − t)`（`t ≤ L/2` 与 `t ≥ L/2` 两段）。
  两条都是初等计算，但本轮时间盒（50 分钟）内未写完；
  按「不改弱、不改名」的要求，没有给出任何弱化版本。 -/

#print axioms z074_mod_lt
#print axioms z074Sectors_eq_startSet_card
#print axioms z074Blk_eq_true
#print axioms z074Blk_eq_false
#print axioms z074_bool_eq_false_of_ne_true
#print axioms z074_bool_eq_true_of_ne_false
#print axioms z074_mod_add_one
#print axioms z074Succ_iter_val
#print axioms z074_mod_eq_of_lt
#print axioms z074_iter_val_eq
#print axioms z074_iter_injOn
#print axioms z074_arc_card
#print axioms z074_mem_arc_self
#print axioms z074_mem_arc_succ
#print axioms z074_iter_surj
#print axioms z074_pred_succ
#print axioms z074_pred_succ_iter
#print axioms z074_iter_last
#print axioms z074_sum_eq
#print axioms z074_card_true_of_sum_zero
#print axioms z074_two_mul_div_two
#print axioms z074_startSet_blk
#print axioms z074_blk_sum_zero
#print axioms z074_blk_of_sectors_one
#print axioms z074_sector_one_card
#print axioms z074_L4_sectors_indistinguishable

end Zero