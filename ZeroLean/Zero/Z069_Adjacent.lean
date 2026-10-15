/-
  Z069（出口 1：$L$ 增长的**组合核**）的机器检查版本。

  Z069 §3 的论证链是：
    (i)   自然嵌入 `w ↦ w_{<p}·(+1,−1)·w_{≥p}`：`𝒲_L → 𝒲_{L+2}`（插入一对自消对）；
    (ii)  该嵌入**满射**（每个长 `L+2` 的平衡词都由长 `L` 的词插一对自消对而来）；
    (iii) 但它**非单射**（同一个长 `L+2` 的词有多个长 `L` 的前身，Z069 (3.3)–(3.5)）；
    (iv)  故**没有典范的逆** ⟹ **无典范转移算符** ⟹ "`L` 增长"是族的参数、不是时间演化。

  本文件机器核验的是 (ii) 的**组合核**：

      $\Sigma(w)=0$ 且 $L\ge 2$ $\Longrightarrow$ $w$ 必含**相邻的异号对**。

  证明链条：和的绝对值 $\le L$，等号只在全同号时取到；故平衡词必同时含 $\pm1$
  （定理 2），于是若不存在相邻异号对，沿相邻关系归纳传递就得到"全词常数"，
  与"常数词的和不为零"（定理 1）矛盾。

  这正是 Z069 §3.2 那一步的依据：一个平衡词里总能**删掉**一对相邻的 `(+1,−1)`（或 `(−1,+1)`）
  而仍得平衡词 —— 即"满射"不是数值巧合，而是组合恒真式；
  而 (iii) 的非单射说明这个删除**没有典范的选法**（删法不唯一），故没有典范转移算符。

  复用：`Zero/Z047_Involution.lean` 的词步与闭合量（`wordStep`、`wordSum`），不重新定义结构。
-/

import Mathlib
import Zero.Z047_Involution
import Zero.Z072_Johnson

namespace Zero

open Finset

/-- **Z069 · 常数词的和不为零**（`L > 0`） -/
theorem z069_const_word_sum_ne_zero {L : ℕ} (hL : 0 < L) (b : Bool) :
    wordSum (fun _ : Fin L => b) ≠ 0 := by
  have hLZ : (L : ℤ) ≠ 0 := by exact_mod_cast hL.ne'
  unfold wordSum
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  cases b
  · simp only [wordStep, Bool.false_eq_true, ite_false, mul_neg, mul_one]
    omega
  · simp only [wordStep, ite_true, mul_one]
    exact hLZ

/-- **Z069 · 平衡词必同时含两种符号**（`L ≥ 2`） -/
theorem z069_has_both_signs {L : ℕ} (hL : 2 ≤ L) (w : Fin L → Bool) (hsum : wordSum w = 0) :
    ∃ i j : Fin L, w i ≠ w j := by
  by_contra h
  simp only [not_exists, not_not] at h
  have hL0 : 0 < L := by omega
  -- 全词同号 ⟹ 和等于常数词的和 = 0，与定理 1 矛盾
  have hconst : wordSum w = wordSum (fun _ : Fin L => w ⟨0, hL0⟩) := by
    unfold wordSum
    exact Finset.sum_congr rfl fun i _ => by rw [h i ⟨0, hL0⟩]
  have hzero : wordSum (fun _ : Fin L => w ⟨0, hL0⟩) = 0 := by rw [← hconst, hsum]
  exact z069_const_word_sum_ne_zero hL0 (w ⟨0, hL0⟩) hzero

/-- **Z069 · 相邻阶梯的传递引理**：若 `v : ℕ → Bool` 在每个相邻对上都不换号，
    则 `v` 在 `[0, L)` 上是常数（对右端点归纳）。这是定理 3 的归纳核，
    陈述在 `ℕ → Bool` 上以避免 `Fin` 证明项差异。 -/
lemma z069_nat_constant_of_no_change {L : ℕ} (v : ℕ → Bool)
    (h : ∀ i, i + 1 < L → v i = v (i + 1)) :
    ∀ j, j < L → v j = v 0 := by
  intro j
  induction j with
  | zero => intro _; rfl
  | succ k ih =>
    intro hk
    have hk' : k < L := by omega
    exact (h k hk).symm.trans (ih hk')

/-- **Z069 · 平衡词必有相邻异号对**（`L ≥ 2`）—— 这就是"满射平凡"的依据。

    陈述注记：原规格写作 `∃ i : ℕ, i + 1 < L ∧ w ⟨i, by omega⟩ ≠ w ⟨i + 1, by omega⟩`，
    但该写法**不可解析**：`∧` 左右两个合取项在同一个局部上下文里分别 elaboration，
    `by omega` 看不到另一个合取项提供的 `i + 1 < L`（实测报 `omega could not prove`）。
    这里改成等价的嵌套存在 `∃ i, ∃ _h : i + 1 < L, …`（命题内容不变，只是把合取换成柯里化），
    使 `by omega` 能看到界 `_h`。**未削弱、未改名**。 -/
theorem z069_exists_adjacent_opposite {L : ℕ} (hL : 2 ≤ L) (w : Fin L → Bool)
    (hsum : wordSum w = 0) :
    ∃ i : ℕ, ∃ _h : i + 1 < L, w ⟨i, by omega⟩ ≠ w ⟨i + 1, by omega⟩ := by
  by_contra h
  simp only [not_exists, not_not] at h
  have hL0 : 0 < L := by omega
  -- 把 `Fin L → Bool` 提升为 `ℕ → Bool`（定义域外填 `false`），避免 `Fin` 证明项差异
  let v : ℕ → Bool := fun j => if hj : j < L then w ⟨j, hj⟩ else false
  have hv : ∀ j (hj : j < L), v j = w ⟨j, hj⟩ := by
    intro j hj
    show (if hj' : j < L then w ⟨j, hj'⟩ else false) = w ⟨j, hj⟩
    rw [dite_eq_left hj]
  have hstep : ∀ i, i + 1 < L → v i = v (i + 1) := by
    intro i hi
    rw [hv i (by omega), hv (i + 1) hi]
    exact h i hi
  have key : ∀ j, j < L → v j = v 0 := z069_nat_constant_of_no_change v hstep
  -- 于是 `w` 是常数词，与定理 1 矛盾
  have hconst : wordSum w = wordSum (fun _ : Fin L => w ⟨0, hL0⟩) := by
    unfold wordSum
    refine Finset.sum_congr rfl fun i _ => ?_
    have hi : (i : ℕ) < L := i.isLt
    rw [← hv (i : ℕ) hi, key (i : ℕ) hi, hv 0 hL0]
  have hzero : wordSum (fun _ : Fin L => w ⟨0, hL0⟩) = 0 := by rw [← hconst, hsum]
  exact z069_const_word_sum_ne_zero hL0 (w ⟨0, hL0⟩) hzero

/-- **Z069 · 与补偿移动图挂钩**：平衡词（`L ≥ 2`）在 Z072 的补偿移动图里**无孤立点**
    （用 Z072 的 `z072Neighbors` 表述：存在 `T`，`T` 与词的 `+1` 位置集互为邻居或相等——
    这里给出更弱但无歧义的版本：词的 `+1` 位置集非空且非全集）

    审计整改（2026-10-15）：本条已由同文件更强的
    `z069_plus_set_neighbors_nonempty` / `z069_plus_set_not_isolated` 直接可得
    （两者都用 `z072Neighbors` 给出"非孤立"的更强结论），故按审计改为 `private`：
    陈述与证明逐字未动，文件内 `z069_plus_set_neighbors_nonempty` 仍引用它，
    `#print axioms` 行一并删除（不再进 `Zero.*` 审计空间）。 -/
private theorem z069_plus_set_nontrivial {L : ℕ} (hL : 2 ≤ L) (w : Fin L → Bool)
    (hsum : wordSum w = 0) :
    (Finset.univ.filter (fun i => w i = true)).card ≠ 0 ∧
      (Finset.univ.filter (fun i => w i = true)).card ≠ L := by
  obtain ⟨i, j, hij⟩ := z069_has_both_signs hL w hsum
  -- 平衡词同时含 `+1` 与 `−1` 的位置
  have hex_t : ∃ i, w i = true := by
    by_cases hi : w i = true
    · exact ⟨i, hi⟩
    · exact ⟨j, by cases hwi : w i <;> cases hwj : w j <;> simp_all⟩
  have hex_f : ∃ j, w j = false := by
    by_cases hj : w j = false
    · exact ⟨j, hj⟩
    · exact ⟨i, by cases hwi : w i <;> cases hwj : w j <;> simp_all⟩
  constructor
  · -- 有 `+1` 位置 ⟹ `+1` 位置集非空
    intro h0
    have hemp : Finset.univ.filter (fun i => w i = true) = ∅ := Finset.card_eq_zero.mp h0
    obtain ⟨k, hk⟩ := hex_t
    have hmem : k ∈ Finset.univ.filter (fun i => w i = true) := by simp [hk]
    rw [hemp] at hmem
    simp at hmem
  · -- 有 `−1` 位置 ⟹ `+1` 位置集不等于全集
    intro hcard
    obtain ⟨k, hk⟩ := hex_f
    have hnot : k ∉ Finset.univ.filter (fun i => w i = true) := by simp [hk]
    have hssub : (Finset.univ.filter (fun i => w i = true)) ⊂ (Finset.univ : Finset (Fin L)) := by
      rw [Finset.ssubset_iff_subset_ne]
      exact ⟨Finset.filter_subset _ _, fun heq => hnot (by rw [heq]; exact Finset.mem_univ k)⟩
    have hlt := Finset.card_lt_card hssub
    rw [Finset.card_univ, Fintype.card_fin] at hlt
    omega

/-- **Z069 · 和的计数刻画**：`Σ(w) = 2·#{i | w i = true} − L`
    （每个 `+1` 位置贡献 2、每个位置再扣 1）。 -/
lemma z069_wordSum_eq_twice_card_sub {L : ℕ} (w : Fin L → Bool) :
    wordSum w = 2 * (((Finset.univ.filter (fun i : Fin L => w i = true)).card : ℤ)) - (L : ℤ) := by
  have hleib : ∀ i : Fin L,
      wordStep (w i) = 2 * (if w i = true then (1 : ℤ) else 0) - 1 := by
    intro i
    cases hwi : w i <;> simp [wordStep]
  unfold wordSum
  rw [Finset.sum_congr rfl (fun i _ => hleib i), Finset.sum_sub_distrib]
  congr 1
  · rw [← Finset.mul_sum]
    congr 1
    exact Finset.sum_boole (R := ℤ) (fun i : Fin L => w i = true) Finset.univ
  · rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one]

/-- **Z069 · 平衡词的权**：`Σ(w) = 0` ⟹ `2·#{i | w i = true} = L`
    （等价地 `+1` 位置集的权是 `L/2`）。 -/
lemma z069_balanced_card {L : ℕ} (w : Fin L → Bool) (hsum : wordSum w = 0) :
    2 * (Finset.univ.filter (fun i : Fin L => w i = true)).card = L := by
  have h := z069_wordSum_eq_twice_card_sub w
  rw [hsum] at h
  omega

/-- **Z069 · 补偿移动图（Z072）里真的没有孤立点**：平衡词的 `+1` 位置集 `S` 在
    Z072 的 `z072Neighbors` 下的邻居集**非空** —— 即 `S` 在补偿移动图里不是孤立点：
    存在把某个 `+1` 位置换成一个 `−1` 位置得到的同规格位置集（Z072 §"补偿移动"）。
    这条直接用了 `Zero/Z072_Johnson.lean` 的 `z072Neighbors` 与 `z072_mem_neighbors`，
    是第 4 条（`z069_plus_set_nontrivial`）想挂钩的图的**真版本**；第 4 条陈述与证明按要求逐字保留，
    但已按审计改为 `private`（见其 docstring）。 -/
theorem z069_plus_set_neighbors_nonempty {L : ℕ} (hL : 2 ≤ L) (w : Fin L → Bool)
    (hsum : wordSum w = 0) :
    (z072Neighbors (Finset.univ.filter (fun i => w i = true))).Nonempty := by
  obtain ⟨hne0, hneL⟩ := z069_plus_set_nontrivial hL w hsum
  -- `S` 非空 ⟹ 取 `b ∈ S`（可删的位置）
  obtain ⟨b, hb⟩ := Finset.card_pos.mp (Nat.pos_of_ne_zero hne0)
  -- `S ≠ univ` ⟹ 取 `a ∉ S`（可添的位置）
  have hproper : (Finset.univ.filter (fun i => w i = true)) ≠ Finset.univ := by
    intro h
    exact hneL (by rw [h, Finset.card_univ, Fintype.card_fin])
  have hex_not : ∃ a, a ∉ Finset.univ.filter (fun i => w i = true) := by
    by_contra h
    simp only [not_exists, not_not] at h
    exact hproper (Finset.eq_univ_of_forall h)
  obtain ⟨a, ha⟩ := hex_not
  exact ⟨insert a ((Finset.univ.filter (fun i => w i = true)).erase b),
    z072_mem_neighbors.mpr ⟨a, ha, b, hb, rfl⟩⟩

/-- **Z069 · 无孤立点（合成版）**：平衡词的 `+1` 位置集 `S` 的**权为 `L/2`**（`2·|S| = L`），
    且 `S` 在补偿移动图里有邻居 —— 即 `S` 是补偿移动图的**非孤立顶点**。 -/
private theorem z069_plus_set_not_isolated {L : ℕ} (hL : 2 ≤ L) (w : Fin L → Bool)
    (hsum : wordSum w = 0) :
    2 * (Finset.univ.filter (fun i => w i = true)).card = L ∧
      (z072Neighbors (Finset.univ.filter (fun i => w i = true))).Nonempty :=
  ⟨z069_balanced_card w hsum, z069_plus_set_neighbors_nonempty hL w hsum⟩

#print axioms z069_const_word_sum_ne_zero
#print axioms z069_has_both_signs
#print axioms z069_nat_constant_of_no_change
#print axioms z069_exists_adjacent_opposite
#print axioms z069_wordSum_eq_twice_card_sub
#print axioms z069_balanced_card
#print axioms z069_plus_set_neighbors_nonempty

end Zero