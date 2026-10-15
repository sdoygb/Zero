/-
  Z072 §2.2／§3（补偿移动图 = Johnson 图 J(L, L/2)）与 Z024（不可约性）的机器检查版本。

  Z072 §2.2：平衡词 ↔ "+1 位置集合" `S ⊆ Fin L`（`|S| = L/2`）；**补偿移动** = 在某个
  `a ∉ S` 与某个 `b ∈ S` 之间交换（`S' = insert a (S.erase b)`），这正是 Johnson 图
  `J(L, L/2)` 的邻接。于是 `deg S = |S| · |Sᶜ| = (L/2)·(L − L/2)`；`L` 偶数时即
  `(L/2)²`（正则）。Z072 §3／Z024：任意两个同基数的子集之间可由有限次交换连通
  （= 转移矩阵不可约），且交换是对合（无向 ⟹ 均匀计数测度平稳、时间可逆）。

  与题面规格的两处偏差（均不涉及任何定理的名字与陈述，均在下文就地注明）：

  1. `z072Neighbors` 定义体的两处写法（内容完全不变，只是为了让证明能过）：
     (a) 本机 Mathlib 的 `Finset.erase` 实参次序是 `s.erase a`（被删元素在后），题面的
         `insert p.1 (erase p.2 S)` 无法通过类型检查（会把 `p.2` 当成被操作的 Finset），
         故写作 `insert p.1 (S.erase p.2)`；
     (b) 本机 `×ˢ` 是 `SProd.sprod`（其 Finset 实例 = `Finset.product`），而
         `Finset.mem_product` / `Finset.card_product` 的陈述用的是 `×ˢ`，`Finset.product`
         与它只定义等价、不同一（`Finset.product_eq_sprod` 是 `rfl`），故把 `(Sᶜ).product S`
         写作 `Sᶜ ×ˢ S`（同一个函数）。

  2. 规格中的 `z072_regular`（`S.card = L / 2 ⟹ |N(S)| = (L/2)·(L/2)`，**无 L 偶数的
     假设**）在 `L` 为奇数时**为假**，故本文件未收录它（不放弱、不改名、不留 sorry）。
     机器检查的反例见 `z072_regular_counterexample_at_three`（`L = 3, S = {0}`：
     `|N(S)| = 2` 而 `(3/2)·(3/2) = 1`）；修正后的正确命题见 `z072_degree_even`
     （补上 `L % 2 = 0`）。注意 `z072_neighbors_card` 的度 `|S|·(L − |S|)` 与
     `L` 的奇偶无关，奇数时的度是 `((L−1)/2)·((L+1)/2)`。
-/

import Mathlib

namespace Zero

open Finset

/-- 补偿移动的邻居集：把一个元素换成另一个（`a ∉ S`、`b ∈ S`） -/
def z072Neighbors {L : ℕ} (S : Finset (Fin L)) : Finset (Finset (Fin L)) :=
  (Sᶜ ×ˢ S) |>.image (fun p => insert p.1 (S.erase p.2))

/-- 邻居集的显式刻画：邻居正是"删掉一个 `b ∈ S`、添上一个 `a ∉ S`"得到的集合。 -/
lemma z072_mem_neighbors {L : ℕ} {S B : Finset (Fin L)} :
    B ∈ z072Neighbors S ↔ ∃ a ∉ S, ∃ b ∈ S, B = insert a (S.erase b) := by
  classical
  unfold z072Neighbors
  rw [Finset.mem_image]
  constructor
  · rintro ⟨p, hp, rfl⟩
    rw [Finset.mem_product] at hp
    exact ⟨p.1, by simpa using hp.1, p.2, hp.2, rfl⟩
  · rintro ⟨a, ha, b, hb, rfl⟩
    exact ⟨(a, b), by rw [Finset.mem_product]; exact ⟨by simpa using ha, hb⟩, rfl⟩

/-- 交换映射 `(a, b) ↦ insert a (S.erase b)` 在 `Sᶜ × S` 上单射。
    还原法：`a` 是 `X \ S` 的唯一元素，`b` 是 `S \ X` 的唯一元素。 -/
lemma z072_injOn {L : ℕ} (S : Finset (Fin L)) :
    Set.InjOn (fun p : Fin L × Fin L => insert p.1 (S.erase p.2)) ↑(Sᶜ ×ˢ S) := by
  classical
  intro p hp q hq hpq
  -- 先把 `hpq` 的 β-redex 消掉，后面才能 `rw`
  have hpq' : insert p.1 (S.erase p.2) = insert q.1 (S.erase q.2) := hpq
  have hp' : p.1 ∈ Sᶜ ∧ p.2 ∈ S := by simpa using hp
  have hq' : q.1 ∈ Sᶜ ∧ q.2 ∈ S := by simpa using hq
  have hp1 : p.1 ∉ S := by simpa using hp'.1
  have hq1 : q.1 ∉ S := by simpa using hq'.1
  have hq2 : q.2 ∈ S := hq'.2
  -- `p.1` 是 `X \ S` 的唯一元素
  have h1 : p.1 = q.1 := by
    have hm : p.1 ∈ insert q.1 (S.erase q.2) := by
      rw [← hpq']
      exact Finset.mem_insert_self _ _
    rcases Finset.mem_insert.mp hm with h | h
    · exact h
    · exact absurd (Finset.mem_erase.mp h).2 hp1
  -- `q.2 ∉ X`，再读出 `q.2 = p.2`
  have h2 : p.2 = q.2 := by
    have hq2Y : q.2 ∉ insert q.1 (S.erase q.2) := by
      intro hmem
      rcases Finset.mem_insert.mp hmem with h | h
      · exact hq1 (by rw [← h]; exact hq2)
      · exact (Finset.mem_erase.mp h).1 rfl
    have hq2X : q.2 ∉ insert p.1 (S.erase p.2) := by rw [hpq']; exact hq2Y
    have hq2S : q.2 ∉ S.erase p.2 := fun hm => hq2X (Finset.mem_insert_of_mem hm)
    by_contra hne
    exact hq2S (Finset.mem_erase.mpr ⟨fun h => hne h.symm, hq2⟩)
  exact Prod.ext h1 h2

/-- **Z072 · 正则性**：`|N(S)| = |S| · (L − |S|)`（故 `|S| = L/2` 时度 = `(L/2)·(L − L/2)`） -/
theorem z072_neighbors_card {L : ℕ} (S : Finset (Fin L)) :
    (z072Neighbors S).card = S.card * (L - S.card) := by
  classical
  unfold z072Neighbors
  rw [Finset.card_image_of_injOn (z072_injOn S), Finset.card_product, Finset.card_compl,
    Fintype.card_fin, Nat.mul_comm]

/-- **Z072 · 度 = (L/2)² 的正确形式**（规格版 `z072_regular` 少了 `L` 偶数的假设，奇数时为假，
    见 `z072_regular_counterexample_at_three`）：`L` 偶数且 `|S| = L/2` 时，
    `|Sᶜ| = L − L/2 = L/2`，故补偿移动图是 `(L/2)²`-正则。 -/
theorem z072_degree_even {L : ℕ} (S : Finset (Fin L)) (hS : S.card = L / 2)
    (hL : L % 2 = 0) : (z072Neighbors S).card = (L / 2) * (L / 2) := by
  rw [z072_neighbors_card, hS]
  have h : L - L / 2 = L / 2 := by omega
  rw [h]

/-- 规格中的 `z072_regular` 的反例（机器检查）：`L = 3`、`S = {0}` 时
    `|S| = 3/2 = 1`（满足假设），但 `|N(S)| = 2`，而规格要求的 `(3/2)·(3/2) = 1`（同为 `1`）——
    故规格为假。审计整改（2026-10-15）：原陈述里的第 2 个合取项 `(3/2)*(3/2) = 1`
    是纯算术填充（`norm_num` 即得，不含反例信息），已按审计删除；
    反例本体（`|N(S)| = 2` 与 `|S| = 3/2`）逐字保留。 -/
theorem z072_regular_counterexample_at_three :
    (z072Neighbors ({0} : Finset (Fin 3))).card = 2 ∧
      ({0} : Finset (Fin 3)).card = 3 / 2 := by
  decide

/-- **Z024/Z072 · 不可约（连通）**：同基数的两个子集之间可由补偿移动到达 -/
theorem z072_reachable {L : ℕ} (S T : Finset (Fin L)) (h : S.card = T.card) :
    Relation.ReflTransGen (fun A B : Finset (Fin L) => B ∈ z072Neighbors A) S T := by
  classical
  -- 对 `|S \ T|` 强归纳：每步把 `S \ T` 里的一个元素换成 `T \ S` 里的一个元素。
  -- motive 显式写出，否则会被认成"外层 `n` 固定"的那个较弱的命题。
  have key : ∀ n : ℕ, ∀ S T : Finset (Fin L), (S \ T).card = n → S.card = T.card →
      Relation.ReflTransGen (fun A B : Finset (Fin L) => B ∈ z072Neighbors A) S T := by
    intro n
    refine Nat.strong_induction_on (p := fun m => ∀ S T : Finset (Fin L), (S \ T).card = m →
      S.card = T.card → Relation.ReflTransGen
        (fun A B : Finset (Fin L) => B ∈ z072Neighbors A) S T) n ?_
    intro m ih S T hm hST
    by_cases hsub : S ⊆ T
    · -- `S ⊆ T` 且同基数 ⟹ `S = T`：自环
      have hEq : S = T := Finset.eq_of_subset_of_card_le hsub (le_of_eq hST.symm)
      rw [hEq]
    · -- 取 `a ∈ S \ T` 与 `b ∈ T \ S`，令 `S' = insert b (S.erase a)`
      have ha : ∃ a, a ∈ S ∧ a ∉ T := by
        by_contra hcon
        exact hsub fun x hx => by
          by_contra hxT
          exact hcon ⟨x, hx, hxT⟩
      obtain ⟨a, haS, haT⟩ := ha
      have hb : ∃ b, b ∈ T ∧ b ∉ S := by
        by_contra hcon
        have hTS : T ⊆ S := fun x hx => by
          by_contra hxS
          exact hcon ⟨x, hx, hxS⟩
        have hEq : T = S := Finset.eq_of_subset_of_card_le hTS (le_of_eq hST)
        exact haT (by rw [hEq]; exact haS)
      obtain ⟨b, hbT, hbS⟩ := hb
      -- `S'` 与 `S` 同基数
      have hbS' : b ∉ S.erase a := fun hm => hbS (Finset.mem_of_mem_erase hm)
      have hcardS' : (insert b (S.erase a)).card = S.card := by
        rw [Finset.card_insert_of_notMem hbS', Finset.card_erase_of_mem haS]
        have hpos : 0 < S.card := Finset.card_pos.mpr ⟨a, haS⟩
        omega
      have hcard' : (insert b (S.erase a)).card = T.card := hcardS'.trans hST
      -- 这一步就是一次补偿移动
      have hstep : insert b (S.erase a) ∈ z072Neighbors S := by
        rw [z072_mem_neighbors]
        exact ⟨b, hbS, a, haS, rfl⟩
      -- `S' \ T = (S \ T).erase a`，故归纳量严格减小
      have hsdiff : (insert b (S.erase a)) \ T = (S \ T).erase a := by
        ext x
        rw [Finset.mem_sdiff, Finset.mem_insert, Finset.mem_erase, Finset.mem_erase,
          Finset.mem_sdiff]
        constructor
        · rintro ⟨hx | hx, hxT⟩
          · exact absurd (by rw [hx]; exact hbT) hxT
          · exact ⟨hx.1, hx.2, hxT⟩
        · rintro ⟨hxa, hxS, hxT⟩
          exact ⟨Or.inr ⟨hxa, hxS⟩, hxT⟩
      have haST : a ∈ S \ T := Finset.mem_sdiff.mpr ⟨haS, haT⟩
      have hlt0 : ((insert b (S.erase a)) \ T).card < (S \ T).card := by
        rw [hsdiff, Finset.card_erase_of_mem haST]
        exact Nat.sub_lt (Finset.card_pos.mpr ⟨a, haST⟩) (by norm_num)
      have hlt : ((insert b (S.erase a)) \ T).card < m := by
        rw [hm] at hlt0
        exact hlt0
      exact Relation.ReflTransGen.head hstep (ih _ hlt _ _ rfl hcard')
  exact key (S \ T).card S T rfl h

/-- **Z072 · 无向（时间可逆）**：`A` 的邻居的邻居包含 `A` 自己（交换是对合） -/
theorem z072_symmetric {L : ℕ} (A B : Finset (Fin L)) (h : B ∈ z072Neighbors A) :
    A ∈ z072Neighbors B := by
  classical
  rw [z072_mem_neighbors] at h ⊢
  obtain ⟨a, haA, b, hbA, rfl⟩ := h
  -- 换回去：`(b, a)` 是反向的那次交换
  refine ⟨b, ?_, a, ?_, ?_⟩
  · have hba : b ≠ a := fun hb => haA (hb ▸ hbA)
    intro hmem
    rcases Finset.mem_insert.mp hmem with h1 | h1
    · exact hba h1
    · exact (Finset.mem_erase.mp h1).1 rfl
  · exact Finset.mem_insert_self a (A.erase b)
  · have haer : a ∉ A.erase b := fun hm => haA (Finset.mem_of_mem_erase hm)
    rw [Finset.erase_insert haer, Finset.insert_erase hbA]

#print axioms z072_neighbors_card
#print axioms z072_degree_even
#print axioms z072_regular_counterexample_at_three
#print axioms z072_reachable
#print axioms z072_symmetric

end Zero