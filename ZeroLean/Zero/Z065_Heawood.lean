/-
  Z065 / Z021 §8：Heawood 图（Fano 平面的关联图）的显式构造与机器核验。

  舞台：Z021 §8「规模与层对齐：升级为 Heawood 图」把全局图 Γ 判定为 Heawood 图
  （|V| = 14、3-正则、围长 6、二部）。这里把它显式造出来并逐条钉死其结构性质。
  Z065 用到的 |E| = 21、圈秩 21 − 14 + 1 = 8 也在此核验（圈秩公式本身见
  Zero/Z065_CycleRank.lean）。

  顶点 `Fin 14` 的编号约定：0..6 是 Fano 平面的**点** p，7..13 是**线** j（线 j 记为顶点 7+j）。
  Fano 线用循环差集：线 j 含点 j、j+1、j+3（模 7）；故 p 在线 j 上 ⟺ (p − j) mod 7 ∈ {0,1,3}。
  Heawood 图 = 点-线关联图：只有"点 ⟷ 线"之间有边。

  证明策略：邻接写成可计算的 `Bool` 判定 `z065HeawoodAdjB`，于是"全体顶点上"的命题
  可以整体交给 `decide`（内核直接归约，不动用任何外部求值器）；含参数的定理由对应的
  "∀ 版本"（同样由 `decide` 一次算完）直接得到。边数既由 `decide` 直接数出，
  也用握手引理 `SimpleGraph.sum_degrees_eq_twice_card_edges` 配 3-正则独立交叉核验一次。
-/

import Mathlib

namespace Zero

/-! ## 一、构造（显式的 Fano 关联图） -/

/-- 点 p 是否在线 j 上：(p − j) mod 7 ∈ {0,1,3}（Fano 线的循环差集）。 -/
def z065OnFanoLine (p j : ℕ) : Bool :=
  let d := (p + 7 - j) % 7
  d == 0 || d == 1 || d == 3

/-- Heawood 邻接（可计算版本）：点（0..6）与线（7..13）之间按 Fano 关联相连；
    点-点、线-线不相邻。 -/
def z065HeawoodAdjB (u v : Fin 14) : Bool :=
  if u.val < 7 ∧ 7 ≤ v.val then z065OnFanoLine u.val (v.val - 7)
  else if v.val < 7 ∧ 7 ≤ u.val then z065OnFanoLine v.val (u.val - 7)
  else false

/-- Heawood 图：14 顶点（7 个点 + 7 条线）的 Fano 关联图。
    `symm` / `loopless` 由 `decide` 在 14×14 上直接算掉（无自环 ⟹ 简单图）。
    注：本工具链里 `Symmetric` / `Irreflexive` 是类 `Std.Symm` / `Std.Irrefl`
    （字段为 `symm` / `irrefl`），故用匿名构造子 `⟨by decide⟩` 交证据。 -/
def z065HeawoodGraph : SimpleGraph (Fin 14) where
  Adj u v := z065HeawoodAdjB u v = true
  symm := ⟨by decide⟩
  loopless := ⟨by decide⟩

/-- 显式登记 `Adj` 的可判定性。`z065HeawoodGraph` 是普通 `def`，实例搜索（在 instance 透明度
    下归约）不会把它展开成 Bool 等式，于是 `DecidableRel z065HeawoodGraph.Adj` 找不到 ——
    这里把它直接登记为实例（值是内核可算的 `Bool` 判定）。 -/
instance z065HeawoodDecidableRel : DecidableRel z065HeawoodGraph.Adj :=
  fun u v => inferInstanceAs (Decidable (z065HeawoodAdjB u v = true))

/-! ## 二、内部辅助引理（"∀ 版本"，由 `decide` 一次算完） -/

/-- 两个 `Bool` 值只取 {true,false}：若 b、x 都与 a 不同，则 b = x。（8 种情形，`decide`。） -/
private theorem z065_bool_eq_of_ne_of_ne :
    ∀ a b x : Bool, a ≠ b → a ≠ x → b = x := by
  decide

/-- 3-正则的"∀ 版本"：14 个顶点的邻域基数全为 3。 -/
private theorem z065_heawood_regular_all :
    ∀ v : Fin 14, (z065HeawoodGraph.neighborFinset v).card = 3 := by
  decide

/-- 无三角形的"∀ 版本"（14³ = 2744 组，`decide`）。 -/
private theorem z065_heawood_triangle_free_all :
    ∀ u v w : Fin 14,
      z065HeawoodGraph.Adj u v → z065HeawoodGraph.Adj v w → ¬ z065HeawoodGraph.Adj u w := by
  decide

/-- 无 4-环的"∀ 版本"：不同顶点至多共享一个邻居（14² 组，`decide`）。 -/
private theorem z065_heawood_no_four_cycle_all :
    ∀ u v : Fin 14, u ≠ v →
      ((z065HeawoodGraph.neighborFinset u) ∩ (z065HeawoodGraph.neighborFinset v)).card ≤ 1 := by
  decide

/-! ## 三、结构性质（声明与 Z065 / Z021 §8 的判据一一对应） -/

-- 审计修订（2026-10-14 夜）：删去 `z065_heawood_adj_irrefl` —— 它是构造里 `loopless` 字段的
-- 投影（零新增信息），属空心；无自环由 `SimpleGraph.loopless` 本身保证。

/-- **Z065 · 3-正则**：每个顶点恰有 3 个邻居（Z021 §8 的 deg = 3）。 -/
theorem z065_heawood_regular (v : Fin 14) : (z065HeawoodGraph.neighborFinset v).card = 3 :=
  z065_heawood_regular_all v

/-- **Z065 · 二部**：点的部分 ⟷ 线的部分（给出显式染色 `v ↦ decide (7 ≤ v.val)`）。 -/
theorem z065_heawood_bipartite :
    ∃ c : Fin 14 → Bool, ∀ u v : Fin 14, z065HeawoodGraph.Adj u v → c u ≠ c v :=
  ⟨fun v => decide (7 ≤ v.val), by decide⟩

/-- **Z065 · 无三角形**：围长 > 3（其实由二部性立即得到，这里独立算一遍）。 -/
theorem z065_heawood_triangle_free (u v w : Fin 14) :
    z065HeawoodGraph.Adj u v → z065HeawoodGraph.Adj v w → ¬ z065HeawoodGraph.Adj u w :=
  z065_heawood_triangle_free_all u v w

/-- **Z065 · 无 4-环**：不同顶点至多共享一个邻居（二部图上这等价于无 C₄）；围长 > 4。 -/
theorem z065_heawood_no_four_cycle (u v : Fin 14) (huv : u ≠ v) :
    ((z065HeawoodGraph.neighborFinset u) ∩ (z065HeawoodGraph.neighborFinset v)).card ≤ 1 :=
  z065_heawood_no_four_cycle_all u v huv

/-- **Z065 · 存在 6-环**：显式给出点-线交替的 6-环
    p₀ — L₀ — p₁ — L₁ — p₂ — L₆ — p₀，即顶点序列 `[0, 7, 1, 8, 2, 13]`
    （L₀ = {0,1,3}、L₁ = {1,2,4}、L₆ = {6,0,2}），故围长 ≤ 6；
    ★ 推理链须写全（审计修订）：无 4-环（上一条）＋ 二部 ⟹ 无奇环（见 §四 `no_odd_closed`，
    它同时排除 3-环与 5-环），两条合起来才得围长 = 6；只写「无三角形 + 无 4-环」漏了 5-环那一环。 -/
theorem z065_heawood_six_cycle :
    ∃ c : Fin 6 → Fin 14,
      (∀ i j : Fin 6, i ≠ j → c i ≠ c j) ∧
      (∀ i : Fin 6, z065HeawoodGraph.Adj (c i) (c ⟨(i.val + 1) % 6, Nat.mod_lt _ (by omega)⟩)) :=
  ⟨![0, 7, 1, 8, 2, 13], by decide, by decide⟩

/-- **Z065 · 边数 21**：直接 `decide` 数出来（`Sym2 (Fin 14)` 只有 105 个元素，内核可算）。
    下一条 private 引理给同一数字的独立交叉核验：3-正则 + 握手引理 ⟹ 42 = 2·21。
    （§ 手算：7 个点 × 3 条关联线 = 7 条线 × 3 个点 = 21 条边。） -/
theorem z065_heawood_edge_count : z065HeawoodGraph.edgeFinset.card = 21 := by
  decide

/-- 交叉核验（不走 `decide` 数边，而走握手引理
    `SimpleGraph.sum_degrees_eq_twice_card_edges` 配 3-正则）：
    2·|E| = Σ_v deg v = 14·3 = 42。 -/
private theorem z065_heawood_two_mul_edge_count :
    2 * z065HeawoodGraph.edgeFinset.card = 42 := by
  have hdeg : ∀ v : Fin 14, z065HeawoodGraph.degree v = 3 := z065_heawood_regular_all
  have hsum : (∑ v : Fin 14, z065HeawoodGraph.degree v) =
      2 * z065HeawoodGraph.edgeFinset.card :=
    SimpleGraph.sum_degrees_eq_twice_card_edges z065HeawoodGraph
  have h42 : (∑ v : Fin 14, z065HeawoodGraph.degree v) = 42 := by
    rw [show (∑ v : Fin 14, z065HeawoodGraph.degree v) = (∑ _v : Fin 14, (3 : ℕ)) from
      Finset.sum_congr rfl fun v _ => hdeg v]
    decide
  rw [← hsum]
  exact h42

/-- **Z065 · 圈秩 = 8**：|E| − |V| + 1 = 21 − 14 + 1 = 8（接 Z065_CycleRank.lean 的圈秩公式）。 -/
theorem z065_heawood_cycle_rank :
    (z065HeawoodGraph.edgeFinset.card : ℤ) - (Fintype.card (Fin 14) : ℤ) + 1 = 8 := by
  have h14 : Fintype.card (Fin 14) = 14 := by decide
  rw [z065_heawood_edge_count, h14]
  norm_num

/-! ## 四、二部 ⟹ 无奇闭道（围长 = 6 的最后一环） -/

/-- **Z065 · 没有奇闭道**：若图有二部染色 `c`，则任何闭道（`f 0 → f 1 → … → f n → f 0`）的长度 n 必为偶数。
    证明：染色沿边走一步必翻转，故 `c (f i) = c (f 0)` 与 i 的奇偶性同步；
    `f n = f 0` 迫使 n 为偶。于是 3-环、5-环等奇环都不存在；结合上面的 6-环即得围长 = 6。 -/
theorem z065_heawood_no_odd_closed (n : ℕ) (f : ℕ → Fin 14)
    (hadj : ∀ i < n, z065HeawoodGraph.Adj (f i) (f (i + 1)))
    (hclosed : f n = f 0) : Even n := by
  obtain ⟨c, hc⟩ := z065_heawood_bipartite
  -- 沿道逐步推进：偶步回到 `c (f 0)` 的色，奇步必与它不同。
  have key : ∀ i : ℕ, i ≤ n →
      (Even i → c (f i) = c (f 0)) ∧ (Odd i → c (f i) ≠ c (f 0)) := by
    intro i hi
    revert hi
    induction i with
    | zero =>
      intro _
      exact ⟨fun _ => rfl, fun h => absurd h Nat.not_odd_zero⟩
    | succ k ih =>
      intro hk
      obtain ⟨hek, hok⟩ := ih (Nat.le_of_succ_le hk)
      have hne : c (f k) ≠ c (f (k + 1)) := hc _ _ (hadj k (Nat.lt_of_succ_le hk))
      exact ⟨fun hev => z065_bool_eq_of_ne_of_ne _ _ _ hne
                (hok (Nat.not_even_iff_odd.mp (Nat.even_add_one.mp hev))),
             fun hod hcon =>
                hne ((hek (Nat.not_odd_iff_even.mp (Nat.odd_add_one.mp hod))).trans hcon.symm)⟩
  obtain ⟨-, hodd⟩ := key n le_rfl
  exact Nat.not_odd_iff_even.mp fun hn => hodd hn (by rw [hclosed])

#print axioms z065_heawood_regular
#print axioms z065_heawood_bipartite
#print axioms z065_heawood_triangle_free
#print axioms z065_heawood_no_four_cycle
#print axioms z065_heawood_six_cycle
#print axioms z065_heawood_edge_count
#print axioms z065_heawood_cycle_rank
#print axioms z065_heawood_no_odd_closed

end Zero