/-
  Z089（粒子谱承载 · K₄ 三频率槽与层间忠实性）的机器检查版本 —— 图论那一半。

  Z089 的判定：耦合计数 `C(D,2) = D` 唯一定出 `D = 3` ⟹ 顶点数 `m = 4` ⟹ 局部舞台
  `Γ = K₄`（**不是具名输入**）。前一半已由 `z021_dim_forced_three` 覆盖；这里补上后一半：
  **4 个顶点、每点度 3 的简单图必是完全图 K₄** —— 于是"局部舞台是 K₄"确实被计数强制。

  用显式邻接关系而不是 `SimpleGraph`：证明完全在 `Finset` 层面（子集 + 基数相等），
  不依赖图论 API 的细节。
-/

import Mathlib

namespace Zero

/-- **Z089 · 4 顶点 3-正则 ⟹ K₄**：若 4 个顶点上的无自环对称关系每点恰连 3 个，
    则任意两个不同顶点都相连（即图是完全图 K₄）。
    推理：每点的邻域是"去掉自己后的 3 个顶点"的子集，而两者基数都是 3 ⟹ 相等。 -/
theorem z089_four_vertex_three_regular_is_complete
    (Adj : Fin 4 → Fin 4 → Prop) [DecidableRel Adj]
    (hirrefl : ∀ v, ¬ Adj v v)
    (hsymm : ∀ u v, Adj u v → Adj v u)
    (hdeg : ∀ v, (Finset.univ.filter (fun u => Adj v u)).card = 3) :
    ∀ u v, u ≠ v → Adj u v := by
  intro u v huv
  -- v 的邻域 ⊆ univ.erase v
  have hsub : (Finset.univ.filter (fun w => Adj v w)) ⊆ Finset.univ.erase v := by
    intro w hw
    rw [Finset.mem_filter] at hw
    rw [Finset.mem_erase]
    exact ⟨fun h => hirrefl v (h ▸ hw.2), Finset.mem_univ w⟩
  -- univ.erase v 的基数是 3
  have hcard : (Finset.univ.erase v : Finset (Fin 4)).card = 3 := by
    rw [Finset.card_erase_of_mem (Finset.mem_univ v), Finset.card_univ, Fintype.card_fin]
  -- 子集且基数相等 ⟹ 相等
  have heq : (Finset.univ.filter (fun w => Adj v w)) = Finset.univ.erase v := by
    refine Finset.eq_of_subset_of_card_le hsub ?_
    rw [hcard, hdeg v]
  -- u 在 v 的邻域里
  have hu : u ∈ Finset.univ.erase v := by
    rw [Finset.mem_erase]
    exact ⟨huv, Finset.mem_univ u⟩
  rw [← heq] at hu
  rw [Finset.mem_filter] at hu
  exact hsymm v u hu.2

/-- **Z089 · K₄ 的计数**：`C(4,2) = 6` 条边、每点度 `3`（同时给出 Z021 的 `C(D,2)=D` 解
    `D = 3` 对应的顶点数 `m = D + 1 = 4`）。 -/
theorem z089_K4_counts :
    Nat.choose 4 2 = 6 ∧
      (Finset.univ.filter (fun u : Fin 4 => u ≠ (0 : Fin 4))).card = 3 := by
  refine ⟨?_, ?_⟩ <;> decide

#print axioms z089_four_vertex_three_regular_is_complete
#print axioms z089_K4_counts

end Zero
