/-
  Zb048 §2 的**无序**计数：K₄ 恰有 3 个完美匹配。

  思路：先证「所有二元边」就是 6 条显式边；完美匹配的每条边都是二元边，故它是
  `zb048Edges` 的一个 2 元子集；于是只需在 15 个候选（6 条边选 2）上枚举。
  这样把枚举从「16 个顶点子集」降到「15 个边对」，避开 `Finset (Fin 4)` 的递归展开。
-/

import Zero.Zb048_K4

namespace Zero

open Finset

/-- K₄ 的 6 条边，显式列出。 -/
def zb048Edges : Finset (Finset (Fin 4)) :=
  {zb048Edge 0 1, zb048Edge 0 2, zb048Edge 0 3,
   zb048Edge 1 2, zb048Edge 1 3, zb048Edge 2 3}

/-- 所有二元子集恰好是这 6 条边。 -/
theorem zb048_all_edges :
    (Finset.univ.filter (fun e : Finset (Fin 4) => e.card = 2)) = zb048Edges := by
  decide

/-- 完美匹配的每条边都是二元边，故整个匹配落在 `zb048Edges` 里。 -/
theorem zb048_perfect_subset_edges {M : Finset (Finset (Fin 4))}
    (hM : zb048IsPerfectMatching M) : M ⊆ zb048Edges := by
  intro e he
  have hedge : e.card = 2 := hM.2.1 e he
  have hmem : e ∈ Finset.univ.filter (fun e : Finset (Fin 4) => e.card = 2) := by
    simp [hedge]
  rwa [zb048_all_edges] at hmem

/-- 完美匹配的枚举可以限制在「6 条边选 2」的 15 个候选上。 -/
theorem zb048_candidates :
    (Finset.univ.filter zb048IsPerfectMatching)
      = (zb048Edges.powerset.filter (fun M => M.card = 2)).filter zb048IsPerfectMatching := by
  ext M
  constructor
  · intro hM
    have hM' : zb048IsPerfectMatching M := (Finset.mem_filter.mp hM).2
    have hsub : M ⊆ zb048Edges := zb048_perfect_subset_edges hM'
    have hmem : M ∈ zb048Edges.powerset := Finset.mem_powerset.mpr hsub
    exact Finset.mem_filter.mpr
      ⟨Finset.mem_filter.mpr ⟨hmem, hM'.1⟩, hM'⟩
  · intro hM
    have hM' : zb048IsPerfectMatching M := (Finset.mem_filter.mp hM).2
    exact Finset.mem_filter.mpr ⟨Finset.mem_univ M, hM'⟩

/-- **Zb048 · 无序恰有 3 个**：15 个候选里恰有 3 个完美匹配，即 `M1, M2, M3`。 -/
theorem zb048_three_matchings_unordered :
    (zb048Edges.powerset.filter (fun M => M.card = 2)).filter zb048IsPerfectMatching
      = {zb048M1, zb048M2, zb048M3} := by
  decide

#print axioms zb048_all_edges
#print axioms zb048_perfect_subset_edges
#print axioms zb048_candidates
#print axioms zb048_three_matchings_unordered

end Zero
