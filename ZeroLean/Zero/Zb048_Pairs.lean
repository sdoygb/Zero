/-
  Zb048 §2 的**有序对**计数：6 条边里恰有 6 个有序边对构成完美匹配。

  一个完美匹配是两条互不相交、覆盖四点的边；把两条边按顺序写成 `(e1,e2)`，
  每个匹配贡献 2 个有序对，所以 6 个有序对 = 3 个无序匹配。
  这里只把「有序对 = 6」送进内核；「无序 = 3」需要再证一次除以 2 的双射，
  留给下一批（见 Zb048_K4.lean 的说明）。
-/

import Mathlib

namespace Zero

/-- 六条边的显式编号。 -/
def zb048E : Fin 6 → Finset (Fin 4) :=
  ![({0, 1} : Finset (Fin 4)), {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}]

/-- 有序边对 `(i,j)` 是完美匹配：两条边不同、不相交、并集覆盖四点。 -/
def zb048IsPerfectPair (i j : Fin 6) : Prop :=
  i ≠ j ∧ Disjoint (zb048E i) (zb048E j) ∧
    (∀ v : Fin 4, v ∈ zb048E i ∨ v ∈ zb048E j)

instance zb048PairDecidable (i j : Fin 6) : Decidable (zb048IsPerfectPair i j) := by
  unfold zb048IsPerfectPair
  infer_instance

/-- **Zb048 · 有序完美边对 = 6**：36 个有序边对里恰有 6 个构成完美匹配。 -/
theorem zb048_ordered_pairs :
    (Finset.univ.filter (fun p : Fin 6 × Fin 6 => zb048IsPerfectPair p.1 p.2)).card = 6 := by
  decide

#print axioms zb048_ordered_pairs

end Zero
