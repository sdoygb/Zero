import Mathlib
open scoped Matrix
namespace Zero
variable {V : Type*} [Fintype V] [DecidableEq V] [LinearOrder V]
variable {G : SimpleGraph V} [DecidableRel G.Adj]

noncomputable def pB (G : SimpleGraph V) [DecidableRel G.Adj] :
    Matrix V {e // e ∈ G.edgeSet} ℝ :=
  Matrix.of fun v e =>
    (if (e.val : Sym2 V).inf = v then (-1 : ℝ) else 0) +
      (if (e.val : Sym2 V).sup = v then (1 : ℝ) else 0)

/-- 顶点处的"边值"：`B v S` 在 `v < w` 时是 −1（v 是下端点）、在 `w < v` 时是 +1。 -/
private lemma pB_vertex (v w : V) (hwv : G.Adj v w) :
    pB G v ⟨s(v, w), by simpa using hwv⟩
      = if v < w then -1 else if w < v then 1 else 0 := by
  rcases lt_trichotomy v w with hlt | heq | hgt
  · rw [if_pos hlt, pB, Matrix.of_apply]
    have hinf : ((s(v, w) : Sym2 V).inf : V) = v := by rw [Sym2.inf_mk, inf_eq_left.mpr hlt.le]
    have hsup : ((s(v, w) : Sym2 V).sup : V) = w := by rw [Sym2.sup_mk, sup_eq_right.mpr hlt.le]
    rw [hinf, hsup, if_pos rfl, if_neg (ne_of_lt hlt).symm]
    norm_num
  · exact absurd heq (G.ne_of_adj hwv)
  · rw [if_neg (not_lt.mpr hgt.le), if_pos hgt, pB, Matrix.of_apply]
    have hinf : ((s(v, w) : Sym2 V).inf : V) = w := by rw [Sym2.inf_mk, inf_eq_right.mpr hgt.le]
    have hsup : ((s(v, w) : Sym2 V).sup : V) = v := by rw [Sym2.sup_mk, sup_eq_left.mpr hgt.le]
    rw [hinf, hsup, if_neg (ne_of_lt hgt), if_pos rfl]
    norm_num

/-- 1 的写法。 -/
example (v w : V) (hwv : G.Adj v w) (hlt : v < w) :
    pB G w ⟨s(v, w), by simpa using hwv⟩ = 1 := by
  rw [pB_vertex w v (G.symm hwv)]
  rw [if_neg (not_lt.mpr hlt.le), if_pos hlt]

/-- −1 的写法。 -/
example (v w : V) (hwv : G.Adj v w) (hlt : v < w) :
    pB G v ⟨s(v, w), by simpa using hwv⟩ = -1 := by
  rw [pB_vertex v w hwv, if_pos hlt]
end Zero
