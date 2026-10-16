import Mathlib
open scoped Matrix
namespace Zero
noncomputable def tB (G : SimpleGraph (Fin 3)) [DecidableRel G.Adj] :
    Matrix {e // e ∈ G.edgeSet} (Fin 3) ℝ := Matrix.of fun _ _ => 0
-- A: mulVecLin 的取向（用 mulVec_apply_eq_sum 直接看函数值）
example (G : SimpleGraph (Fin 3)) [DecidableRel G.Adj] (v : Fin 3 → ℝ) (e : {e // e ∈ G.edgeSet}) :
    ((tB G).mulVecLin v) e = ∑ i, tB G e i * v i := by
  rw [Matrix.mulVecLin_apply, Matrix.mulVec_apply_eq_sum]
-- B: 如果它真是 (V→ℝ) → (E→ℝ)，那么下面这条也该过
example (G : SimpleGraph (Fin 3)) [DecidableRel G.Adj] (v : Fin 3 → ℝ) :
    (tB G).mulVecLin v = tB G *ᵥ v := rfl
end Zero
