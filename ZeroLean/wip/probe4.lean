import Mathlib
open scoped Matrix
namespace Zero
variable {V : Type*} [Fintype V] [DecidableEq V]
noncomputable def tB1 (G : SimpleGraph V) [DecidableRel G.Adj] :
    Matrix {e // e ∈ G.edgeSet} V ℝ := Matrix.of fun _ _ => 0
noncomputable def tB2 (G : SimpleGraph V) [DecidableRel G.Adj] :
    Matrix V {e // e ∈ G.edgeSet} ℝ := Matrix.of fun _ _ => 0
#check (tB1 (⊥ : SimpleGraph V))
#check (tB2 (⊥ : SimpleGraph V))
#check (tB1 (⊥ : SimpleGraph V)).mulVecLin
#check (tB2 (⊥ : SimpleGraph V)).mulVecLin
end Zero
