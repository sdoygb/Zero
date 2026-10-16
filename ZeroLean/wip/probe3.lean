import Mathlib
open scoped Matrix
namespace Zero
variable {V : Type*} [Fintype V] [DecidableEq V]
noncomputable def testB (G : SimpleGraph V) [DecidableRel G.Adj] :
    Matrix {e // e ∈ G.edgeSet} V ℝ := Matrix.of fun _ _ => 0
noncomputable def testT (G : SimpleGraph V) [DecidableRel G.Adj] :
    ({e // e ∈ G.edgeSet} → ℝ) →ₗ[ℝ] (V → ℝ) := (testB G).mulVecLin
#check @testT
end Zero
