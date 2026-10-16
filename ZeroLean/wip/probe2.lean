import Mathlib
open scoped Matrix
namespace Zero
variable {V : Type*} [Fintype V] [DecidableEq V]
noncomputable def testB (G : SimpleGraph V) [DecidableRel G.Adj] : Matrix V {e // e ∈ G.edgeSet} ℝ :=
  Matrix.of fun _ _ => 0
-- 检查各记法的实际类型
#check @Matrix.vecMulLinear
#check (testB (⊥ : SimpleGraph V)).vecMulLinear
#check (testB (⊥ : SimpleGraph V))ᵀ.mulVecLin
#check ((testB (⊥ : SimpleGraph V)).row)
#check ((testB (⊥ : SimpleGraph V)).col)
end Zero
