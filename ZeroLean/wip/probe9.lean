import Mathlib
open scoped Matrix
namespace Zero
noncomputable def tA : Matrix (Fin 2) (Fin 3) ℝ := Matrix.of fun _ _ => 0
noncomputable def tT : (Fin 3 → ℝ) →ₗ[ℝ] (Fin 2 → ℝ) := tA.mulVecLin
#check tT
noncomputable def tT2 : (Fin 3 → ℝ) →ₗ[ℝ] (Fin 2 → ℝ) := tA.vecMulLinear
noncomputable def tT3 : (Fin 3 → ℝ) →ₗ[ℝ] (Fin 2 → ℝ) := (Matrix.transpose tA).mulVecLin
#check tT2
end Zero
