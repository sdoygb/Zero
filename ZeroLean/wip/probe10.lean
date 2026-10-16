import Mathlib
open scoped Matrix
namespace Zero
noncomputable def tA : Matrix (Fin 2) (Fin 3) ℝ := Matrix.of fun _ _ => 0
#check (tA.mulVecLin : (Fin 3 → ℝ) →ₗ[ℝ] (Fin 2 → ℝ))   -- 期望：成立
#check (tA.mulVecLin : (Fin 2 → ℝ) →ₗ[ℝ] (Fin 3 → ℝ))   -- 期望：失败
example : (tA.mulVecLin : (Fin 3 → ℝ) →ₗ[ℝ] (Fin 2 → ℝ)) = tA.mulVecLin := rfl
end Zero
