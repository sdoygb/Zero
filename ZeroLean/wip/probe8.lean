import Mathlib
open scoped Matrix
namespace Zero
variable {Vι : Type*} [Fintype Vι]
variable {Eι : Type*}
-- 用具体的有限类型，看 mulVecLin 的真实取向
example (M : Matrix Vι Eι ℝ) : (M.mulVecLin : (Eι → ℝ) →ₗ[ℝ] (Vι → ℝ)) := M.mulVecLin
example (M : Matrix Eι Vι ℝ) : (M.mulVecLin : (Vι → ℝ) →ₗ[ℝ] (Eι → ℝ)) := M.mulVecLin
end Zero
