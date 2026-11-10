/-
  Zb080 · 4d Hodge 对偶：2-形式星算子的代数核

  对应篇：Zb080 §2–§4。
    · 欧氏 4d：star^2 = 1（六维 2-形式空间）；
    · 洛伦兹型：一个 6×6 实复结构模型 J^2 = -1；
      这不是完整 4-立方体离散复形的 Hodge 算子。
  自对偶/反自对偶分解与 F∧F 恒等式由脚本 `ZeroGPU/src/zb080_hodge_duality.py` 核验。
-/

import Mathlib

namespace Zero

open Matrix

/-- 欧氏 4d 的 Hodge 星算子，基 (12,13,14,23,24,34)：
`*(12)=34, *(13)=-24, *(14)=23, *(23)=14, *(24)=-13, *(34)=12`。 -/
def zb080Star : Matrix (Fin 6) (Fin 6) ℤ :=
  !![0, 0, 0, 0, 0, 1;
     0, 0, 0, 0, -1, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, -1, 0, 0, 0, 0;
     1, 0, 0, 0, 0, 0]

/-- **Zb080 定理 2.1（欧氏）**：`star^2 = 1`。 -/
theorem zb080_star_sq : zb080Star * zb080Star = 1 := by
  ext i j <;> fin_cases i <;> fin_cases j <;> decide

/-- 洛伦兹型的复结构模型：`J^2 = -1`（6×6 实矩阵）。 -/
def zb080J : Matrix (Fin 6) (Fin 6) ℤ :=
  !![0, -1, 0, 0, 0, 0;
     1, 0, 0, 0, 0, 0;
     0, 0, 0, -1, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, 0, 0, 0, 0, -1;
     0, 0, 0, 0, 1, 0]

/-- **Zb080 定理 4.1（洛伦兹型）**：`J^2 = -1`。 -/
theorem zb080_J_sq : zb080J * zb080J = -1 := by
  ext i j <;> fin_cases i <;> fin_cases j <;> decide

#print axioms zb080_star_sq
#print axioms zb080_J_sq

end Zero
