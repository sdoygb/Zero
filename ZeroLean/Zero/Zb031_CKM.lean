/-
  Zb031 §3 的机器检查版本：CKM 矩阵的实参数计数。

  3×3 复矩阵 = 18 个实参数；幺正性 `U†U = 1` 给 9 个独立实约束，
  剩 9 个；再除去 5 个不可观测的相位（3 个夸克场 + 2 个相对相位），
  物理参数 = 4 = 3 个混合角 + 1 个 CP 相位。

  这是纯计数（正文把它记作「导出」）；本文件只把这些整数关系送进内核。
-/

import Mathlib

namespace Zero

/-- 3×3 复矩阵的实参数个数：`2 × 9 = 18`。 -/
theorem zb031_complex_matrix_real_dim : (2 * 9 : ℕ) = 18 := by norm_num

/-- 幺正性给出的独立实约束：`3 × 3 = 9`。 -/
theorem zb031_unitary_constraints : (3 * 3 : ℕ) = 9 := by norm_num

/-- 扣掉幺正约束后的实参数：`18 − 9 = 9`。 -/
theorem zb031_after_unitarity : (18 - 9 : ℕ) = 9 := by norm_num

/-- 再扣掉 5 个不可观测相位：`9 − 5 = 4`。 -/
theorem zb031_physical_parameters : (9 - 5 : ℕ) = 4 := by norm_num

/-- 4 个物理参数 = 3 个混合角 + 1 个 CP 相位。 -/
theorem zb031_three_angles_one_phase : (3 + 1 : ℕ) = 4 := by norm_num

#print axioms zb031_complex_matrix_real_dim
#print axioms zb031_unitary_constraints
#print axioms zb031_after_unitarity
#print axioms zb031_physical_parameters
#print axioms zb031_three_angles_one_phase

end Zero
