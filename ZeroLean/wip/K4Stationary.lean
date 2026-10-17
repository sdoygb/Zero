/-
  K4 的均匀测度是平稳的：$K_4$ 是 3-正则图，故每个顶点的入权重都是 $3/4$。

  这条把"四态等权（每态权重 $1/4$）"从**统计读法**升级为**定理**：
  均匀分布是 $K_4$ 上（无向、正则）行走的平稳测度，权重 $1/4$ 由正则性给出。
-/

import Mathlib

namespace Zero

/-- $K_4$ 的邻接关系：$v\neq w$ 时相连。 -/
def k4Adj (v w : Fin 4) : ℝ := if v ≠ w then 1 else 0

/-- **均匀测度平稳**：每个顶点从邻居收到的总权重 $=3/4$。 -/
theorem k4_uniform_stationary (v : Fin 4) :
    (∑ w : Fin 4, k4Adj v w * (1 / 4)) = 3 / 4 := by
  fin_cases v <;>
    simp [k4Adj, Fin.sum_univ_four] <;>
    norm_num

/-- **总权重为 1**：$4\times\frac14=1$。 -/
theorem k4_uniform_total : (∑ _v : Fin 4, (1 / 4 : ℝ)) = 1 := by
  simp [Fin.sum_univ_four]

#print axioms k4_uniform_stationary
#print axioms k4_uniform_total

end Zero
