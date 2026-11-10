/-
  Zb078 · K4 的线图 = 八面体：线图邻接恒等式

  对应篇：Zb078 §3（定理 3.1）。
  只核验代数核：`A_L = Bᵀ B - 2 I`，其中
    · B 是 K4 的无向关联矩阵（4×6）；
    · A_L 是线图邻接矩阵（6×6，边相邻当且仅当共享顶点）。
  置换同构 L(K4) ≅ 八面体由脚本 `ZeroGPU/src/zb078_k4_line_graph.py` 核验。
-/

import Mathlib

namespace Zero

open Matrix

/-- K4 的 6 条边：01, 02, 03, 12, 13, 23。 -/
def zb078Edge : Fin 6 → Fin 4 × Fin 4
  | 0 => (0, 1)
  | 1 => (0, 2)
  | 2 => (0, 3)
  | 3 => (1, 2)
  | 4 => (1, 3)
  | 5 => (2, 3)

/-- K4 的无向关联矩阵 B（4×6）。 -/
def zb078Inc : Matrix (Fin 4) (Fin 6) ℤ :=
  fun v e => if (zb078Edge e).1 = v ∨ (zb078Edge e).2 = v then 1 else 0

/-- K4 的线图邻接矩阵（6×6）。 -/
def zb078LineAdj : Matrix (Fin 6) (Fin 6) ℤ :=
  fun e f =>
    if e = f then 0
    else if (zb078Edge e).1 = (zb078Edge f).1 ∨ (zb078Edge e).1 = (zb078Edge f).2 ∨
            (zb078Edge e).2 = (zb078Edge f).1 ∨ (zb078Edge e).2 = (zb078Edge f).2 then 1
    else 0

/-- **Zb078 定理 3.1（代数核）**：`A_L = Bᵀ B - 2 I`。 -/
theorem zb078_line_adj :
    zb078LineAdj = zb078Incᵀ * zb078Inc - (2 : ℤ) • (1 : Matrix (Fin 6) (Fin 6) ℤ) := by
  ext e f
  fin_cases e <;> fin_cases f <;> decide

#print axioms zb078_line_adj

end Zero
