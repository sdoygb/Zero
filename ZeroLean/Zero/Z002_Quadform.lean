/-
  Z002 的机器检查版本：**图 Laplacian 的二次型**三条基本性质 ——
  非负性、零点刻画（连通图上 ⟺ 常向量）、像落在零和超平面里（行和为 0）。

  对应 Z002 里的哪几句（Z002《从 Zero 条件退化出广义相对论》§3.4 引理 3，原文写的是
  `E(f) = Σ_e K_e (f_i − f_j)² ≥ 0`、`E(f) = 0 ⟺ f 为常量`、`E = fᵀ L f`、
  `ker L = span{1}`、`L(H_Q) ⊆ H_Q`、`L|_{H_Q} ≻ 0`；本文件是**单位导纳** `K_e ≡ 1`
  的特例，即 Mathlib 的无权图 Laplacian `L = D − A`）：
  * `z002_energy_nonneg` ↔ 引理 3 的 `E(f) ≥ 0`（`E = fᵀ L f` 就是这里的 `f ⬝ᵥ (L f)`）。
  * `z002_energy_eq_zero_iff` ↔ 引理 3 的 `E(f) = 0 ⟺ f 为常量`（等价写法 `ker L = span{1}`）。
    **连通性在这里用上**：先把"平方和为零"压到每条边上相等，再沿 Walk 传播到全图；
    不连通时结论为假（每个连通分支各取一个常数即可），故 `hconn` 是承重假设
    （已用机器反例核验：无边图 `⊥` 上取非常数 `f`，能量恒为 0 而 `f` 非常数）。
  * `z002_mulVec_mem_zeroSum` ↔ 引理 3 的 `L(H_Q) ⊆ H_Q`，`H_Q = {x : Σ_v x_v = 0}`：
    机器版本证的是更强的"`L` 的每一行和为 0"（由 `L` 对称即得列和也为 0），
    取 `Σ_v x_v = 0` 即得 `L(H_Q) ⊆ H_Q`。
    （`im L = H_Q` 这种满射性在 Z002 引理 2 是对关联矩阵 `B` 说的 `im B = H_Q`，
    本文件不涉及；`ker L = 常量` 与"行和为 0"合起来即 Z064 的横截结论
    `z063_eq_zero_of_sum_zero`。）

  与 Z064 的关系：`Zero.Z063_Laplacian` 的 `z063_quadform` 把同一个二次型写成
  "逐顶点邻居和"的形式 `Σ_v f v · (Σ_{w ∈ N(v)} (f v − f w)) = Σ_v Σ_{w∈N(v)} (f v − f w)²/2`；
  本文件用 Mathlib 的 `SimpleGraph.lapMatrix_toLinearMap₂'`（"逐边差平方和 / 2"，
  同上式的右手边）与 `SimpleGraph.lapMatrix_toLinearMap₂'_apply'_eq_zero_iff_forall_reachable`
  给出同一条恒等式的矩阵写法，两处结论因此可以互相对照。

  用到的 Mathlib 引理（注释中逐条标注）：
  * `Matrix.toLinearMap₂'_apply'`：`toLinearMap₂' ℝ L f f = f ⬝ᵥ (L *ᵥ f)`
    —— 把题设里的"点积"与双线性型对接起来。
  * `SimpleGraph.lapMatrix_toLinearMap₂'`：`f ⬝ᵥ (L f) = (Σ_{i,j, i~j} (f i − f j)²)/2`
    —— 非负性只是"平方和 / 2"。
  * `SimpleGraph.lapMatrix_toLinearMap₂'_apply'_eq_zero_iff_forall_reachable`：
    `f ⬝ᵥ (L f) = 0 ↔ ∀ i j, i 可达 j → f i = f j`
    —— 内部正是"`sum_eq_zero_iff_of_nonneg` ⟹ 每项为零"＋Walk 归纳。
  * `SimpleGraph.lapMatrix_mulVec_one_eq_zero`：`L *ᵥ 1 = 0`（行和为 0）。
  * `SimpleGraph.isSymm_lapMatrix`：`Lᵀ = L`（把 `⬝ᵥ` 挪到左边的关键）。
  * `Matrix.dotProduct_mulVec`、`Matrix.vecMul_transpose`：`1 ⬝ᵥ (L x) = (Lᵀ *ᵥ 1) ⬝ᵥ x`。
-/

import Mathlib
import Zero.Z063_Laplacian

open scoped Matrix

namespace Zero

/-- **Z002 · 二次型非负**：`f ⬝ᵥ (L f) ≥ 0`（逐边差平方和的一半） -/
theorem z002_energy_nonneg {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] (f : V → ℝ) :
    0 ≤ f ⬝ᵥ (G.lapMatrix ℝ *ᵥ f) := by
  -- 先把点积认成双线性型 `toLinearMap₂' ℝ L f f`，再换成"平方和 / 2"
  rw [← Matrix.toLinearMap₂'_apply', SimpleGraph.lapMatrix_toLinearMap₂' ℝ G f]
  -- 剩下的就是"平方和的一半 ≥ 0"：外层和内层和式逐项非负
  refine div_nonneg (Finset.sum_nonneg fun i _ => Finset.sum_nonneg fun j _ => ?_) (by norm_num)
  by_cases h : G.Adj i j
  · -- 相邻项：`(f i − f j)² ≥ 0`
    simpa [h] using sq_nonneg (f i - f j)
  · -- 非相邻项：就是 `0 ≥ 0`
    simp [h]

/-- **Z002 · 零点刻画**：连通图上二次型为零 ⟺ `f` 是常数 -/
theorem z002_energy_eq_zero_iff {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] (hconn : G.Connected) (f : V → ℝ) :
    f ⬝ᵥ (G.lapMatrix ℝ *ᵥ f) = 0 ↔ ∀ u v : V, f u = f v := by
  -- 同一条恒等式把左边认成 `toLinearMap₂'`，再用 Mathlib 已证的
  -- "零点 ⟺ 可达处相等"（其内部即"平方和为零 ⟹ 每边相等"＋Walk 归纳）
  rw [← Matrix.toLinearMap₂'_apply',
    SimpleGraph.lapMatrix_toLinearMap₂'_apply'_eq_zero_iff_forall_reachable ℝ G f]
  constructor
  · -- 正向：可达处相等 ＋ 连通性（每对顶点都可达）⟹ 全域常数
    intro h u v
    exact h u v (hconn u v)
  · -- 反向：常向量 ⟹ 任意两顶点（特别地，可达的两点）相等
    intro h i j _
    exact h i j

/-- **Z002 · 像落在零和超平面里**：`Σ_v (L x) v = 0`（行和为 0） -/
theorem z002_mulVec_mem_zeroSum {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] (x : V → ℝ) :
    ∑ v, (G.lapMatrix ℝ *ᵥ x) v = 0 := by
  -- 求和 = 与常向量 `1` 的点积（`one_dotProduct`：`1 ⬝ᵥ y = Σ y`）
  have h1 : (∑ v, (G.lapMatrix ℝ *ᵥ x) v) = (1 : V → ℝ) ⬝ᵥ (G.lapMatrix ℝ *ᵥ x) :=
    (one_dotProduct _).symm
  -- 把 `⬝ᵥ` 挪到左边：`1 ⬝ᵥ (L x) = (1 ᵥ* L) ⬝ᵥ x = (1 ᵥ* Lᵀ) ⬝ᵥ x = (L *ᵥ 1) ⬝ᵥ x`（用了对称性 `Lᵀ = L`）
  rw [h1, Matrix.dotProduct_mulVec, ← SimpleGraph.isSymm_lapMatrix ℝ G,
    Matrix.vecMul_transpose (G.lapMatrix ℝ) (1 : V → ℝ),
    SimpleGraph.lapMatrix_mulVec_one_eq_zero ℝ]
  -- 剩下 `0 ⬝ᵥ x = 0`（`zero_dotProduct`）
  simp

#print axioms z002_energy_nonneg
#print axioms z002_energy_eq_zero_iff
#print axioms z002_mulVec_mem_zeroSum

end Zero