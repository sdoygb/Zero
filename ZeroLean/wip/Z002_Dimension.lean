/-
  Z002 §6.7 维数账的第一段（机器核验版，设计稿）。

  文档 §6.7 把"闭合＝运动"写成计数读法 `dim ker B = dim im B`，并核验"两者相等仅当 m = 4"。
  本文件先把其中可机器核验的**骨架**钉住：

  * `z002B`：关联矩阵（顶点 × 边），边 `{a<b}` 在 a 处 `−1`、在 b 处 `+1`；
  * `z002T`：关联映射（边函数 ↦ 顶点函数）；
  * `z002ZeroSum` = 运动空间 `H_Q = {Σ x = 0}`；
  * `z002B_sum_col`：**列和为零**（每条边的两端各 ∓1）；
  * `z002T_range_le_zeroSum`：**像 ⊆ H_Q**（列和为零 + 求和交换）。

  ~~待补~~：`z002ZeroSum_finrank`（H_Q 维数 = card V − 1）、
  `z002T_ker_eq_zero`（连通图上 T 单射）、`z002_dim_eq_iff`（主定理）。
  这些的完整草稿见 `wip/Z002_Dimension_probe.lean`（尚欠 Sym2 端点的簿记）。
-/

import Mathlib
import Zero.Z064_Laplacian

open scoped Matrix

namespace Zero

variable {V : Type*} [Fintype V] [DecidableEq V] [LinearOrder V]
variable {G : SimpleGraph V} [DecidableRel G.Adj]

/-- 关联矩阵（顶点 × 边）：`B v e = ±1`，边 `{a<b}` 在 `a` 处 `−1`、在 `b` 处 `+1`。 -/
noncomputable def z002B (G : SimpleGraph V) [DecidableRel G.Adj] :
    Matrix V {e // e ∈ G.edgeSet} ℝ :=
  Matrix.of fun v e =>
    (if (e.val : Sym2 V).inf = v then (-1 : ℝ) else 0) +
      (if (e.val : Sym2 V).sup = v then (1 : ℝ) else 0)

/-- 关联映射 `T`（边函数 ↦ 顶点函数）：`(T y)_v = Σ_e B v e · y e`（流量取边界）。 -/
noncomputable def z002T (G : SimpleGraph V) [DecidableRel G.Adj] :
    ({e // e ∈ G.edgeSet} → ℝ) →ₗ[ℝ] (V → ℝ) where
  toFun y := fun v => ∑ e, z002B G v e * y e
  map_add' x y := by
    ext v
    simp only [Pi.add_apply]
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun e _ => by ring
  map_smul' c x := by
    ext v
    simp only [Pi.smul_apply, smul_eq_mul, RingHom.id_apply]
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun e _ => by ring

/-- 圈空间（闭合回路空间）`= ker T`。 -/
noncomputable def z002Closure (G : SimpleGraph V) [DecidableRel G.Adj] :
    Submodule ℝ ({e // e ∈ G.edgeSet} → ℝ) :=
  LinearMap.ker (z002T G)

/-- 运动空间（零和超平面）`H_Q = {x : Σ_v x_v = 0}`。 -/
noncomputable def z002ZeroSum : Submodule ℝ (V → ℝ) :=
  LinearMap.ker
    { toFun := fun x : V → ℝ => ∑ v, x v
      map_add' := fun x y => Finset.sum_add_distrib
      map_smul' := fun c x => by
        show (∑ v : V, (c • x) v) = (RingHom.id ℝ) c • ∑ v, x v
        rw [Finset.smul_sum]
        exact Finset.sum_congr rfl fun v _ => rfl }

private lemma ite_eq_comm (a v : V) (c : ℝ) :
    (if a = v then c else 0) = if v = a then c else 0 := by
  by_cases h : a = v
  · rw [ite_eq_left h, ite_eq_left h.symm]
  · rw [ite_eq_right h, ite_eq_right (Ne.symm h)]

/-- **列和为零**：每一列（一条边）的两端各贡献 `−1` 与 `+1`。 -/
theorem z002B_sum_col (G' : SimpleGraph V) [DecidableRel G'.Adj]
    (e : {e // e ∈ G'.edgeSet}) :
    ∑ v, z002B G' v e = 0 := by
  classical
  simp only [z002B, Matrix.of_apply, Finset.sum_add_distrib]
  rw [Finset.sum_congr rfl fun v _ => ite_eq_comm _ _ (-1 : ℝ),
    Finset.sum_congr rfl fun v _ => ite_eq_comm _ _ (1 : ℝ)]
  rw [show (∑ v : V, (if v = ((e.val : Sym2 V).inf : V) then (-1 : ℝ) else 0))
      = -1 from by
        rw [Finset.sum_ite_eq']
        simp]
  rw [show (∑ v : V, (if v = ((e.val : Sym2 V).sup : V) then (1 : ℝ) else 0))
      = 1 from by
        rw [Finset.sum_ite_eq']
        simp]
  norm_num

/-- **像落在运动空间里**：`T` 的像 ⊆ `H_Q`。 -/
theorem z002T_range_le_zeroSum :
    LinearMap.range (z002T G) ≤ z002ZeroSum := by
  rintro x ⟨y, rfl⟩
  rw [z002ZeroSum, LinearMap.mem_ker]
  show (∑ v : V, (∑ e, z002B G v e * y e)) = 0
  rw [Finset.sum_comm]
  refine Finset.sum_eq_zero fun e _ => ?_
  rw [← Finset.sum_mul, z002B_sum_col G e, zero_mul]

end Zero
