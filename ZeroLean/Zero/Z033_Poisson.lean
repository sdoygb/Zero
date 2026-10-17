/-
  Z033（《引力的观测对应：离散泊松方程成立，而牛顿极限与引力波无载体》）的机器检查版本：
  **连通图上 Poisson 方程 `L φ = −x` 的可解性 ⟺ 零和，且解模常数唯一**（Z033 §2 定理 16）。

  对应 Z033 的哪一句：
  * `z033_sum_zero_of_solvable` ↔ Z033 §2 的**必要性**：`L = D − A` 的行和恒为 0，
    故 `Σ_v (L φ)_v = 0`（本文 `z033_sum_mulVec`），于是 `L φ = −x` 逼出 `Σ x = 0`：
    零和是**可解的必要条件**（Z033 式 (1.1) 的符号约定 `L = −Δ`、`Lφ = −x`）。
  * `z033_diff_constant` ↔ Z033 §2 的**唯一性**：两解之差 `d = φ − ψ` 满足 `L d = 0`，
    连通图上 `ker L = span{1}`（Z064 已机器核验），故 `d` 为常数：解在**模常数**下唯一。
  * `z033_solvable_of_sum_zero` ↔ Z033 §2 的**充分性（Fredholm 条件自动满足）**：
    原文写 `im L = (ker L)^⊥ = {1}^⊥ = H_Q`，`x ∈ H_Q ⟹ x ∈ im L`。
    这里用等价而更便于机器化的形式：`L` 保持零和超平面 `H_Q`（本文 `z033_sum_mulVec`），
    且 `L|_{H_Q}` **单射**（`ker L ∩ H_Q = {0}`，即 Z064 的 `z063_eq_zero_of_sum_zero`），
    有限维下 `H_Q → H_Q` 的单射自映射满射，故 `−x ∈ H_Q = im L|_{H_Q}`。

  用到的 Mathlib 引理（注释中逐条标注）：
  * `SimpleGraph.lapMatrix_mulVec_one_eq_zero`：`L · 1 = 0`，即 **L 的行和为 0**。
  * `Matrix.dotProduct_mulVec`、`Matrix.mulVec_transpose`：
    `1 ⬝ (L f) = (1 ᵥ* L) ⬝ f = (Lᵀ · 1) ⬝ f`，配 `G.isSymm_lapMatrix`（`Lᵀ = L`）。
  * `SimpleGraph.lapMatrix_mulVec_apply'`：`(L f) v = Σ_{w ∈ N(v)} (f v − f w)`。
  * `SimpleGraph.lapMatrix_mulVec_eq_zero_iff_forall_reachable`：
    `L f = 0 ↔ ∀ i j, G.Reachable i j → f i = f j`（连通性由此进入）。
  * `LinearMap.injective_iff_surjective`：有限维空间的自映射单射 ⟺ 满射
    （这就是"Fredholm 条件"在这里的机器形式）。
  * `LinearMap.mem_ker`：`f ∈ ker g ↔ g f = 0`。
  * 复用的既有机器核验：`Zero.Z063_Laplacian` 的 `z063_eq_zero_of_sum_zero`
    （`ker L ∩ H_Q = {0}`，即 `L|_{H_Q} ≻ 0`）。
-/

import Mathlib
import Zero.Z063_Laplacian

open scoped Matrix

namespace Zero

/-- **辅助（Z033 §2）· 行和为 0 ⟹ 分量的总和为 0**：
`Σ_v (L f)_v = 0`，即 `1 ⬝ (L f) = 0`。
这是 Z033 §2 "可解性 ⟹ 零和" 的全部内容（`L φ = −x` 时左边为 `Σ(−x)`）。 -/
private theorem z033_sum_mulVec {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] (f : V → ℝ) :
    ∑ v, (G.lapMatrix ℝ *ᵥ f) v = 0 := by
  -- `Σ = 1 ⬝ (L f)`，再用 `dotProduct_mulVec` 把左边的 1 移过去：`1 ⬝ (L f) = (1 ᵥ* L) ⬝ f`
  have h1 : (∑ v, (G.lapMatrix ℝ *ᵥ f) v) = (1 : V → ℝ) ⬝ᵥ (G.lapMatrix ℝ *ᵥ f) := by
    simp [dotProduct]
  -- `1 ᵥ* L = Lᵀ · 1 = L · 1 = 0`（`L` 对称 ＋ 行和为 0）
  rw [h1, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose,
    Matrix.IsSymm.eq (SimpleGraph.isSymm_lapMatrix ℝ G),
    SimpleGraph.lapMatrix_mulVec_one_eq_zero ℝ G]
  -- 剩下 `(0 : V → ℝ) ⬝ f = 0`
  simp [dotProduct]

/-- **辅助（Z033 §2）· 零和线性泛函** `f ↦ Σ_v f v`，其核即 Z033 的零和超平面 `H_Q`。 -/
def z033SumLM {V : Type*} [Fintype V] : (V → ℝ) →ₗ[ℝ] ℝ where
  toFun f := ∑ v, f v
  map_add' f g := by simp [Finset.sum_add_distrib]
  map_smul' c f := by simp [Finset.mul_sum]

/-- **辅助（Z033 §2）· 零和超平面** `H_Q = {x | Σ_v x v = 0}`（Z033 的"零和约束"）。 -/
def z033ZeroSum {V : Type*} [Fintype V] : Submodule ℝ (V → ℝ) :=
  LinearMap.ker (z033SumLM (V := V))

private lemma z033_mem_zeroSum {V : Type*} [Fintype V] (f : V → ℝ) :
    f ∈ z033ZeroSum (V := V) ↔ ∑ v, f v = 0 :=
  LinearMap.mem_ker

/-- **辅助（Z033 §2）· `L` 保持零和超平面**：`x ∈ H_Q ⟹ L x ∈ H_Q`。 -/
lemma z033_mulVec_mem_zeroSum {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] (f : V → ℝ) :
    G.lapMatrix ℝ *ᵥ f ∈ z033ZeroSum (V := V) := by
  rw [z033_mem_zeroSum]
  exact z033_sum_mulVec G f

/-- **辅助（Z033 §2）· 零和超平面上的 Laplacian 自映射** `L|_{H_Q} : H_Q → H_Q`。 -/
def z033LapOnZeroSum {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] : z033ZeroSum (V := V) →ₗ[ℝ] z033ZeroSum (V := V) where
  toFun f := ⟨G.lapMatrix ℝ *ᵥ (f : V → ℝ), z033_mulVec_mem_zeroSum G (f : V → ℝ)⟩
  map_add' f g := by
    apply Subtype.ext
    simp only [Submodule.coe_add, Matrix.mulVec_add]
  map_smul' c f := by
    apply Subtype.ext
    simp only [Submodule.coe_smul, Matrix.mulVec_smul, RingHom.id_apply]

@[simp] private lemma z033LapOnZeroSum_apply {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] (f : z033ZeroSum (V := V)) :
    ((z033LapOnZeroSum G f : z033ZeroSum (V := V)) : V → ℝ)
      = G.lapMatrix ℝ *ᵥ (f : V → ℝ) := rfl

/-- **Z033 · 可解 ⟹ 零和**：行和为 0 ⟹ `Σx = 0` -/
theorem z033_sum_zero_of_solvable {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] (x : V → ℝ)
    (h : ∃ φ : V → ℝ, G.lapMatrix ℝ *ᵥ φ = -x) : ∑ v, x v = 0 := by
  obtain ⟨φ, hφ⟩ := h
  -- 对 `L φ = −x` 两边求和，左边由 `z033_sum_mulVec` 为零
  have h0 : ∑ v, (G.lapMatrix ℝ *ᵥ φ) v = 0 := z033_sum_mulVec G φ
  rw [hφ] at h0
  simpa using h0

/-- **Z033 · 解模常数唯一**：连通图上两个解的差是常数 -/
theorem z033_diff_constant {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] (hconn : G.Connected) (x φ ψ : V → ℝ)
    (hφ : G.lapMatrix ℝ *ᵥ φ = -x) (hψ : G.lapMatrix ℝ *ᵥ ψ = -x) :
    ∀ u v : V, φ u - ψ u = φ v - ψ v := by
  -- 两解之差落在 `ker L` 里（线性性）
  have hzero : G.lapMatrix ℝ *ᵥ (φ - ψ) = 0 := by
    rw [Matrix.mulVec_sub, hφ, hψ, sub_self]
  intro u v
  -- 连通图：`ker L = 常数`（Z064 定理 28）
  have hconst :=
    (SimpleGraph.lapMatrix_mulVec_eq_zero_iff_forall_reachable (G := G)).mp hzero u v (hconn u v)
  simpa [Pi.sub_apply] using hconst

/-- **Z033 · 连通 + 零和 ⟹ 可解** -/
theorem z033_solvable_of_sum_zero {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] (hconn : G.Connected) (x : V → ℝ) (hx : ∑ v, x v = 0) :
    ∃ φ : V → ℝ, G.lapMatrix ℝ *ᵥ φ = -x := by
  -- ① `L|_{H_Q}` 单射：`ker L ∩ H_Q = {0}`（Z064 的 `z063_eq_zero_of_sum_zero`）。
  --    这一步正是 Z033 §2 的 `L|_{H_Q} ≻ 0`。
  have hinj : Function.Injective (z033LapOnZeroSum G) := by
    intro f g hfg
    -- 由 `L f = L g` 得 `L (f − g) = 0`
    have h1 : G.lapMatrix ℝ *ᵥ ((f : V → ℝ) - (g : V → ℝ)) = 0 := by
      rw [Matrix.mulVec_sub, ← z033LapOnZeroSum_apply G f, ← z033LapOnZeroSum_apply G g,
        sub_eq_zero]
      exact congrArg Subtype.val hfg
    -- `f, g` 都在 `H_Q` 里，故 `Σ (f − g) = 0`；再把 `L (f − g) = 0` 换成逐顶点形式
    have hsum : ∑ v, ((f : V → ℝ) - (g : V → ℝ)) v = 0 := by
      have hf : ∑ v, (f : V → ℝ) v = 0 := (z033_mem_zeroSum (f : V → ℝ)).mp f.2
      have hg : ∑ v, (g : V → ℝ) v = 0 := (z033_mem_zeroSum (g : V → ℝ)).mp g.2
      simp only [Pi.sub_apply, Finset.sum_sub_distrib, hf, hg, sub_zero]
    have hf0 : ∀ v, ∑ w ∈ G.neighborFinset v,
        (((f : V → ℝ) - (g : V → ℝ)) v - ((f : V → ℝ) - (g : V → ℝ)) w) = 0 := by
      intro v
      rw [← SimpleGraph.lapMatrix_mulVec_apply' G v ((f : V → ℝ) - (g : V → ℝ))]
      exact congrFun h1 v
    exact Subtype.ext (sub_eq_zero.mp (z063_eq_zero_of_sum_zero G hconn _ hf0 hsum))
  -- ② 有限维下 `H_Q → H_Q` 的单射自映射满射 —— 这就是"Fredholm 条件自动满足"
  have hsurj : Function.Surjective (z033LapOnZeroSum G) :=
    LinearMap.injective_iff_surjective.mp hinj
  -- ③ 零和的 `−x` 在 `H_Q` 里，故是 `L|_{H_Q}` 的像
  have hxS : (-x) ∈ z033ZeroSum (V := V) := by
    rw [z033_mem_zeroSum]
    simpa using hx
  obtain ⟨φ, hφ⟩ := hsurj ⟨-x, hxS⟩
  refine ⟨(φ : V → ℝ), ?_⟩
  have hcoe : ((z033LapOnZeroSum G φ : z033ZeroSum (V := V)) : V → ℝ) = -x := by
    rw [hφ]
  simpa using hcoe

#print axioms z033_mulVec_mem_zeroSum
#print axioms z033_sum_zero_of_solvable
#print axioms z033_diff_constant
#print axioms z033_solvable_of_sum_zero

end Zero