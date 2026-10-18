/-
  Z002 §6.7（并 §6.8 的实例）：**维数账**。

  文档 §6.7 把"闭合＝运动"写成计数读法
    `dim ker B = dim im B`（闭合回路自由度 = 运动自由度），并核验"两者相等**仅当** m = 4"。
  §6.8 又用同一读法核了五个实例（$K_4$、$C_4$、$4\times4$ 方格、Heawood、$K_4\square C_4$）。

  记号钉死（本文件落实 §6.7 的那句限定）：`B` 取 **1-骨架的关联映射**（边流量 ↦ 顶点边界），于是
    `ker B = Z_1 ⊂ ℝ^E`（圈空间）、`im B = H_Q ⊂ ℝ^V`（零和超平面）：**两端不在同一空间**，
  读法比的是这两个数。本文件给出：

  * `z002B`：关联矩阵（顶点 × 边），边 `{a<b}` 在 a 处 `−1`、在 b 处 `+1`；
  * `z002T = B.mulVecLin`：关联映射（边函数 ↦ 顶点函数），`z002Closure = ker T` = 圈空间；
  * `z002ZeroSum` = `H_Q = {Σ x = 0}`；
  * `z002B_sum_col`：列和为零；`z002T_range_le_zeroSum`：像 ⊆ H_Q；
  * `z002ZeroSum_finrank`：`finrank H_Q = card V − 1`；
  * `z002B_sup_eq_one`／`z002B_inf_eq_neg_one`：端点处列系数；
  * `z002B_mul_eq_indicator_diff`：逐边分解 `B v e · y e = (sup=v) − (inf=v)`；
  * `z002T_sup_eq_inf`：核条件 `T y = 0` ⟹ `Σ_{sup=v} y = Σ_{inf=v} y`（端点流量守恒）。

  说明（诚实边界）：本篇是**维数账的输入层**——定义与端点分解已证；主定理
  `finrank (圈空间) = finrank (range T) ⟺ m = 4` 还需要"关联映射的秩 = card V − 1"
  （连通图），Mathlib 的 `IncMatrix` 未提供该秩结论，故尚未在此形式化。

  口径注（承重）："闭合回路空间 = 关联映射的核"只在**实数系数**下成立（整数系数下还有挠元），
  与 Zero 全部取实系数一致。
-/

import Mathlib
import Zero.Z063_Laplacian

open scoped Matrix

namespace Zero

variable {V : Type*} [Fintype V] [DecidableEq V] [LinearOrder V]
variable {G : SimpleGraph V} [DecidableRel G.Adj]

/-! ## 一、关联矩阵与关联映射 -/

/-- 关联矩阵（顶点 × 边）：`B v e = ±1`，边 `{a<b}` 在 `a` 处 `−1`、在 `b` 处 `+1`。 -/
noncomputable def z002B (G : SimpleGraph V) [DecidableRel G.Adj] :
    Matrix V {e // e ∈ G.edgeSet} ℝ :=
  Matrix.of fun v e =>
    (if (e.val : Sym2 V).inf = v then (-1 : ℝ) else 0) +
      (if (e.val : Sym2 V).sup = v then (1 : ℝ) else 0)

/-- 关联映射 `T = B ·`（边函数 ↦ 顶点函数）：`(T y)_v = Σ_e B v e · y e`（流量的边界）。 -/
noncomputable def z002T (G : SimpleGraph V) [DecidableRel G.Adj] :
    ({e // e ∈ G.edgeSet} → ℝ) →ₗ[ℝ] (V → ℝ) :=
  (z002B G).mulVecLin

/-- 圈空间（闭合回路空间）`= ker T`。 -/
noncomputable def z002Closure (G : SimpleGraph V) [DecidableRel G.Adj] :
    Submodule ℝ ({e // e ∈ G.edgeSet} → ℝ) :=
  LinearMap.ker (z002T G)

/-- 和泛函 `Σ : x ↦ Σ_v x_v`。 -/
noncomputable def z002Sum : (V → ℝ) →ₗ[ℝ] ℝ :=
  { toFun := fun x => ∑ v, x v
    map_add' := fun x y => Finset.sum_add_distrib
    map_smul' := fun c x => by
      show (∑ v : V, (c • x) v) = (RingHom.id ℝ) c • ∑ v, x v
      rw [Finset.smul_sum]
      exact Finset.sum_congr rfl fun v _ => rfl }

/-- 运动空间（零和超平面）`H_Q = ker Σ = {x : Σ_v x_v = 0}`。 -/
noncomputable def z002ZeroSum : Submodule ℝ (V → ℝ) := LinearMap.ker (z002Sum (V := V))

/-- `Σ` 的逐点写法。 -/
theorem z002Sum_apply (x : V → ℝ) : z002Sum x = ∑ v, x v := rfl

/-! ## 二、基本引理 -/

private lemma ite_eq_comm (a v : V) (c : ℝ) :
    (if a = v then c else 0) = if v = a then c else 0 := by
  by_cases h : a = v
  · rw [ite_eq_left h, ite_eq_left h.symm]
  · rw [ite_eq_right h, ite_eq_right (Ne.symm h)]

/-- `inf = sup` 的对称二元组必是对角元（边集不含对角元，用得上）。 -/
private lemma isDiag_of_inf_eq_sup (z : Sym2 V)
    (h : ((z.inf : V) : V) = ((z.sup : V) : V)) : z.IsDiag := by
  induction z using Sym2.ind with
  | _ a b =>
    rw [Sym2.mk_isDiag_iff]
    rw [Sym2.inf_mk, Sym2.sup_mk] at h
    rcases le_total a b with hle | hle
    · rw [inf_eq_left.mpr hle, sup_eq_right.mpr hle] at h
      exact h
    · rw [inf_eq_right.mpr hle, sup_eq_left.mpr hle] at h
      exact h.symm

/-- `Sym2` 的 `sup` 落在边上。 -/
private lemma sup_mem_edge (e : {e // e ∈ G.edgeSet}) :
    ((e.val : Sym2 V).sup : V) ∈ (e.val : Sym2 V) := by
  obtain ⟨z, hz⟩ := e
  induction z using Sym2.ind with
  | _ a b =>
    rw [Sym2.sup_mk]
    rcases le_total a b with h | h
    · rw [sup_eq_right.mpr h]
      exact Sym2.mem_iff'.mpr (Or.inr rfl)
    · rw [sup_eq_left.mpr h]
      exact Sym2.mem_iff'.mpr (Or.inl rfl)

/-- 边 `s(a,b) ∈ edgeSet` 给出 `G.Adj a b`。 -/
private lemma adj_of_mem {a b : V} (hz : (s(a, b) : Sym2 V) ∈ G.edgeSet) : G.Adj a b := by
  simpa using hz

/-- 不是端点的顶点处系数为零。 -/
private lemma z002B_eq_zero (e : {e // e ∈ G.edgeSet}) {v : V}
    (h1 : v ≠ ((e.val : Sym2 V).inf : V)) (h2 : v ≠ ((e.val : Sym2 V).sup : V)) :
    z002B G v e = 0 := by
  rw [z002B, Matrix.of_apply, ite_eq_right (fun hc => h1 hc.symm),
    ite_eq_right (fun hc => h2 hc.symm)]
  norm_num

/-- 单边系数：边 `{a,b}`（`a < b`）在 `a` 处 `−1`、在 `b` 处 `+1`。 -/
private lemma z002B_inf_sup (a b : V) (hab : G.Adj a b) (hlt : a < b) :
    z002B G a ⟨s(a, b), by simpa using hab⟩ = -1 ∧
      z002B G b ⟨s(a, b), by simpa using hab⟩ = 1 := by
  have hinf : ((s(a, b) : Sym2 V).inf : V) = a := by rw [Sym2.inf_mk, inf_eq_left.mpr hlt.le]
  have hsup : ((s(a, b) : Sym2 V).sup : V) = b := by rw [Sym2.sup_mk, sup_eq_right.mpr hlt.le]
  constructor
  · rw [z002B, Matrix.of_apply, hinf, hsup, ite_eq_left rfl,
      ite_eq_right (ne_of_lt hlt).symm]
    norm_num
  · rw [z002B, Matrix.of_apply, hinf, hsup, ite_eq_right (ne_of_lt hlt), ite_eq_left rfl]
    norm_num

/-! ## 三、列和为 0 ⟹ 像落在运动空间里 -/

/-- **列和为零**：每一列（一条边）的两端各贡献 `−1` 与 `+1`。 -/
theorem z002B_sum_col (G' : SimpleGraph V) [DecidableRel G'.Adj]
    (e : {e // e ∈ G'.edgeSet}) :
    ∑ v, z002B G' v e = 0 := by
  classical
  simp only [z002B, Matrix.of_apply, Finset.sum_add_distrib]
  rw [Finset.sum_congr rfl fun v _ => ite_eq_comm _ _ (-1 : ℝ),
    Finset.sum_congr rfl fun v _ => ite_eq_comm _ _ (1 : ℝ)]
  rw [show (∑ v : V, (if v = ((e.val : Sym2 V).inf : V) then (-1 : ℝ) else 0)) = -1 from by
        rw [Finset.sum_ite_eq']; simp,
      show (∑ v : V, (if v = ((e.val : Sym2 V).sup : V) then (1 : ℝ) else 0)) = 1 from by
        rw [Finset.sum_ite_eq']; simp]
  norm_num

/-- **像落在运动空间里**：`range T ≤ H_Q`。 -/
theorem z002T_range_le_zeroSum :
    LinearMap.range (z002T G) ≤ z002ZeroSum := by
  rintro x ⟨y, rfl⟩
  rw [z002ZeroSum, LinearMap.mem_ker, z002Sum_apply]
  show (∑ v : V, (∑ e, z002B G v e * y e)) = 0
  rw [Finset.sum_comm]
  refine Finset.sum_eq_zero fun e _ => ?_
  rw [← Finset.sum_mul, z002B_sum_col G e, zero_mul]

/-! ## 四、运动空间的维数 -/

/-- **运动空间维数**：`finrank H_Q = card V − 1`。 -/
theorem z002ZeroSum_finrank (hV : 1 ≤ Nat.card V) :
    Module.finrank ℝ (z002ZeroSum (V := V)) = Nat.card V - 1 := by
  have hcard : Nat.card V = Fintype.card V := Nat.card_eq_fintype_card
  have hne : z002Sum (V := V) ≠ 0 := by
    intro hzero
    have hpos : 0 < Nat.card V := by omega
    haveI : Nonempty V := by
      haveI := Fintype.ofFinite V
      rw [Nat.card_eq_fintype_card] at hpos
      exact Fintype.card_pos_iff.mp hpos
    obtain ⟨i⟩ := ‹Nonempty V›
    have h3 : (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v) = 0 := by
      have h4 := DFunLike.congr_fun hzero (Pi.single i (1 : ℝ) : V → ℝ)
      rwa [z002Sum_apply] at h4
    have h5 : (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v) = 1 := by
      rw [show (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v)
          = ∑ v : V, (if v = i then (1 : ℝ) else 0) from
          Finset.sum_congr rfl fun v _ => by simp [Pi.single_apply, eq_comm]]
      rw [Finset.sum_ite_eq']
      simp
    rw [h5] at h3
    norm_num at h3
  have hdim := Module.Dual.finrank_ker_add_one_of_ne_zero hne
  rw [Module.finrank_fintype_fun_eq_card] at hdim
  rw [z002ZeroSum]
  omega

/-! ## 五、核条件的端点分解 -/

/-- `sup = v` 时该列系数为 `+1`。 -/
private lemma z002B_sup_eq_one (e : {e // e ∈ G.edgeSet}) {v : V}
    (hv : ((e.val : Sym2 V).sup : V) = v) : z002B G v e = 1 := by
  obtain ⟨e, hz⟩ := e
  induction e using Sym2.ind with
  | _ a b =>
    have hab : G.Adj a b := adj_of_mem hz
    rcases lt_trichotomy a b with hlt | heq | hgt
    · have hs : ((s(a, b) : Sym2 V).sup : V) = b := by
        rw [Sym2.sup_mk, sup_eq_right.mpr hlt.le]
      rw [hs] at hv
      rw [hv.symm]
      exact (z002B_inf_sup a b hab hlt).2
    · exact absurd heq (G.ne_of_adj hab)
    · have hs : ((s(a, b) : Sym2 V).sup : V) = a := by
        rw [Sym2.sup_mk, sup_eq_left.mpr hgt.le]
      rw [hs] at hv
      rw [hv.symm]
      have hba : (⟨s(b, a), by simpa using hab.symm⟩ : {e // e ∈ G.edgeSet}) =
          ⟨s(a, b), hz⟩ := Subtype.ext (Sym2.eq_swap)
      rw [← hba]
      exact (z002B_inf_sup b a hab.symm hgt).2

/-- `inf = v` 时该列系数为 `-1`。 -/
private lemma z002B_inf_eq_neg_one (e : {e // e ∈ G.edgeSet}) {v : V}
    (hv : ((e.val : Sym2 V).inf : V) = v) : z002B G v e = -1 := by
  obtain ⟨e, hz⟩ := e
  induction e using Sym2.ind with
  | _ a b =>
    have hab : G.Adj a b := adj_of_mem hz
    rcases lt_trichotomy a b with hlt | heq | hgt
    · have hi : ((s(a, b) : Sym2 V).inf : V) = a := by
        rw [Sym2.inf_mk, inf_eq_left.mpr hlt.le]
      rw [hi] at hv
      rw [hv.symm]
      exact (z002B_inf_sup a b hab hlt).1
    · exact absurd heq (G.ne_of_adj hab)
    · have hi : ((s(a, b) : Sym2 V).inf : V) = b := by
        rw [Sym2.inf_mk, inf_eq_right.mpr hgt.le]
      rw [hi] at hv
      rw [hv.symm]
      have hba : (⟨s(b, a), by simpa using hab.symm⟩ : {e // e ∈ G.edgeSet}) =
          ⟨s(a, b), hz⟩ := Subtype.ext (Sym2.eq_swap)
      rw [← hba]
      exact (z002B_inf_sup b a hab.symm hgt).1

/-- **逐边分解**：`B v e · y e = (sup = v 的指示) − (inf = v 的指示)`。 -/
private lemma z002B_mul_eq_indicator_diff (y : {e // e ∈ G.edgeSet} → ℝ) (v : V)
    (e : {e // e ∈ G.edgeSet}) :
    z002B G v e * y e =
      (if ((e.val : Sym2 V).sup : V) = v then y e else 0)
        - (if ((e.val : Sym2 V).inf : V) = v then y e else 0) := by
  by_cases hsup : ((e.val : Sym2 V).sup : V) = v
  · have hcoe : z002B G v e = 1 := z002B_sup_eq_one e hsup
    have hinf : ((e.val : Sym2 V).inf : V) ≠ v := by
      intro hc
      exact G.not_isDiag_of_mem_edgeSet e.prop (isDiag_of_inf_eq_sup _ (hc.trans hsup.symm))
    rw [hcoe, one_mul, if_pos hsup, if_neg hinf, sub_zero]
  · by_cases hinf : ((e.val : Sym2 V).inf : V) = v
    · have hcoe : z002B G v e = -1 := z002B_inf_eq_neg_one e hinf
      rw [hcoe, neg_one_mul, if_neg hsup, if_pos hinf, zero_sub]
    · have hcoe : z002B G v e = 0 := z002B_eq_zero e (fun hc => hinf hc.symm) (fun hc => hsup hc.symm)
      rw [hcoe, zero_mul, if_neg hsup, if_neg hinf, sub_self]

/-- **核条件的端点分解**：`(T y)_v = 0` ⟹ `Σ_{sup=v} y = Σ_{inf=v} y`。 -/
theorem z002T_sup_eq_inf (y : {e // e ∈ G.edgeSet} → ℝ) (h : z002T G y = 0) (v : V) :
    (∑ e : {e // e ∈ G.edgeSet}, (if ((e.val : Sym2 V).sup : V) = v then y e else 0))
      = ∑ e : {e // e ∈ G.edgeSet}, (if ((e.val : Sym2 V).inf : V) = v then y e else 0) := by
  have h1 : (∑ e : {e // e ∈ G.edgeSet}, z002B G v e * y e) = 0 := by
    have h2 := congrFun h v
    rw [show z002T G y = (z002B G) *ᵥ y from rfl, Matrix.mulVec_apply_eq_sum] at h2
    exact h2
  have hsplit : (∑ e : {e // e ∈ G.edgeSet}, z002B G v e * y e)
      = (∑ e : {e // e ∈ G.edgeSet}, (if ((e.val : Sym2 V).sup : V) = v then y e else 0))
        - ∑ e : {e // e ∈ G.edgeSet}, (if ((e.val : Sym2 V).inf : V) = v then y e else 0) := by
    rw [← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun e _ => z002B_mul_eq_indicator_diff y v e
  rw [hsplit] at h1
  linarith

end Zero