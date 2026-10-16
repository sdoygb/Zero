/-
  Z002 §6.7 的**维数账**：关联矩阵的核与像 —— `dim (圈空间) = dim (运动空间) ⟺ m = 4`。
-/

import Mathlib
import Zero.Z064_Laplacian

open scoped Matrix

namespace Zero

variable {V : Type*} [Fintype V] [DecidableEq V] [LinearOrder V]
variable {G : SimpleGraph V} [DecidableRel G.Adj]

/-! ## 一、关联矩阵与关联映射 -/

/-- 关联矩阵（顶点 × 边）：`B v e = ±1`，边 `e = {a<b}` 在 `a` 处取 `−1`、在 `b` 处取 `+1`。 -/
noncomputable def z002B (G : SimpleGraph V) [DecidableRel G.Adj] :
    Matrix V {e // e ∈ G.edgeSet} ℝ :=
  Matrix.of fun v e =>
    (if (e.val : Sym2 V).inf = v then (-1 : ℝ) else 0) +
      (if (e.val : Sym2 V).sup = v then (1 : ℝ) else 0)

/-- 关联映射 `T`（边函数 ↦ 顶点函数）：`(T y)_v = Σ_e B v e · y e`，即**流量的边界**。 -/
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

/-- 圈空间（闭合回路空间）`= ker T`：边界为零的流量。 -/
noncomputable def z002Closure (G : SimpleGraph V) [DecidableRel G.Adj] :
    Submodule ℝ ({e // e ∈ G.edgeSet} → ℝ) :=
  LinearMap.ker (z002T G)

/-- 运动空间（零和超平面）`H_Q = {x : Σ_v x_v = 0}`。 -/
noncomputable def z002Sum : (V → ℝ) →ₗ[ℝ] ℝ :=
  { toFun := fun x => ∑ v, x v
    map_add' := fun x y => Finset.sum_add_distrib
    map_smul' := fun c x => by
      show (∑ v : V, (c • x) v) = (RingHom.id ℝ) c • ∑ v, x v
      rw [Finset.smul_sum]
      exact Finset.sum_congr rfl fun v _ => rfl }

noncomputable def z002ZeroSum : Submodule ℝ (V → ℝ) := LinearMap.ker (z002Sum (V := V))

/-! ## 二、基本引理 -/

private lemma ite_eq_comm (a v : V) (c : ℝ) :
    (if a = v then c else 0) = if v = a then c else 0 := by
  by_cases h : a = v
  · rw [if_pos h, if_pos h.symm]
  · rw [if_neg h, if_neg (Ne.symm h)]

private lemma sum_univ_ite_eq (x : V) (c : ℝ) :
    ∑ v : V, (if v = x then c else 0) = c := by
  rw [Finset.sum_ite_eq' Finset.univ x (fun _ => c)]
  simp

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

private lemma inf_ne_sup (e : {e // e ∈ G.edgeSet}) :
    ((e.val : Sym2 V).inf : V) ≠ ((e.val : Sym2 V).sup : V) :=
  fun h => G.not_isDiag_of_mem_edgeSet e.prop (isDiag_of_inf_eq_sup e.val h)

/-- 边 `s(a,b) ∈ edgeSet` 给出 `G.Adj a b`。 -/
private lemma adj_of_mem {a b : V} (hz : (s(a, b) : Sym2 V) ∈ G.edgeSet) : G.Adj a b := by
  simpa using hz

/-- 单边系数：边 `{a,b}`（`a < b`）在 `a` 处取 `−1`、在 `b` 处取 `+1`。 -/
private lemma z002B_inf_sup (a b : V) (hab : G.Adj a b) (hlt : a < b) :
    z002B G a ⟨s(a, b), by simpa using hab⟩ = -1 ∧
      z002B G b ⟨s(a, b), by simpa using hab⟩ = 1 := by
  have hinf : ((s(a, b) : Sym2 V).inf : V) = a := by rw [Sym2.inf_mk, inf_eq_left.mpr hlt.le]
  have hsup : ((s(a, b) : Sym2 V).sup : V) = b := by rw [Sym2.sup_mk, sup_eq_right.mpr hlt.le]
  constructor
  · rw [z002B, Matrix.of_apply, hinf, hsup, if_pos rfl, if_neg (ne_of_lt hlt).symm]
    norm_num
  · rw [z002B, Matrix.of_apply, hinf, hsup, if_neg (ne_of_lt hlt), if_pos rfl]
    norm_num

/-- 只有两端非零：`v` 既不是 `inf` 也不是 `sup` 时 `B v e = 0`。 -/
private lemma z002B_eq_zero (e : {e // e ∈ G.edgeSet}) {v : V}
    (h1 : v ≠ ((e.val : Sym2 V).inf : V)) (h2 : v ≠ ((e.val : Sym2 V).sup : V)) :
    z002B G v e = 0 := by
  rw [z002B, Matrix.of_apply, if_neg (Ne.symm h1), if_neg (Ne.symm h2)]
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
  rw [sum_univ_ite_eq, sum_univ_ite_eq]
  norm_num

/-- 每一列都在运动空间里。 -/
theorem z002B_col_mem_zeroSum (G' : SimpleGraph V) [DecidableRel G'.Adj]
    (e : {e // e ∈ G'.edgeSet}) :
    (z002B G').col e ∈ z002ZeroSum := by
  have hcol : (z002B G').col e = fun v => z002B G' v e := rfl
  rw [z002ZeroSum, LinearMap.mem_ker, hcol]
  show (∑ v, z002B G' v e) = 0
  exact z002B_sum_col G' e

/-- **像落在运动空间里**：`range T ≤ H_Q`。 -/
theorem z002T_range_le_zeroSum :
    LinearMap.range (z002T G) ≤ z002ZeroSum := by
  rintro x ⟨y, rfl⟩
  rw [z002ZeroSum, LinearMap.mem_ker]
  show (∑ v : V, (∑ e, z002B G v e * y e)) = 0
  change (∑ v : V, (∑ e, z002B G v e * y e)) = 0
  show ∑ v : V, (∑ e, z002B G v e * y e) = 0
  rw [Finset.sum_comm]
  refine Finset.sum_eq_zero fun e _ => ?_
  rw [← Finset.sum_mul, z002B_sum_col G e, zero_mul]

/-! ## 四、运动空间的维数 -/

/-- 运动空间的维数：`finrank H_Q = card V − 1`。 -/
theorem z002ZeroSum_finrank (hV : 2 ≤ Nat.card V) :
    Module.finrank ℝ (z002ZeroSum (V := V)) = Nat.card V - 1 := by
  have hpos : 0 < Nat.card V := by omega
  haveI : Nonempty V := by
    haveI := Fintype.ofFinite V
    rw [Nat.card_eq_fintype_card] at hpos
    exact Fintype.card_pos_iff.mp hpos
  obtain ⟨i⟩ := ‹Nonempty V›
  have hsingle : (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v) = 1 := by
    rw [show (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v)
        = ∑ v : V, (if v = i then (1 : ℝ) else 0) from
        Finset.sum_congr rfl fun v _ => by rw [Pi.single_apply]; by_cases h : v = i <;> simp [h]]
    rw [Finset.sum_ite_eq']
    simp
  have hne : z002Sum (V := V) ≠ 0 := by
    intro hzero
    have h1 : (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v) = 0 := by
      have h3 : (∑ v : V, (Pi.single i (1 : ℝ) : V → ℝ) v) = 0 :=
        DFunLike.congr_fun hzero (Pi.single i (1 : ℝ) : V → ℝ)
      simpa using h3
    rw [hsingle] at h1
    norm_num at h1
  have hdim := Module.Dual.finrank_ker_add_one_of_ne_zero hne
  rw [Module.finrank_fintype_fun_eq_card] at hdim
  have hcard : Nat.card V = Fintype.card V := Nat.card_eq_fintype_card
  rw [z002ZeroSum]
  omega

/-! ## 五、核（圈空间）：连通图上边界为零的流量为零 -/

/-- `sup` 端点处的核条件：以 `v` 为超端点的那些边上流量之和为零。 -/
private lemma z002T_sup_eq_zero (y : {e // e ∈ G.edgeSet} → ℝ) (h : z002T G y = 0) (v : V) :
    ∑ e : {e // e ∈ G.edgeSet}, (if ((e.val : Sym2 V).sup : V) = v then y e else 0) = 0 := by
  have hh := congrFun h v
  simp only [z002T, LinearMap.coe_mk, AddHom.coe_mk, Pi.zero_apply] at hh
  rw [← hh]
  refine Finset.sum_congr rfl fun e _ => ?_
  by_cases hv : ((e.val : Sym2 V).sup : V) = v
  · rw [ite_eq_left hv]
    have hcoe : z002B G v e = 1 := by
      obtain ⟨e, hz⟩ := e
      induction e using Sym2.ind with
      | _ a b =>
        have hab : G.Adj a b := adj_of_mem hz
        rcases lt_trichotomy a b with hlt | heq | hgt
        · have hvb : v = b := by
            have hs : ((s(a, b) : Sym2 V).sup : V) = b := by
              rw [Sym2.sup_mk, sup_eq_right.mpr hlt.le]
            rw [hs] at hv; exact hv.symm
          subst hvb
          exact (z002B_inf_sup a b hab hlt).2
        · exact absurd heq (G.ne_of_adj hab)
        · have hva : v = a := by
            have hs : ((s(a, b) : Sym2 V).sup : V) = a := by
              rw [Sym2.sup_mk, sup_eq_left.mpr hgt.le]
            rw [hs] at hv; exact hv.symm
          subst hva
          exact (z002B_inf_sup b a hab.symm hgt).2
    rw [hcoe, one_mul]
  · rw [ite_eq_right hv]
    have hcoe : z002B G v e = 0 := by
      by_contra hne
      obtain ⟨e, hz⟩ := e
      induction e using Sym2.ind with
      | _ a b =>
        have hab : G.Adj a b := adj_of_mem hz
        have hmem_sup : ((s(a, b) : Sym2 V).sup : V) ∈ (s(a, b) : Sym2 V) :=
          sup_mem_edge ⟨s(a, b), hz⟩
        rw [Sym2.mem_iff'] at hmem_sup
        rcases hmem_sup with h2 | h2
        · rw [h2] at hv
          exact hne (by
            rw [← h2]
            rcases lt_trichotomy a b with hlt | heq | hgt
            · exact (z002B_inf_sup a b hab hlt).1
            · exact absurd heq (G.ne_of_adj hab)
            · exact (z002B_inf_sup b a hab.symm hgt).2)
        · rw [h2] at hv
          exact hne (by
            rw [← h2]
            rcases lt_trichotomy a b with hlt | heq | hgt
            · exact (z002B_inf_sup a b hab hlt).2
            · exact absurd heq (G.ne_of_adj hab)
            · exact (z002B_inf_sup b a hab.symm hgt).1)
    rw [hcoe, zero_mul]

/-- **每条边上的流量为零**：`T y = 0` 且 `G` 连通 ⟹ `y e = 0`。 -/
theorem z002T_edge_eq_zero (hconn : G.Connected) (y : {e // e ∈ G.edgeSet} → ℝ)
    (h : z002T G y = 0) (e : {e // e ∈ G.edgeSet}) : y e = 0 := by
  -- 取该边的一个端点作为"超端点"
  set a : V := (e.val : Sym2 V).out.1 with ha
  set b : V := (e.val : Sym2 V).out.2 with hb
  have heq : (s(a, b) : Sym2 V) = (e.val : Sym2 V) := (e.val).out_eq
  rcases lt_trichotomy a b with hlt | heq' | hgt
  · have hsum := z002T_sup_eq_zero (G := G) y h b
    have hsub : (⟨s(a, b), by rw [heq]; exact e.prop⟩ : {e // e ∈ G.edgeSet}) = e := by
      exact Subtype.ext heq
    rw [Finset.sum_eq_single (⟨s(a, b), by rw [heq]; exact e.prop⟩ :
      {e // e ∈ G.edgeSet})] at hsum
    · rw [ite_eq_left (by rw [Sym2.sup_mk, sup_eq_right.mpr hlt.le])] at hsum
      rw [hsub] at hsum
      exact hsum
    · intro e' _ he'
      rw [ite_eq_right]
      intro hc
      -- 若 e' 的超端点也是 b，则 e' 与 {a,b} 是同一个无向对
      exact he' (Subtype.ext (by
        rw [← Sym2.mem_iff']
        rw [← hc]
        exact sup_mem_edge e'))
    · intro hmem; exact absurd (Finset.mem_univ _) hmem
  · exact absurd heq' (G.ne_of_adj (adj_of_mem (by rw [heq]; exact e.prop)))
  · have hsum := z002T_sup_eq_zero (G := G) y h a
    have hsub : (⟨s(a, b), by rw [heq]; exact e.prop⟩ : {e // e ∈ G.edgeSet}) = e := by
      exact Subtype.ext heq
    rw [Finset.sum_eq_single (⟨s(a, b), by rw [heq]; exact e.prop⟩ :
      {e // e ∈ G.edgeSet})] at hsum
    · rw [ite_eq_left (by rw [Sym2.sup_mk, sup_eq_left.mpr hgt.le])] at hsum
      rw [hsub] at hsum
      exact hsum
    · intro e' _ he'
      rw [ite_eq_right]
      intro hc
      exact he' (Subtype.ext (by
        rw [← Sym2.mem_iff']
        rw [← hc]
        exact sup_mem_edge e'))
    · intro hmem; exact absurd (Finset.mem_univ _) hmem

/-- 连通图上关联映射的单射性。 -/
theorem z002T_injective (hconn : G.Connected) : Function.Injective (z002T G) := by
  intro x y hxy
  have hz : z002T G (x - y) = 0 := by rw [map_sub, hxy, sub_self]
  have h0 : x - y = 0 := z002T_ker_eq_zero hconn (x - y) hz
  exact sub_eq_zero.mp h0

/-! ## 六、维数账 -/

/-- **像维数**：连通图上 `finrank (range T) = card V − 1`。 -/
theorem z002_finrank_range (hconn : G.Connected) (hV : 2 ≤ Nat.card V) :
    Module.finrank ℝ (LinearMap.range (z002T G)) = Nat.card V - 1 := by
  have hmono := Submodule.finrank_mono (z002T_range_le_zeroSum (G := G))
  rw [z002ZeroSum_finrank (V := V) hV] at hmono
  have hker : LinearMap.ker (z002T G) = ⊥ :=
    LinearMap.ker_eq_bot.mpr (z002T_injective hconn)
  have hrn := LinearMap.finrank_range_add_finrank_ker (z002T G)
  rw [hker, finrank_bot, Nat.add_zero] at hrn
  have hcard : Module.finrank ℝ ({e // e ∈ G.edgeSet} → ℝ)
      = Nat.card {e // e ∈ G.edgeSet} := by
    rw [Nat.card_eq_fintype_card, Module.finrank_fintype_fun_eq_card]
  rw [hcard] at hrn
  omega

/-- **核维数**：连通图上 `finrank (ker T) = card E − (card V − 1)`。 -/
theorem z002_finrank_closure (hconn : G.Connected) (hV : 2 ≤ Nat.card V) :
    Module.finrank ℝ (z002Closure G)
      = Nat.card {e // e ∈ G.edgeSet} - (Nat.card V - 1) := by
  have hrn := LinearMap.finrank_range_add_finrank_ker (z002T G)
  have hcard : Module.finrank ℝ ({e // e ∈ G.edgeSet} → ℝ)
      = Nat.card {e // e ∈ G.edgeSet} := by
    rw [Nat.card_eq_fintype_card, Module.finrank_fintype_fun_eq_card]
  rw [hcard, z002_finrank_range hconn hV, add_comm] at hrn
  have h := Nat.add_left_cancel hrn
  rw [z002Closure]
  exact h.symm

/-! ## 七、完全图：主定理 -/

/-- 完全图的边数：`C(m,2)`。 -/
theorem z002_card_edge (m : ℕ) :
    Nat.card {e // e ∈ (completeGraph (Fin m)).edgeSet} = m * (m - 1) / 2 := by
  rw [Nat.card_eq_fintype_card, Fintype.card_coe,
    ← SimpleGraph.card_edgeFinset_top_eq_card_choose_two (V := Fin m), Fintype.card_fin,
    Nat.choose_two_right]

/-- **主定理（Z002 §6.7 的那条读法）**：`dim (圈空间) = dim (运动空间) ⟺ m = 4`。 -/
theorem z002_dim_eq_iff (m : ℕ) (hm : 2 ≤ m) :
    Module.finrank ℝ (z002Closure (completeGraph (Fin m)))
        = Module.finrank ℝ (LinearMap.range (z002T (completeGraph (Fin m))))
      ↔ m = 4 := by
  have hconn : (completeGraph (Fin m)).Connected := by
    haveI : Nonempty (Fin m) := ⟨⟨0, by omega⟩⟩
    exact SimpleGraph.connected_top
  have hV : 2 ≤ Nat.card (Fin m) := by rw [Nat.card_fin]; exact hm
  have hker := z002_finrank_closure (G := completeGraph (Fin m)) hconn hV
  have him := z002_finrank_range (G := completeGraph (Fin m)) hconn hV
  rw [z002_card_edge m, Nat.card_fin] at hker
  rw [Nat.card_fin] at him
  rw [hker, him]
  have h2 : ∀ k : ℕ, 0 < k → k * 2 / 2 = k := fun k hk => Nat.mul_div_cancel _ (by omega)
  constructor
  · intro h
    have h4 : m * (m - 1) / 2 = m - 1 := by omega
    have h5 : ((m - 1) * m / 2) * 2 = (m - 1) * 2 := by
      rw [Nat.mul_comm (m - 1) m, h4]
    rw [h2 (m - 1) (by omega), h2 m (by omega)] at h5
    exact Nat.mul_left_cancel (by omega : 0 < m - 1) h5
  · intro h
    subst h
    norm_num

end Zero
