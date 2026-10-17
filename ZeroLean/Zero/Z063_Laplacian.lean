/-
  Z063（"局域 Gauss 不可能"）的机器检查版本：连通图上 **散度为零 ⟹ 沿路为常数**，
  等价地 `ker L = 常量`（`L = D − A` 为图 Laplacian）。

  对应 Z063《攻 8b：加"局域 Gauss 定律" —— 不可能成功》的哪一句：
  * `z063_quadform` ↔ Z063 §2 定理 28 的**能量恒等式**来源：
    `⟨f, L f⟩ = Σ_{e={u,v}∈E}(f u − f v)² ≥ 0`（原文 "由连通性 + 线性代数给出"，
    即先有二次型恒等式，才有下一条的"每条边上相等"）。
  * `z063_div_zero_adj_eq` ↔ Z063 §2 证明的第一步：由 `L f = 0`
    得**每条边上** `f u = f v`（此处尚未用连通性）。
  * `z063_connected_constant` ↔ Z063 式 (0.1) 定理 28：
    `ker Bᵀ = span{(1,…,1)}` 即"连通图上散度为零的解只有常量"（连通性在此用上：
    沿 Walk 传播边上相等）。
  * `z063_eq_zero_of_sum_zero` ↔ Z063 §2 的**横截**结论：
    `ker Bᵀ ∩ H_Q = {0}`，`H_Q = {x : Σ_i x_i = 0}`；
    它是 Z063 结论"局域 Gauss ＋ 全局零和 什么都不许（不是禁闭，是真空）"的机器版本。
  * 与 Z024（离散连续性方程 ⟹ 总荷守恒）的关系：Z024 引理 1 `x(τ) − x(0) = B f`
    （"荷由流决定"）给的是**全局零和** `Σ x = 0` 这一条；
    `z063_eq_zero_of_sum_zero` 正是"这条全局零和"与"局域 Gauss（散度为零）"合起来
    只剩平凡解的那一步（Z063 §3.2：Zero 的零和 = "荷由流决定" 的陈述）。

  用到的 Mathlib 引理（注释中逐条标注）：
  * `SimpleGraph.lapMatrix_mulVec_apply'`：`(L f) v = Σ_{w ∈ N(v)} (f v − f w)`
    —— 把题设的"逐顶点零和"认成 `L f = 0`。
  * `SimpleGraph.lapMatrix_mulVec_eq_zero_iff_forall_adj`：`L f = 0 ↔ ∀ u v, u ~ v → f u = f v`。
  * `SimpleGraph.lapMatrix_toLinearMap₂'`：`f ⬝ᵥ (L f) = (Σ_{i,j, i~j} (f i − f j)²)/2`
    —— Mathlib 已证的**二次型恒等式**（本条即 Z063 的 `⟨f, L f⟩ = Σ_{edges}(Δf)²`）。
  * `Matrix.toLinearMap₂'_apply'`：`toLinearMap₂' ℝ L f f = f ⬝ᵥ (L *ᵥ f)`。
  * `SimpleGraph.neighborFinset_eq_filter` + `Finset.sum_filter`
    （`∑ w ∈ G.neighborFinset v, g w = ∑ w, if G.Adj v w then g w else 0`）。
  * `Finset.sum_div`（`(Σ g)/2 = Σ g/2`）、`Finset.sum_const`、`nsmul_eq_mul`。
  * 连通性：`SimpleGraph.Connected`（structure，带 `CoeFun` 实例，`hconn u v : G.Reachable u v`，
    即 `Nonempty (G.Walk u v)`），对其 `Walk` 做归纳传播边上相等。

  ⚠ 审计整改注明（2026-10-15）：本文件的 `z063_quadform`、`z063_div_zero_adj_eq`、
  `z063_connected_constant` 三条**不是新数学**：本文件是这些 Mathlib 引理在 Zero 记法下的复述（供 Z002/Z033 引用）——
  依次换记法复述 `SimpleGraph.lapMatrix_toLinearMap₂'`、
  `SimpleGraph.lapMatrix_mulVec_eq_zero_iff_forall_adj`、
  `SimpleGraph.lapMatrix_mulVec_eq_zero_iff_forall_reachable`（第三条多用一次连通性沿 Walk 归纳）。
  按审计要求**不删、不改陈述**：它们保留公开，仅在此注明其"复述"性质，
  以免被误当成 Zero 的新命题。本文件真正承重的结论是最后一条
  `z063_eq_zero_of_sum_zero`（局域 Gauss ＋ 全局零和 = 真空），保留公开。
-/

import Mathlib

open scoped Matrix

namespace Zero

/-- **Z063 · 二次型恒等式**：`⟨f, L f⟩ = Σ_{e ∈ E} (f u − f v)²`（右边按每个顶点算邻居和） -/
theorem z063_quadform {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] (f : V → ℝ) :
    (∑ v, f v * (∑ w ∈ G.neighborFinset v, (f v - f w)))
      = ∑ v, ∑ w ∈ G.neighborFinset v, (f v - f w) ^ 2 / 2 := by
  -- 左边 = `f ⬝ᵥ (L f)`：逐顶点用 `lapMatrix_mulVec_apply'` 认出 `(L f) v = Σ_{w∈N(v)} (f v − f w)`
  have hL : (∑ v, f v * (∑ w ∈ G.neighborFinset v, (f v - f w)))
      = Matrix.toLinearMap₂' ℝ (G.lapMatrix ℝ) f f := by
    rw [Matrix.toLinearMap₂'_apply', dotProduct]
    refine Finset.sum_congr rfl fun v _ => ?_
    rw [SimpleGraph.lapMatrix_mulVec_apply']
  -- 右边 = （按"每顶点数一遍邻居"的双重和）/2：邻居和就是 `univ.filter (G.Adj v)`
  have key : ∀ v : V, (∑ w ∈ G.neighborFinset v, (f v - f w) ^ 2 / 2)
      = (∑ w : V, if G.Adj v w then (f v - f w) ^ 2 else 0) / 2 := by
    intro v
    rw [← Finset.sum_div, SimpleGraph.neighborFinset_eq_filter, Finset.sum_filter]
  -- 中间项由 Mathlib 的 `lapMatrix_toLinearMap₂'` 给出（同一个二次型的"边上差平方和"写法）
  rw [hL, SimpleGraph.lapMatrix_toLinearMap₂' ℝ G f, Finset.sum_div,
    Finset.sum_congr rfl (fun v _ => (key v).symm)]

/-- **Z063 · 散度为零 ⟹ 每条边上相等**（连通性还没用上） -/
theorem z063_div_zero_adj_eq {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] (f : V → ℝ)
    (hf : ∀ v, ∑ w ∈ G.neighborFinset v, (f v - f w) = 0) :
    ∀ u v : V, G.Adj u v → f u = f v := by
  -- 题设逐顶点零和 = `L f = 0`
  have hmv : G.lapMatrix ℝ *ᵥ f = 0 := by
    funext v
    rw [SimpleGraph.lapMatrix_mulVec_apply', hf v]
    rfl
  -- 再用 Mathlib 的 ker 刻画 `L f = 0 ↔ 每条边上 f 相等`
  exact fun u v huv =>
    (SimpleGraph.lapMatrix_mulVec_eq_zero_iff_forall_adj (G := G)).mp hmv u v huv

/-- **Z063 · 连通 ⟹ 全域为常数**（ker L = 常量） -/
theorem z063_connected_constant {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] (hconn : G.Connected) (f : V → ℝ)
    (hf : ∀ v, ∑ w ∈ G.neighborFinset v, (f v - f w) = 0) :
    ∀ u v : V, f u = f v := by
  have hadj : ∀ u v : V, G.Adj u v → f u = f v := z063_div_zero_adj_eq G f hf
  intro u v
  -- 连通性：`hconn u v : G.Reachable u v`（= `Nonempty (G.Walk u v)`），沿 Walk 归纳
  obtain ⟨p⟩ := hconn u v
  induction p with
  | nil => rfl
  | cons hAdj _ ih => exact (hadj _ _ hAdj).trans ih

/-- **Z063 · 与零和超平面横截**：`f` 若同时"散度为零"且 `Σ f = 0`，则 `f = 0` -/
theorem z063_eq_zero_of_sum_zero {V : Type*} [Fintype V] [DecidableEq V] (G : SimpleGraph V)
    [DecidableRel G.Adj] (hconn : G.Connected) (f : V → ℝ)
    (hf : ∀ v, ∑ w ∈ G.neighborFinset v, (f v - f w) = 0) (hsum : ∑ v, f v = 0) :
    f = 0 := by
  -- 先由连通性得 `f ≡ f v`（Z063 定理 28：ker Bᵀ = span{1}）
  have hconst : ∀ u v : V, f u = f v := z063_connected_constant G hconn f hf
  funext v
  -- 常量函数的总和 = |V| · (f v)；`G.Connected` 自带 `Nonempty V`，故 |V| ≠ 0
  have h1 : (∑ u, f u) = (Finset.univ : Finset V).card • f v := by
    rw [Finset.sum_congr rfl (fun u _ => hconst u v), Finset.sum_const]
  rw [hsum, Finset.card_univ, nsmul_eq_mul] at h1
  have hcard : (Fintype.card V : ℝ) ≠ 0 := by
    have hpos : 0 < Fintype.card V := Fintype.card_pos_iff.mpr ⟨v⟩
    exact_mod_cast Nat.pos_iff_ne_zero.mp hpos
  -- 于是 |V| · f v = 0 逼出 f v = 0：这就是"局域 Gauss ＋ 全局零和 = 真空"
  exact (mul_eq_zero.mp h1.symm).resolve_left hcard

#print axioms z063_quadform
#print axioms z063_div_zero_adj_eq
#print axioms z063_connected_constant
#print axioms z063_eq_zero_of_sum_zero

end Zero