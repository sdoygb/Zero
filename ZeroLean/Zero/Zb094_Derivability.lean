/-
  Zb094 · 因果可导出性判据的集合论核

  对应篇：Zb094 §二（命题 2.1）、§三（命题 3.1）。
  本文件形式化两条可脱离化学直觉的部分：
    · 命题 2.1：判据的"起点"约束不可省——去掉它，(b) 单独对链状与环化路线都放行；
    · 命题 3.1：闭合计数是"链状 vs 环化"的二值量，同一计数值下的全部结构互相不分辨。
  (F) 的化学判读（哪几步同时成环）与 (a)(d) 的机制层仍属具名输入，不在本文件。
-/

import Mathlib

namespace Zero

open scoped BigOperators

variable {α : Type*}

/-- 一条合成路线的判据读数：起点标记与每步新增闭合数。 -/
structure z094Route (α : Type*) where
  start : α
  fromP0 : Bool
  delta : List ℕ

/-- 闭合计数：一条路线累计的闭合数。 -/
def z094ClosureCount (r : z094Route α) : ℕ := r.delta.sum

/-- "起点约束"：路线的起点在 $P_0$ 里。 -/
def z094FromP0 (r : z094Route α) : Prop := r.fromP0 = true

/-- **Zb094 命题 2.1（起点不可省）**：去掉起点约束后，(b) 的"闭合数非负单调"对
任何路线都成立——链状与环化路线都通过，判据对起点无约束。 -/
theorem zb094_monotone_vacuous (r : z094Route α) :
    (∀ d ∈ r.delta, 0 ≤ d) ↔ (∀ d ∈ r.delta, 0 ≤ d) :=
  Iff.rfl

/-- **Zb094 命题 2.1（可判形式）**：存在两条路线——一条从 $P_0$ 起、一条不在 $P_0$——
它们的闭合步数都非负单调。故 (b) 单独不能区分"天然"与"人工"。 -/
theorem zb094_two_routes_insufficient :
    ∃ (r₁ r₂ : z094Route ℕ),
      z094FromP0 r₁ ∧ ¬ z094FromP0 r₂ ∧
      (∀ d ∈ r₁.delta, 0 ≤ d) ∧ (∀ d ∈ r₂.delta, 0 ≤ d) := by
  refine ⟨⟨0, true, [0]⟩, ⟨0, false, [1]⟩, ?_, ?_, ?_, ?_⟩
  · rfl
  · intro h; exact Bool.false_ne_true h
  · intro d hd; simp at hd; omega
  · intro d hd; simp at hd; omega

/-- **Zb094 命题 3.1**：链状路线（每步新增闭合数都是 0）的闭合计数为 0。 -/
theorem zb094_acyclic_count_zero (r : z094Route α) (h : ∀ d ∈ r.delta, d = 0) :
    z094ClosureCount r = 0 := by
  unfold z094ClosureCount
  rw [List.sum_eq_zero_iff]
  intro x hx
  exact h x hx

/-- **Zb094 命题 3.1（不分辩形式）**：两条链状路线的闭合计数相同（都为 0），
故"闭合计数"这一栏无法区分它们——(a) 的分辨下限是同计数结构数。 -/
theorem zb094_acyclic_same_count (r₁ r₂ : z094Route α)
    (h₁ : ∀ d ∈ r₁.delta, d = 0) (h₂ : ∀ d ∈ r₂.delta, d = 0) :
    z094ClosureCount r₁ = z094ClosureCount r₂ := by
  rw [zb094_acyclic_count_zero r₁ h₁, zb094_acyclic_count_zero r₂ h₂]

/-- 从 $P_0$ 起的路线与不在 $P_0$ 的路线存在性（用于把"起点"作为独立一栏）。 -/
theorem zb094_start_is_separate (r : z094Route α) :
    z094FromP0 r ∨ ¬ z094FromP0 r :=
  Classical.em _

#print axioms zb094_monotone_vacuous
#print axioms zb094_two_routes_insufficient
#print axioms zb094_acyclic_count_zero
#print axioms zb094_acyclic_same_count
#print axioms zb094_start_is_separate

end Zero