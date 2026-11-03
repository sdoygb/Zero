/-
  共扼检验（原型）：把「共扼不变」变成每条定理的**证明义务**。

  用户构想：
    · 假设共扼是普遍属性；
    · 每条定理登记一个 ConjEntry，proof 字段就是「该定理在共扼下不变」的证明；
    · 证明失败/未登记 ⟹ 缺少一条共扼引理/定理，正是「缺的定理」；
    · 普通 Lean 能过、共扼检验不过的定理会立刻暴露。

  本文件：框架 + 最小玩具实例（S₃ = Perm (Fin 3) 作用在 Fin 3 上）。
  真正的库定理登记由后续逐条补。
-/

import Mathlib

namespace Zero

/-- 共扼不变：谓词 `P` 在群 `G` 对 `X` 的作用下不变。 -/
def ConjInvariant (G X : Type*) [Group G] [MulAction G X] (P : X → Prop) : Prop :=
  ∀ (g : G) (x : X), P (g • x) ↔ P x


/-- 共扼不变（显式作用版）：不依赖 typeclass，直接用作用函数。 -/
def ConjInvariantFun {G X : Type*} (act : G → X → X) (P : X → Prop) : Prop :=
  ∀ g x, P (act g x) ↔ P x


/-- 共扼陈述：把 `P` 沿作用 `act g` 推过去。 -/
def conjStmt {G X : Type*} (act : G → X → X) (g : G) (P : X → Prop) : Prop :=
  ∀ x, P (act g x)

/-- **共扼迁移**：若 `P` 在该作用下不变，则 `P x` 给出共扼陈述 `P (act g x)`。 -/
theorem conj_of_invariant {G X : Type*} (act : G → X → X) (g : G) (P : X → Prop)
    (hinv : ∀ x, P (act g x) ↔ P x) {x : X} (h : P x) : P (act g x) :=
  (hinv x).mpr h

/-- 一条「待补的共扼定理」：source 是已有定理，candidate 是它的共扼陈述，
missing 是推出 candidate 所需的共扼不变引理（若 missing 假，candidate 是独立新定理）。 -/
structure ConjCandidate where
  name : String
  source : String
  candidate : Prop
  missing : Prop

/-- 登记的共扼条目；`proof` 就是该条的共扼检验。 -/
structure ConjEntry where
  name : String
  claim : Prop
  proof : claim

/-- 正例：恒真谓词在 S₃ 下不变。 -/
example : ConjInvariant (Equiv.Perm (Fin 3)) (Fin 3) (fun _ => True) := by
  intro g x
  simp

/-- 负例：`x = 0` 在 S₃ 下**不**不变；检验失败，暴露「缺共扼引理」。 -/
example : ¬ ConjInvariant (Equiv.Perm (Fin 3)) (Fin 3) (fun x : Fin 3 => x = 0) := by
  intro h
  have hh := h (Equiv.swap (0 : Fin 3) (1 : Fin 3)) 0
  simp at hh

/-- 共扼登记表：每个条目都带证明；证明失败则 build 失败。 -/
def conjRegistry : List ConjEntry := [
  ⟨"toy_true_invariant",
    ConjInvariant (Equiv.Perm (Fin 3)) (Fin 3) (fun _ => True),
    by
      intro g x
      simp⟩,
  ⟨"toy_point_not_invariant",
    ¬ ConjInvariant (Equiv.Perm (Fin 3)) (Fin 3) (fun x : Fin 3 => x = 0),
    by
      intro h
      have hh := h (Equiv.swap (0 : Fin 3) (1 : Fin 3)) 0
      simp at hh⟩
]

end Zero
