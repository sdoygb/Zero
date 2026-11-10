/-
  共扼检验的**真实例子**：Zb067 的 `wordSum` 在位置置换下不变。

  这是「普通 Lean 能过 + 共扼检验」的第一条真定理：
  `wordSum (fun i => w (g.symm i)) = wordSum w`，故 `wordSum w = 0` 在位置置换下不变。
-/

import Zero.Conjugation
import Zero.Zb067_AuxWords

namespace Zero

/-- **Zb067 · 真实共扼检验**：`wordSum w = 0` 在位置置换下不变。 -/
theorem zb067_wordSum_conj_invariant (L : ℕ) :
    ConjInvariantFun
      (fun (g : Equiv.Perm (Fin L)) (w : Fin L → Bool) => fun i => w (g.symm i))
      (fun w : Fin L → Bool => wordSum w = 0) := by
  intro g w
  have h : wordSum (fun i => w (g.symm i)) = wordSum w := by
    simpa only [wordSum] using Equiv.sum_comp g.symm (fun i => wordStep (w i))
  change (wordSum (fun i => w (g.symm i)) = 0) ↔ (wordSum w = 0)
  rw [h]

/-- 真实登记表（与玩具表分开）。 -/
def conjRegistryReal : List ConjEntry := [
  ⟨"zb067_wordSum_conj_invariant",
    ∀ L : ℕ, ConjInvariantFun
      (fun (g : Equiv.Perm (Fin L)) (w : Fin L → Bool) => fun i => w (g.symm i))
      (fun w : Fin L → Bool => wordSum w = 0),
    by
      intro L
      exact zb067_wordSum_conj_invariant L⟩
]

end Zero

open Zero

/-- 正例迁移：已有 `wordSum w = 0`，沿位置置换推得 `wordSum (g • w) = 0`。 -/
example (L : ℕ) (g : Equiv.Perm (Fin L)) (w : Fin L → Bool) (h : wordSum w = 0) :
    wordSum (fun i => w (g.symm i)) = 0 := by
  exact conj_of_invariant
    (fun (g : Equiv.Perm (Fin L)) (w : Fin L → Bool) => fun i => w (g.symm i))
    g (fun w => wordSum w = 0) (zb067_wordSum_conj_invariant L g) h

/-- 负例：`w 0 = true` 的共扼陈述是 `w 1 = true`；所需的 missing 是交换不变性，它是**假**的。
    所以 candidate 不是 source 的推论，而是一条**独立的新定理**（共扼猜测的形状）。 -/
def swapWord : Fin 2 → Bool := ![true, false]

def candidateSwap0 : ConjCandidate :=
  ⟨"swap_w0",
   "swapWord 0 = true",
   swapWord 1 = true,
   (∀ w : Fin 2 → Bool, (w 1 = true) ↔ (w 0 = true))⟩

/-- 上面的 missing 确实是假：candidate 不能由 source 推出。 -/
theorem candidateSwap0_missing_false :
    ¬ (∀ w : Fin 2 → Bool, (w 1 = true) ↔ (w 0 = true)) := by
  intro h
  have := h (fun i => decide (i = 0))
  simp at this
