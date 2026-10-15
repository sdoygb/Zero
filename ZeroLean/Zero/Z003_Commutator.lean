/-
  Z003（量子力学导出纲领）与 Z060（有限忠实迹下不存在 [X,P]=i·1）的机器检查版本。

  核心命题：在有限维矩阵代数中，迹是线性的且 tr(XP)=tr(PX)，所以
    tr([X,P]) = tr(XP−PX) = 0
  但 tr(i·I) = i·n ≠ 0（n = 维数 > 0）。矛盾。
  因此 [X,P]=i·1 在任何有限维表示中无解 —— 这是"位置/动量共轭变量需要无穷维
  （即 Wey1 群表示论）"的算术内核。
-/

import Mathlib

namespace Zero

open Matrix

variable {n : ℕ} [NeZero n]

/-- **Z003/Z060 · 有限维无共轭对**：对任意两个 n×n 方阵 X、P，
    `tr(XP − PX) = 0`，因此 `XP − PX` 不可能等于 `i • 1`（其迹为 `i•n ≠ 0`）。 -/
theorem z003_no_fin_dim_commutator (X P : Matrix (Fin n) (Fin n) ℂ) :
    trace (X * P - P * X) = 0 := by
  rw [sub_eq_add_neg, trace_add]
  have h : trace (X * P) = trace (P * X) := Matrix.trace_mul_comm X P
  rw [h, trace_neg, add_neg_cancel]

/-- 具体推论：`XP − PX = i • 1` 在 n > 0 时不可能（左边迹 0，右边迹 i·n）。 -/
theorem z003_commutator_ne_scalar (X P : Matrix (Fin n) (Fin n) ℂ) :
    X * P - P * X ≠ (Complex.I : ℂ) • (1 : Matrix (Fin n) (Fin n) ℂ) := by
  intro h
  have htr := congrArg trace h
  rw [z003_no_fin_dim_commutator, trace_smul, trace_one] at htr
  -- 左边 = 0，右边 = i • (n : ℂ)
  have hn : (n : ℂ) ≠ 0 := by exact_mod_cast NeZero.pos n |>.ne'
  simpa using hn (by simpa using htr)

#print axioms z003_no_fin_dim_commutator
#print axioms z003_commutator_ne_scalar

end Zero
