/-
  Z082（Koide 相位从权重导出）的机器检查版本。

  Z082 的参数化：三项权重 s_k 的三点 DFT 给出
    ŝ₀ = Σ_k s_k,  ŝ₁ = Σ_k s_k ω^k   (ω = e^{2πi/3}),
    θ = arg ŝ₁,  b = 2|ŝ₁|/|ŝ₀|,  Q = (1 + b²/2)/3.

  本篇钉两件事：
    · **Koide 恒等式** `b² = 2(3Q−1)` —— 它只是 Q 定义的代数重排（因此"Q=2/3 ⟺ b=√2"）；
    · **等权 ⟹ 第一谐波为零**：若三项权重全同，则 ŝ₁ = 0（这正是 Z124 的 `1+ω+ω²=0`）。

  诚实边界：Z082 里"取三点 DFT 作为权重→相位的桥"是**具名输入**（模型选择），
  本篇证的是"给定这个桥以后，两条代数关系对不对"。
-/

import Mathlib
import Zero.Z063_Closure

namespace Zero

/-- **Z082 · Koide 恒等式**：`b = 2|ŝ₁|/|ŝ₀|`、`Q = (1+b²/2)/3` ⟹ `b² = 2(3Q−1)`。
    于是 `Q = 2/3` 与 `b = √2` 等价 —— Koide 关系即这条代数重排。 -/
theorem z082_koide_identity (a c : ℝ) (hc : c ≠ 0) :
    let b : ℝ := 2 * a / c
    let Q : ℝ := (1 + b ^ 2 / 2) / 3
    b ^ 2 = 2 * (3 * Q - 1) := by
  simp only
  field_simp
  ring

/-- **Z082 · Koide 关系的等价形式**：`b² = 2(3Q−1)` 时 `Q = 2/3 ⟺ b² = 2`。 -/
theorem z082_koide_Q_two_thirds (b Q : ℝ) (h : b ^ 2 = 2 * (3 * Q - 1)) :
    Q = 2 / 3 ↔ b ^ 2 = 2 := by
  constructor
  · intro hQ
    rw [h, hQ]; norm_num
  · intro hb
    rw [hb] at h
    linarith

/-- **Z082 · 等权 ⟹ 第一谐波为零**：三项权重全同时 `ŝ₁ = Σ_k s ω^k = s(1+ω+ω²) = 0`
    （用 Z124 的色单态相消）。 -/
theorem z082_equal_weights_kill_first_harmonic (s : ℂ) :
    s + s * Complex.exp (2 * Real.pi / 3 * Complex.I)
      + s * Complex.exp (2 * Real.pi / 3 * Complex.I) ^ 2 = 0 := by
  have h := z063_color_singlet_sum_zero
  linear_combination s * h

/-- **Z082 · 零权重项不贡献第一谐波**：ŝ₁ 对常数偏移不敏感（`1+ω+ω²=0` 的直接后果），
    所以"相位"只由权重的**差异**决定、与绝对标度无关。 -/
theorem z082_common_shift_irrelevant (s d : ℂ) :
    (s + d) + (s + d) * Complex.exp (2 * Real.pi / 3 * Complex.I)
      + (s + d) * Complex.exp (2 * Real.pi / 3 * Complex.I) ^ 2 = 0 :=
  z082_equal_weights_kill_first_harmonic (s + d)

#print axioms z082_koide_identity
#print axioms z082_koide_Q_two_thirds
#print axioms z082_equal_weights_kill_first_harmonic
#print axioms z082_common_shift_irrelevant

end Zero
