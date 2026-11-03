/-
  Zb014 §1 的机器检查版本：受迫稳态的频率锁定与幅度线性。

  正文定理 1.2 的核心是标量方程
      `(c²λ − f² + 2iγf) U = g`
  在 `γ f ≠ 0` 时有唯一解；频率 `f` 是方程的参数，故稳态频率恒为 `f`，
  与驱动幅度 `g` 无关。命题 1.3 再把 `U = g U₀` 的线性性写出来。

  未覆盖：矩阵版 `M = c²L − f² + 2iγf I` 的可逆性（需要谱/正定性），
  以及推论 1.4 的波数反解（已在 Zb005 证过 `c=1` 的版本）。
-/

import Mathlib

namespace Zero

/-- **Zb014 · 频率锁定（标量核）**：`z = a + i b`、`b ≠ 0` 时，
    方程 `z U = g` 有唯一解 `U = g / z`。 -/
theorem zb014_frequency_locked (a b g : ℝ) (hb : b ≠ 0) :
    ∃! U : ℂ, ((a : ℂ) + (b : ℂ) * Complex.I) * U = (g : ℂ) := by
  have hz : ((a : ℂ) + (b : ℂ) * Complex.I) ≠ 0 := by
    intro h
    have him := congrArg Complex.im h
    simp at him
    exact hb him
  refine ⟨(g : ℂ) / ((a : ℂ) + (b : ℂ) * Complex.I), ?_, ?_⟩
  · field_simp [hz]
  · intro U hU
    calc
      U = (((a : ℂ) + (b : ℂ) * Complex.I) * U) / ((a : ℂ) + (b : ℂ) * Complex.I) := by
        field_simp [hz]
      _ = (g : ℂ) / ((a : ℂ) + (b : ℂ) * Complex.I) := by rw [hU]

/-- **Zb014 · 幅度线性**：若 `U₀` 是单位源响应（`z U₀ = 1`），
    则 `g U₀` 是源 `g` 的响应（`z (g U₀) = g`）。 -/
theorem zb014_amplitude_linear (z : ℂ) (g : ℝ) (U₀ : ℂ) (h : z * U₀ = 1) :
    z * ((g : ℂ) * U₀) = (g : ℂ) := by
  calc
    z * ((g : ℂ) * U₀) = (g : ℂ) * (z * U₀) := by ring
    _ = (g : ℂ) := by rw [h, mul_one]

/-- **Zb014 · 频率不进入幅度比**：两个驱动幅度 `g₁, g₂` 的响应之比是 `g₁/g₂`，
    与频率参数 `f` 无关（这里用 `z` 代表 `c²λ − f² + 2iγf`）。 -/
theorem zb014_amplitude_ratio (z : ℂ) (g₁ g₂ : ℝ) (U₀ : ℂ)
    (h : z * U₀ = 1) (hg₂ : g₂ ≠ 0) :
    ((g₁ : ℂ) * U₀) / ((g₂ : ℂ) * U₀) = (g₁ : ℂ) / (g₂ : ℂ) := by
  have hU₀ : U₀ ≠ 0 := by
    intro h0
    rw [h0, mul_zero] at h
    exact zero_ne_one h
  field_simp [hU₀, hg₂]

#print axioms zb014_frequency_locked
#print axioms zb014_amplitude_linear
#print axioms zb014_amplitude_ratio

end Zero
