/-
  Zb063 §2.5 的机器检查版本：平均场块的 Bogoliubov 对角化与自洽代数。

  对应正文：
    · 定理 2.4：块 `[[0, −Δ], [−Δ, 2ε]]` 的本征值是 `ε ± E`（`E² = ε² + Δ²`），
      以及 BCS 系数的归一 `u² + v² = 1`；
    · 定理 2.5：单能级形态的自洽 `Δ = V Δ/(2E)`、`Δ ≠ 0` ⇒ `1 = V/(2E)`；
    · 命题 2.7：`K = 1` 时基态能量的两种写法一致。

  说明：本文件不引入任何 Zero 侧公理；依赖清单应只有
  `propext` / `Quot.sound` / `Classical.choice`（见文件末 `#print axioms`）。
-/

import Mathlib

namespace Zero

/-- **Zb063 · BCS 系数归一**：`u² + v² = 1`，其中 `u² = (1 + ε/E)/2`、`v² = (1 − ε/E)/2`。 -/
theorem zb063_bcs_coeff_sq (ε E : ℝ) (hE : E ≠ 0) :
    (1 + ε / E) / 2 + (1 - ε / E) / 2 = 1 := by
  field_simp [hE]
  ring

/-- **Zb063 · 平均场块的本征值（上支）**：`[[0, −Δ], [−Δ, 2ε]]` 的特征多项式在 `ε + E` 处为零。 -/
theorem zb063_block_eigen_plus (ε Δ E : ℝ) (hE : E ^ 2 = ε ^ 2 + Δ ^ 2) :
    (0 - (ε + E)) * (2 * ε - (ε + E)) - (-Δ) * (-Δ) = 0 := by
  nlinarith [hE]

/-- **Zb063 · 平均场块的本征值（下支）**：特征多项式在 `ε − E` 处为零。 -/
theorem zb063_block_eigen_minus (ε Δ E : ℝ) (hE : E ^ 2 = ε ^ 2 + Δ ^ 2) :
    (0 - (ε - E)) * (2 * ε - (ε - E)) - (-Δ) * (-Δ) = 0 := by
  nlinarith [hE]

/-- **Zb063 · 自洽 ⇒ 能隙方程（单能级形态）**：`Δ = V Δ/(2E)`、`Δ ≠ 0` ⇒ `1 = V/(2E)`。 -/
theorem zb063_gap_selfconsistent (V Δ E : ℝ) (hE : E ≠ 0) (hΔ : Δ ≠ 0)
    (h : Δ = V * Δ / (2 * E)) : 1 = V / (2 * E) := by
  have hfac : Δ * (1 - V / (2 * E)) = 0 := by
    have h1 : Δ * (1 - V / (2 * E)) = Δ - V * Δ / (2 * E) := by ring
    rw [h1]
    linarith
  have hzero : 1 - V / (2 * E) = 0 := (mul_eq_zero.mp hfac).resolve_left hΔ
  linarith

/-- **Zb063 · 基态能量的两种写法（`K = 1`）**：`ε − E + Δ²/V = 2ε v² − Δ²/V`。 -/
theorem zb063_mf_energy_single (ε Δ E V : ℝ) (hE : E ≠ 0)
    (hE2 : E ^ 2 = ε ^ 2 + Δ ^ 2) (hgap : 1 = V / (2 * E)) :
    ε - E + Δ ^ 2 / V = 2 * ε * ((1 - ε / E) / 2) - Δ ^ 2 / V := by
  have hV : V = 2 * E := by
    field_simp [hE] at hgap
    linarith
  rw [hV]
  field_simp [hE]
  nlinarith [hE2]

#print axioms zb063_bcs_coeff_sq

#print axioms zb063_block_eigen_plus

#print axioms zb063_block_eigen_minus

#print axioms zb063_gap_selfconsistent

/-- **Zb063 · 单能级阈值**：`K = 1` 时 `1 = V/(2√(ε²+Δ²))` 有正解 `Δ` 当且仅当 `V > 2|ε|`
    （§2 定理 2.2 的 `K = 1` 形态）。 -/
theorem zb063_threshold_single (ε V : ℝ) :
    (∃ D : ℝ, 0 < D ∧ 1 = V / (2 * Real.sqrt (ε ^ 2 + D ^ 2))) ↔ 2 * |ε| < V := by
  constructor
  · rintro ⟨D, hD, h⟩
    have hEpos : 0 < Real.sqrt (ε ^ 2 + D ^ 2) := Real.sqrt_pos.mpr (by positivity)
    have hV : V = 2 * Real.sqrt (ε ^ 2 + D ^ 2) := by
      field_simp [hEpos.ne'] at h
      linarith
    have hEgt : |ε| < Real.sqrt (ε ^ 2 + D ^ 2) := by
      rw [← Real.sqrt_sq_eq_abs ε]
      exact Real.sqrt_lt_sqrt (by positivity) (by nlinarith)
    nlinarith
  · intro hV
    have hε : |ε| < V / 2 := by linarith
    have hVpos : 0 < V := by nlinarith [abs_nonneg ε]
    have hrad : 0 < (V / 2) ^ 2 - ε ^ 2 := by nlinarith [hε, abs_nonneg ε, sq_abs ε]
    refine ⟨Real.sqrt ((V / 2) ^ 2 - ε ^ 2), Real.sqrt_pos.mpr hrad, ?_⟩
    have hsq : ε ^ 2 + (Real.sqrt ((V / 2) ^ 2 - ε ^ 2)) ^ 2 = (V / 2) ^ 2 := by
      rw [Real.sq_sqrt hrad.le]
      ring
    have hE : Real.sqrt (ε ^ 2 + (Real.sqrt ((V / 2) ^ 2 - ε ^ 2)) ^ 2) = V / 2 := by
      rw [hsq, Real.sqrt_sq_eq_abs, abs_of_pos (by linarith : (0:ℝ) < V / 2)]
    rw [hE]
    field_simp [hVpos.ne']

/-- **Zb063 · $K=2$ 的 $P=1$ 块**：矩阵 `[[2ε₁−V, −V], [−V, 2ε₂−V]]` 的特征多项式在
    `ε₁+ε₂−V ± √((ε₁−ε₂)²+V²)` 处为零。 -/
theorem zb063_pair_block_k2 (e1 e2 V E : ℝ) (hE : E ^ 2 = (e1 - e2) ^ 2 + V ^ 2) :
    (2 * e1 - V - (e1 + e2 - V + E)) * (2 * e2 - V - (e1 + e2 - V + E)) - V ^ 2 = 0 ∧
    (2 * e1 - V - (e1 + e2 - V - E)) * (2 * e2 - V - (e1 + e2 - V - E)) - V ^ 2 = 0 := by
  constructor <;> nlinarith [hE]

/-- **Zb063 · $P=1$ 块的一般形态**（命题 2.12）：若 $1=V\sum_k(2\varepsilon_k-\lambda)^{-1}$，
    则 $v_k=(2\varepsilon_k-\lambda)^{-1}$ 是 $P=1$ 块的本征向量、本征值 $\lambda$。 -/
theorem zb063_pair_block_general {n : ℕ} (ε : Fin n → ℝ) (V lam : ℝ)
    (hden : ∀ k, 2 * ε k - lam ≠ 0)
    (h : 1 = V * ∑ k, (2 * ε k - lam)⁻¹) :
    ∀ k, 2 * ε k * (2 * ε k - lam)⁻¹ - V * (∑ l, (2 * ε l - lam)⁻¹)
      = lam * (2 * ε k - lam)⁻¹ := by
  intro k
  have hv : V * (∑ l, (2 * ε l - lam)⁻¹) = 1 := by linarith
  rw [hv]
  field_simp [hden k]
  ring

#print axioms zb063_pair_block_general

#print axioms zb063_pair_block_k2

#print axioms zb063_threshold_single

#print axioms zb063_mf_energy_single

end Zero
