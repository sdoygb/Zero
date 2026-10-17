/-
  Z044（能谱之桥能否推掉 · §4「模流在账本态下为何平凡」/ §5「去掉块内标量会怎样」）的机器检查版本。

  要验的命题（有限维 2×2 显式模型）：密度矩阵 `ρ = diag(p, q)`（`p, q > 0`，`p_0 = p`、`p_1 = q`），
  模流在对角态上是相位缩放
    `σ_t(A)_{ij} = exp(i · t · (log p_i − log p_j)) · A_{ij}`。
  取非对角基元 `A = E_{01}`：`σ_t(E_{01}) = e^{i t (log p − log q)} · E_{01}`。于是
    · `p = q`（**中心态**／块标量）⟹ 模流**平凡**：`z044_flow_trivial_of_eq`；
    · `p ≠ q`（**非块标量**）⟹ 存在 `t` 与 `A` 使 `σ_t(A) ≠ A`：`z044_nontrivial`——
      见证取 `t (log p − log q) = π`，此时 `σ_t(E_{01}) = −E_{01} ≠ E_{01}`。

  这是 Z044 §5.3 对 **Z040 定理 34** 的**更正**：「非平凡内禀模流 ⟺ A 是因子」**过强**——
  在 2×2（块内）上，只要态**不是块标量**，模流就已非平凡，**不需要**换代数型（换成因子）。
  即：使 `Δ = 1` 的是**块内标量性**，不是块对角性；去掉「块内标量」这一**额外假设**即可。

  ⚠ 口径（机器核验的边界）：本文件只验 2×2 有限维**显式模型**里的相位公式与两条推论。
  「`σ_t` 是 von Neumann 代数上的模流（Tomita–Takesaki）」「KMS 条件」「`Δ` 的算子定义」
  这些**不在**本文件的形式化范围内——这里形式化的是对角态的相位公式本身。
  数值旁证见 Z044 正文与 `_staging/Z044_energy_bridge_check.py`（§C：块内不等的 6 个判据）。
-/

import Mathlib

-- 声明级关闭 unusedVariables：逐字陈述里的假设一律不删（任务纪律），故就地静音该 linter。
set_option linter.unusedVariables false

namespace Zero

/-- 非对角基元 `E_{01}`（在 `Fin 2` 上） -/
def z045OffDiag : Matrix (Fin 2) (Fin 2) ℂ :=
  Matrix.of (fun i j => if i = 0 ∧ j = 1 then 1 else 0)

/-- 模流在 2×2 对角态上：`σ_t(A)_{ij} = e^{i t (log p_i − log p_j)} A_{ij}`（`p_0 = p`、`p_1 = q`） -/
noncomputable def z045Flow (p q t : ℝ) (A : Matrix (Fin 2) (Fin 2) ℂ) : Matrix (Fin 2) (Fin 2) ℂ :=
  Matrix.of (fun i j => Complex.exp (Complex.I * (t * (Real.log (if i = 0 then p else q)
    - Real.log (if j = 0 then p else q)))) * A i j)

/-- **Z044 · 中心态 ⟹ 模流平凡**：`p = q` 时 `σ_t = id` -/
theorem z044_flow_trivial_of_eq (p t : ℝ) (A : Matrix (Fin 2) (Fin 2) ℂ) :
    z045Flow p p t A = A := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [z045Flow]

/-- **Z044 · 非对角元的相位**：`σ_t(E_{01}) = e^{i t (log p − log q)} • E_{01}` -/
theorem z044_flow_offdiag (p q t : ℝ) :
    z045Flow p q t z045OffDiag
      = (Complex.exp (Complex.I * (t * (Real.log p - Real.log q)))) • z045OffDiag := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [z045Flow, z045OffDiag]

/-- **Z044 · 非块标量 ⟹ 模流非平凡**：`p ≠ q`、`p,q > 0` ⟹ 存在 `t` 与 `A`，`σ_t(A) ≠ A` -/
theorem z044_nontrivial (p q : ℝ) (hp : 0 < p) (hq : 0 < q) (hpq : p ≠ q) :
    ∃ t : ℝ, z045Flow p q t z045OffDiag ≠ z045OffDiag := by
  -- 相位系数取 `−1` 的充分条件：`t (log p − log q) = π`
  have hlog : Real.log p - Real.log q ≠ 0 := by
    intro h
    exact hpq (Real.log_injOn_pos hp hq (sub_eq_zero.mp h))
  -- 分母在 ℂ 里非零（`↑(log p) − ↑(log q) = ↑(log p − log q)`，由 `Complex.ofReal_sub`）
  have hden : ((Real.log p : ℂ) - (Real.log q : ℂ)) ≠ 0 := by
    intro h
    have h1 : (Real.log p : ℂ) = (Real.log q : ℂ) := sub_eq_zero.mp h
    exact hlog (sub_eq_zero.mpr (Complex.ofReal_inj.mp h1))
  -- 相位系数：`I · (t (log p − log q)) = ↑π · I`（`t = π/(log p − log q)` 消去分母），实部为 `exp(πi) = −1`
  have hcoeff :
      Complex.exp (Complex.I * ((Real.pi / (Real.log p - Real.log q))
        * (Real.log p - Real.log q))) = -1 := by
    have harg : Complex.I * ((Real.pi / (Real.log p - Real.log q))
        * (Real.log p - Real.log q)) = ↑Real.pi * Complex.I := by
      rw [div_mul_cancel₀ (↑Real.pi) hden]
      exact mul_comm _ _
    rw [harg]
    exact Complex.exp_pi_mul_I
  refine ⟨Real.pi / (Real.log p - Real.log q), ?_⟩
  -- 系数为 `−1` 时 `σ_t(E_{01}) = −E_{01} ≠ E_{01}`（在 `(0,1)` 分量上收尾）
  intro hcontra
  rw [z044_flow_offdiag] at hcontra
  have h01 := congrFun (congrFun hcontra 0) 1
  simp [Matrix.smul_apply, z045OffDiag] at h01
  rw [hcoeff] at h01
  exact absurd h01 (by norm_num)

/-! ### 要素来源表（`#print axioms`） -/

#print axioms z044_flow_trivial_of_eq
#print axioms z044_flow_offdiag
#print axioms z044_nontrivial

end Zero