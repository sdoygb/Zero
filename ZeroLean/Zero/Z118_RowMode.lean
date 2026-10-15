/-
  Z118（能量标度桥 · 行模反演）的机器检查版本。

  要验的命题：1D 色散关系 `λ(q) = 2 − 2 cos q = 4 sin²(q/2)`，以及它的**逆** ——
  给定 `λ = f²`（`q ∈ [0, π]`、`0 ≤ f ≤ 2`），行模反演解出 `q = 2 · arcsin(f/2)`；
  反向亦然。这三条是 Z118 §2 桥选项 A 里 "ω = f → q = 2 arcsin(f/2) → 定态散射"
  这一步的纯数学核（数值侧见 `ZeroGPU/src/k4_quantum_doubleslit_A.py`）。

  白话线路：
    · 恒等式：`cos q = cos(2·(q/2)) = cos²(q/2) − sin²(q/2) = 1 − 2 sin²(q/2)`；
    · 正向：`4 sin²(q/2) = f²` ⟹ `2 sin(q/2) = ±f`；但 `q/2 ∈ [0, π/2]` 给 `sin(q/2) ≥ 0`、
      又 `f ≥ 0`，两个非负量之积仅当都为 0 时为 0 ⟹ `sin(q/2) = f/2`；
      再由 `q/2 ∈ [−π/2, π/2]` 用 `Real.arcsin_sin` 取逆得 `arcsin(f/2) = q/2`；
    · 逆向：`q/2 = arcsin(f/2)`，由 `|f/2| ≤ 1` 用 `Real.sin_arcsin` 得 `sin(q/2) = f/2`，
      回代恒等式即得 `2 − 2 cos q = 4·(f/2)² = f²`。

  ⚠ 实现口径（本机 Lean 4.35.0 + Mathlib 实测）：
    任务书建议的两个引理名 `Real.cos_sq_half` / `Real.sin_sq_half` **在本版 Mathlib 里不存在**
    （`Unknown constant`）。本文件改用等价且在库的两条：`Real.cos_two_mul'`
    （`cos (2x) = cos²x − sin²x`）与 `Real.cos_sq'`（`cos²x = 1 − sin²x`）。
    三条定理的名字与陈述与任务书逐字一致，未做任何削弱；无 `sorry`/`admit`/`native_decide`。

  ⚠ 审计整改（2026-10-15）：原陈述里对证明并非必需的假设已按审计**删除**
    （删前提＝加强命题）：正向 `z118_inverse_of_dispersion` 的 `hf2 : f ≤ 2`；
    逆向 `z118_dispersion_of_inverse` 的 `hq0 : 0 ≤ q`、`hqpi : q ≤ Real.pi`。
    删后两条仍证得出（正向的乘积分解只用 `0 ≤ f` 与 `sin(q/2) ≥ 0`，
    `f ≤ 2` 本可由 `2 − 2cos q = f² ≤ 4` 推出；逆向由 `h` 直接读出 `q/2 = arcsin(f/2)`，
    `q ∈ [0,π]` 从未参与），故不再需要原先为"不删假设"而设的
    `set_option linter.unusedVariables false`，该选项已删。

  登记状态（2026-10-15 补注）：本文件**已接入** `Zero.lean`，进 `lake build` 与 CI 审计空间（原注释写于接入前，已更正）。单文件验收：
      `cd ZeroLean && lake env lean Zero/Z118_RowMode.lean`。
-/

import Mathlib

namespace Zero

/-- **Z118 · 色散恒等式**：`2 − 2cos q = 4 sin²(q/2)` -/
theorem z118_dispersion_identity (q : ℝ) :
    2 - 2 * Real.cos q = 4 * Real.sin (q / 2) ^ 2 := by
  -- 倍角 + 平方关系：cos q = 1 − 2 sin²(q/2)
  have h2 : Real.cos q = 1 - 2 * Real.sin (q / 2) ^ 2 := by
    have h := Real.cos_two_mul' (q / 2)
    rw [show 2 * (q / 2) = q by ring, Real.cos_sq'] at h
    linarith
  linarith

/-- **Z118 · 反演（正向）**：`q ∈ [0, π]`、`0 ≤ f`、`2 − 2cos q = f²` ⟹ `q = 2·arcsin(f/2)`
    （原陈述带的 `f ≤ 2` 已被审计删去：它在本证明中未被用到） -/
theorem z118_inverse_of_dispersion (q f : ℝ) (hq0 : 0 ≤ q) (hqpi : q ≤ Real.pi)
    (hf0 : 0 ≤ f)
    (h : 2 - 2 * Real.cos q = f ^ 2) :
    q = 2 * Real.arcsin (f / 2) := by
  -- q/2 ∈ [0, π/2]，故 sin(q/2) ≥ 0
  have hq2_nonneg : 0 ≤ q / 2 := by linarith
  have hq2_le : q / 2 ≤ Real.pi / 2 := by linarith
  have hsin_nonneg : 0 ≤ Real.sin (q / 2) :=
    Real.sin_nonneg_of_nonneg_of_le_pi hq2_nonneg (by linarith)
  -- 恒等式把 `2 − 2cos q = f²` 变成 `(2 sin(q/2))² = f²`
  have hdisp := z118_dispersion_identity q
  have hsq : (2 * Real.sin (q / 2)) ^ 2 = f ^ 2 := by nlinarith [hdisp, h]
  -- 两个非负量：平方相等 ⟹ 相等
  have hsin_eq : Real.sin (q / 2) = f / 2 := by
    have hprod : (2 * Real.sin (q / 2) - f) * (2 * Real.sin (q / 2) + f) = 0 := by
      nlinarith [hsq]
    rcases mul_eq_zero.mp hprod with hcase | hcase
    · linarith
    · linarith
  -- q/2 ∈ [−π/2, π/2]，可对 sin 取逆
  have harc : Real.arcsin (Real.sin (q / 2)) = q / 2 := Real.arcsin_sin (by linarith) hq2_le
  rw [hsin_eq] at harc
  linarith

/-- **Z118 · 反演（逆向）**：`0 ≤ f ≤ 2`、`q = 2·arcsin(f/2)` ⟹ `2 − 2cos q = f²`
    （原陈述带的 `hq0 : 0 ≤ q`、`hqpi : q ≤ Real.pi` 已被审计删去：本证明中未用到） -/
theorem z118_dispersion_of_inverse (q f : ℝ)
    (hf0 : 0 ≤ f) (hf2 : f ≤ 2)
    (h : q = 2 * Real.arcsin (f / 2)) :
    2 - 2 * Real.cos q = f ^ 2 := by
  have hq2 : q / 2 = Real.arcsin (f / 2) := by linarith
  -- f/2 ∈ [−1, 1]（由 0 ≤ f ≤ 2），故 sin(arcsin(f/2)) = f/2
  have hsin : Real.sin (q / 2) = f / 2 := by
    rw [hq2, Real.sin_arcsin (by linarith) (by linarith)]
  have hdisp := z118_dispersion_identity q
  nlinarith [hdisp, hsin]

/-! ### 要素来源表（`#print axioms`） -/

#print axioms z118_dispersion_identity
#print axioms z118_inverse_of_dispersion
#print axioms z118_dispersion_of_inverse

end Zero