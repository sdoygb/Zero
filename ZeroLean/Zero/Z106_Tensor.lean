/-
  Z106（张量层）的机器检查版本：每顶点对称张量场的**无迹投影**与**逐分量色散**。

  两件事（都是有限维线性代数，不含物理）：
    · **TT 无迹投影**：`P h = h − (tr h / n) • 1` 满足 `tr (P h) = 0`、`P (P h) = P h`
      （幂等投影），且对已无迹的 `h` 是恒等 —— 即"无迹守恒"是投影的结构性质；
    · **逐分量色散**：张量场按 `L ⊗ I` 演化（图 Laplacian 逐分量作用），
      若标量场有色散 `L v = λ v`，则张量场 `v ⊗ w` 同满足 `ω² = c²λ`（谱不因张量层而变）。

  数值互证：`ZeroGPU/src/k4_grid_wave.py` / `k4_dispersion.py`（标量版）与 Z106 正文的
  grid6×6 实测（max|tr(h)| = 0.0、色散 ω² = c²λ_k）。
-/

import Mathlib

open scoped Matrix

namespace Zero

/-! ## 一、TT 无迹投影 -/

/-- **Z106 · 无迹投影**：`P h = h − (tr h / n) • 1`（逐顶点作用；`n` 是每顶点张量的阶）。
    TT 规范要求"剪切无体积变化"，其代数实现就是这条投影。 -/
noncomputable def z106TracelessProj (n : ℕ) (h : Matrix (Fin n) (Fin n) ℝ) :
    Matrix (Fin n) (Fin n) ℝ :=
  h - (Matrix.trace h / (n : ℝ)) • 1

/-- **Z106 · 投影结果无迹**：`tr (P h) = 0`。 -/
theorem z106_proj_traceless {n : ℕ} (hn : 0 < n) (h : Matrix (Fin n) (Fin n) ℝ) :
    Matrix.trace (z106TracelessProj n h) = 0 := by
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  rw [z106TracelessProj, Matrix.trace_sub, Matrix.trace_smul, Matrix.trace_one, Fintype.card_fin,
    smul_eq_mul, div_mul_cancel₀ _ hn', sub_self]

/-- **Z106 · 投影幂等**：`P (P h) = P h`。 -/
theorem z106_proj_idempotent {n : ℕ} (hn : 0 < n) (h : Matrix (Fin n) (Fin n) ℝ) :
    z106TracelessProj n (z106TracelessProj n h) = z106TracelessProj n h := by
  have htr : Matrix.trace (z106TracelessProj n h) = 0 := z106_proj_traceless hn h
  rw [z106TracelessProj, htr]
  simp

/-- **Z106 · 已无迹者不动**：`tr h = 0 ⟹ P h = h`（故无迹性**严格守恒**，而不是逐代衰减）。 -/
theorem z106_proj_fixes_traceless {n : ℕ} (h : Matrix (Fin n) (Fin n) ℝ)
    (hh : Matrix.trace h = 0) :
    z106TracelessProj n h = h := by
  rw [z106TracelessProj, hh]
  simp

/-! ## 二、张量层的逐分量色散 -/

/-- **Z106 · 张量场算子 `L ⊗ I`**：对每顶点的 `m` 分量张量场 `H`，`L` 只作用在顶点指标上，
    分量指标不动 —— 这就是正文的"图 Laplacian 逐分量作用"。 -/
noncomputable def z106TensorAct {n m : ℕ} (L : Matrix (Fin n) (Fin n) ℝ)
    (H : Fin n → Fin m → ℝ) : Fin n → Fin m → ℝ :=
  fun i a => ∑ j, L i j * H j a

/-- **Z106 · 逐分量色散**：标量本征态 `L v = λ v` 与任意分量图案 `w` 的张量积 `v ⊗ w`
    是 `L ⊗ I` 的本征态，本征值**仍是** `λ` —— 张量层不改变色散 `ω² = c²λ`。 -/
theorem z106_tensor_eigen {n m : ℕ} (L : Matrix (Fin n) (Fin n) ℝ) (v : Fin n → ℝ)
    (lam : ℝ) (hv : L *ᵥ v = lam • v) (w : Fin m → ℝ) :
    z106TensorAct L (fun i a => v i * w a) = fun i a => lam * (v i * w a) := by
  funext i a
  simp only [z106TensorAct]
  have hsum : (∑ j, L i j * (v j * w a)) = (∑ j, L i j * v j) * w a := by
    rw [Finset.sum_mul]
    exact Finset.sum_congr rfl fun j _ => by ring
  rw [hsum]
  have hvi : (∑ j, L i j * v j) = lam * v i := by
    have := congrFun hv i
    simpa [Matrix.mulVec, dotProduct] using this
  rw [hvi]
  ring

#print axioms z106_proj_traceless
#print axioms z106_proj_idempotent
#print axioms z106_proj_fixes_traceless
#print axioms z106_tensor_eigen

end Zero