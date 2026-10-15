/-
  Z021 §3 的**秩唯一性裁决**（ZeroGPU/doc/coupling_verdict.md §3）的机器检查版本。

  文档 §3.1 的论证（数目级）：
    等長等角族的 Gram 矩阵 `G = (1−ρ)I + ρJ`（`J` = 全 1 矩阵，`m` = 通道数 = `D+1`），
    其特征值为 `1−ρ`（重 `m−1`）与 `1+(m−1)ρ`（重 `1`），于是
      · `ρ = 1`        ⟹ `G = J`，秩 `1`；
      · `ρ = −1/(m−1)` ⟹ 秩 `m−1 = D`（全 1 向量落核）；
      · 其余           ⟹ 秩 `m`（超出口空间 `H_Q` 的维数 `D`）。
    要求秩恰为 `D` 就把 `ρ = −1/D` **唯一**逼出来 —— 在等長等角族内，正则单纯形是唯一的。

  本文件把"其余 ⟹ 秩 m"这一步做成精确的机器检查（这是文档里唯一没有构造性核验的一步）：
    · `z021_gram_mulVec_one`：行和公式 `G·1 = (1+(m−1)ρ)·1`（文档 §3.1 的 λ₂ 特征向量）；
    · `z021_gram_rank_forces_rho`：**主定理**。`ρ ≠ 1` 且秩 `= m−1` ⟹ `ρ = −1/(m−1)`；
    · `z021_gram_kernel_one`：`ρ = −1/(m−1)` 处全 1 向量落核（故秩 `≤ m−1`，上界一侧）。

  证明路线（与文档 §3.1 的秩表一致，但用核-秩定理代替特征值分解）：
    `1+(m−1)ρ ≠ 0` 且 `ρ ≠ 1` 时，`G x = 0` 逐分量给出 `(1−ρ) x i + ρ·S = 0`（`S = ∑ x`），
    故 `x` 是常向量 `x = c·1`；再用行和公式 `c·(1+(m−1)ρ) = 0` 得 `c = 0`，即 `ker G = ⊥`。
    秩-零化度（`LinearMap.finrank_range_add_finrank_ker`）于是给 `rank G = m`，
    与假设 `rank G = m−1` 矛盾（`m ≥ 3`）。故必有 `1+(m−1)ρ = 0`，解出 `ρ = −1/(m−1)`。

  数值互证：`ZeroGPU/doc/coupling_verdict.md` §3.1（三情形秩表）。

  工程注记：
    · `Matrix.rank A` 的定义就是 `finrank R (LinearMap.range A.mulVecLin)`，
      所以核-秩定理可以直接接上，不需要行列式或特征值分解；
    · 全 1 矩阵另取名 `z021Ones`（`z021Gram` 的定义逐字不动）。
-/

import Mathlib
import Zero.Z021_Dict

namespace Zero

open scoped Matrix

/-! ## 〇、Gram 矩阵（等长等角族：对角 `1`、非对角 `ρ`） -/

/-- Gram 矩阵：等长等角族的 G = (1−ρ)I + ρJ（J = 全 1 矩阵），m = 通道数 = D+1 -/
noncomputable def z021Gram (m : ℕ) (ρ : ℝ) : Matrix (Fin m) (Fin m) ℝ :=
  (1 - ρ) • (1 : Matrix (Fin m) (Fin m) ℝ) + ρ • Matrix.of (fun _ _ => (1 : ℝ))

/-! ## 一、全 1 矩阵 `J` 及其作用（`z021Gram` 的 `ρ` 项） -/

/-- 全 1 矩阵 `J`（`z021Gram` 定义里的 `ρ` 项；只为叙述方便另起名字，不改 `z021Gram` 本身）。 -/
noncomputable def z021Ones (m : ℕ) : Matrix (Fin m) (Fin m) ℝ := Matrix.of (fun _ _ => (1 : ℝ))

/-- `z021Gram` 的等价写法 `G = (1−ρ)·I + ρ·J`（按定义即是；文档 §3.1 的 (3.1) 式）。 -/
private theorem z021Gram_eq (m : ℕ) (ρ : ℝ) :
    z021Gram m ρ = (1 - ρ) • (1 : Matrix (Fin m) (Fin m) ℝ) + ρ • z021Ones m := rfl

/-- `J` 的作用：`J·x = (∑_i x_i)·1`（全 1 矩阵把向量映到其分量的常向量）。 -/
theorem z021Ones_mulVec (m : ℕ) (x : Fin m → ℝ) :
    z021Ones m *ᵥ x = (∑ i, x i) • (fun _ : Fin m => (1 : ℝ)) := by
  ext i
  simp only [z021Ones, Matrix.of_apply, Matrix.mulVec, dotProduct, one_mul, Pi.smul_apply,
    smul_eq_mul, mul_one]

/-- `J·1 = m·1`（`J` 作用在全 1 向量上）。 -/
theorem z021Ones_mulVec_one (m : ℕ) :
    z021Ones m *ᵥ (fun _ : Fin m => (1 : ℝ)) = (m : ℝ) • (fun _ : Fin m => (1 : ℝ)) := by
  rw [z021Ones_mulVec]
  congr 1
  simp [Finset.sum_const, Finset.card_univ, Fintype.card_fin]

/-- **Z021 · §3.1 的行和公式（λ₂ 特征向量）**：`G` 作用在全 1 向量上等于 `(1+(m−1)ρ)` 倍的全 1 向量。
    右端的系数就是文档 §3.1 里的特征值 `λ₂ = 1+(m−1)ρ`。 -/
theorem z021_gram_mulVec_one (m : ℕ) (ρ : ℝ) :
    z021Gram m ρ *ᵥ (fun _ : Fin m => (1 : ℝ))
      = (1 + ((m : ℝ) - 1) * ρ) • (fun _ : Fin m => (1 : ℝ)) := by
  rw [z021Gram_eq, Matrix.add_mulVec, Matrix.smul_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec,
    z021Ones_mulVec_one]
  ext i
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, mul_one]
  ring

/-- **Z021 · §3.1 的核向量**：`1+(m−1)ρ = 0` 时全 1 向量在 `G` 的核里，
    故 `rank G ≤ m−1`（文档 §3.1 秩表里 `ρ = −1/(m−1)` 那一行的上界一侧）。 -/
theorem z021_gram_kernel_one (m : ℕ) (ρ : ℝ) (h : 1 + ((m : ℝ) - 1) * ρ = 0) :
    z021Gram m ρ *ᵥ (fun _ : Fin m => (1 : ℝ)) = 0 := by
  rw [z021_gram_mulVec_one, h, zero_smul]

/-! ## 二、核-秩路线：`ρ ≠ 1` 且 `1+(m−1)ρ ≠ 0` 时 `G` 满秩 -/

/-- `G` 的逐分量作用：`(G x)_i = (1−ρ)·x_i + ρ·(∑_j x_j)`。
    这是核向量分析的基本式（把 `G = (1−ρ)I + ρJ` 与 `z021Ones_mulVec` 合起来）。 -/
theorem z021Gram_mulVec (m : ℕ) (ρ : ℝ) (x : Fin m → ℝ) :
    z021Gram m ρ *ᵥ x
      = (1 - ρ) • x + (ρ * ∑ i, x i) • (fun _ : Fin m => (1 : ℝ)) := by
  rw [z021Gram_eq, Matrix.add_mulVec, Matrix.smul_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec,
    z021Ones_mulVec, smul_smul]

/-- **核平凡**：`ρ ≠ 1` 且 `1+(m−1)ρ ≠ 0` 时，`G` 的核只有零向量。
    论证：`G x = 0` ⟹ 逐分量 `(1−ρ)x_i + ρS = 0`（`S = ∑ x`）⟹ `x = c·1`（`ρ ≠ 1`）；
    再用行和公式 `0 = G x = c·(1+(m−1)ρ)·1` 得 `c = 0`（`1+(m−1)ρ ≠ 0`）。 -/
theorem z021Gram_ker_eq_bot (m : ℕ) (hm : 0 < m) (ρ : ℝ) (hρ : ρ ≠ 1)
    (h : 1 + ((m : ℝ) - 1) * ρ ≠ 0) :
    LinearMap.ker (z021Gram m ρ).mulVecLin = ⊥ := by
  rw [Submodule.eq_bot_iff]
  intro x hx
  rw [LinearMap.mem_ker, Matrix.mulVecLin_apply] at hx
  -- 展开形式（逐分量分析用）与折叠形式 `hx` 各留一份
  have hx' : (1 - ρ) • x + (ρ * ∑ i, x i) • (fun _ : Fin m => (1 : ℝ)) = 0 := by
    rw [← z021Gram_mulVec]
    exact hx
  have h1ρ : (1 : ℝ) - ρ ≠ 0 := by
    intro hh
    exact hρ (by linarith)
  -- 逐分量：`(1−ρ)·x_i + ρ·S = 0`
  have hcomp : ∀ i, (1 - ρ) * x i + ρ * (∑ j, x j) = 0 := by
    intro i
    have hi := congrFun hx' i
    simpa only [Pi.add_apply, Pi.smul_apply, Pi.zero_apply, smul_eq_mul, mul_one] using hi
  -- 于是 `x` 是常向量
  set c : ℝ := -(ρ * ∑ j, x j) / (1 - ρ) with hc
  have hxc : ∀ i, x i = c := by
    intro i
    rw [hc, eq_div_iff h1ρ]
    linarith [hcomp i]
  have hxfun : x = c • (fun _ : Fin m => (1 : ℝ)) := by
    funext i
    simpa only [Pi.smul_apply, smul_eq_mul, mul_one] using hxc i
  -- 行和公式把常数逼成 0
  have hzero : c * (1 + ((m : ℝ) - 1) * ρ) = 0 := by
    have h1 : (c * (1 + ((m : ℝ) - 1) * ρ)) • (fun _ : Fin m => (1 : ℝ)) = 0 := by
      rw [← smul_smul, ← z021_gram_mulVec_one m ρ, ← Matrix.mulVec_smul, ← hxfun]
      exact hx
    have h2 := congrFun h1 ⟨0, hm⟩
    simpa only [Pi.smul_apply, Pi.zero_apply, smul_eq_mul, mul_one] using h2
  rcases mul_eq_zero.mp hzero with hc0 | hA
  · rw [hxfun, hc0, zero_smul]
  · exact absurd hA h

/-- **满秩**：`ρ ≠ 1` 且 `1+(m−1)ρ ≠ 0` 时 `rank G = m`（文档 §3.1 秩表"其余 ⟹ 秩 m"一行）。
    由核平凡 + 秩-零化度定理 `finrank range + finrank ker = finrank (Fin m → ℝ) = m` 得到。 -/
theorem z021Gram_rank_eq_card (m : ℕ) (hm : 0 < m) (ρ : ℝ) (hρ : ρ ≠ 1)
    (h : 1 + ((m : ℝ) - 1) * ρ ≠ 0) : (z021Gram m ρ).rank = m := by
  have hker := z021Gram_ker_eq_bot m hm ρ hρ h
  have hrn := LinearMap.finrank_range_add_finrank_ker (z021Gram m ρ).mulVecLin
  rw [hker, finrank_bot, add_zero, Module.finrank_fintype_fun_eq_card, Fintype.card_fin] at hrn
  show Module.finrank ℝ (LinearMap.range (z021Gram m ρ).mulVecLin) = m
  exact hrn

/-! ## 三、唯一性：秩 `m−1` 逼出 `ρ = −1/(m−1)`（文档 §3.1 的方框） -/

/-- **唯一性（coupling_verdict §3）**：`m ≥ 3`、`ρ ≠ 1`，若 Gram 矩阵的秩恰为 `m−1`，
    则 `ρ = −1/(m−1)`（即 `D = m−1` 时的 `−1/D`）。 -/
theorem z021_gram_rank_forces_rho (m : ℕ) (hm : 3 ≤ m) (ρ : ℝ) (hρ : ρ ≠ 1)
    (hrank : (z021Gram m ρ).rank = m - 1) :
    ρ = -1 / ((m : ℝ) - 1) := by
  have hm0 : 0 < m := by omega
  have hm1 : (m : ℝ) - 1 ≠ 0 := by
    have h2 : (1 : ℝ) < (m : ℝ) := by exact_mod_cast (by omega : 1 < m)
    linarith
  by_cases h : 1 + ((m : ℝ) - 1) * ρ = 0
  · rw [eq_div_iff hm1]
    linarith
  · exfalso
    have hcard := z021Gram_rank_eq_card m hm0 ρ hρ h
    rw [hrank] at hcard
    omega

#print axioms z021Ones_mulVec
#print axioms z021Ones_mulVec_one
#print axioms z021_gram_mulVec_one
#print axioms z021_gram_kernel_one
#print axioms z021Gram_mulVec
#print axioms z021Gram_ker_eq_bot
#print axioms z021Gram_rank_eq_card
#print axioms z021_gram_rank_forces_rho

end Zero