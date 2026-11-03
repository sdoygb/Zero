/-
  Zb041（量子特征锚定与超荷内生 · 统一群边界）第三节「K₄ 对称性分解 → 粒子种类」的机器检查版本。

  要验的命题：K₄ 邻接矩阵 `A`（对角 0、其余 1）的谱为 `{3, −1×3}`，且 −1 特征空间（3 维）
  在三循环 `C₃`（顶点 0 固定，1→2→3→1，即 Zb041 说的"3 轴旋转"）下分解为 `1 ⊕ ω ⊕ ω²`，
  其中 `ω = exp(2πi/3)`；`(0,1,ω,ω²)` 是 `A` 的 −1 特征向量，并按 `ω` 变换。

  数值核验（本机实跑，numpy）：
    · `ZeroGPU/src/zero_sim.py` 结构层的 `k4` / `adj_matrix` / `spectrum`：K₄ 谱 = `{3, −1, −1, −1}`；
    · `ZeroGPU/src/k4_anchors.py`（Zb041 关联脚本）：−1 特征空间维数 3，`C₃` 作用其上的
      表示特征值 = `{1, ω, ω²}` = `{1, −1/2+0.866i, −1/2−0.866i}`；
    · `A *ᵥ (0,1,ω,ω²) + (0,1,ω,ω²) = 0`（残差 `0`），且 `C₃` 作用后 = `ω · (0,1,ω,ω²)`；
      `(3,−1,−1,−1)` 同为 −1 特征向量但在 `C₃` 下不动（平凡分量 `1`）。

  ⚠ 审计修订（2026-10-15，本次机器核验，**未改动任何给定定理的名字与陈述**）：
    任务书原本列出第 4 条 `zb041_minus_one_shape`：
      `(v : Fin 4 → ℂ) (hv : zb041K4Adj *ᵥ v = -v) : v 0 = 0 ∧ v 1 + v 2 + v 3 = 0`
    这条**陈述为假**，本文件用 `zb041_minus_one_shape_refuted` 给出机器核验的反例
    （`v = (−1,1,1,−1)`：`A *ᵥ v = −v` 成立，但 `v 0 = −1 ≠ 0`，且 `v 1+v 2+v 3 = 1 ≠ 0`）。
    成因：`(A *ᵥ v) i = (∑ j, v j) − v i`，故 `A *ᵥ v = −v ⟺ ∑ j, v j = 0`（`zb041_minus_one_space`）——
    −1 特征空间是"分量和为零"的 3 维超平面，**不是** `{(0,a,b,c) : a+b+c=0}`（后者只有 2 维，
    只是其中 `ω/ω²` 那两个一维分量）。Zb041 正文第三节的原话是"−1 特征空间**在 C₃ 下分解**为
    `{1, ω, ω²}`"，与本文件一致；`v 0 = 0` 只对 `ω`、`ω²` 分量成立，见 `zb041_minus_one_shape_of_C3`
    （新增条目，非任务书条目）。因此该条不写入本文件（写入必然报错），以反例 + 正确刻画代替；
    其余各条名字、陈述与任务书逐字一致。

  ⚠ 二次审计整改（2026-10-15）：`zb041_minus_one_span` 名实不符（名字说 "span"，
    内容实为 `LinearIndependent`，且第三向量 `(0,1,1,1)` 不在 −1 空间里）⇒
    按审计**改名为** `zb041_minus_one_linearIndependent`：**只改名字，陈述逐字未动**
    （文件内引用与 `#print axioms` 行同步更新）。其余各条名字与陈述均未改动。
-/

import Mathlib

open scoped Matrix

namespace Zero

/-- `ω = exp(2πi/3)` -/
noncomputable def zb041Omega : ℂ := Complex.exp (2 * Real.pi * Complex.I / 3)

/-! ### 技术引理（private，不进 `Zero.*` 审计空间） -/

/-- `ω` 是 3 次本原单位根：`Complex.isPrimitiveRoot_exp` 直接给出（`exp(2πi/n)` 是 `n` 次本原根）。 -/
private lemma zb041Omega_isPrimitiveRoot : IsPrimitiveRoot zb041Omega 3 := by
  have h := Complex.isPrimitiveRoot_exp 3 (by norm_num)
  simpa [zb041Omega] using h

/-- `1 + ω + ω² = 0`：由 `(ω−1)(ω²+ω+1) = ω³−1 = 0` 与 `ω ≠ 1` 得到。 -/
private lemma zb041_omega_sum_aux : 1 + zb041Omega + zb041Omega ^ 2 = 0 := by
  have hfac : (zb041Omega - 1) * (zb041Omega ^ 2 + zb041Omega + 1) = 0 := by
    linear_combination zb041Omega_isPrimitiveRoot.pow_eq_one
  have hsub : zb041Omega - 1 ≠ 0 :=
    sub_ne_zero.mpr (zb041Omega_isPrimitiveRoot.ne_one (by norm_num))
  have hq : zb041Omega ^ 2 + zb041Omega + 1 = 0 := by
    rcases mul_eq_zero.mp hfac with h | h
    · exact absurd h hsub
    · exact h
  linear_combination hq

/-- `ω² ≠ ω`（供线性无关证明用）：否则 `ω(ω−1) = 0`，而 `ω ≠ 0`、`ω ≠ 1`。 -/
private lemma zb041Omega_ne_sq : zb041Omega ≠ zb041Omega ^ 2 := by
  intro h
  have hzero : zb041Omega * (zb041Omega - 1) = 0 := by linear_combination -h
  rcases mul_eq_zero.mp hzero with h1 | h1
  · exact zb041Omega_isPrimitiveRoot.ne_zero (by norm_num) h1
  · exact zb041Omega_isPrimitiveRoot.ne_one (by norm_num) (sub_eq_zero.mp h1)

/-- `ω³ = 1`（`ω` 是 3 次本原根 ⟹ 是 `X³ − 1` 的根）。 -/
theorem zb041_omega_cube : zb041Omega ^ 3 = 1 :=
  zb041Omega_isPrimitiveRoot.pow_eq_one

/-- `ω ≠ 1`：3 次本原根的定义正是"阶恰为 3"，故不是 1。 -/
theorem zb041_omega_ne_one : zb041Omega ≠ 1 :=
  zb041Omega_isPrimitiveRoot.ne_one (by norm_num)

/-- `1 + ω + ω² = 0`：`Z₃` 特征标的正交关系（`1 ⊕ ω ⊕ ω²` 里平凡表示的和为零）。 -/
theorem zb041_omega_sum : 1 + zb041Omega + zb041Omega ^ 2 = 0 :=
  zb041_omega_sum_aux

/-- K₄ 邻接矩阵 -/
def zb041K4Adj : Matrix (Fin 4) (Fin 4) ℂ := Matrix.of (fun i j => if i = j then 0 else 1)

/-- **(0,1,ω,ω²) 是 −1 特征向量** —
    逐分量算：`(A *ᵥ v) i = (∑ j, v j) − v i`，而 `∑ j, v j = 1 + ω + ω² = 0`。 -/
theorem zb041_minus_one_eigen :
    zb041K4Adj *ᵥ ![0, 1, zb041Omega, zb041Omega ^ 2] = -![0, 1, zb041Omega, zb041Omega ^ 2] := by
  funext i
  fin_cases i <;>
    simp [zb041K4Adj, Matrix.mulVec, dotProduct, Fin.sum_univ_four] <;>
    linear_combination zb041_omega_sum

/-- **分量式**（新增条目，非任务书条目）：`(A *ᵥ v) i = (∑ j, v j) − v i` ——
    这是 −1 特征空间之所以是"分量和为零"超平面的**全部原因**，也是本文件后续各条的承重步骤。 -/
theorem zb041_mulVec_apply (v : Fin 4 → ℂ) (i : Fin 4) :
    (zb041K4Adj *ᵥ v) i = (v 0 + v 1 + v 2 + v 3) - v i := by
  fin_cases i <;>
    simp [zb041K4Adj, Matrix.mulVec, dotProduct, Fin.sum_univ_four] <;>
    ring

/-- **−1 特征空间的正确刻画**（新增条目，非任务书条目）：
    `A *ᵥ v = −v ⟺ v 0 + v 1 + v 2 + v 3 = 0`。
    由分量式 `zb041_mulVec_apply`：`A *ᵥ v = −v ⟺ (∑ v) − v i = −v i (∀i) ⟺ ∑ v = 0`。 -/
theorem zb041_minus_one_space (v : Fin 4 → ℂ) :
    zb041K4Adj *ᵥ v = -v ↔ v 0 + v 1 + v 2 + v 3 = 0 := by
  constructor
  · intro hv
    have h := congrFun hv 0
    rw [zb041_mulVec_apply] at h
    simp at h
    linear_combination h
  · intro hv
    funext i
    rw [zb041_mulVec_apply, hv]
    simp

/-- **反例：任务书第 4 条 `zb041_minus_one_shape` 的陈述为假**（新增条目，非任务书条目）。
    `v = (−1,1,1,−1)` 满足 `A *ᵥ v = −v`（分量和 = 0），但 `v 0 = −1 ≠ 0`
    （且 `v 1 + v 2 + v 3 = 1 ≠ 0`）——故 `A *ᵥ v = −v → v 0 = 0 ∧ v 1+v 2+v 3 = 0` 不成立。
    这正是"−1 空间 ⊋ ω/ω² 分量"的显式见证。 -/
theorem zb041_minus_one_shape_refuted :
    ∃ v : Fin 4 → ℂ, zb041K4Adj *ᵥ v = -v ∧ v 0 ≠ 0 := by
  refine ⟨![-1, 1, 1, -1], ?_, ?_⟩
  · funext i
    fin_cases i <;> simp [zb041K4Adj, Matrix.mulVec, dotProduct, Fin.sum_univ_four]
  · norm_num

/-- 三循环 `C₃`：顶点 0 固定，`1 → 2 → 3 → 1`（Zb041 的"3 轴旋转"）。 -/
def zb041Cycle3 : Fin 4 → Fin 4 := ![0, 2, 3, 1]

/-- **修正版形状结论**（新增条目，非任务书条目；**不是**原陈述的替身，原陈述 `zb041_minus_one_shape` 已由
    `zb041_minus_one_shape_refuted` 否决，该名字**不复用**）：
    加 `C₃` 假设后形状结论成立 —— 若 `v` 既在 −1 特征空间里、又按 `ω` 变换
    （`v (c i) = ω * v i`，即落在 `ω` 分量上），则 `v 0 = 0 ∧ v 1 + v 2 + v 3 = 0`。
    证明：分量 0 处 `v 0 = ω * v 0` ⟹ `(1−ω) v 0 = 0` ⟹ `v 0 = 0`（因 `ω ≠ 1`）；
    再由 `zb041_minus_one_space` 得 `v 1 + v 2 + v 3 = 0`。 -/
theorem zb041_minus_one_shape_of_C3 (v : Fin 4 → ℂ)
    (hv : zb041K4Adj *ᵥ v = -v) (hc : ∀ i, v (zb041Cycle3 i) = zb041Omega * v i) :
    v 0 = 0 ∧ v 1 + v 2 + v 3 = 0 := by
  have hsum : v 0 + v 1 + v 2 + v 3 = 0 := (zb041_minus_one_space v).mp hv
  have h0 : v 0 = 0 := by
    have hc0 : v 0 = zb041Omega * v 0 := by
      simpa [zb041Cycle3] using hc 0
    have hz : (1 - zb041Omega) * v 0 = 0 := by linear_combination hc0
    rcases mul_eq_zero.mp hz with h | h
    · exact absurd (sub_eq_zero.mp h).symm zb041_omega_ne_one
    · exact h
  exact ⟨h0, by linear_combination hsum - h0⟩

/-- **平凡分量 `1`**（新增条目，非任务书条目）：`(3,−1,−1,−1)` 是 −1 特征向量，且在 `C₃` 下不动 ——
    这就是 Zb041 第三节里"`C₃` 平凡表示 → 色单态"的那个一维分量。 -/
theorem zb041_trivial_component :
    zb041K4Adj *ᵥ ![3, -1, -1, -1] = -![3, -1, -1, -1] ∧
      (∀ i, (![3, -1, -1, -1] : Fin 4 → ℂ) (zb041Cycle3 i) = ![3, -1, -1, -1] i) := by
  constructor
  · funext i
    fin_cases i <;> simp [zb041K4Adj, Matrix.mulVec, dotProduct, Fin.sum_univ_four] <;> norm_num
  · intro i
    fin_cases i <;> simp [zb041Cycle3]

/-- **Z₃ 分解**：(0,1,ω,ω²) 与它的共轭、以及 (0,1,1,1) 三者**线性无关**（任务书原文即 `LinearIndependent` 形式）。
    证明用 `Fintype.linearIndependent_iff` 取系数向量 `g`：
    分量 1/2/3 给出 `g 0+g 1+g 2 = 0`、`g 0·ω + g 1·ω² + g 2 = 0`、`g 0·ω² + g 1·ω + g 2 = 0`；
    三式相加并减去 `g 0·(1+ω+ω²)`、`g 1·(1+ω+ω²)` 得 `3 g 2 = 0` ⟹ `g 2 = 0`；
    再由 `g 0 = g 1`（第 2 式减第 3 式 + `ω ≠ ω²`）与第 1 式得 `2 g 0 = 0` ⟹ `g 0 = g 1 = 0`。

    ⚠ 口径澄清：这条证的是"三者**线性无关**"，**不是**"三者**张成** −1 空间" ——
    第三向量 `(0,1,1,1)` 压根不在 −1 空间里（分量和 = `3 ≠ 0`，见
    `zb041_third_vector_not_in_minus_one`），它落在 `3` 特征方向上；前两个才在 `ω/ω²` 分量里。
    真正张成 −1 空间的是 `(0,1,ω,ω²)`、`(0,1,ω²,ω)` 加**平凡分量** `(3,−1,−1,−1)`
    （`zb041_trivial_component`），三条各自独立。

    ⚠ 改名（2026-10-15 二次审计整改）：本条原名 `zb041_minus_one_span` **名实不符**
    （名字说 span，内容却是 `LinearIndependent`），已按审计改名为
    `zb041_minus_one_linearIndependent` —— **陈述逐字未动**，只改名字。 -/
theorem zb041_minus_one_linearIndependent :
    LinearIndependent ℂ ![![0, 1, zb041Omega, zb041Omega ^ 2],
                         ![0, 1, zb041Omega ^ 2, zb041Omega],
                         ![0, 1, 1, 1]] := by
  rw [Fintype.linearIndependent_iff]
  intro g hg
  have e1 := congrFun hg 1
  have e2 := congrFun hg 2
  have e3 := congrFun hg 3
  simp only [Fin.sum_univ_three] at e1 e2 e3
  simp at e1 e2 e3
  -- e1 : g 0 + g 1 + g 2 = 0
  -- e2 : g 0 * ω + g 1 * ω² + g 2 = 0
  -- e3 : g 0 * ω² + g 1 * ω + g 2 = 0
  have hg2 : g 2 = 0 := by
    have h3 : 3 * g 2 = 0 := by
      linear_combination e1 + e2 + e3 - (g 0) * zb041_omega_sum - (g 1) * zb041_omega_sum
    rcases mul_eq_zero.mp h3 with h | h
    · norm_num at h
    · exact h
  have h01 : g 0 = g 1 := by
    have hd : (g 0 - g 1) * (zb041Omega - zb041Omega ^ 2) = 0 := by
      linear_combination e2 - e3
    rcases mul_eq_zero.mp hd with h | h
    · linear_combination h
    · exact absurd (sub_eq_zero.mp h) zb041Omega_ne_sq
  have hg0 : g 0 = 0 := by
    have h2 : 2 * g 0 = 0 := by linear_combination e1 + h01 - hg2
    rcases mul_eq_zero.mp h2 with h | h
    · norm_num at h
    · exact h
  have hg1 : g 1 = 0 := by linear_combination hg0 - h01
  intro i
  fin_cases i
  · exact hg0
  · exact hg1
  · exact hg2

/-- **任务书第三个向量不在 −1 空间里**（新增条目，非任务书条目）：`(0,1,1,1)` 的分量和为 `3`，
    故 `A *ᵥ (0,1,1,1) = (3,2,2,2) ≠ −(0,1,1,1)` ——
    所以 `zb041_minus_one_linearIndependent` 证的是"三者线性无关"，不是"三者张成 −1 空间"。 -/
theorem zb041_third_vector_not_in_minus_one :
    zb041K4Adj *ᵥ ![0, 1, 1, 1] ≠ -![0, 1, 1, 1] := by
  intro h
  have h0 := congrFun h 0
  norm_num [zb041_mulVec_apply] at h0

/-! ### 要素来源表（`#print axioms`） -/

#print axioms zb041_omega_cube
#print axioms zb041_omega_ne_one
#print axioms zb041_omega_sum
#print axioms zb041_minus_one_eigen
#print axioms zb041_mulVec_apply

/-! ## 共扼闭包（复共扼 ω ↔ ω²） -/

/-- B：`ω²` 也是三次单位根。 -/
theorem zb041_omega_sq_cube : (zb041Omega ^ 2) ^ 3 = 1 := by
  have h : (zb041Omega ^ 2) ^ 3 = (zb041Omega ^ 3) ^ 2 := by ring
  rw [h, zb041_omega_cube]
  norm_num

/-- B：`ω² ≠ 1`。 -/
theorem zb041_omega_sq_ne_one : zb041Omega ^ 2 ≠ 1 := by
  intro h
  have h3 : zb041Omega = 1 := by
    calc zb041Omega = zb041Omega * 1 := by ring
      _ = zb041Omega * zb041Omega ^ 2 := by rw [h]
      _ = zb041Omega ^ 3 := by ring
      _ = 1 := zb041_omega_cube
  exact zb041_omega_ne_one h3

/-- B：`1 + ω² + (ω²)² = 0`（`ω²` 版正交关系）。 -/
theorem zb041_omega_sum_conj : 1 + zb041Omega ^ 2 + (zb041Omega ^ 2) ^ 2 = 0 := by
  have h4 : (zb041Omega ^ 2) ^ 2 = zb041Omega := by
    calc (zb041Omega ^ 2) ^ 2 = zb041Omega ^ 4 := by ring
      _ = zb041Omega ^ 3 * zb041Omega := by ring
      _ = 1 * zb041Omega := by rw [zb041_omega_cube]
      _ = zb041Omega := by ring
  rw [h4]
  simpa [add_comm, add_left_comm, add_assoc] using zb041_omega_sum

/-- B：复共扼把 `-1` 特征向量 `![0,1,ω,ω²]` 送到 `![0,1,ω²,ω]`。 -/
theorem zb041_minus_one_eigen_conj :
    zb041K4Adj.mulVec ![0, 1, zb041Omega ^ 2, zb041Omega] =
      -![0, 1, zb041Omega ^ 2, zb041Omega] := by
  funext i
  fin_cases i <;>
    simp [zb041K4Adj, Matrix.mulVec, dotProduct, Fin.sum_univ_four] <;>
    linear_combination zb041_omega_sum

/-- A：复共扼把线性无关三元组 `(v₊,v₋,t)` 换成 `(v₋,v₊,t)`，仍线性无关。 -/
theorem zb041_linearIndependent_conj :
    LinearIndependent ℂ
      ![![0, 1, zb041Omega ^ 2, zb041Omega],
        ![0, 1, zb041Omega, zb041Omega ^ 2],
        ![0, 1, 1, 1]] := by
  have h := zb041_minus_one_linearIndependent
  have h' := h.comp (Equiv.swap (0 : Fin 3) 1) (Equiv.injective _)
  simpa using h'

#print axioms zb041_omega_sq_cube
#print axioms zb041_omega_sq_ne_one
#print axioms zb041_omega_sum_conj
#print axioms zb041_minus_one_eigen_conj
#print axioms zb041_linearIndependent_conj


/-- B：`star ω = ω²`（复共扼把 ω 送到 ω²）。 -/
theorem zb041_star_omega : star zb041Omega = zb041Omega ^ 2 := by
  have h1 : star zb041Omega = Complex.exp (-(2 * Real.pi * Complex.I / 3)) := by
    rw [zb041Omega]
    change starRingEnd ℂ (Complex.exp (2 * Real.pi * Complex.I / 3)) =
      Complex.exp (-(2 * Real.pi * Complex.I / 3))
    have h : starRingEnd ℂ (Complex.exp (2 * Real.pi * Complex.I / 3)) =
        Complex.exp (starRingEnd ℂ (2 * Real.pi * Complex.I / 3)) := by
      simpa using (Complex.exp_conj (2 * Real.pi * Complex.I / 3)).symm
    rw [h]
    congr 1
    simp only [map_div₀, map_mul, Complex.conj_I, map_ofNat, Complex.conj_ofReal]
    ring
  have h2 : zb041Omega ^ 2 = Complex.exp (4 * Real.pi * Complex.I / 3) := by
    rw [zb041Omega, ← Complex.exp_nat_mul]
    congr 1
    ring
  have h2pi : Complex.exp (2 * Real.pi * Complex.I) = 1 := by
    rw [Complex.exp_eq_one_iff]
    exact ⟨1, by norm_num⟩
  rw [h1, h2]
  rw [show (4 * Real.pi * Complex.I / 3) =
      -(2 * Real.pi * Complex.I / 3) + 2 * Real.pi * Complex.I by ring]
  rw [Complex.exp_add, h2pi, mul_one]

/-- B：`star (ω²) = ω`。 -/
theorem zb041_star_omega_sq : star (zb041Omega ^ 2) = zb041Omega := by
  rw [star_pow, zb041_star_omega]
  calc (zb041Omega ^ 2) ^ 2 = zb041Omega ^ 4 := by ring
    _ = zb041Omega ^ 3 * zb041Omega := by ring
    _ = 1 * zb041Omega := by rw [zb041_omega_cube]
    _ = zb041Omega := by ring

/-- A：`star` 与 K₄ 邻接矩阵交换（A 实）。 -/
theorem zb041_star_mulVec (v : Fin 4 → ℂ) :
    star (zb041K4Adj.mulVec v) = zb041K4Adj.mulVec (fun i => star (v i)) := by
  funext i
  fin_cases i <;>
    simp [zb041K4Adj, Matrix.mulVec, dotProduct, Fin.sum_univ_four, star_sum, star_mul,
      mul_comm, mul_left_comm, mul_assoc]

/-- C：`ω²` 条件分支——复共扼把 `ω²` 条件搬成 `ω` 条件，再用原定理。 -/
theorem zb041_minus_one_shape_of_C3_conj (v : Fin 4 → ℂ)
    (hv : zb041K4Adj.mulVec v = -v)
    (hc : ∀ i, v (zb041Cycle3 i) = zb041Omega ^ 2 * v i) :
    v 0 = 0 ∧ v 1 + v 2 + v 3 = 0 := by
  let w : Fin 4 → ℂ := fun i => star (v i)
  have hw_eigen : zb041K4Adj.mulVec w = -w := by
    have h := congrArg star hv
    rw [zb041_star_mulVec, star_neg] at h
    change zb041K4Adj.mulVec (fun i => star (v i)) = -(fun i => star (v i)) at h
    simpa [w] using h
  have hsq : (starRingEnd ℂ zb041Omega) ^ 2 = zb041Omega := by
    change (star zb041Omega) ^ 2 = zb041Omega
    rw [← star_pow, zb041_star_omega_sq]
  have hw_C3 : ∀ i, w (zb041Cycle3 i) = zb041Omega * w i := by
    intro i
    have h := hc i
    simp [w, h, star_mul, hsq, mul_comm]
  have h := zb041_minus_one_shape_of_C3 w hw_eigen hw_C3
  constructor
  · have := congrArg star h.1
    simpa [w] using this
  · have := congrArg star h.2
    simpa [w, map_add] using this

#print axioms zb041_minus_one_space
#print axioms zb041_minus_one_shape_refuted
#print axioms zb041_minus_one_shape_of_C3
#print axioms zb041_third_vector_not_in_minus_one
#print axioms zb041_trivial_component
#print axioms zb041_minus_one_linearIndependent
#print axioms zb041_star_omega
#print axioms zb041_omega_sq_cube
#print axioms zb041_omega_sum_conj
#print axioms zb041_minus_one_eigen_conj
#print axioms zb041_linearIndependent_conj
#print axioms zb041_minus_one_shape_of_C3_conj

end Zero