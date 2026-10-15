/-
  Z105（"Zero 内生标量引力波：色散普适性与辐射功率 P ∝ A²"）的机器检查版本：
  **图上波的色散关系 ω² = c²λ_k**。

  对应 Z105 的哪一句：Z105 把 L1 历史记录层里的势取为按 **L 本征模式**演化的量
  （`L v = λ_k v`，`L` 为图 Laplacian），演化方程是波动方程 `∂ₜ²φ = −c² L φ`；
  Z105 §2 的结论"ω² = c²λ_k 色散在 2D 格点图普适成立"在正文里被点为**数学定理**：
  本征模式解满足 `ω² = c²λ_k`，数值（蛙跳）只是忠实实现它。本文件就是这条数学定理的
  机器核验，且两条方向都证：

  * `z105_wave_solution` ↔ **充分性**：`L v = λ v`、`ω² = c²λ` ⟹ `u(t) = cos(ωt)•v`
    确实满足 `ü = −c² L u`（Z105 正文"本征模式 e^{iωt}φ_k 满足 ω²=c²λ_k"的实余弦版；
    L 实线性，故实版蕴含复指数版）。
  * `z105_dispersion_relation` ↔ **必要性（色散关系本身，即"色散普适性"）**：
    `v ≠ 0` 且驻波拟设**确实是**解 ⟹ 必有 `ω² = c²λ`。
    取 `t = 0`：左端 `= (−ω²) • v`，右端 `= (−c²) • (L v) = (−c²λ) • v`；
    由 `v ≠ 0` 逐坐标消去标量即得。**全程只用到"`v` 是 `L` 的本征向量"这一条**——
    不涉及图的形状、顶点数、紫外截断或 `c²` 的取值，这正是 Z105 说的"不是 path16 特例，
    是 L 驱动动力学的普适性质"在推导层面的确切含义。

  两条导数引理（`z105_deriv_wave` / `z105_deriv2_wave`）把 `u` 的一阶、二阶导写成显式函数，
  供上面两条主定理共用；它们也是"ω² 从时间二阶导数来"这一步的机器证据。

  用到的 Mathlib 引理（注释中逐条标注）：
  * `Real.hasDerivAt_cos` / `Real.hasDerivAt_sin`：`cos`、`sin` 在一点可导及其导数。
  * `hasDerivAt_id` + `HasDerivAt.const_mul`：内层 `t ↦ ω * t` 的导数为 `ω`。
  * `HasDerivAt.comp`：链式法则，得 `d/dt cos(ω t) = −ω sin(ω t)`、`d/dt sin(ω t) = ω cos(ω t)`。
  * `HasDerivAt.smul_const`：**向量值**情形 `c(t) • v` 的导数是 `c'(t) • v`；
    这是把 `deriv` 直接用在 `ℝⁿ`-值函数 `u` 上的关键（无需绕道 Fréchet 导数）。
  * `HasDerivAt.deriv`：`HasDerivAt f f' x → deriv f x = f'`；配合 `funext` 得函数相等。
  * `HasDerivAt.neg` / `HasDerivAt.congr_deriv`：二阶导数化形（`−(ω(ωcos)) = −(ω²cos)`）。
  * `Matrix.mulVec_smul`：`M *ᵥ (b • v) = b • (M *ᵥ v)`（把 `cos(ωt)` 提到 `L` 外面）。
  * `smul_smul`（`a • b • v = (a*b) • v`）、`Pi.smul_apply`（逐坐标展开 `•`）。
  * `mul_right_cancel₀`：由坐标 `v i ≠ 0` 从 `a * v i = c * v i` 得 `a = c`。
  * `Function.ne_iff`：`v ≠ 0` 给出坐标 `i` 使 `v i ≠ 0`；`neg_injective` 收尾。

  注：文件里 `L` 是**任意**实矩阵（不限于图 Laplacian），结论因此更强——
  色散关系只用"`v` 是 `L` 的本征向量"这一条，与 Z105 的"普适性"主张一致。
-/

import Mathlib

open scoped Matrix

namespace Zero

/-! ### 辅助：标量链式法则的两条导数（`cos(ω·t)` 与其负 `ω sin(ω·t)`）

本文件自用的技术引理一律写 `private lemma`（名字变 `_private.…`）：
不进 `Zero.*` 审计空间，也不用在 `#print axioms` 里交代。 -/

/-- 辅助引理：`t ↦ cos (ω * t)` 在 `s` 处可导，导数为 `−(ω sin (ω s))`。 -/
private lemma deriv_cos_mul (ω s : ℝ) :
    HasDerivAt (fun t : ℝ => Real.cos (ω * t)) (-(ω * Real.sin (ω * s))) s := by
  -- 链式法则：`cos` 在 `ω * s` 处可导，内层 `t ↦ ω * t` 的导数为 `ω`
  have h1 := (Real.hasDerivAt_cos (ω * s)).comp s ((hasDerivAt_id s).const_mul ω)
  -- `comp` 给出 `−sin(ωs) * ω`，与 `−(ω sin(ωs))` 用交换律对齐
  exact h1.congr_deriv (by ring)

/-- 辅助引理：`t ↦ −(ω sin (ω * t))` 在 `s` 处可导，导数为 `−(ω² cos (ω s))`。 -/
private lemma deriv_neg_sin_mul (ω s : ℝ) :
    HasDerivAt (fun t : ℝ => -(ω * Real.sin (ω * t))) (-(ω ^ 2 * Real.cos (ω * s))) s := by
  -- `d/dt sin(ω t) = ω cos(ω t)`（链式法则）
  have h1 := (Real.hasDerivAt_sin (ω * s)).comp s ((hasDerivAt_id s).const_mul ω)
  have h1' : HasDerivAt (fun t : ℝ => Real.sin (ω * t)) (ω * Real.cos (ω * s)) s :=
    h1.congr_deriv (by ring)
  -- 乘常数 `ω`、取负，得 `d/dt (−(ω sin(ω t))) = −(ω·(ω cos(ω t)))`
  exact (h1'.const_mul ω).neg.congr_deriv (by ring)

/-! ### Z105 的两条导数引理（把 `deriv` 写成显式函数） -/

/-- 一阶导函数（把 `deriv` 写成显式函数，供两条主定理共用） -/
theorem z105_deriv_wave {n : ℕ} (v : Fin n → ℝ) (ω : ℝ) :
    deriv (fun s : ℝ => Real.cos (ω * s) • v)
      = fun s : ℝ => (-(ω * Real.sin (ω * s))) • v := by
  funext s
  -- `HasDerivAt.smul_const`：向量值 `c(t) • v` 的导数逐分量就是 `c'(t) • v`
  exact (deriv_cos_mul ω s).smul_const v |>.deriv

/-- 二阶导函数 -/
theorem z105_deriv2_wave {n : ℕ} (v : Fin n → ℝ) (ω : ℝ) :
    deriv (fun s : ℝ => (-(ω * Real.sin (ω * s))) • v)
      = fun s : ℝ => (-(ω ^ 2 * Real.cos (ω * s))) • v := by
  funext s
  exact (deriv_neg_sin_mul ω s).smul_const v |>.deriv

/-! ### Z105 主定理：色散关系的充分性与必要性 -/

/-- **Z105 · 充分性**：若 `L v = lam • v` 且 `ω² = c²·lam`，则 `u(t) = cos(ω t) • v`
    满足波动方程 `ü = −c²·L u` -/
theorem z105_wave_solution {n : ℕ} (L : Matrix (Fin n) (Fin n) ℝ) (v : Fin n → ℝ)
    (lam c ω : ℝ) (hv : L *ᵥ v = lam • v) (hω : ω ^ 2 = c ^ 2 * lam) :
    deriv (deriv (fun s : ℝ => Real.cos (ω * s) • v))
      = fun t : ℝ => (-(c ^ 2)) • (L *ᵥ (Real.cos (ω * t) • v)) := by
  -- 两次求导：`ü = −ω² cos(ωt) • v`
  rw [z105_deriv_wave, z105_deriv2_wave]
  funext t
  -- 右端：`L *ᵥ (cos(ωt) • v) = cos(ωt) • (L *ᵥ v) = cos(ωt) • (lam • v)`
  simp only [Matrix.mulVec_smul, hv, smul_smul]
  -- 两边都是 `v` 的标量倍，只需标量相等；用本征方程给的 `ω² = c²lam`
  congr 1
  rw [hω]
  ring

/-- **Z105 · 必要性（色散关系）**：若 `v ≠ 0`、`L v = lam • v`，且驻波**确实是**解
    （即上式对一切 `t` 成立），则必有 `ω² = c²·lam` -/
theorem z105_dispersion_relation {n : ℕ} (L : Matrix (Fin n) (Fin n) ℝ) (v : Fin n → ℝ)
    (lam c ω : ℝ) (hv : L *ᵥ v = lam • v) (hv0 : v ≠ 0)
    (hsol : deriv (deriv (fun s : ℝ => Real.cos (ω * s) • v))
      = fun t : ℝ => (-(c ^ 2)) • (L *ᵥ (Real.cos (ω * t) • v))) :
    ω ^ 2 = c ^ 2 * lam := by
  -- 在 `t = 0` 处取值（`cos 0 = 1`）
  have h0 := congrFun hsol 0
  rw [z105_deriv_wave, z105_deriv2_wave] at h0
  simp only [mul_zero, Real.cos_zero, mul_one, one_smul, hv, smul_smul] at h0
  -- 现在 `h0 : (−ω²) • v = ((−c²) * lam) • v`；由 `v ≠ 0` 取一个非零坐标
  obtain ⟨i, hi⟩ := Function.ne_iff.mp hv0
  simp only [Pi.zero_apply] at hi
  have hscalar : -(ω ^ 2) = -(c ^ 2 * lam) := by
    have h := congrFun h0 i
    simp only [Pi.smul_apply] at h
    -- 标量在非零坐标上可右消去
    simpa only [neg_mul] using mul_right_cancel₀ hi h
  -- 去掉两边的负号，得色散关系 `ω² = c²·lam`
  exact neg_injective hscalar

#print axioms z105_deriv_wave
#print axioms z105_deriv2_wave
#print axioms z105_wave_solution
#print axioms z105_dispersion_relation

end Zero
