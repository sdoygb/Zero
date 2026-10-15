/-
  Z109（"能量守恒耦合：从作用量设计，长期稳定"）的机器检查版本：
  **对称耦合势 ⟹ 能量守恒**（图 Laplacian 耦合的波动方程，二次能量的时间导数）。

  对应 Z109 的哪一句：Z109 §2 说耦合从势能泛函变分派生、"L 自伴 → 变分派生 →
  能量守恒"（守恒的对偶是"耦合项对称、无正反馈"）。本文件把这个机制放在最小模型上
  机器核验：位移场 `u : ℝ → ℝⁿ` 演化满足 `ü = −c² L u`（`L` 为对称矩阵，图上
  Laplacian 即其特例），二次能量
  `E(t) = ½⟨u̇, u̇⟩ + ½⟨u, L u⟩`（内积 `x ⬝ᵥ y = ∑ i, x i * y i`）。

  三条目标定理逐条交代：

  * `z109_deriv_kinetic`（动能项的导数）：`d/dt ½⟨u̇,u̇⟩ = ⟨ü, u̇⟩`——**已证**。
    逐坐标用乘积法则（`HasDerivAt.mul`）再求和（`HasDerivAt.sum`）。
    （审计整改 2026-10-15：原陈述带的 `h1`（`u` 的一阶导数）在证明中确实未被用到，
    已按审计**删除该假设** —— 删前提＝加强命题；原先为它保留的声明级
    `set_option linter.unusedVariables false` 已一并删除。）
  * `z109_deriv_potential`（势能项的导数）：`d/dt ½⟨u,Lu⟩ = ⟨u̇, L u⟩`——**已证**。
    关键两步：`L *ᵥ u s` 的导数是 `L *ᵥ u̇ s`（逐坐标 `∑ j, L i j * u s j` + 求和法则），
    以及用对称性 `Lᵀ = L` 把交叉项 `⟨u, L u̇⟩` 换成 `⟨u̇, L u⟩`
    （`Matrix.dotProduct_mulVec` + `Matrix.mulVec_transpose` + `dotProduct_comm`），
    于是 `d/dt ½⟨u,Lu⟩ = ½(⟨u̇,Lu⟩ + ⟨u̇,Lu⟩) = ⟨u̇,Lu⟩`。
  * `z109_energy_conserved`（能量守恒本身）：**这一条的陈述是假的，未证，且本文件给出
    机器核验的反例**——见文件末尾"反例"一节。原因一句话：
    `dE/dt = ⟨ü,u̇⟩ + ⟨u̇,Lu⟩`（两条导数引理相加），代入 `ü = −c²Lu` 得
    `dE/dt = (1 − c²)·⟨Lu, u̇⟩`，**只有 `c² = 1` 时才恒为零**；
    陈述里对 `c` 无任何约束，故 `c = 2` 即为反例。
    守恒的二次型应当是 `½⟨u̇,u̇⟩ + (c²/2)⟨u,Lu⟩`（把 `c²` 吸进势能，或等价地令
    `L` 已经含 `c²`）——这条**修正版**在下面证成 `z109_energy_conserved_corrected`；
    本文件还把原陈述的**一般导数公式**证成 `z109_energy_deriv`，
    再给出具体反例 `z109_energy_conserved_counterexample`（`n = 1`、`L = 1`、`c = 2`、
    `u(t) = cos 2t`：满足全部假设，而 `dE/dt` 在 `t = π/8` 处为正，解析值 3 ≠ 0）。

  Z109 正文对应的数值结论（能量漂移 ≈ 0、CFL 稳定性）不在机器检查范围内；本文件只覆盖
  "L 自伴 ⟹ 二次能量守恒"这一步的解析机制，以及它在陈述层面的确切适用范围。

  用到的 Mathlib 引理（注释中逐条标注）：
  * `hasDerivAt_pi`：向量值函数可导 ⟺ 逐坐标可导——把 `Fin n → ℝ` 值函数拆成坐标。
  * `HasDerivAt.mul` / `HasDerivAt.sum` / `HasDerivAt.const_mul` / `HasDerivAt.add`：
    逐坐标乘积法则、有限和与常数倍。
  * `hasDerivAt_const`：常函数的导数为 0（用于 `L` 的行系数）。
  * `dotProduct`（记号 `⬝ᵥ`，定义即 `∑ i, x i * y i`，rfl 可展开）、
    `Matrix.mulVec`（定义即 `fun i => ∑ j, M i j * v j`，rfl 可展开）。
  * `Matrix.dotProduct_mulVec`：`v ⬝ᵥ (A *ᵥ w) = (v ᵥ* A) ⬝ᵥ w`；
    `Matrix.mulVec_transpose`：`Aᵀ *ᵥ x = x ᵥ* A`——合起来把 `L` 从内积右边换到左边。
  * `dotProduct_comm`：`x ⬝ᵥ y = y ⬝ᵥ x`（实内积对称）。
  * `smul_dotProduct`：`(a • v) ⬝ᵥ w = a * (v ⬝ᵥ w)`（代入 `ü` 后把标量提出来）。
  * `HasDerivAt.deriv`：`HasDerivAt f f' x → deriv f x = f'`。
  * 反例一节另用：`Real.hasDerivAt_cos` / `Real.hasDerivAt_sin`、`hasDerivAt_id`、
    `HasDerivAt.comp`（链式法则）、`Matrix.one_mulVec`、`Real.cos_pos_of_mem_Ioo`、
    `Real.sin_pos_of_pos_of_lt_pi`（在 `t = π/8` 处判定导数非零）。
-/

import Mathlib

open scoped Matrix

namespace Zero

/-- 速度的平方模：`⟨u̇, u̇⟩` -/
noncomputable def z109SpeedSq {n : ℕ} (u : ℝ → Fin n → ℝ) (t : ℝ) : ℝ :=
  deriv u t ⬝ᵥ deriv u t

/-- 势能：`½⟨u, L u⟩` -/
noncomputable def z109Pot {n : ℕ} (L : Matrix (Fin n) (Fin n) ℝ) (u : ℝ → Fin n → ℝ) (t : ℝ) : ℝ :=
  (1 / 2) * (u t ⬝ᵥ (L *ᵥ u t))

/-! ### 技术引理（本文件自用，一律 `private`） -/

/-- 坐标面：`Fin n → ℝ` 值函数可导 ⟹ 每个坐标函数可导。 -/
private lemma z109HasDerivAt_apply {n : ℕ} {f : ℝ → Fin n → ℝ} {f' : Fin n → ℝ} {t : ℝ}
    (h : HasDerivAt f f' t) (i : Fin n) : HasDerivAt (fun s => f s i) (f' i) t :=
  (hasDerivAt_pi).mp h i

/-- 实内积的乘积法则：`d/dt ⟨f, g⟩ = ⟨f', g⟩ + ⟨f, g'⟩`（逐坐标 + 求和法则）。 -/
private lemma z109HasDerivAt_dotProduct {n : ℕ} (f g : ℝ → Fin n → ℝ) (f' g' : Fin n → ℝ)
    (t : ℝ) (hf : ∀ i, HasDerivAt (fun s => f s i) (f' i) t)
    (hg : ∀ i, HasDerivAt (fun s => g s i) (g' i) t) :
    HasDerivAt (fun s => f s ⬝ᵥ g s) (f' ⬝ᵥ g t + f t ⬝ᵥ g') t := by
  -- 逐坐标：`d/dt (f s i * g s i) = f' i * g t i + f t i * g' i`
  have hterm : ∀ i : Fin n, HasDerivAt (fun s => f s i * g s i) (f' i * g t i + f t i * g' i) t :=
    fun i => (hf i).mul (hg i)
  -- 有限和：`∑` 的导数是导数的 `∑`（`HasDerivAt.sum` 给的是"函数的和"）
  have hsum := HasDerivAt.sum (u := (Finset.univ : Finset (Fin n)))
    (A := fun i s => f s i * g s i) (A' := fun i => f' i * g t i + f t i * g' i)
    (fun i _ => hterm i)
  -- 左端："函数的和"写成"lambda 的和"，再展开成 `⬝ᵥ`；右端：和的分配律 + `⬝ᵥ` 展开
  convert hsum using 1 <;>
    first
      | (funext s; simp only [Finset.sum_apply, dotProduct])
      | simp only [dotProduct, Finset.sum_add_distrib]

/-- 矩阵—向量乘法的导数：`d/dt (L *ᵥ u s) = L *ᵥ u̇ s`（逐行 + 求和法则）。 -/
private lemma z109HasDerivAt_mulVec {n : ℕ} (L : Matrix (Fin n) (Fin n) ℝ) (u : ℝ → Fin n → ℝ)
    (u' : Fin n → ℝ) (t : ℝ) (hu : ∀ j, HasDerivAt (fun s => u s j) (u' j) t) :
    ∀ i, HasDerivAt (fun s => (L *ᵥ u s) i) ((L *ᵥ u') i) t := by
  intro i
  -- 第 `i` 行：`(L *ᵥ u s) i = ∑ j, L i j * u s j`，逐项乘积法则（`L i j` 是常数）
  have hterm : ∀ j : Fin n, HasDerivAt (fun s => L i j * u s j) (L i j * u' j) t := fun j =>
    (by
      -- 常函数 `L i j` 乘 `u · j`：先把函数的乘积改写成 lambda 形式，再化简 `0 * _ + _`
      have h := (hasDerivAt_const (x := t) (c := L i j)).mul (hu j)
      have hfun : ((fun _ : ℝ => L i j) * fun s => u s j) = fun s => L i j * u s j := by
        funext s; simp only [Pi.mul_apply]
      rw [hfun] at h
      simpa only [zero_mul, zero_add] using h)
  have hsum := HasDerivAt.sum (u := (Finset.univ : Finset (Fin n)))
    (A := fun j s => L i j * u s j) (A' := fun j => L i j * u' j) (fun j _ => hterm j)
  convert hsum using 1 <;>
    first
      | (funext s; simp only [Finset.sum_apply, Matrix.mulVec, dotProduct])
      | rfl

/-! ### Z109 的两条导数引理 -/

/-- **Z109 · 动能项的导数**：`d/dt ½⟨u̇,u̇⟩ = ⟨ü, u̇⟩`
    （只需二阶导可导；原陈述里"`u` 一阶可导"的假设在证明中未被用到，已删） -/
theorem z109_deriv_kinetic {n : ℕ} (u : ℝ → Fin n → ℝ) (t : ℝ)
    (h2 : HasDerivAt (deriv u) (deriv (deriv u) t) t) :
    HasDerivAt (fun s : ℝ => (1 / 2) * (deriv u s ⬝ᵥ deriv u s))
      ((deriv (deriv u) t) ⬝ᵥ (deriv u t)) t := by
  -- `d/dt ⟨u̇,u̇⟩ = ⟨ü,u̇⟩ + ⟨u̇,ü⟩ = 2⟨ü,u̇⟩`（实内积对称）；乘 `½` 后只剩一项
  have hdot := z109HasDerivAt_dotProduct (deriv u) (deriv u) (deriv (deriv u) t)
    (deriv (deriv u) t) t (fun i => z109HasDerivAt_apply h2 i) (fun i => z109HasDerivAt_apply h2 i)
  have hhalf := hdot.const_mul (1 / 2 : ℝ)
  have hcomm : (deriv u t) ⬝ᵥ (deriv (deriv u) t) = (deriv (deriv u) t) ⬝ᵥ (deriv u t) :=
    dotProduct_comm _ _
  convert hhalf using 1
  rw [hcomm]; ring

/-- **Z109 · 势能项的导数（对称 L）**：`d/dt ½⟨u,Lu⟩ = ⟨u̇, L u⟩` -/
theorem z109_deriv_potential {n : ℕ} (L : Matrix (Fin n) (Fin n) ℝ) (hsym : Lᵀ = L)
    (u : ℝ → Fin n → ℝ) (t : ℝ) (h1 : HasDerivAt u (deriv u t) t) :
    HasDerivAt (fun s : ℝ => (1 / 2) * (u s ⬝ᵥ (L *ᵥ u s)))
      ((deriv u t) ⬝ᵥ (L *ᵥ u t)) t := by
  have hu : ∀ i, HasDerivAt (fun s => u s i) (deriv u t i) t :=
    fun i => z109HasDerivAt_apply h1 i
  -- `d/dt ⟨u, Lu⟩ = ⟨u̇, Lu⟩ + ⟨u, L u̇⟩`；对称性把第二项换成第一项
  have hdot := z109HasDerivAt_dotProduct u (fun s => L *ᵥ u s) (deriv u t) (L *ᵥ deriv u t) t hu
    (z109HasDerivAt_mulVec L u (deriv u t) t hu)
  have hsymm : u t ⬝ᵥ (L *ᵥ deriv u t) = (deriv u t) ⬝ᵥ (L *ᵥ u t) := by
    rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose L (u t), hsym, dotProduct_comm]
  have hhalf := hdot.const_mul (1 / 2 : ℝ)
  convert hhalf using 1
  rw [hsymm]; ring

/-! ### `z109_energy_conserved` 的机器核验：陈述不成立

把上面两条导数引理相加即得能量的一般导数公式
`dE/dt = ⟨ü,u̇⟩ + ⟨u̇,Lu⟩ = (1 − c²)·⟨Lu,u̇⟩`（第二式代入 `ü = −c²Lu` 并用实内积对称）。
`(1 − c²)` 一般不为零，故"能量守恒"的结论不成立；下面先机器核验这条一般公式，
再给具体反例。 -/

/-- **Z109 · 反例的核心（一般导数公式）**：在 `z109_energy_conserved` 的全部假设下，
    `dE/dt = (1 − c²)·⟨L u, u̇⟩`——一般不为零，故该陈述为假；
    只有 `c² = 1`（或 `⟨Lu,u̇⟩ = 0`）时它才退化为守恒。 -/
theorem z109_energy_deriv {n : ℕ} (L : Matrix (Fin n) (Fin n) ℝ) (hsym : Lᵀ = L) (c : ℝ)
    (u : ℝ → Fin n → ℝ) (hdiff : ∀ t, HasDerivAt u (deriv u t) t)
    (hdiff2 : ∀ t, HasDerivAt (deriv u) (deriv (deriv u) t) t)
    (hu : ∀ t, deriv (deriv u) t = (-(c ^ 2)) • (L *ᵥ u t)) (t : ℝ) :
    deriv (fun s : ℝ => (1 / 2) * (deriv u s ⬝ᵥ deriv u s) + z109Pot L u s) t
      = (1 - c ^ 2) * (L *ᵥ u t ⬝ᵥ deriv u t) := by
  -- 两项各自的导数相加
  have hsum : HasDerivAt (fun s : ℝ => (1 / 2) * (deriv u s ⬝ᵥ deriv u s) + z109Pot L u s)
      ((deriv (deriv u) t) ⬝ᵥ (deriv u t) + (deriv u t) ⬝ᵥ (L *ᵥ u t)) t :=
    (z109_deriv_kinetic u t (hdiff2 t)).add (z109_deriv_potential L hsym u t (hdiff t))
  rw [hsum.deriv]
  -- 代入波动方程 `ü = −c²Lu`，把标量提出来，并用实内积对称对齐两项
  have hcomm : (deriv u t) ⬝ᵥ (L *ᵥ u t) = (L *ᵥ u t) ⬝ᵥ (deriv u t) := dotProduct_comm _ _
  rw [hu t, smul_dotProduct, hcomm]
  ring

/-! ### 修正后的能量守恒（势能项带 `c²` 权重）——机器核验

把 `c²` 吸进势能：`E_c(t) = ½⟨u̇,u̇⟩ + (c²/2)⟨u,Lu⟩`，则
`dE_c/dt = ⟨ü,u̇⟩ + c²⟨u̇,Lu⟩ = −c²⟨Lu,u̇⟩ + c²⟨u̇,Lu⟩ = 0`。
这就是波动方程 `ü = −c²Lu` 下真正守恒的二次型（等价说法：`L` 已含 `c²`，或 `c = ±1`）。
Z109 的"从作用量设计 ⟹ 能量守恒"在解析层面的确切内容是这一条。 -/

/-- 带 `c²` 权重的势能：`(c²/2)⟨u, L u⟩ = c²·z109Pot L u` -/
noncomputable def z109PotC {n : ℕ} (L : Matrix (Fin n) (Fin n) ℝ) (c : ℝ) (u : ℝ → Fin n → ℝ)
    (t : ℝ) : ℝ :=
  (c ^ 2) * z109Pot L u t

/-- **Z109 · 修正后的能量守恒（机器核验）**：`Lᵀ = L` 且 `ü = −c²Lu` ⟹ 修正能量
    `½⟨u̇,u̇⟩ + (c²/2)⟨u,Lu⟩` 的导数恒为 0（对**任意** `c`）。
    与 `z109_energy_conserved` 的差别只在势能项的 `c²` 权重——原陈述缺这一权重，故为假。 -/
theorem z109_energy_conserved_corrected {n : ℕ} (L : Matrix (Fin n) (Fin n) ℝ) (hsym : Lᵀ = L)
    (c : ℝ) (u : ℝ → Fin n → ℝ) (hdiff : ∀ t, HasDerivAt u (deriv u t) t)
    (hdiff2 : ∀ t, HasDerivAt (deriv u) (deriv (deriv u) t) t)
    (hu : ∀ t, deriv (deriv u) t = (-(c ^ 2)) • (L *ᵥ u t)) :
    ∀ t, deriv (fun s : ℝ => (1 / 2) * (deriv u s ⬝ᵥ deriv u s) + z109PotC L c u s) t = 0 := by
  intro t
  -- 动能项导数 + `c²` 倍的势能项导数
  have hsum : HasDerivAt
      (fun s : ℝ => (1 / 2) * (deriv u s ⬝ᵥ deriv u s) + z109PotC L c u s)
      ((deriv (deriv u) t) ⬝ᵥ (deriv u t) + (c ^ 2) * ((deriv u t) ⬝ᵥ (L *ᵥ u t))) t :=
    (z109_deriv_kinetic u t (hdiff2 t)).add
      ((z109_deriv_potential L hsym u t (hdiff t)).const_mul (c ^ 2))
  rw [hsum.deriv, hu t, smul_dotProduct]
  have hcomm : (deriv u t) ⬝ᵥ (L *ᵥ u t) = (L *ᵥ u t) ⬝ᵥ (deriv u t) := dotProduct_comm _ _
  rw [hcomm]
  ring

/-! ### 反例：`z109_energy_conserved` 的陈述不成立（机器核验）

具体数据取 `n = 1`、`L = 1`（对称）、`c = 2`、`u(t) = cos(2t)·e₁`：
* `u̇(t) = −2 sin(2t)·e₁`，`ü(t) = −4 cos(2t)·e₁ = −c²·(L u(t))`，假设全满足；
* 但 `E(t) = ½⟨u̇,u̇⟩ + ½⟨u,Lu⟩ = 2sin²(2t) + ½cos²(2t)`，由 `z109_energy_deriv`
  `E'(t) = (1 − c²)·⟨Lu,u̇⟩ = −3·(cos(2t)·(−2 sin(2t))) = 6 sin(2t)cos(2t)`，
  在 `t = π/8` 处 `= 6·(√2/2)·(√2/2) = 3 ≠ 0`。故"能量守恒"的陈述为假。 -/

/-- 反例场的一阶导：`d/dt cos(2t) = −2 sin(2t)`（链式法则 + `hasDerivAt_pi` 拆坐标） -/
private lemma z109Cex_hasDeriv (t : ℝ) :
    HasDerivAt (fun s : ℝ => fun _ : Fin 1 => Real.cos (2 * s))
      (fun _ : Fin 1 => (-2) * Real.sin (2 * t)) t := by
  rw [hasDerivAt_pi]
  intro i
  exact ((Real.hasDerivAt_cos (2 * t)).comp t ((hasDerivAt_id t).const_mul 2)).congr_deriv (by ring)

/-- 反例场的二阶导：`d/dt (−2 sin 2t) = −4 cos 2t` -/
private lemma z109Cex_hasDeriv2 (t : ℝ) :
    HasDerivAt (fun s : ℝ => fun _ : Fin 1 => (-2) * Real.sin (2 * s))
      (fun _ : Fin 1 => (-4) * Real.cos (2 * t)) t := by
  rw [hasDerivAt_pi]
  intro i
  have h : HasDerivAt (fun s : ℝ => Real.sin (2 * s)) (Real.cos (2 * t) * 2) t :=
    ((Real.hasDerivAt_sin (2 * t)).comp t ((hasDerivAt_id t).const_mul 2)).congr_deriv (by ring)
  exact (h.const_mul (-2)).congr_deriv (by ring)

/-- 反例场的一阶导函数 -/
private lemma z109Cex_deriv :
    deriv (fun s : ℝ => fun _ : Fin 1 => Real.cos (2 * s))
      = fun t : ℝ => fun _ : Fin 1 => (-2) * Real.sin (2 * t) :=
  funext fun t => (z109Cex_hasDeriv t).deriv

/-- 反例场的二阶导函数 -/
private lemma z109Cex_deriv2 :
    deriv (deriv (fun s : ℝ => fun _ : Fin 1 => Real.cos (2 * s)))
      = fun t : ℝ => fun _ : Fin 1 => (-4) * Real.cos (2 * t) := by
  rw [z109Cex_deriv]
  exact funext fun t => (z109Cex_hasDeriv2 t).deriv

/-- 反例数据满足"一阶可导"假设 -/
private lemma z109Cex_hdiff : ∀ t, HasDerivAt (fun s : ℝ => fun _ : Fin 1 => Real.cos (2 * s))
    (deriv (fun s : ℝ => fun _ : Fin 1 => Real.cos (2 * s)) t) t := fun t => by
  rw [z109Cex_deriv]
  exact z109Cex_hasDeriv t

/-- 反例数据满足"二阶可导"假设 -/
private lemma z109Cex_hdiff2 : ∀ t,
    HasDerivAt (deriv (fun s : ℝ => fun _ : Fin 1 => Real.cos (2 * s)))
      (deriv (deriv (fun s : ℝ => fun _ : Fin 1 => Real.cos (2 * s))) t) t := fun t => by
  rw [z109Cex_deriv2, z109Cex_deriv]
  exact z109Cex_hasDeriv2 t

/-- 反例数据满足波动方程 `ü = −c²·L u`（`c = 2`、`L = 1`：两边都是 `−4 cos 2t`） -/
private lemma z109Cex_hPDE : ∀ t,
    deriv (deriv (fun s : ℝ => fun _ : Fin 1 => Real.cos (2 * s))) t
      = (-((2 : ℝ) ^ 2)) • ((1 : Matrix (Fin 1) (Fin 1) ℝ) *ᵥ
          (fun s : ℝ => fun _ : Fin 1 => Real.cos (2 * s)) t) := by
  intro t
  rw [z109Cex_deriv2, Matrix.one_mulVec]
  funext i
  simp only [Pi.smul_apply, smul_eq_mul]
  ring

/-- **Z109 · 反例（机器核验）**：`z109_energy_conserved` 的陈述**为假**。
    取 `n = 1`、`L = 1`、`c = 2`、`u(t) = cos(2t)`：该组数据满足陈述里的**全部**假设
    （`Lᵀ = L`、一阶/二阶可导、`ü = −c²Lu`），但能量在 `t = π/8` 处的导数**为正**
    （文件里机器核验的是该量的正性，故 ≠ 0；解析值 = `(1−4)·(cos(π/4)·(−2 sin(π/4)))`
    = `6·(√2/2)²` = `3`）。
    故"对任意 `c` 能量守恒"不成立；守恒要求势能项带 `c²`（即 `½⟨u̇,u̇⟩ + (c²/2)⟨u,Lu⟩`），
    或等价地 `c = ±1` / `L` 已含 `c²`。 -/
theorem z109_energy_conserved_counterexample :
    ∃ (L : Matrix (Fin 1) (Fin 1) ℝ) (c : ℝ) (u : ℝ → Fin 1 → ℝ),
      Lᵀ = L ∧ (∀ t, HasDerivAt u (deriv u t) t) ∧
      (∀ t, HasDerivAt (deriv u) (deriv (deriv u) t) t) ∧
      (∀ t, deriv (deriv u) t = (-(c ^ 2)) • (L *ᵥ u t)) ∧
      ∃ t, deriv (fun s : ℝ => (1 / 2) * (deriv u s ⬝ᵥ deriv u s) + z109Pot L u s) t ≠ 0 := by
  refine ⟨1, 2, (fun s : ℝ => fun _ : Fin 1 => Real.cos (2 * s)), Matrix.transpose_one,
    z109Cex_hdiff, z109Cex_hdiff2, z109Cex_hPDE, ?_⟩
  -- 取 `t = π/8`：用一般公式（`z109_energy_deriv`）算出导数，再判定它非零
  refine ⟨Real.pi / 8, ?_⟩
  have hgen := z109_energy_deriv (1 : Matrix (Fin 1) (Fin 1) ℝ) Matrix.transpose_one 2
    (fun s : ℝ => fun _ : Fin 1 => Real.cos (2 * s)) z109Cex_hdiff z109Cex_hdiff2
    z109Cex_hPDE (Real.pi / 8)
  rw [hgen, z109Cex_deriv, Matrix.one_mulVec]
  simp only [dotProduct, Fin.sum_univ_one]
  rw [show (2 : ℝ) * (Real.pi / 8) = Real.pi / 4 by ring]
  -- 目标：`(1 − 2²)·(cos(π/4)·(−2 sin(π/4))) ≠ 0`——负负得正
  have hc : 0 < Real.cos (Real.pi / 4) :=
    Real.cos_pos_of_mem_Ioo ⟨by linarith [Real.pi_pos], by linarith [Real.pi_pos]⟩
  have hs : 0 < Real.sin (Real.pi / 4) :=
    Real.sin_pos_of_pos_of_lt_pi (by linarith [Real.pi_pos]) (by linarith [Real.pi_pos])
  have hneg : Real.cos (Real.pi / 4) * ((-2) * Real.sin (Real.pi / 4)) < 0 :=
    mul_neg_of_pos_of_neg hc (by linarith)
  have hpos : 0 < (1 - 2 ^ 2) * (Real.cos (Real.pi / 4) * ((-2) * Real.sin (Real.pi / 4))) :=
    mul_pos_of_neg_of_neg (by norm_num) hneg
  exact ne_of_gt hpos

#print axioms z109_deriv_kinetic
#print axioms z109_deriv_potential
#print axioms z109_energy_deriv
#print axioms z109_energy_conserved_corrected
#print axioms z109_energy_conserved_counterexample

end Zero