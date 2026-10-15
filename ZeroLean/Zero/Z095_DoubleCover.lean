/-
  Z095（无质量规范玻色子承载 · 胶子 = K₄ 双覆盖 8 顶点）与
  Z106（张量层 · 每顶点对称张量度规扰动 / 无迹偏振色散）的机器检查版本。

  要验的命题：K₄□K₂（K₄ 的双重覆盖，8 顶点 = `Fin 4 × Fin 2`）的邻接谱为
  `{4, 2, 0³, (−2)³}`，且图 4-正则（Z095 第三节把它的 8 个顶点用作 SU(3) 伴随表示、
  即"胶子色八重"的承载；Z106 张量层用的是同一条笛卡尔积分解）。

  谱为何是这四个数（Z106 张量层的结构来源）：
    K₄ 的邻接谱 `{3, (−1)³}`，K₂ 的谱 `{1, −1}`；**笛卡尔积的谱 = 和集**：
    `{3+1, 3−1, (−1)+1 (重 3), (−1)−1 (重 3)} = {4, 2, 0³, (−2)³}` ✓
    即"层内"自由度与"层间"自由度的本征值**相加**——这正是 Z106 把每顶点扰动写成
    层内 ⊕ 层间之和的那条分解。数值交叉核对：`ZeroGPU/src/zero_sim.py` 的结构层
    `double_cover` + `spectrum` 自检判据里有 "K₄□K₂ 谱 {4,2,0³,(−2)³}"（现跑全过）。

  本文件钉住的（全部机器核验，无 `sorry`）：
    · 行和 = 4（4-正则）、迹 = 0、`∑∑ = 32`（= 2|E|，|E| = 16）；
    · **张量本征关系**（承重）：`A₄ v = λ₄ v ∧ A₂ w = λ₂ w ⟹ M (v⊗w) = (λ₄+λ₂)(v⊗w)`
      —— 这就是"笛卡尔积谱 = 和集"本身，且证明不依赖 K₄/K₂ 的具体值
      （只用两块的**零对角**邻接结构，见 `z095_mulVec_apply` 的分量式）；
    · **零化多项式** `(M−4·1)(M−2·1)M(M+2·1) = 0` ⟹ 谱 ⊆ {4,2,0,−2}（谱的**上界**）；
    · 4 / 2 / 0 / −2 各有**显式本征向量见证**（谱的**下界**）。

  重数的旁注（非本文件定理，只作口径说明）：把谱写成重数 `m₄+m₂+m₀+m₋₂ = 8`，
  迹给 `4m₄+2m₂−2m₋₂ = 0`，而对称 0/1 矩阵有 `trace(M²) = ∑_{p,q} M p q`
  （= 本文件第 3 条的 32）给 `16m₄+4m₂+4m₋₂ = 32`；加上"四个值都出现"（见证），
  唯一解 `m₄ = m₂ = 1`、`m₀ = m₋₂ = 3` ✓ —— 故文件名里的 `{4, 2, 0³, (−2)³}` 与
  见证条一致（`trace(M²) = ∑∑M` 用到 M 对称，由定义即见）。

  ⚠ 对任务书第 6 条（`z095_four_eigenvalues`）的修订 —— **定理名与"四个本征值各有一个
  显式见证"的内容不变，只换了两个见证向量**：
    ① 任务书原第 3 个见证 `fun p => if p.1 = 0 then -1 else 1` **不是**本征向量。
       设 `x p = if p.1 = 0 then −1 else 1`（只依赖第一坐标 f = (−1,1,1,1)），则
       `(M x)(v,a) = ∑_{w≠v} f w + f v = (∑_w f w) − f v + f v = ∑_w f w = 2`（与 (v,a) 无关），
       即 `M x` 恒为 2 的常函数，而 `x` 取值 ±1 —— 没有任何 λ 使 `M x = λ x`。
    ② 任务书原第 4 个见证 `fun p => if p.2 = 0 then 1 else -1` 与**第 2 个见证逐字相同**，
       它见证的是本征值 **2**，不是 −2。
    按任务书指示换成「K₄ 的 −1 特征向量 ⊗ K₂ 的对应向量」，取 Z090 用过的 −1 特征向量
    `f = (3,−1,−1,−1)`（Z090：−1 特征空间 = 分量和为零的超平面）：
      第 3 个：`f ⊗ 1`     ⟹ 本征值 (−1) + (+1) = 0   ✓（0³ 的第一个）
      第 4 个：`f ⊗ (1,−1)` ⟹ 本征值 (−1) + (−1) = −2 ✓
    这两个替换后的见证在本文件中机器核验通过。
    （任务书陈述里未标注类型的字面量 `1`/`-1` 在本文件写成 `(1 : ℝ)`：否则 `2 • (fun p => …)`
      的字面量会被推成 ℤ，合成 `HSMul ℕ (→ ℤ) (→ ℝ)` 直接失败。）

  机器核验（2026-10-15）：`lake env lean Zero/Z095_DoubleCover.lean` 零 error 零 warning
  （约 1m45s，其中 `z095_zeroing` 的 64 个分量式穷举占大头）；五条 `#print axioms` 全为
  `[propext, Classical.choice, Quot.sound]`，无 `sorryAx`、无 `native_decide`。

  ⚠ 审计整改（2026-10-15）：`z095_trace_zero` 内容极薄（就是"邻接矩阵无自环"，
  `simp` 一步），已按审计改为 `private` —— 不再进 `Zero.*` 审计空间，
  其 `#print axioms` 行一并删除；陈述与证明逐字未动。

  环境坑（本文件实测，供后人复用）：
    ① `simp` 会把 0/1 指示函数求和改写成 `Finset.card`（`Finset.sum_boole`），使 `norm_num`
       无从下手 —— 本文件一律写 `-Finset.sum_boole` 关掉它；
    ② 大矩阵恒等式的分量式穷举会造出**矩阵数词** `((4 : Matrix …) i j)`，它不是 `Pi` 数词
       （`Matrix.ofNat_apply` **没有** `@[simp]`），默认 simp 化简不掉，残目标停在纯数词却
       关不掉 —— 必须把 `Matrix.ofNat_apply` 显式写进 simp 清单；
    ③ 提升上限要用**声明级** `set_option maxHeartbeats N`（独立命令）：写成证明内的
       `set_option maxHeartbeats N in <tactic>` 不抬高声明累计上限，仍会在 200000 处超时。

  要素来源表见文末 `#print axioms`（预期全为 propext / Classical.choice / Quot.sound）。
-/

import Mathlib

open scoped Matrix

namespace Zero

/-- K₄□K₂ 的邻接矩阵（顶点 = Fin 4 × Fin 2） -/
def z095Adj : Matrix (Fin 4 × Fin 2) (Fin 4 × Fin 2) ℝ :=
  Matrix.of (fun p q => if (p.1 = q.1 ∧ p.2 ≠ q.2) ∨ (p.2 = q.2 ∧ p.1 ≠ q.1) then 1 else 0)

/-! ### 技术引理（private，不进 `Zero.*` 审计空间） -/

/-- **笛卡尔积的邻接分解**：`K₄□K₂` 的邻接 = "同层（a 相同）的 K₄ 部分" + "同顶点（v 相同）的 K₂ 部分"；
    两块的对角都空（`if i₀ = i then 0 else 1`），故不会重复计数。 -/
private lemma z095_decomp (i₀ : Fin 4) (a₀ : Fin 2) (i : Fin 4) (a : Fin 2) :
    z095Adj (i₀, a₀) (i, a) =
      (if a₀ = a then (if i₀ = i then 0 else 1) else 0) +
        (if i₀ = i then (if a₀ = a then 0 else 1) else 0) := by
  by_cases h1 : i₀ = i <;> by_cases h2 : a₀ = a <;> simp [z095Adj, h1, h2]

/-- `∑ a : Fin 2, (if c = a then t else 0) = t`（只有一项存活）。 -/
private lemma z095_sum_fin2_const (c : Fin 2) (t : ℝ) :
    (∑ a : Fin 2, (if c = a then t else 0)) = t := by
  fin_cases c <;> simp

/-- `∑ i : Fin 4, (if c = i then t else 0) = t`（只有一项存活）。 -/
private lemma z095_sum_fin4_const (c : Fin 4) (t : ℝ) :
    (∑ i : Fin 4, (if c = i then t else 0)) = t := by
  fin_cases c <;> simp

/-- `∑ i : Fin 4, (if c = i then f i else 0) = f c`（逐项版，`i` 处留着求和变量）。 -/
private lemma z095_sum_fin4_ite (c : Fin 4) (f : Fin 4 → ℝ) :
    (∑ i : Fin 4, (if c = i then f i else 0)) = f c := by
  fin_cases c <;> simp

/-- 命题 `c` 下的整体常数分解：`∑ x, (if c then f x else 0) = if c then ∑ x, f x else 0`。 -/
private lemma z095_sum_ite_const {α : Type*} [Fintype α] (c : Prop) [Decidable c]
    (f : α → ℝ) :
    (∑ x : α, (if c then f x else 0)) = if c then (∑ x : α, f x) else 0 := by
  by_cases h : c <;> simp [h]

/-- **矩阵—向量乘法的分量式**（承重引理）：在顶点 `(i₀, a₀)` 处，像 = 同层 3 个邻居
    （`i ≠ i₀`，同 `a₀`）+ 同顶点 1 个邻居（`a ≠ a₀`）。
    这就是 Z106 张量层"层内 ⊕ 层间"两个自由度的分解；也说明谱的和集性质与具体图无关。 -/
private lemma z095_mulVec_apply (x : Fin 4 × Fin 2 → ℝ) (i₀ : Fin 4) (a₀ : Fin 2) :
    (z095Adj *ᵥ x) (i₀, a₀) =
      (∑ i : Fin 4, (if i₀ = i then (0 : ℝ) else 1) * x (i, a₀)) +
        (∑ a : Fin 2, (if a₀ = a then (0 : ℝ) else 1) * x (i₀, a)) := by
  rw [Matrix.mulVec, dotProduct, Fintype.sum_prod_type]
  have hterm : ∀ (i : Fin 4) (a : Fin 2), z095Adj (i₀, a₀) (i, a) * x (i, a) =
      (if a₀ = a then (if i₀ = i then (0 : ℝ) else 1) * x (i, a₀) else 0) +
        (if i₀ = i then (if a₀ = a then (0 : ℝ) else 1) * x (i, a) else 0) := by
    intro i a
    rw [z095_decomp, add_mul]
    by_cases h1 : i₀ = i <;> by_cases h2 : a₀ = a <;> simp [h1, h2]
  simp only [hterm, Finset.sum_add_distrib]
  congr 1
  · refine Finset.sum_congr rfl (fun i _ => ?_)
    exact z095_sum_fin2_const a₀ ((if i₀ = i then (0 : ℝ) else 1) * x (i, a₀))
  · rw [Finset.sum_congr rfl (fun i _ =>
      z095_sum_ite_const (i₀ = i) (fun a => (if a₀ = a then (0 : ℝ) else 1) * x (i, a)))]
    exact z095_sum_fin4_ite i₀ (fun i => ∑ a : Fin 2, (if a₀ = a then (0 : ℝ) else 1) * x (i, a))

/-! ### 任务书要求的六条 -/

/-- **Z095 · 4-正则**：每行和为 4 -/
theorem z095_row_sum (p : Fin 4 × Fin 2) : ∑ q, z095Adj p q = 4 := by
  obtain ⟨i, a⟩ := p
  fin_cases i <;> fin_cases a <;>
    simp [z095Adj, Fintype.sum_prod_type, Fin.sum_univ_four, Fin.sum_univ_two, -Finset.sum_boole] <;>
    norm_num

/-- **Z095 · 迹为 0**（技术条目，已按审计私有化）：无自环（`(v,a)` 与自身既不"同层异 a"也不"同 a 异 v"）。 -/
private theorem z095_trace_zero : Matrix.trace z095Adj = 0 := by
  simp [Matrix.trace, Matrix.diag, z095Adj]

/-- **Z095 · 边数**：`∑ p q, z095Adj p q = 2 * 16`（= 2|E|，|E| = 16 = 8·4/2） -/
theorem z095_two_mul_edges : ∑ p : Fin 4 × Fin 2, ∑ q : Fin 4 × Fin 2, z095Adj p q = 32 := by
  simp [z095Adj, Fintype.sum_prod_type, Fin.sum_univ_four, Fin.sum_univ_two, -Finset.sum_boole]
  norm_num

/-- **Z095 · 张量本征关系**：若 `A₄ v = lam₄ • v` 且 `A₂ w = lam₂ • w` 则积向量是本征向量，
    本征值相加（这就是"笛卡尔积谱 = 和集"）。

    证明走 `z095_mulVec_apply` 的分量式：像 = `（层内 K₄ 部分） + （层间 K₂ 部分）`，
    两个和分别提出标量 `w a₀` / `v i₀` 后用 `hv`/`hw` 换成本征值，最后 `ring`。 -/
theorem z095_eigen_product (v : Fin 4 → ℝ) (w : Fin 2 → ℝ) (lam4 lam2 : ℝ)
    (hv : (Matrix.of (fun i j : Fin 4 => if i = j then 0 else 1)) *ᵥ v = lam4 • v)
    (hw : (Matrix.of (fun i j : Fin 2 => if i = j then 0 else 1)) *ᵥ w = lam2 • w) :
    z095Adj *ᵥ (fun p : Fin 4 × Fin 2 => v p.1 * w p.2)
      = (lam4 + lam2) • (fun p : Fin 4 × Fin 2 => v p.1 * w p.2) := by
  have hv' : ∀ i : Fin 4,
      (∑ j : Fin 4, (if i = j then (0 : ℝ) else 1) * v j) = lam4 * v i := by
    intro i
    have h := congrFun hv i
    simpa [Matrix.mulVec, dotProduct] using h
  have hw' : ∀ a : Fin 2,
      (∑ b : Fin 2, (if a = b then (0 : ℝ) else 1) * w b) = lam2 * w a := by
    intro a
    have h := congrFun hw a
    simpa [Matrix.mulVec, dotProduct] using h
  funext p
  obtain ⟨i₀, a₀⟩ := p
  rw [z095_mulVec_apply]
  simp only [Pi.smul_apply, smul_eq_mul]
  have h1 : (∑ i : Fin 4, (if i₀ = i then (0 : ℝ) else 1) * (v i * w a₀)) =
      lam4 * v i₀ * w a₀ := by
    have hsplit : (∑ i : Fin 4, (if i₀ = i then (0 : ℝ) else 1) * (v i * w a₀)) =
        (∑ i : Fin 4, (if i₀ = i then (0 : ℝ) else 1) * v i) * w a₀ := by
      rw [Finset.sum_mul]
      exact Finset.sum_congr rfl (fun i _ => by ring)
    rw [hsplit, hv' i₀]
  have h2 : (∑ a : Fin 2, (if a₀ = a then (0 : ℝ) else 1) * (v i₀ * w a)) =
      v i₀ * (lam2 * w a₀) := by
    have hsplit : (∑ a : Fin 2, (if a₀ = a then (0 : ℝ) else 1) * (v i₀ * w a)) =
        v i₀ * (∑ a : Fin 2, (if a₀ = a then (0 : ℝ) else 1) * w a) := by
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl (fun a _ => by ring)
    rw [hsplit, hw' a₀]
  rw [h1, h2]
  ring

set_option maxHeartbeats 8000000

/-- **Z095 · 零化多项式**：`(M − 4·1)(M − 2·1)M(M + 2·1) = 0` ⟹ 谱 ⊆ {4, 2, 0, −2} -/
theorem z095_zeroing :
    (z095Adj - 4 • (1 : Matrix (Fin 4 × Fin 2) (Fin 4 × Fin 2) ℝ))
      * (z095Adj - 2 • (1 : Matrix (Fin 4 × Fin 2) (Fin 4 × Fin 2) ℝ))
      * z095Adj
      * (z095Adj + 2 • (1 : Matrix (Fin 4 × Fin 2) (Fin 4 × Fin 2) ℝ)) = 0 := by
  ext p q
  obtain ⟨i, a⟩ := p
  obtain ⟨j, b⟩ := q
  fin_cases i <;> fin_cases a <;> fin_cases j <;> fin_cases b <;>
    simp [z095Adj, Matrix.mul_apply, Matrix.sub_apply, Matrix.add_apply, Matrix.zero_apply,
      Matrix.ofNat_apply, Fintype.sum_prod_type, Fin.sum_univ_four, Fin.sum_univ_two,
      -Finset.sum_boole] <;>
    norm_num

/-- **Z095 · 四个本征值都有显式见证**（重数由零化多项式 + 迹 + 平方迹 32 钉住；
    第 3、4 个见证向量已按任务书指示换成 `f ⊗ 1` 与 `f ⊗ (1,−1)`，`f = (3,−1,−1,−1)`
    是 K₄ 的 −1 特征向量，见文件头注释的修订说明）。 -/
theorem z095_four_eigenvalues :
    z095Adj *ᵥ (fun _ : Fin 4 × Fin 2 => (1 : ℝ)) = 4 • (fun _ : Fin 4 × Fin 2 => (1 : ℝ)) ∧
    z095Adj *ᵥ (fun p : Fin 4 × Fin 2 => if p.2 = 0 then (1 : ℝ) else -1)
      = 2 • (fun p : Fin 4 × Fin 2 => if p.2 = 0 then (1 : ℝ) else -1) ∧
    z095Adj *ᵥ (fun p : Fin 4 × Fin 2 => if p.1 = 0 then (3 : ℝ) else -1) = 0 ∧
    z095Adj *ᵥ (fun p : Fin 4 × Fin 2 =>
        (if p.1 = 0 then (3 : ℝ) else -1) * (if p.2 = 0 then (1 : ℝ) else -1))
      = -2 • (fun p : Fin 4 × Fin 2 =>
        (if p.1 = 0 then (3 : ℝ) else -1) * (if p.2 = 0 then (1 : ℝ) else -1)) := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;>
    · funext p
      obtain ⟨i, a⟩ := p
      fin_cases i <;> fin_cases a <;>
        simp [z095Adj, Matrix.mulVec, dotProduct, Fintype.sum_prod_type, Fin.sum_univ_four,
          Fin.sum_univ_two, -Finset.sum_boole] <;>
        norm_num

/-! ### 要素来源表（`#print axioms`） -/

#print axioms z095_row_sum
#print axioms z095_two_mul_edges
#print axioms z095_eigen_product
#print axioms z095_zeroing
#print axioms z095_four_eigenvalues

end Zero