/-
  Z021 · 正则 D-单纯形的规范构造与等角性 —— 机器检查版

  对应文档：`ZeroGPU/doc/coupling_verdict.md`（耦合方式裁决：正则单纯形在"无偏好"下唯一）
  对应篇：Z021（维数字典；通道数 m = D + 1，维数 D = m − 1）

  文档 §2.2/§2.3 的核验表是**数值**核验（零和残差 ~10⁻¹⁶、夹角"精确"）；
  本文件把同一批断言升级成精确的机器检查定理，逐条对应关系写在每条注释里。

  · §2.1 构造：u_a = e_a − (1/m)·1
  · §1.3 约束二（零和）：∑_a u_a = 0，即通道向量落在超平面 H_Q = {x : ∑_a x_a = 0} 内
  · §1.3 约束三（通道无偏好）→ §2.2 等长：⟨u_a,u_a⟩ = 1 − 1/m，与 a 无关
  · §2.2 等角：⟨u_a,u_b⟩ = −1/m（a ≠ b）；归一化后夹角余弦 ρ = −1/(m−1) = −1/D
  · §2.3 逐维核验表：D = 4（m = 5）给 ρ = −1/4；D = 3（m = 4）给 ρ = −1/3（四面体角）
  · §1.1 表 / §2.3 表 D = 3 行：C(4,2) = 6 对关系、C(3,2) = 3 独立身份

  注：文档 §2.2 表里"等长 |u_a| = 1"一行指的是**归一化后**的向量；
  本文件用未归一化的 u_a = e_a − (1/m)·1，故长度为 1 − 1/m（正值，m ≥ 1），
  等角性以 ⟨u_a,u_b⟩/⟨u_a,u_a⟩ = −1/D 的形式给出，与归一化后的夹角一致。
-/

import Mathlib

namespace Zero

noncomputable section

/-- 中心化标准基：u_a = e_a − (1/m)·1（对应 coupling_verdict.md §2.1 的构造式）。 -/
def z021SimplexVec (m : ℕ) (a : Fin m) : Fin m → ℝ :=
  fun i => (if i = a then (1:ℝ) else 0) - 1/(m:ℝ)

/-- 标准内积（自给自足，不用 Mathlib 的内积实例）。 -/
def z021Dot {m : ℕ} (x y : Fin m → ℝ) : ℝ := ∑ i, x i * y i

/-! ### 内部辅助引理

三个只在本文件内使用的技术引理：把"指示函数"的和算掉，
从而把内积 ⟨u_a, u_b⟩ 化成单项式表达式（见 coupling_verdict.md §2.2）。 -/

/-- 单点指示函数：`i = a` 取 1，否则取 0。 -/
private lemma z021_sum_ind {m : ℕ} (a : Fin m) :
    (∑ i : Fin m, (if i = a then (1:ℝ) else 0)) = 1 := by
  rw [Fintype.sum_ite_eq']

/-- 两个指示函数乘积的和：`∑_i [i=a]·[i=b] = [a=b]`（即 Gram 矩阵的对角性）。 -/
private lemma z021_sum_ind_mul {m : ℕ} (a b : Fin m) :
    (∑ i : Fin m, (if i = a then (1:ℝ) else 0) * (if i = b then (1:ℝ) else 0))
      = if a = b then 1 else 0 := by
  have h : (∑ i : Fin m, (if i = a then (1:ℝ) else 0) * (if i = b then (1:ℝ) else 0))
      = ∑ i : Fin m, (if i = a then (if i = b then (1:ℝ) else 0) else 0) := by
    refine Finset.sum_congr rfl (fun i _ => ?_)
    by_cases hia : i = a <;> simp [hia]
  rw [h, Fintype.sum_ite_eq']

/-- 关键恒等式（coupling_verdict.md §2.2 的全部内容）：
    `⟨u_a, u_b⟩ = δ_ab − 1/m`。展开后用上面两条指示函数引理收尾。 -/
private lemma z021Dot_simplexVec (m : ℕ) (hm : 0 < m) (a b : Fin m) :
    z021Dot (z021SimplexVec m a) (z021SimplexVec m b)
      = (if a = b then (1:ℝ) else 0) - 1/(m:ℝ) := by
  have expand : ∀ i : Fin m, z021SimplexVec m a i * z021SimplexVec m b i
      = (if i = a then (1:ℝ) else 0) * (if i = b then (1:ℝ) else 0)
        - (1/(m:ℝ)) * ((if i = a then (1:ℝ) else 0) + (if i = b then (1:ℝ) else 0)
            - 1/(m:ℝ)) := by
    intro i
    simp only [z021SimplexVec]
    ring
  rw [z021Dot]
  simp only [expand]
  rw [Finset.sum_sub_distrib, ← Finset.mul_sum]
  rw [Finset.sum_sub_distrib, Finset.sum_add_distrib]
  rw [z021_sum_ind, z021_sum_ind, z021_sum_ind_mul]
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin]
  simp only [nsmul_eq_mul]
  have hm' : (m:ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hm.ne'
  field_simp
  ring

/-! ### 主要定理（与 coupling_verdict.md 逐条对应） -/

/-- **§1.3 约束二（零和）**：`∑_a u_a = 0`。
    即 m 个通道向量都落在超平面 `H_Q = {x : ∑_a x_a = 0}` 内（`dim H_Q = m − 1 = D`）。
    文档 §2.2 表里这一行记的是"残差 ~10⁻¹⁶"（数值）；这里是精确等式。 -/
theorem z021_simplex_sum_zero (m : ℕ) (hm : 0 < m) :
    (∑ a : Fin m, z021SimplexVec m a) = 0 := by
  funext i
  rw [Finset.sum_apply, Pi.zero_apply]
  simp only [z021SimplexVec]
  rw [Finset.sum_sub_distrib, Fintype.sum_ite_eq]
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin]
  simp only [nsmul_eq_mul]
  have hm' : (m:ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hm.ne'
  field_simp
  ring

/-- **§2.2 等长（自内积）**：`⟨u_a, u_a⟩ = 1 − 1/m`。
    对应文档 §2.2 表"等长"一行（那里写的是归一化后的 `|u_a| = 1`）。 -/
theorem z021_simplex_self (m : ℕ) (hm : 0 < m) (a : Fin m) :
    z021Dot (z021SimplexVec m a) (z021SimplexVec m a) = 1 - 1/(m:ℝ) := by
  rw [z021Dot_simplexVec m hm a a]
  simp

/-- **§2.2 等角（互内积）**：`a ≠ b ⟹ ⟨u_a, u_b⟩ = −1/m`。
    对应文档 §2.2 表"等角：u_a·u_b = −1/D"一行（未归一化时右端是 −1/m）。 -/
theorem z021_simplex_pair (m : ℕ) (hm : 1 < m) (a b : Fin m) (hab : a ≠ b) :
    z021Dot (z021SimplexVec m a) (z021SimplexVec m b) = -1/(m:ℝ) := by
  rw [z021Dot_simplexVec m (by omega) a b, ite_eq_right hab]
  ring

/-- **§3.1 归一化夹角余弦**：`ρ = ⟨u_a,u_b⟩/⟨u_a,u_a⟩ = −1/(m−1) = −1/D`。
    这就是文档里"夹角 arccos(−1/D)"的代数内容（`D = m − 1`）。 -/
theorem z021_simplex_rho (m : ℕ) (hm : 1 < m) (a b : Fin m) (hab : a ≠ b) :
    z021Dot (z021SimplexVec m a) (z021SimplexVec m b) / z021Dot (z021SimplexVec m a) (z021SimplexVec m a)
      = -1/((m:ℝ) - 1) := by
  rw [z021_simplex_pair m hm a b hab, z021_simplex_self m (by omega) a]
  have hm0 : (m:ℝ) ≠ 0 := by positivity
  have hm1 : (m:ℝ) - 1 ≠ 0 := by
    have h : (1:ℝ) < (m:ℝ) := by exact_mod_cast hm
    linarith
  have h1 : (1:ℝ) - 1/(m:ℝ) = ((m:ℝ) - 1)/(m:ℝ) := by field_simp
  rw [h1, div_eq_iff (div_ne_zero hm1 hm0)]
  field_simp

/-- **§1.3 约束三（通道无偏好）→ §2.2 等长**：所有通道向量自内积相同（"通道无偏好"的一半）。 -/
theorem z021_simplex_equilateral (m : ℕ) (hm : 1 < m) (a b : Fin m) :
    z021Dot (z021SimplexVec m a) (z021SimplexVec m a)
      = z021Dot (z021SimplexVec m b) (z021SimplexVec m b) := by
  rw [z021_simplex_self m (by omega) a, z021_simplex_self m (by omega) b]

/-- **§2.3 逐维核验表，D = 4 行（m = 5）**：夹角余弦 = −1/4（104.4775°）。 -/
theorem z021_simplex_rho_five (a b : Fin 5) (hab : a ≠ b) :
    z021Dot (z021SimplexVec 5 a) (z021SimplexVec 5 b) / z021Dot (z021SimplexVec 5 a) (z021SimplexVec 5 a)
      = -1/4 := by
  rw [z021_simplex_rho 5 (by norm_num) a b hab]
  norm_num

/-- **§2.3 逐维核验表，D = 3 行（m = 4，四面体）**：夹角余弦 = −1/3（109.4712°）。 -/
theorem z021_simplex_rho_tetra (a b : Fin 4) (hab : a ≠ b) :
    z021Dot (z021SimplexVec 4 a) (z021SimplexVec 4 b) / z021Dot (z021SimplexVec 4 a) (z021SimplexVec 4 a)
      = -1/3 := by
  rw [z021_simplex_rho 4 (by norm_num) a b hab]
  norm_num

-- 审计修订（2026-10-14 夜）：删去 `z021_tetra_angle_cos` —— 它是 `cos∘arccos = id` 的往返
-- 恒等式，与四面体角无独立内容（真内容是 `z021_simplex_rho_tetra` 的 ρ = −1/3）。

-- 同时删去 `z021_simplex_counts`（纯数字算术 C(4,2)=6 ∧ C(3,2)=3）。
-- 判据：只读定理类型、不读注释；依赖清单干净不等于有内容。

#print axioms z021_simplex_sum_zero
#print axioms z021_simplex_self
#print axioms z021_simplex_pair
#print axioms z021_simplex_rho
#print axioms z021_simplex_equilateral
#print axioms z021_simplex_rho_five
#print axioms z021_simplex_rho_tetra

end

end Zero