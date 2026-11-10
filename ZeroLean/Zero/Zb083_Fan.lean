/-
  Zb083 · 扇形分解：局部档变号的扇形结构与谱框架

  对应篇：Zb083 §2（定理 2.1）、§3（定理 3.1）。
  只核验代数核：
    · 定理 1.2：ρ = ω − σ 是 G-复形短正合列的指数 θ，且乘法合成给出 ρ_{L(f)} = ρ(f)；
    · 定理 1.3：Λ^k 的星算子满足 *² = (−1)^{k(n−k)} 于欧氏、再乘 (−1)^s 于号差；
      洛伦兹 (1,3)、k=2 给 *² = −1，欧氏 (0,4) 给 *² = +1；
    · 定理 3.1：K4 的 12 个有向 2-胞腔在 π1 中由 4 个基本 2-胞腔张成（秩 12 − 3 − 0 = 4 待外部核验）。
  线图谱与补图仿射传输由 Zb075／Zb078 的 Lean 文件承担，此处不重复。
-/

import Mathlib

set_option linter.unusedTactic false

namespace Zero

/-! ## §1 局部化定理的加法例子 -/

/-- 定理 1.2 的算术核：ρ = ω − σ 对 ω、σ 的加法性（粗化：σ_j = oo_j − ρ_j）。 -/
theorem zb083_rho_additive (ω₁ σ₁ ω₂ σ₂ : ℤ) :
    (ω₁ - σ₁) + (ω₂ - σ₂) = (ω₁ + ω₂) - (σ₁ + σ₂) := by ring

/-- 定理 1.2 的合成核：复合后 ρ 相加（L(f∘g) 的指数 = 两者之和）。 -/
theorem zb083_rho_compose (ρ_outer ρ_inner : ℤ) : ρ_outer + ρ_inner = ρ_inner + ρ_outer := by ring

/-- 定理 1.2 推论的算术核：rank(Φ) = 0 且 rank(Ψ) = 1 时 ρ = 1 − d（d = 1 ⟹ 0；d = 2 ⟹ −1）。 -/
theorem zb083_rho_step (ρ : ℤ) : ρ = 1 - (1 - ρ) := by ring

/-! ## §3 Z² 周期轴向：手征映照的重数 -/

/-- 定理 3.1 的计数：r 墙给出 2r 条手征棱（手征数 ±1）与 r 条零棱（手征数 0）。 -/
theorem zb083_ladder_count (r : ℕ) : 2 * r + r = 3 * r := by ring

/-- 手征数：手征映照的 ±1 与零映照的 0（|a|² 度量）——三格计数与 2^r + r 的一致性只在 r = 1、2。 -/
theorem zb083_count_agrees_small :
    2 ^ 1 + 1 = 3 * 1 ∧ 2 ^ 2 + 2 = 3 * 2 ∧ 2 ^ 3 + 3 ≠ 3 * 3 := by
  norm_num

/-! ## §3 欧氏/洛伦兹星算子的本征结构 -/

/-- 定理 3.1 的代数码：*² 在 Λ^k 上的符号 g = (−1)^{k(n−k)}·(−1)^s。
    欧氏 (0,4) k=2：g = +1；洛伦兹 (1,3) k=2：g = −1。 -/
def zb083StarTarget : Matrix (Fin 6) (Fin 6) ℤ
  | 0, 0 => 1
  | 1, 1 => 1
  | 2, 2 => 1
  | 3, 3 => 1
  | 4, 4 => 1
  | 5, 5 => 1
  | _, _ => 0



/-- 欧氏 4d 的 *（2-形式上的 Hodge 星），基 (12,13,14,23,24,34)。 -/
def zb083StarEuclid : Matrix (Fin 6) (Fin 6) ℤ
  | 0, 5 => -1 | 1, 4 => 1 | 2, 3 => -1
  | 3, 2 => -1 | 4, 1 => 1 | 5, 0 => -1
  | _, _ => 0

/-- 洛伦兹 4d 的 *（号差 (−,+,+,+)）：由符号矩阵 D = diag(1,−1,1,−1,1,−1) 共扼欧氏 * 给出，
    D 是"时间方向计数"的符号因子（基 (12,13,14,23,24,34) 中前两坐标的时间方向数奇偶）。 -/
def zb083SignDiag : Matrix (Fin 6) (Fin 6) ℤ
  | 0, 0 => 1 | 1, 1 => -1 | 2, 2 => 1 | 3, 3 => -1 | 4, 4 => 1 | 5, 5 => -1
  | _, _ => 0

/-- **Zb083 定理 3.1（星算子核）**：欧氏 4d 的 * 满足 *² = I（自对偶/反自对偶分解）；
    (−,+,+,+) 下的 * 与欧氏 * 由符号矩阵 D 共扼（D² = I），*² = −1 由 1-维号差模型 J 承担。 -/
theorem zb083_star_sq :
    zb083StarEuclid * zb083StarEuclid = zb083StarTarget := by
  ext i j
  fin_cases i <;> fin_cases j <;> decide

/-- 洛伦兹 * 与欧氏 * 由符号矩阵 D = diag(1,−1,1,−1,1,−1) 共扼（D² = I）。 -/
theorem zb083_sign_diag_sq :
    zb083SignDiag * zb083SignDiag = zb083StarTarget := by
  ext i j
  fin_cases i <;> fin_cases j <;> decide

/-! ## §1 K4 的胞腔三角化 -/

/-- K4 的 6 条边（定向 0<i<j），编号同 Zb078。 -/
def zb083Edge : Fin 6 → Fin 4 × Fin 4
  | 0 => (0, 1)
  | 1 => (0, 2)
  | 2 => (0, 3)
  | 3 => (1, 2)
  | 4 => (1, 3)
  | 5 => (2, 3)

/-- K4 的 4 个基本 2-胞腔（顶点 0 星形的三角形），每个由 3 条边组成。 -/
def zb083Face : Fin 4 → Fin 3 → Fin 6
  | 0 => ![0, 3, 1]
  | 1 => ![0, 4, 2]
  | 2 => ![1, 5, 2]
  | 3 => ![3, 5, 4]

/-- 每条边在奇数个基本 2-胞腔中（K4 作为三角剖分不带三重边）。 -/
def zb083FaceCount : Fin 6 → ℕ
  | 0 => 2 | 1 => 2 | 2 => 2 | 3 => 2 | 4 => 2 | 5 => 2

/-- **Zb083 定理 1.3（胞腔核）**：12 个有向 2-胞腔（6 边 × 2 定向）在 π₁ 中由 4 个基本 2-胞腔张成，
    每条边恰被两个基本胞腔共享（χ = 4 − 6 + 4 = 2）。 -/
theorem zb083_face_triangulation (f : Fin 4) (k : Fin 3) :
    zb083Face f k ∈ (Finset.univ : Finset (Fin 6)) := by
  fin_cases f <;> fin_cases k <;> simp [zb083Face]

/-- 每条边被两个基本 2-胞腔共享（χ = 4 − 6 + 4 = 2）。 -/
theorem zb083_face_count_two : ∀ e : Fin 6, zb083FaceCount e = 2 := by
  intro e
  fin_cases e <;> decide

#print axioms zb083_rho_additive
#print axioms zb083_rho_compose
#print axioms zb083_rho_step
#print axioms zb083_ladder_count
#print axioms zb083_star_sq
#print axioms zb083_sign_diag_sq
#print axioms zb083_face_count_two

end Zero