/-
  Zb076 · 结构共扼第一批：Lean 机检的七条代数核

  对应篇：Zb076 §1–§3。
    · C02 转置：rank Aᵀ = rank A；
    · C09 Schur 补：det (fromBlocks A B C D) = det A * det (D - C A⁻¹ B)；
    · C04 二部符号：D A D = -A（块形状）；
    · C05 Seidel 切换：S(A^σ) = D S(A) D（A 对角为零）；
    · C06 面积反序：area(w.reverse) = -area(w)（平衡词）；
    · C07 面积取反：area(-w) = -area(w)；
    · C10 子群共扼：Stab(g • x) = g Stab(x) g⁻¹。
  其余条目（线图谱、约束转置）由脚本
  `ZeroGPU/src/zb076_conj_batch.py` 数值核验；Lean 化列入下一批。
-/

import Mathlib

namespace Zero

open Matrix

/-- **Zb076 C02**：转置不改变秩。 -/
theorem zb076_rank_transpose {m n : Type*} [Fintype m] [Fintype n]
    {R : Type*} [Field R] (A : Matrix m n R) : Aᵀ.rank = A.rank :=
  Matrix.rank_transpose A

/-- **Zb076 C09**：分块行列式 = `det A * det (D - C A⁻¹ B)`（Schur 补）。 -/
theorem zb076_det_schur {m n : Type*} [Fintype m] [Fintype n]
    [DecidableEq m] [DecidableEq n]
    (A : Matrix m m ℝ) (B : Matrix m n ℝ) (C : Matrix n m ℝ) (D : Matrix n n ℝ) [Invertible A] :
    (Matrix.fromBlocks A B C D).det = A.det * (D - C * ⅟A * B).det :=
  Matrix.det_fromBlocks₁₁ A B C D

/-- **Zb076 C10**：稳定子在共扼下的变换 `Stab(g • x) = g Stab(x) g⁻¹`。 -/
theorem zb076_stabilizer_conj {G X : Type*} [Group G] [MulAction G X] (g : G) (x : X) :
    MulAction.stabilizer G (g • x) =
      (MulAction.stabilizer G x).map (MulAut.conj g).toMonoidHom :=
  MulAction.stabilizer_smul_eq_stabilizer_map_conj g x

/-- 二部邻接的块形状：`[[0, B], [Bᵀ, 0]]`（C04）。 -/
def zb076BipAdj {m n : ℕ} (B : Matrix (Fin m) (Fin n) ℝ) :
    Matrix (Sum (Fin m) (Fin n)) (Sum (Fin m) (Fin n)) ℝ :=
  Matrix.fromBlocks 0 B Bᵀ 0

/-- 部符号：左 +1、右 −1（C04）。 -/
def zb076BipSign (m n : ℕ) : Matrix (Sum (Fin m) (Fin n)) (Sum (Fin m) (Fin n)) ℝ :=
  Matrix.diagonal (Sum.elim (fun _ => (1 : ℝ)) (fun _ => (-1 : ℝ)))

/-- **Zb076 C04**：二部符号共扼 `D A D = -A`。 -/
theorem zb076_bipartite_sign {m n : ℕ} (B : Matrix (Fin m) (Fin n) ℝ) :
    zb076BipSign m n * zb076BipAdj B * zb076BipSign m n = -zb076BipAdj B := by
  ext i j
  cases i <;> cases j <;>
    simp [zb076BipAdj, zb076BipSign, Matrix.fromBlocks]

/-- Seidel 矩阵 `S(A) = J - I - 2A`（C05）。 -/
def zb076Seidel {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  fun i j => 1 - (if i = j then 1 else 0) - 2 * A i j

/-- Seidel 切换图（C05）。 -/
noncomputable def zb076Switch {n : ℕ} (σ : Fin n → ℝ) (A : Matrix (Fin n) (Fin n) ℝ) :
    Matrix (Fin n) (Fin n) ℝ :=
  fun i j => if i = j then 0 else if σ i * σ j = 1 then A i j else 1 - A i j

/-- **Zb076 C05**：Seidel 切换等价 `S(A^σ) = D S(A) D`（`A` 的对角为零）。 -/
theorem zb076_seidel_switch {n : ℕ} (σ : Fin n → ℝ) (A : Matrix (Fin n) (Fin n) ℝ)
    (hσ : ∀ i, σ i = 1 ∨ σ i = -1) (hA : ∀ i, A i i = 0) :
    zb076Seidel (zb076Switch σ A) =
      Matrix.diagonal σ * zb076Seidel A * Matrix.diagonal σ := by
  have hmain : ∀ (a s t : ℝ), (s = 1 ∨ s = -1) → (t = 1 ∨ t = -1) →
      (1 - (if s * t = 1 then 2 * a else 2 * (1 - a))) = s * (1 - 2 * a) * t := by
    intro a s t hs ht
    rcases hs with rfl | rfl <;> rcases ht with rfl | rfl <;> norm_num <;> ring
  ext i j
  by_cases hij : i = j
  · subst hij
    simp [zb076Seidel, zb076Switch, hA i]
  · simpa [zb076Seidel, zb076Switch, hij] using hmain (A i j) (σ i) (σ j) (hσ i) (hσ j)

/-- 词的前缀和序列（Zb076 §2）：从初值 b 出发，逐步加 w_i。 -/
def zb076Scan : List ℤ → ℤ → List ℤ
  | [], b => [b]
  | x :: xs, b => b :: zb076Scan xs (b + x)

/-- 词的前缀和之和（面积）。 -/
def zb076Area (w : List ℤ) : ℤ := (zb076Scan w 0).sum

/-- 取反把前缀和整体取负（带任意初值）。 -/
theorem zb076_scan_neg (w : List ℤ) (b : ℤ) :
    zb076Scan (w.map (fun x => -x)) (-b) = (zb076Scan w b).map (fun x => -x) := by
  induction w generalizing b with
  | nil => simp [zb076Scan]
  | cons x xs ih =>
      simp [zb076Scan]
      have h : -b + -x = -(b + x) := by ring
      rw [h]
      exact ih (b + x)

/-- 列表求和取反。 -/
theorem zb076_sum_neg (l : List ℤ) :
    (l.map (fun x => -x)).sum = -l.sum := by
  induction l with
  | nil => simp
  | cons x xs ih =>
      simp [ih]
      ring

/-- **Zb076 C07**：面积在词取反下变号。 -/
theorem zb076_area_neg (w : List ℤ) :
    zb076Area (w.map (fun x => -x)) = -zb076Area w := by
  unfold zb076Area
  rw [show zb076Scan (w.map (fun x => -x)) 0 = (zb076Scan w 0).map (fun x => -x) by
    simpa using zb076_scan_neg w 0]
  exact zb076_sum_neg (zb076Scan w 0)

/-- 前缀和序列的和：`sum (scan w b) = (|w|+1) b + sum (scan w 0)`。 -/
theorem zb076Scan_sum (w : List ℤ) (b : ℤ) :
    (zb076Scan w b).sum = (w.length + 1) * b + (zb076Scan w 0).sum := by
  induction w generalizing b with
  | nil => simp [zb076Scan]
  | cons x xs ih =>
      rw [zb076Scan.eq_2, zb076Scan.eq_2]
      simp only [List.sum_cons, List.length_cons, Nat.cast_add, Nat.cast_one, zero_add]
      rw [ih (b + x), ih x]
      ring

/-- 面积的前插公式：`area (x :: xs) = (|xs|+1) x + area xs`。 -/
theorem zb076Area_cons (x : ℤ) (xs : List ℤ) :
    zb076Area (x :: xs) = (xs.length + 1) * x + zb076Area xs := by
  unfold zb076Area
  rw [zb076Scan.eq_2]
  simp only [List.sum_cons, zero_add]
  rw [zb076Scan_sum xs x]

/-- 面积的尾接公式：`area (xs ++ [x]) = area xs + sum xs + x`。 -/
theorem zb076Area_append_singleton (xs : List ℤ) (x : ℤ) :
    zb076Area (xs ++ [x]) = zb076Area xs + xs.sum + x := by
  induction xs with
  | nil => simp [zb076Area, zb076Scan]
  | cons y ys ih =>
      rw [List.cons_append, zb076Area_cons, ih, zb076Area_cons]
      simp only [List.sum_cons, List.length_append, List.length_singleton,
        Nat.cast_add, Nat.cast_one]
      ring

/-- **Zb076 C06**：`area(w^R) + area(w) = (|w|+1) sum(w)`。 -/
theorem zb076Area_reverse_add (w : List ℤ) :
    zb076Area w.reverse + zb076Area w = (w.length + 1) * w.sum := by
  induction w with
  | nil => simp [zb076Area, zb076Scan]
  | cons x xs ih =>
      rw [List.reverse_cons, zb076Area_append_singleton, zb076Area_cons]
      rw [List.sum_reverse]
      simp only [List.sum_cons, List.length_cons, Nat.cast_add, Nat.cast_one] at *
      linarith

/-- **Zb076 C06**：平衡词上面积反号。 -/
theorem zb076_area_reverse (w : List ℤ) (h : w.sum = 0) :
    zb076Area w.reverse = -zb076Area w := by
  have hh := zb076Area_reverse_add w
  rw [h] at hh
  linarith

#print axioms zb076_rank_transpose
#print axioms zb076_det_schur
#print axioms zb076_stabilizer_conj
#print axioms zb076_bipartite_sign
#print axioms zb076_seidel_switch
#print axioms zb076_area_neg
#print axioms zb076Area_reverse_add
#print axioms zb076_area_reverse

end Zero
