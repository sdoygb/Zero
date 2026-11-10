/-
  共扼检验第一条真实登记：Zb001 `zb001_kernel_eq_const`。

  作用：S₄ 顶点置换的一个生成元——交换 a、b 两个分量（`V4.swapAB`）。
  要证：`P001 v := (K4L v = 0 ↔ 四分量全等)` 在 swap 下不变。
  这条 missing 证出来后，`zb001_kernel_eq_const` 的共扼陈述即由 `conj_of_invariant` 得到。
-/

import Zero.Conjugation
import Zero.Zb001_Zero

namespace Zero

/-- 交换 V4 的 a、b 两个分量（S₄ 的一个换位生成元）。 -/
def V4.swapAB (v : V4) : V4 := ⟨v.b, v.a, v.c, v.d⟩

/-- `K4L` 与 swap 交换（K₄ Laplacian 的顶点置换等变性）。 -/
theorem K4L_swapAB (v : V4) : K4L (V4.swapAB v) = V4.swapAB (K4L v) := by
  cases v with
  | mk a b c d =>
    simp [K4L, V4.swapAB]
    constructor <;> omega

/-- 「四分量全等」在 swap 下不变。 -/
theorem swapAB_eq_const_iff (v : V4) :
    (V4.swapAB v).a = (V4.swapAB v).b ∧
      (V4.swapAB v).b = (V4.swapAB v).c ∧
      (V4.swapAB v).c = (V4.swapAB v).d ↔
    (v.a = v.b ∧ v.b = v.c ∧ v.c = v.d) := by
  simp only [V4.swapAB]
  constructor <;> intro h <;> omega

/-- swap 保持零向量。 -/
theorem swapAB_zero : V4.swapAB (V4.zero) = V4.zero := by
  rfl

/-- swap 可逆，故 `swapAB w = 0 ↔ w = 0`。 -/
theorem swapAB_eq_zero_iff (w : V4) : V4.swapAB w = V4.zero ↔ w = V4.zero := by
  constructor
  · intro h
    have hh := congrArg V4.swapAB h
    simp [V4.swapAB, V4.zero] at hh
    cases w with
    | mk a b c d =>
      simp_all [V4.zero]
  · intro h
    rw [h]
    exact swapAB_zero

/-- Zb001 内核定理的谓词。 -/
def P001 (v : V4) : Prop :=
  K4L v = V4.zero ↔ (v.a = v.b ∧ v.b = v.c ∧ v.c = v.d)

/-- **第一条共扼检验**：`P001` 在 swap 下不变。 -/
theorem P001_swap (v : V4) : P001 (V4.swapAB v) ↔ P001 v := by
  unfold P001
  rw [K4L_swapAB, swapAB_eq_const_iff, swapAB_eq_zero_iff]

/-- swap 是单射（其实是 involution）。 -/
theorem swapAB_injective : Function.Injective V4.swapAB := by
  intro x y h
  calc
    x = V4.swapAB (V4.swapAB x) := by rfl
    _ = V4.swapAB (V4.swapAB y) := by rw [h]
    _ = y := by rfl

/-- `smul 4` 与 swap 交换。 -/
theorem smul_swapAB (v : V4) : V4.smul 4 (V4.swapAB v) = V4.swapAB (V4.smul 4 v) := by
  rfl

/-- 四分量之和在 swap 下不变。 -/
theorem sum_swapAB (v : V4) :
    (V4.swapAB v).a + (V4.swapAB v).b + (V4.swapAB v).c + (V4.swapAB v).d
      = v.a + v.b + v.c + v.d := by
  simp [V4.swapAB]
  ring

/-- λ = 4 特征方程的 swap 等变。 -/
theorem K4L_smul_swap_iff (v : V4) :
    (K4L (V4.swapAB v) = V4.smul 4 (V4.swapAB v)) ↔
      (K4L v = V4.smul 4 v) := by
  rw [K4L_swapAB, smul_swapAB]
  exact ⟨fun h => swapAB_injective h, fun h => congrArg V4.swapAB h⟩

/-- Zb001 第二条的谓词。 -/
def P002 (v : V4) : Prop :=
  K4L v = V4.smul 4 v ↔ v.a + v.b + v.c + v.d = 0

/-- `P002` 在 swap 下不变。 -/
theorem P002_swap (v : V4) : P002 (V4.swapAB v) ↔ P002 v := by
  unfold P002
  rw [K4L_smul_swap_iff, sum_swapAB]

/-- Zb001 第三条的谓词（带假设的定理）。 -/
def P003 (v : V4) : Prop :=
  K4L v = V4.smul 4 v → v.a + v.b + v.c + v.d = 0

/-- `P003` 在 swap 下不变。 -/
theorem P003_swap (v : V4) : P003 (V4.swapAB v) ↔ P003 v := by
  unfold P003
  rw [K4L_smul_swap_iff, sum_swapAB]


/- ================= 完整 S₄ 作用 ================= -/

/-- 四分量向量表示。 -/
def V4.toVec (v : V4) : Fin 4 → Int := ![v.a, v.b, v.c, v.d]
/-- 从四分量向量构造。 -/
def V4.ofVec (f : Fin 4 → Int) : V4 := ⟨f 0, f 1, f 2, f 3⟩
/-- S₄ 作用：按 `σ` 重排四个分量。 -/
def V4.perm (σ : Equiv.Perm (Fin 4)) (v : V4) : V4 :=
  V4.ofVec (fun i => V4.toVec v (σ.symm i))

@[simp] theorem V4.toVec_ofVec (f : Fin 4 → Int) : V4.toVec (V4.ofVec f) = f := by
  funext i; fin_cases i <;> rfl
@[simp] theorem V4.ofVec_toVec (v : V4) : V4.ofVec (V4.toVec v) = v := by
  cases v; rfl
theorem V4.toVec_injective : Function.Injective V4.toVec := by
  intro x y h; have := congrArg V4.ofVec h; simpa using this

@[simp] theorem V4.toVec_zero : V4.toVec (V4.zero) = fun _ => 0 := by
  funext i; fin_cases i <;> rfl

@[simp] theorem V4.toVec_smul (k : Int) (v : V4) :
    V4.toVec (V4.smul k v) = fun i => k * V4.toVec v i := by
  funext i; fin_cases i <;> rfl

@[simp] theorem V4.toVec_perm (σ : Equiv.Perm (Fin 4)) (v : V4) :
    V4.toVec (V4.perm σ v) = fun i => V4.toVec v (σ.symm i) := by
  funext i
  simp [V4.perm]

theorem sum_four (f : Fin 4 → Int) :
    (∑ i : Fin 4, f i) = f 0 + f 1 + f 2 + f 3 := by
  rw [Fin.sum_univ_succ]
  simp [Fin.sum_univ_three]
  ring

theorem sum_four_perm (σ : Equiv.Perm (Fin 4)) (f : Fin 4 → Int) :
    f (σ.symm 0) + f (σ.symm 1) + f (σ.symm 2) + f (σ.symm 3)
      = f 0 + f 1 + f 2 + f 3 := by
  rw [← sum_four (fun i => f (σ.symm i)), ← sum_four f]
  exact Equiv.sum_comp σ.symm f

def K4fun (f : Fin 4 → Int) : Fin 4 → Int := fun i => 4 * f i - (f 0 + f 1 + f 2 + f 3)

theorem V4.toVec_K4L (v : V4) : V4.toVec (K4L v) = K4fun (V4.toVec v) := by
  funext i
  fin_cases i <;> simp [K4L, K4fun, V4.toVec] <;> ring

theorem K4L_perm (σ : Equiv.Perm (Fin 4)) (v : V4) :
    K4L (V4.perm σ v) = V4.perm σ (K4L v) := by
  apply V4.toVec_injective
  rw [V4.toVec_K4L]
  rw [V4.toVec_perm]
  rw [V4.toVec_perm, V4.toVec_K4L]
  funext i
  simp [K4fun, sum_four_perm, Equiv.symm_symm, Equiv.symm_apply_apply]

theorem V4.perm_inv_perm (σ : Equiv.Perm (Fin 4)) (v : V4) :
    V4.perm (Equiv.symm σ) (V4.perm σ v) = v := by
  apply V4.toVec_injective
  simp only [V4.toVec_perm]
  funext i
  simp [Equiv.symm_symm, Equiv.symm_apply_apply]

@[simp] theorem V4.perm_zero (σ : Equiv.Perm (Fin 4)) : V4.perm σ V4.zero = V4.zero := by
  apply V4.toVec_injective
  simp [V4.toVec_perm, V4.toVec_zero]

theorem V4.perm_injective (σ : Equiv.Perm (Fin 4)) : Function.Injective (V4.perm σ) := by
  intro x y h
  have := congrArg (V4.perm (Equiv.symm σ)) h
  simpa [V4.perm_inv_perm] using this

theorem V4.perm_eq_zero_iff (σ : Equiv.Perm (Fin 4)) (v : V4) :
    V4.perm σ v = V4.zero ↔ v = V4.zero := by
  constructor
  · intro h
    have := congrArg (V4.perm (Equiv.symm σ)) h
    rwa [V4.perm_inv_perm, V4.perm_zero] at this
  · intro h
    rw [h]
    exact V4.perm_zero σ

def allEq4 (f : Fin 4 → Int) : Prop := ∀ i j, f i = f j

theorem V4.allEq_iff (v : V4) :
    (v.a = v.b ∧ v.b = v.c ∧ v.c = v.d) ↔ allEq4 (V4.toVec v) := by
  constructor
  · rintro ⟨hab, hbc, hcd⟩ i j
    fin_cases i <;> fin_cases j <;> simp [V4.toVec] <;> omega
  · intro h
    exact ⟨by simpa [V4.toVec] using h 0 1,
           by simpa [V4.toVec] using h 1 2,
           by simpa [V4.toVec] using h 2 3⟩

theorem allEq4_perm (σ : Equiv.Perm (Fin 4)) (f : Fin 4 → Int) :
    allEq4 (fun i => f (σ.symm i)) ↔ allEq4 f := by
  constructor
  · intro h i j
    have := h (σ i) (σ j)
    simpa using this
  · intro h i j
    exact h (σ.symm i) (σ.symm j)

theorem V4.perm_allEq_iff (σ : Equiv.Perm (Fin 4)) (v : V4) :
    ((V4.perm σ v).a = (V4.perm σ v).b ∧
      (V4.perm σ v).b = (V4.perm σ v).c ∧
      (V4.perm σ v).c = (V4.perm σ v).d) ↔
    (v.a = v.b ∧ v.b = v.c ∧ v.c = v.d) := by
  rw [V4.allEq_iff, V4.allEq_iff]
  rw [V4.toVec_perm]
  exact allEq4_perm σ (V4.toVec v)

theorem smul_perm (σ : Equiv.Perm (Fin 4)) (v : V4) :
    V4.smul 4 (V4.perm σ v) = V4.perm σ (V4.smul 4 v) := by
  apply V4.toVec_injective
  simp [V4.toVec_perm, V4.toVec_smul]

theorem V4.sum_eq_toVec_sum (v : V4) :
    v.a + v.b + v.c + v.d = ∑ i : Fin 4, V4.toVec v i := by
  rw [sum_four]
  rfl

theorem sum_perm (σ : Equiv.Perm (Fin 4)) (v : V4) :
    (V4.perm σ v).a + (V4.perm σ v).b + (V4.perm σ v).c + (V4.perm σ v).d
      = v.a + v.b + v.c + v.d := by
  rw [V4.sum_eq_toVec_sum, V4.sum_eq_toVec_sum]
  rw [V4.toVec_perm]
  exact Equiv.sum_comp σ.symm (V4.toVec v)

theorem P001_perm (σ : Equiv.Perm (Fin 4)) (v : V4) :
    P001 (V4.perm σ v) ↔ P001 v := by
  unfold P001
  rw [K4L_perm, V4.perm_eq_zero_iff, V4.perm_allEq_iff]

theorem P002_perm (σ : Equiv.Perm (Fin 4)) (v : V4) :
    P002 (V4.perm σ v) ↔ P002 v := by
  unfold P002
  have hL : (K4L (V4.perm σ v) = V4.smul 4 (V4.perm σ v)) ↔
      (K4L v = V4.smul 4 v) := by
    rw [K4L_perm, smul_perm]
    exact ⟨fun h => V4.perm_injective σ h, fun h => congrArg (V4.perm σ) h⟩
  rw [hL, sum_perm]

theorem P003_perm (σ : Equiv.Perm (Fin 4)) (v : V4) :
    P003 (V4.perm σ v) ↔ P003 v := by
  unfold P003
  have hL : (K4L (V4.perm σ v) = V4.smul 4 (V4.perm σ v)) ↔
      (K4L v = V4.smul 4 v) := by
    rw [K4L_perm, smul_perm]
    exact ⟨fun h => V4.perm_injective σ h, fun h => congrArg (V4.perm σ) h⟩
  rw [hL, sum_perm]

/-- 登记：Zb001 三条定理的共扼检验。 -/
def conjRegistryZb001 : List ConjEntry := [
  ⟨"zb001_kernel_eq_const",
    ∀ (σ : Equiv.Perm (Fin 4)) (v : V4), P001 (V4.perm σ v) ↔ P001 v,
    by
      intro σ v
      exact P001_perm σ v⟩,
  ⟨"zb001_eigen_four_iff",
    ∀ (σ : Equiv.Perm (Fin 4)) (v : V4), P002 (V4.perm σ v) ↔ P002 v,
    by
      intro σ v
      exact P002_perm σ v⟩,
  ⟨"zb001_eigen_four_zero_sum",
    ∀ (σ : Equiv.Perm (Fin 4)) (v : V4), P003 (V4.perm σ v) ↔ P003 v,
    by
      intro σ v
      exact P003_perm σ v⟩,
  ⟨"K4L_eq_zero_iff",
    ∀ (σ : Equiv.Perm (Fin 4)) (v : V4), P001 (V4.perm σ v) ↔ P001 v,
    by
      intro σ v
      exact P001_perm σ v⟩,
  ⟨"K4L_eigen_four_iff",
    ∀ (σ : Equiv.Perm (Fin 4)) (v : V4), P002 (V4.perm σ v) ↔ P002 v,
    by
      intro σ v
      exact P002_perm σ v⟩,
  ⟨"eigen_four_zero_sum",
    ∀ (σ : Equiv.Perm (Fin 4)) (v : V4), P003 (V4.perm σ v) ↔ P003 v,
    by
      intro σ v
      exact P003_perm σ v⟩
]

end Zero
