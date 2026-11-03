/-
  Zb048 §3 的机器检查版本：`S₄` 在三个完美匹配上的作用、核的阶、满射。

  口径：
    · `Equiv.Perm (Fin 4)` 作用在边上：`e ↦ e.image σ`；
    · 作用在匹配上：`M ↦ M.image (actEdge σ)`；
    · 三个匹配 `M1,M2,M3` 编号为 `Fin 3`，作用后的编号由 `idx σ i` 给出；
    · 24 个置换是有限类型，全部用 `decide` 核验：作用良定义、`idx σ` 是双射、
      `idx` 满足乘法、核的阶为 4、像的阶为 6。

  诚实边界：这里证的是「作用在 3 个匹配上给出 S₄→S₃ 的满同态，核阶 4」；
  核**恰为** `V₄` 的四个具体置换的列举留作下一批（阶 4 已足够支撑定理 3.2 的计数）。
-/

import Zero.Zb048_K4
import Zero.Zb048_Count

namespace Zero

open Finset

/-- 置换作用在一条边上：逐点像。 -/
def zb048ActEdge (σ : Equiv.Perm (Fin 4)) (e : Finset (Fin 4)) : Finset (Fin 4) :=
  e.image σ

/-- 置换作用在一个匹配上：逐边像。 -/
def zb048ActM (σ : Equiv.Perm (Fin 4)) (M : Finset (Finset (Fin 4))) :
    Finset (Finset (Fin 4)) :=
  M.image (zb048ActEdge σ)

/-- 三个匹配的编号。 -/
def zb048MatchOf : Fin 3 → Finset (Finset (Fin 4))
  | 0 => zb048M1
  | 1 => zb048M2
  | 2 => zb048M3

/-- 作用后的匹配编号：`σ · M_i = M_{idx σ i}`。 -/
def zb048Idx (σ : Equiv.Perm (Fin 4)) (i : Fin 3) : Fin 3 :=
  if zb048ActM σ (zb048MatchOf i) = zb048M1 then 0
  else if zb048ActM σ (zb048MatchOf i) = zb048M2 then 1
  else 2

/-- **Zb048 · 作用良定义**：每个置换把每个匹配送到三个匹配之一。 -/
theorem zb048_action_maps :
    ∀ (σ : Equiv.Perm (Fin 4)) (i : Fin 3),
      zb048ActM σ (zb048MatchOf i) = zb048MatchOf (zb048Idx σ i) := by
  decide

/-- **Zb048 · 诱导映射是双射**：`idx σ : Fin 3 → Fin 3` 是置换。 -/
theorem zb048_idx_bijective :
    ∀ σ : Equiv.Perm (Fin 4), Function.Bijective (zb048Idx σ) := by
  decide

/-- **Zb048 · 诱导映射是同态**：`idx (σ * τ) = idx σ ∘ idx τ`。 -/
theorem zb048_idx_mul :
    ∀ (σ τ : Equiv.Perm (Fin 4)) (i : Fin 3),
      zb048Idx (σ * τ) i = zb048Idx σ (zb048Idx τ i) := by
  decide

/-- **Zb048 · 核的阶 = 4**：同时固定三个匹配的置换恰有 4 个。 -/
theorem zb048_kernel_card :
    (Finset.univ.filter (fun σ : Equiv.Perm (Fin 4) => ∀ i : Fin 3, zb048Idx σ i = i)).card = 4 := by
  decide

/-- **Zb048 · 像的阶 = 6**：24 个置换在 3 个匹配上给出 6 个不同的置换。 -/
theorem zb048_image_card :
    (Finset.univ.image (fun σ : Equiv.Perm (Fin 4) => (zb048Idx σ : Fin 3 → Fin 3))).card = 6 := by
  decide

/-! ## 共扼闭包 -/

theorem zb048_perfect_actM (σ : Equiv.Perm (Fin 4)) (i : Fin 3) :
    zb048IsPerfectMatching (zb048ActM σ (zb048MatchOf i)) := by
  rw [zb048_action_maps]
  generalize hk : zb048Idx σ i = k
  fin_cases k <;>
    simp [zb048MatchOf, zb048_M1_perfect, zb048_M2_perfect, zb048_M3_perfect]

/-- A：三个匹配的集合在 S4 作用下不变（`three_matchings_unordered` 的共扼闭包）。 -/
theorem zb048_perfect_set_S4 :
    ∀ σ : Equiv.Perm (Fin 4),
      ({zb048M1, zb048M2, zb048M3} : Finset (Finset (Finset (Fin 4)))).image (zb048ActM σ) =
        {zb048M1, zb048M2, zb048M3} := by
  decide

/-- A：核在 S4 内自同构下不变（`kernel_eq_V4` 的共扼闭包）。 -/
theorem zb048_kernel_conj :
    ∀ τ : Equiv.Perm (Fin 4),
      (Finset.univ.filter
          (fun σ : Equiv.Perm (Fin 4) => ∀ i : Fin 3, zb048Idx (τ * σ * τ⁻¹) i = i)) =
        (Finset.univ.filter
          (fun σ : Equiv.Perm (Fin 4) => ∀ i : Fin 3, zb048Idx σ i = i)) := by
  decide

/-- B：`idx` 同态在共扼下的显式公式（`idx_mul` 的共扼形式）。 -/
theorem zb048_idx_conj (τ σ : Equiv.Perm (Fin 4)) (i : Fin 3) :
    zb048Idx (τ * σ * τ⁻¹) i =
      zb048Idx τ (zb048Idx σ (zb048Idx τ⁻¹ i)) := by
  rw [zb048_idx_mul, zb048_idx_mul]

#print axioms zb048_action_maps
#print axioms zb048_idx_bijective
#print axioms zb048_idx_mul
#print axioms zb048_kernel_card

/-- 核的四个元素：`1` 与三个双对换。 -/
def zb048V4 : Finset (Equiv.Perm (Fin 4)) :=
  {1,
   Equiv.swap (0 : Fin 4) 1 * Equiv.swap (2 : Fin 4) 3,
   Equiv.swap (0 : Fin 4) 2 * Equiv.swap (1 : Fin 4) 3,
   Equiv.swap (0 : Fin 4) 3 * Equiv.swap (1 : Fin 4) 2}

/-- **Zb048 · 核恰为 `V₄`**：同时固定三个匹配的置换恰是 `1` 与三个双对换。 -/
theorem zb048_kernel_eq_V4 :
    (Finset.univ.filter (fun σ : Equiv.Perm (Fin 4) => ∀ i : Fin 3, zb048Idx σ i = i))
      = zb048V4 := by
  decide

#print axioms zb048_kernel_eq_V4

#print axioms zb048_image_card

#print axioms zb048_perfect_actM
#print axioms zb048_perfect_set_S4
#print axioms zb048_kernel_conj
#print axioms zb048_idx_conj

end Zero
