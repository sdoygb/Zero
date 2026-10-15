/-
  Z020（确定性＋有限状态 ⟹ 轨道终将周期）的机器检查版本。

  审计修订（2026-10-14）：原版**从未编译通过**（用了不存在的 `Fintype.card_univ`，
  且证明结构有类型错误）——它只是"磁盘上的文件"，没进过 CI。这里重写：
  用鸽巢原理（`Fintype.card_le_of_injective`）：若前 `N+1` 个迭代点两两不同，
  就得到 `N+1 ≤ N` 的矛盾。

  本轮加强（2026-10-14 夜，来自 `ZeroGPU/doc/longterm_eval.md` 的"长时演化否证"）：
    · **终将周期的严格版**：给出预周期 `ν` 与周期 `μ > 0`，并证明从 `ν` 起**逐点**周期；
    · **环上记录层必为常数**：若观测量沿轨道单调不减（时间箭头），而轨道周期为 `μ`，
      则观测量在整条环上取**常值** —— 即"记录层在环上继续增长"与"确定性＋有限状态"不相容。
-/

import Mathlib

namespace Zero

/-- **Z020 · 确定性有限系统终将周期**：有限状态集上的自映射，迭代序列必有重复
    （即轨道最终进入循环）。 -/
theorem z020_iterate_eventually_periodic {α : Type*} [Fintype α] [DecidableEq α]
    (f : α → α) (x : α) : ∃ i j : ℕ, i < j ∧ f^[i] x = f^[j] x := by
  by_contra h
  push Not at h
  have hinj : Function.Injective (fun n : Fin (Fintype.card α + 1) => f^[(n : ℕ)] x) := by
    intro m n hmn
    by_contra hne
    have hne' : (m : ℕ) ≠ (n : ℕ) := fun h => hne (Fin.ext h)
    rcases lt_or_gt_of_ne hne' with hlt | hlt
    · exact h m n hlt hmn
    · exact h n m hlt hmn.symm
  have hle := Fintype.card_le_of_injective _ hinj
  simp only [Fintype.card_fin] at hle
  omega

/-- **Z020 · 终将周期的严格版**（`longterm_eval.md` §3 的定理）：
    存在预周期 `ν` 与周期 `μ > 0`，使轨道从 `ν` 起**逐点**以 `μ` 为周期。 -/
theorem z020_iterate_periodic_eventually {α : Type*} [Fintype α] [DecidableEq α]
    (f : α → α) (x : α) :
    ∃ ν μ : ℕ, 0 < μ ∧ f^[ν + μ] x = f^[ν] x ∧ ∀ k : ℕ, f^[ν + k + μ] x = f^[ν + k] x := by
  obtain ⟨i, j, hij, hfij⟩ := z020_iterate_eventually_periodic f x
  have hbase : f^[i + (j - i)] x = f^[i] x := by
    have h : i + (j - i) = j := by omega
    rw [h]; exact hfij.symm
  refine ⟨i, j - i, by omega, hbase, ?_⟩
  intro k
  have hidx : i + k + (j - i) = k + (i + (j - i)) := by omega
  have hidx2 : i + k = k + i := by omega
  calc f^[i + k + (j - i)] x = f^[k + (i + (j - i))] x := by rw [hidx]
    _ = f^[k] (f^[i + (j - i)] x) := by rw [Function.iterate_add_apply]
    _ = f^[k] (f^[i] x) := by rw [hbase]
    _ = f^[k + i] x := by rw [Function.iterate_add_apply]
    _ = f^[i + k] x := by rw [hidx2]

/-- **Z020 · 环上记录层必为常数**（`longterm_eval.md` §2.2 的严格化，同时是一处**否证**）：
    观测量 `r` 沿轨道单调不减（时间箭头）、轨道自第 `ν` 步起以 `μ` 为周期，则 `r` 在整个环上取常值。
    ★ 故实测数据（活动层在环上有 3 个取值、记录层跨 8.31 个数量级单调增长）按本条**不可能同时成立**，
    除非**记录层不是有限状态的一部分** —— 该数据是对"记录层属于状态"的否证，而不是"与连续性相容"。
    工程注记（审计修订）：原陈述多带了 `hμ : 0 < μ`，证明并未用到，且 `μ = 0` 时结论仍真，故删去。 -/
theorem z020_monotone_on_cycle_constant {α : Type*} (f : α → α) (r : α → ℚ) (x : α)
    (hmono : ∀ y, r y ≤ r (f y)) (ν μ : ℕ) (hper : f^[ν + μ] x = f^[ν] x) :
    ∀ i j : ℕ, i ≤ μ → j ≤ μ → r (f^[ν + i] x) = r (f^[ν + j] x) := by
  -- 单步单调
  have hstep : ∀ n : ℕ, r (f^[n] x) ≤ r (f^[n + 1] x) := by
    intro n
    have h : f^[n + 1] x = f (f^[n] x) := Function.iterate_succ_apply' f n x
    rw [h]
    exact hmono (f^[n] x)
  -- 轨道上单调不减
  have hmono' : ∀ n m : ℕ, n ≤ m → r (f^[n] x) ≤ r (f^[m] x) := by
    intro n m hnm
    obtain ⟨c, rfl⟩ := Nat.exists_eq_add_of_le hnm
    induction c with
    | zero => simp
    | succ c ih =>
      have h1 : r (f^[n + c] x) ≤ r (f^[(n + c) + 1] x) := hstep (n + c)
      have hidx : n + (c + 1) = (n + c) + 1 := by omega
      rw [hidx]
      exact le_trans (ih (by omega)) h1
  -- 状态自身逐点周期（把环上任意点平移回 ν 起的同一相位）
  have hstate : ∀ k : ℕ, f^[ν + k + μ] x = f^[ν + k] x := by
    intro k
    have hidx : k + (ν + μ) = ν + k + μ := by omega
    have hidx2 : k + ν = ν + k := by omega
    calc f^[ν + k + μ] x = f^[k + (ν + μ)] x := by rw [hidx]
      _ = f^[k] (f^[ν + μ] x) := by rw [Function.iterate_add_apply]
      _ = f^[k] (f^[ν] x) := by rw [hper]
      _ = f^[k + ν] x := by rw [Function.iterate_add_apply]
      _ = f^[ν + k] x := by rw [hidx2]
  intro i j hi hj
  rcases le_total i j with hij | hji
  · have h1 : r (f^[ν + i] x) ≤ r (f^[ν + j] x) := hmono' _ _ (by omega)
    have h2 : r (f^[ν + j] x) ≤ r (f^[ν + i + μ] x) := hmono' _ _ (by omega)
    have h3 : f^[ν + i + μ] x = f^[ν + i] x := hstate i
    rw [h3] at h2
    exact le_antisymm h1 h2
  · have h1 : r (f^[ν + j] x) ≤ r (f^[ν + i] x) := hmono' _ _ (by omega)
    have h2 : r (f^[ν + i] x) ≤ r (f^[ν + j + μ] x) := hmono' _ _ (by omega)
    have h3 : f^[ν + j + μ] x = f^[ν + j] x := hstate j
    rw [h3] at h2
    exact le_antisymm h2 h1

#print axioms z020_iterate_eventually_periodic
#print axioms z020_iterate_periodic_eventually
#print axioms z020_monotone_on_cycle_constant

end Zero