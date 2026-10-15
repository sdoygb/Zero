/-
  Z021 §6／§7 的**耦合计数匹配**（ZeroGPU/doc/coupling_map.md）的机器检查版本。

  两条结构不同的假设，各自唯一给出 `D = 3`：
    · (3.1) 独立成对身份 ⟷ 格点轴：`C(D,2) = D`（这一半已在 `Z021_Dimension.lean` 里）；
    · (3.2) 成对关系 ⟷ 邻居方向：`C(D+1,2) = 2D`（本文件新增）。
  并记下两族的"方向集大小不同"（`m = D+1 ≠ 2D`，除 D=1 外），即耦合的通道方向→格点方向
  映射**不能是恒等**。

  数值互证：`ZeroGPU/src/coupling_map.py`（14 条断言）。
-/

import Mathlib
import Zero.Z021_Dimension
import Zero.Z021_Dict

namespace Zero

/-! ## 一、假设 (3.2)：成对关系数 = 邻居方向数 -/

/-- **Z021 · 耦合计数假设 (3.2)**：`C(D+1,2) = 2D ⟺ D = 0 ∨ D = 3`（非平凡解唯一是 `D = 3`）。
    `C(D+1,2) = D(D+1)/2`，方程化为 `D(D−3) = 0`。 -/
theorem z021_relation_dirs (D : ℕ) : (D + 1).choose 2 = 2 * D ↔ D = 0 ∨ D = 3 := by
  constructor
  · intro h
    have h2 : 2 * (D + 1).choose 2 = 2 * (2 * D) := by rw [h]
    have h3 : 2 * (D + 1).choose 2 = (D + 1) * D := by
      rw [z021_two_mul_choose_two, Nat.add_sub_cancel]
    have hmul : (D + 1) * D = 2 * (2 * D) := by omega
    have hmulz : ((D : ℤ) + 1) * (D : ℤ) = 2 * (2 * (D : ℤ)) := by exact_mod_cast hmul
    have hzero : (D : ℤ) * ((D : ℤ) - 3) = 0 := by nlinarith
    rcases mul_eq_zero.mp hzero with hD0 | hD3
    · left; omega
    · right; omega
  · rintro (rfl | rfl) <;> decide

/-- **Z021 · 假设 (3.2) 的非平凡解唯一**：`D ≥ 1` 时 `C(D+1,2) = 2D ⟺ D = 3`。 -/
theorem z021_relation_dirs_forced_three (D : ℕ) (hD : 1 ≤ D) :
    (D + 1).choose 2 = 2 * D ↔ D = 3 := by
  rw [z021_relation_dirs]
  constructor
  · rintro (h | h)
    · omega
    · exact h
  · intro h; exact Or.inr h

/-! ## 二、两条独立假设同解 -/

/-- **Z021 · 两条独立假设同解**：`D ≥ 1` 时
    (3.1) `C(D,2) = D`（独立成对身份 ⟷ 轴）与
    (3.2) `C(D+1,2) = 2D`（成对关系 ⟷ 邻居方向）
    都**且只**在 `D = 3` 成立 —— 单条假设的解可能是巧合，两条结构不同的假设同解使巧合解释力下降。 -/
theorem z021_coupling_two_hypotheses (D : ℕ) (hD : 1 ≤ D) :
    (D.choose 2 = D ↔ D = 3) ∧ ((D + 1).choose 2 = 2 * D ↔ D = 3) :=
  ⟨z021_dim_forced_three D (by omega), z021_relation_dirs_forced_three D hD⟩

/-! ## 三、通道方向与格点方向的数目不匹配 -/

/-- **Z021 · 方向集大小不同**：`D ≥ 2` 时通道数 `m = D+1` 与格点邻居数 `2D` 不相等，
    故"通道方向 → 格点方向"的映射**不能是恒等**（需要一条计数规则）。 -/
theorem z021_channel_ne_neighbours (D : ℕ) (hD : 2 ≤ D) : D + 1 ≠ 2 * D := by
  omega

-- 审计修订（2026-10-14 夜）：删去 `z021_coupling_counts_three`（纯数字合取，含 `2*3=6` 模式）。
-- 「D = 3 时两侧计数同时闭合」的内容由 `z021_coupling_two_hypotheses` 承担。

#print axioms z021_relation_dirs
#print axioms z021_relation_dirs_forced_three
#print axioms z021_coupling_two_hypotheses
#print axioms z021_channel_ne_neighbours

end Zero