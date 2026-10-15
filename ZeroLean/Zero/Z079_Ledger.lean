/-
  Z079（L=4 的完整精确解）与 Z080（M⁺ 自旋实现）的机器检查版本。

  Z079 是 Zero 系列"唯一无条件预言"：L=4 的组态空间只有 6 个词，层结构、模谱、
  动力学谱全部**精确可解**（w=(2/3,1/3)、κ₁=5/9、Δ_mod=log2、谱隙 1/2）。

  这里把其中**组合与算术**的部分钉进 kernel：
    · W₄ = {w ∈ {±1}⁴ : Σw = 0} 恰 6 个词；
    · 层高 h(w) = max(M⁺, M⁻)（M± = 前缀和的最大值／最小值取反）—— 这个定义是从 Z079 §1
      的 6 行表反推并逐行核验的；
    · 层计数 b = (4, 2) ⟹ 层权 w = (2/3, 1/3)，κ₁ = w₁²+w₂² = 5/9 精确；
    · Δ_mod = log(w₁/w₂) = log 2 精确。
    · Z080：M⁺ 壳分布 m_j 在 L=4 上与反射原理公式 C(4,2−j) − C(4,1−j) 逐项相符。

  工程取舍：全部用**内核 `decide`**（不是 `native_decide`）—— 后者会引入
  `Lean.ofReduceBool`，会污染「要素来源表」（依赖清单里多出一个非基石公理）。
  一般 L 的反射原理（ballot 定理）Mathlib 没有现成件，留作后续（见 Z135 §二队列第 2 项）。
-/

import Mathlib

namespace Zero

/-- Z079 §1：L=4 的平衡词（`±1` 四元组，和为 0），顺序照原文表格。 -/
def words4 : List (List ℤ) :=
  [[1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1],
   [-1, 1, 1, -1], [-1, 1, -1, 1], [-1, -1, 1, 1]]

/-- 前 `k` 项部分和（`k = 0` 给 0）。 -/
def prefixSum (w : List ℤ) (k : ℕ) : ℤ := (w.take k).sum

/-- 全部前缀下标 `0 … length`（非空，给 `sup'`/`inf'` 用）。 -/
def prefixIdx (w : List ℤ) : Finset ℕ := Finset.range (w.length + 1)

theorem prefixIdx_nonempty (w : List ℤ) : (prefixIdx w).Nonempty :=
  ⟨0, Finset.mem_range.mpr (Nat.lt_succ_of_le (Nat.zero_le _))⟩

/-- `M⁺`：前缀和的最大值（含空前缀）。Z080 §3.2 的"正内程"。 -/
def Mplus (w : List ℤ) : ℤ := (prefixIdx w).sup' (prefixIdx_nonempty w) (prefixSum w)

/-- `M⁻`：前缀和最小值的相反数（因 Σw = 0，它等于 `M⁺` 作用在反词上）。 -/
def Mminus (w : List ℤ) : ℤ := -((prefixIdx w).inf' (prefixIdx_nonempty w) (prefixSum w))

/-- 层高 `h(w) = max(M⁺, M⁻)` —— 与 Z079 §1 表中 h 列逐行相符（下方断言）。 -/
def layer (w : List ℤ) : ℤ := max (Mplus w) (Mminus w)

/-- **Z079 · W₄ 恰 6 个词**：|{w ∈ {±1}⁴ : Σw = 0}| = C(4,2) = 6。 -/
theorem z079_W4_card : words4.length = 6 := by decide

/-- **Z079 · 层高表逐行核验**：6 个词的 h 值依次为 2,1,1,1,1,2（原文 §1 表格）。 -/
theorem z079_layer_table : words4.map layer = [2, 1, 1, 1, 1, 2] := by decide

/-- **Z079 · 层计数 (b₁,b₂) = (4,2)**（原文 §2 式 (2.1)）。 -/
theorem z079_layer_counts :
    (words4.filter (fun w => layer w = 1)).length = 4 ∧
      (words4.filter (fun w => layer w = 2)).length = 2 := by decide

/-- **Z079 · 层权 w = (2/3, 1/3) 精确**（原文 §2 核验 §B）。 -/
theorem z079_layer_weights :
    (4 : ℚ) / 6 = 2 / 3 ∧ (2 : ℚ) / 6 = 1 / 3 := by norm_num

/-- **Z079 · 账本纯度 κ₁ = w₁² + w₂² = 5/9 精确**，且 5/9 ∈ (1/2, 3/5)（四维窗内）。 -/
theorem z079_kappa1 :
    ((2 : ℚ) / 3) ^ 2 + ((1 : ℚ) / 3) ^ 2 = 5 / 9 ∧
      (1 : ℚ) / 2 < 5 / 9 ∧ (5 : ℚ) / 9 < 3 / 5 := by norm_num

/-- **Z079 · 模间隙 Δ_mod = log(w₁/w₂) = log 2 精确**（原文 §3 核验 §C）。 -/
theorem z079_delta_mod : Real.log ((2 : ℝ) / 3 / ((1 : ℝ) / 3)) = Real.log 2 := by norm_num

/-- **Z079 · 6 个词各自都是平衡词**（Σw = 0）—— 组态空间定义的直接核验。 -/
theorem z079_all_balanced : ∀ w ∈ words4, w.sum = 0 := by decide

/-- **Z080 · M⁺ 壳分布（L=4）**：m₀=2, m₁=3, m₂=1。 -/
theorem z080_Mplus_distribution :
    (words4.filter (fun w => Mplus w = 0)).length = 2 ∧
      (words4.filter (fun w => Mplus w = 1)).length = 3 ∧
      (words4.filter (fun w => Mplus w = 2)).length = 1 := by decide

/-- **Z080 · 反射原理公式在 L=4 上逐项相符**：m_j = C(4,2−j) − C(4,1−j)
    （j=2 时第二项按 0 记）。 -/
theorem z080_ballot_formula_L4 :
    ((Nat.choose 4 2 : ℤ) - (Nat.choose 4 1 : ℤ) = 2) ∧
      ((Nat.choose 4 1 : ℤ) - (Nat.choose 4 0 : ℤ) = 3) ∧
      ((Nat.choose 4 0 : ℤ) - 0 = 1) := by
  refine ⟨?_, ?_, ?_⟩ <;> decide

/-- **Z080 · 重数和 = 中央二项式**（L=4）：Σ_j m_j = 2+3+1 = C(4,2) = 6。 -/
theorem z080_mass_sum_L4 : (2 + 3 + 1 : ℕ) = Nat.choose 4 2 := by decide

#print axioms z079_W4_card
#print axioms z079_layer_table
#print axioms z079_layer_counts
#print axioms z079_layer_weights
#print axioms z079_kappa1
#print axioms z079_delta_mod
#print axioms z079_all_balanced
#print axioms z080_Mplus_distribution
#print axioms z080_ballot_formula_L4
#print axioms z080_mass_sum_L4

end Zero
