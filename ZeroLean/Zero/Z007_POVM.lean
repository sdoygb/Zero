/-
  Z007（概率与测量 · 效应与退相干）的机器检查版本 —— POVM 与无信号条件。

  Z007 的两条 F 类：
    · **POVM 归一 ⟹ 概率归一**：`Σᵢ Mᵢ = 1` 时 `Σᵢ tr(Mᵢρ) = tr ρ`（Born 规则的线性内核）；
    · **无信号条件 ⟹ 迹保持**：Kraus 算子满足 `Σᵢ Kᵢ†Kᵢ = 1` 时 `Σᵢ tr(KᵢρKᵢ†) = tr ρ`
      —— 这就是"信道不丢概率"（无信号）的矩阵形式，也是 Z116 单粒子劈裂的前提。

  两条都是**有限维线性代数**：只用迹的线性、循环性（`tr(AB) = tr(BA)`）与求和分配律。

  诚实边界：Z007 的另一条（"干涉交叉项需非中心态或显式类间相位条件"）是**具名输入**
  （类间相位条件 𝒞_S），不在本篇证明范围。
-/

import Mathlib

namespace Zero

open Matrix

/-- **Z007 · POVM 归一 ⟹ 概率归一**：若 `Σᵢ Mᵢ = 1`，则对任意 ρ 有 `Σᵢ tr(Mᵢρ) = tr ρ`。
    （Born 规则的线性内核：测量结果概率之和恒为 1，与 ρ 无关。） -/
theorem z007_povm_probability_sum {n m : ℕ} (M : Fin m → Matrix (Fin n) (Fin n) ℂ)
    (h : ∑ i, M i = 1) (ρ : Matrix (Fin n) (Fin n) ℂ) :
    ∑ i, (M i * ρ).trace = ρ.trace := by
  rw [← Matrix.trace_sum, ← Finset.sum_mul, h, Matrix.one_mul]

/-- **Z007 · 无信号条件 ⟹ 迹保持**：Kraus 算子 `{Kᵢ}` 满足 `Σᵢ Kᵢ†Kᵢ = 1` 时，
    演化后的总概率守恒：`Σᵢ tr(Kᵢ ρ Kᵢ†) = tr ρ`。
    证明：迹的循环性把 `tr(KρK†)` 换成 `tr(K†K ρ)`，再用求和分配与归一条件。 -/
theorem z007_no_signal_trace {n m : ℕ} (K : Fin m → Matrix (Fin n) (Fin n) ℂ)
    (h : ∑ i, (K i).conjTranspose * K i = 1) (ρ : Matrix (Fin n) (Fin n) ℂ) :
    ∑ i, (K i * ρ * (K i).conjTranspose).trace = ρ.trace := by
  -- 逐项：tr(KρK†) = tr((K†K)ρ)
  have hterm : ∀ i : Fin m,
      (K i * ρ * (K i).conjTranspose).trace = (((K i).conjTranspose * K i) * ρ).trace := by
    intro i
    rw [Matrix.trace_mul_comm (K i * ρ) ((K i).conjTranspose), ← Matrix.mul_assoc]
  calc ∑ i, (K i * ρ * (K i).conjTranspose).trace
      = ∑ i, (((K i).conjTranspose * K i) * ρ).trace := Finset.sum_congr rfl fun i _ => hterm i
    _ = (∑ i, ((K i).conjTranspose * K i) * ρ).trace := by rw [Matrix.trace_sum]
    _ = ((∑ i, (K i).conjTranspose * K i) * ρ).trace := by rw [Finset.sum_mul]
    _ = ((1 : Matrix (Fin n) (Fin n) ℂ) * ρ).trace := by rw [h]
    _ = ρ.trace := by rw [Matrix.one_mul]

#print axioms z007_povm_probability_sum
#print axioms z007_no_signal_trace

end Zero
