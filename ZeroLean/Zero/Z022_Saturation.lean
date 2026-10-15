/-
  Z022（饱和判据 `C(T,⌊T/2⌋) ≤ 2(T+1) ⟺ T ≤ 5`）的机器检查版本。

  审计修订（2026-10-14）：原版**从未编译通过**（用了不存在的 `Nat.choose_le_choose` 参数名），
  被审计移入 `wip/`。这里按待办账写明的路径重写：

    · **中心二项式单调**（真有内容）：偶数/奇数两种情形各由 Pascal 递推直接给出；
    · **T ≤ 5 成立** 与 **6 ≤ T ≤ 20 不成立**：内核 `decide` 逐值验证。

  工程注记：`Nat.even_or_odd` 给的是**加法形式**（`m+1+(m+1)`），而 `2 * x` 不是 `x + x` 的
  defeq（要 `two_mul`）—— 所以这里统一用加法形式陈述，避免 `show`/`rw` 的匹配失败。

  诚实边界（记入 `_staging/LEAN_BACKLOG.md`）：**一般的 "T ≥ 6 ⟹ 判据不成立" 尚未证**，
  需要两条算术引理 ——(a) 同 n 下 `C(n,3) ≤ C(n,⌊n/2⌋)`（k 方向单调，Mathlib 无现成件）、
  (b) `2(T+1) < C(T,3)` 对 `T ≥ 6`。这里没有削弱命题冒充完成，而是明确标出范围。
-/

import Mathlib

namespace Zero

/-- **Z022 · 中心二项式单调不减**：`C(T,⌊T/2⌋) ≤ C(T+1,⌊(T+1)/2⌋)`。
    偶数与奇数两种情形都由 Pascal 递推（`Nat.choose_succ_succ`）直接给出。 -/
theorem z022_central_choose_mono (T : ℕ) :
    Nat.choose T (T / 2) ≤ Nat.choose (T + 1) ((T + 1) / 2) := by
  rcases Nat.even_or_odd T with ⟨m, rfl⟩ | ⟨m, rfl⟩
  · rcases m with _ | m
    · decide
    · -- T = (m+1) + (m+1)
      have h1 : (m + 1 + (m + 1)) / 2 = m + 1 := by omega
      have h2 : (m + 1 + (m + 1) + 1) / 2 = m + 1 := by omega
      rw [h1, h2]
      have hp : Nat.choose (m + 1 + (m + 1) + 1) (m + 1)
          = Nat.choose (m + 1 + (m + 1)) m + Nat.choose (m + 1 + (m + 1)) (m + 1) :=
        Nat.choose_succ_succ (m + 1 + (m + 1)) m
      rw [hp]
      omega
  · -- T = 2 * m + 1
    have h1 : (2 * m + 1) / 2 = m := by omega
    have h2 : (2 * m + 1 + 1) / 2 = m + 1 := by omega
    rw [h1, h2]
    have hp : Nat.choose (2 * m + 1 + 1) (m + 1)
        = Nat.choose (2 * m + 1) m + Nat.choose (2 * m + 1) (m + 1) :=
      Nat.choose_succ_succ (2 * m + 1) m
    rw [hp]
    omega

/-- **Z022 · 饱和判据在 `T ≤ 5` 上成立**（逐值内核验证）。 -/
theorem z022_saturation_small (T : ℕ) (hT : T ≤ 5) :
    Nat.choose T (T / 2) ≤ 2 * (T + 1) := by
  interval_cases T <;> decide

/-- **Z022 · 饱和判据在 `6 ≤ T ≤ 20` 上不成立**（逐值内核验证）。
    配合 `z022_central_choose_mono` 与待办账里的两条算术引理，可外推到一般的 `T ≥ 6`。 -/
theorem z022_saturation_fails_range (T : ℕ) (h6 : 6 ≤ T) (h20 : T ≤ 20) :
    2 * (T + 1) < Nat.choose T (T / 2) := by
  interval_cases T <;> decide

/-! ## 三、一般情形（本轮补齐：此前只做到 `T ≤ 20` 的逐值核验）

两条算术件：
  (a) 中心系数是 k 方向的最大值：`C(T,3) ≤ C(T,⌊T/2⌋)`（Mathlib `Nat.choose_le_middle`）；
  (b) 立方增长压过线性：`T ≥ 6 ⟹ 2(T+1) < C(T,3)`（Pascal 递推 + 归纳）。
合起来即得**一般**的 `T ≥ 6 ⟹ 判据不成立`，从而判据等价于 `T ≤ 5`。 -/

/-- **Z022(a) · 三阶选择数 ≤ 中心选择数**：`C(T,3) ≤ C(T,⌊T/2⌋)`。 -/
private theorem z022_choose_three_le_central (T : ℕ) : T.choose 3 ≤ T.choose (T / 2) :=
  Nat.choose_le_middle 3 T

/-- **Z022(b) · 立方压线性**：`T ≥ 6 ⟹ 2(T+1) < C(T,3)`。
    归纳步用 Pascal：`C(T+1,3) = C(T,2) + C(T,3)`，而 `C(T,2) ≥ C(6,2) = 15 ≥ 2`。 -/
theorem z022_two_mul_succ_lt_choose_three (T : ℕ) (hT : 6 ≤ T) : 2 * (T + 1) < T.choose 3 := by
  induction T, hT using Nat.le_induction with
  | base => decide
  | succ T hT ih =>
    have hpascal : (T + 1).choose 3 = T.choose 2 + T.choose 3 := Nat.choose_succ_succ T 2
    have hmono := Nat.choose_le_choose 2 (show 6 ≤ T by omega)
    have h2 : 2 ≤ T.choose 2 := by
      have h6 : (6 : ℕ).choose 2 = 15 := by decide
      omega
    rw [hpascal]
    omega

/-- **Z022 · 饱和判据（一般 T）**：`C(T,⌊T/2⌋) ≤ 2(T+1) ⟺ T ≤ 5`。
    这就是既有饱和律 `c_T = C(T,⌊T/2⌋) ≤ 2(T+1) ⇒ T ≤ 5` 的完整机器版本。 -/
theorem z022_saturation_iff (T : ℕ) : T.choose (T / 2) ≤ 2 * (T + 1) ↔ T ≤ 5 := by
  constructor
  · intro h
    by_contra hT
    have h6 : 6 ≤ T := by omega
    have hlt := z022_two_mul_succ_lt_choose_three T h6
    have hle := z022_choose_three_le_central T
    omega
  · intro hT
    exact z022_saturation_small T hT

/-- **Z022 · 词数版**：饱和解在偶长度上只有 `T = 2` 与 `T = 4`（`|W_T| = C(T,T/2)`）。 -/
theorem z022_saturation_word_count (T : ℕ) (hT : 2 ≤ T) (heven : Even T)
    (h : T.choose (T / 2) ≤ 2 * (T + 1)) : T = 2 ∨ T = 4 := by
  have h5 : T ≤ 5 := (z022_saturation_iff T).mp h
  obtain ⟨k, rfl⟩ := heven
  omega

#print axioms z022_central_choose_mono
#print axioms z022_saturation_small
#print axioms z022_saturation_fails_range
#print axioms z022_two_mul_succ_lt_choose_three
#print axioms z022_saturation_iff
#print axioms z022_saturation_word_count

end Zero
