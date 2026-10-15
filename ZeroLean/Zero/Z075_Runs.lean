/-
  Z075（段结构：粒子数与种类）§3 闭式的机器检查版本。

  Z075 §3 的闭式（篇内枚举核到 L ≤ 16，**只做枚举核对、未给解析证明**）：

      b_r := #{w ∈ W_L : runs(w) = r}
           = 2·C(n−1, c−1)²            ,  r = 2c   （**偶数段**）
           = 2·C(n−1, c)·C(n−1, c−1)   ,  r = 2c+1 （**奇数段**）
      （L = 2n，C = `Nat.choose`；合计 Σ_r b_r = C(2n, n)。）

  段结构把词分成 **重子型**（`r = 2`：一段全 +、一段全 −）与 **介子型**
  （初稿口径 `r ≥ 4`；Z075 §3.1a 更正为 `r ≥ 3` —— `r = 3` 确实出现，如 L = 4 的 `(+−−+)`）。

  **本轮做到哪、没做到哪（诚实边界，必读）**：
    · 口径：`z075Runs` = 线性极大常值段数（`1 + #转向位置`）；`z075Count` = 恰有 `r` 段的
      **平衡词**（`wordSum = 0`，即零和词全体）个数。题面原写法里的 `wordShell L 0` 是
      "和 0 且层高 0"（几乎只是空词），与"平衡词全体"不是一回事，故按题面授权改为
      `wordSum = 0`（名字与另两条陈述未动）。
    · **一般情形（任意 n、任意 c）本轮未证**。卡点是把
      "词 ↦（正段长度分拆，负段长度分拆）" 这个双射形式化所需的 stars-and-bars：
      把 `n` 拆成 `c` 个正整数共 `C(n−1, c−1)` 种。本版 Mathlib **没有** `Nat.compositions`
      之类的现成件，自建需对 `c` 归纳并用曲棍球棒恒等式（`Nat.sum_range_choose`）。
    · 本轮给出 **L ≤ 10 的全部 r**、以及 **L = 12 的抽样（偶段 c = 3、奇段 c = 3）**
      的机器核对（`decide`，内核归约，不用 `native_decide`）。
    · Z075 §3 的**合计式**按题面陈述（`∀ n : ℕ`，无 `1 ≤ n`）在原样下**为假**：
      `n = 0` 时左端是空和 `0`、右端 `Nat.choose 0 0 = 1`。下面给出机器核验的反例；
      按"不改弱、不改名"的要求，本条**不**给出改名后的版本。
-/

import Mathlib
import Zero.Z047_Involution

namespace Zero

/-- 段数（线性极大常值段数）：`1 + #{转向位置}`。

与题面写法 `1 + (Finset.range (L-1)).filter (fun i => w i ≠ w (i+1)) |>.card` 在 `L ≥ 1` 上一致；
这里用 `Fin m` 的 `castSucc`/`succ` 索引（题面那种 `by omega` 写法在 `def` 里没有
`i < L-1` 的假设可依，不通过）。`L = 0`（空词）取 `0`：空词无段，
且不影响任何一条定理（`n = 0` 时假设自相矛盾，见下）。 -/
def z075Runs {L : ℕ} (w : Fin L → Bool) : ℕ :=
  match L, w with
  | 0, _ => 0
  | m + 1, w => 1 + ((Finset.univ : Finset (Fin m)).filter
      (fun i => w i.castSucc ≠ w i.succ)).card

/-- 恰有 `r` 段的平衡词个数。平衡 = `wordSum = 0`（零和 `±1` 词全体）；`true ↦ +1`、`false ↦ −1`
（Z047 `wordStep`）。 -/
def z075Count (L r : ℕ) : ℕ :=
  ((Finset.univ : Finset (Fin L → Bool)).filter (fun w => wordSum w = 0)).filter
    (fun w => z075Runs w = r) |>.card

/-! ### 一、合计式的边界反例（n = 0） -/

/-- Z075 §3 合计式在 `n = 0` 时的左端：`r` 取遍 `Icc 2 0 = ∅`，是空和 `0`。 -/
private theorem z075_runs_sum_zero_lhs :
    (∑ r ∈ Finset.Icc 2 (2 * 0), z075Count (2 * 0) r) = 0 := by decide

/-- Z075 §3 合计式在 `n = 0` 时的右端：`C(0, 0) = 1`。 -/
private theorem z075_runs_sum_zero_rhs : Nat.choose (2 * 0) 0 = 1 := by decide

/-- **边界反例**：合计式 `∑ r ∈ Icc 2 (2n), b_r = C(2n, n)` 在 `n = 0` 时**为假**
（左端 `0`、右端 `1`）。故题面那一条（无 `1 ≤ n`）不能成立，本轮不改弱、不改名，只报反例。 -/
private theorem z075_runs_sum_zero_counterexample :
    ¬ (∑ r ∈ Finset.Icc 2 (2 * 0), z075Count (2 * 0) r) = Nat.choose (2 * 0) 0 := by decide

/-! ### 二、闭式的机器核对（L ≤ 10 全部 r；L = 12 抽样）

`L = 12` 的核对要枚举 4096 个词，内核归约（`decide`）需要放开递归深度；
本段之后的命令用下面这个选项（只作用于本文件）。 -/

set_option maxRecDepth 100000

/-- **Z075 · 闭式（偶数段，重子型 r = 2 含于 c = 1 这一项）· L ≤ 10 机器核对**：
`b_{2c} = 2·C(n−1, c−1)²`，`L = 2n`。 -/
theorem z075_runs_even_smallL (n c : ℕ) (hn : n ≤ 5) (hc : 1 ≤ c) (hcn : c ≤ n) :
    z075Count (2 * n) (2 * c) = 2 * (Nat.choose (n - 1) (c - 1)) ^ 2 := by
  interval_cases n <;> interval_cases c <;> decide

/-- **Z075 · 闭式（奇数段）· L ≤ 10 机器核对**：`b_{2c+1} = 2·C(n−1, c)·C(n−1, c−1)`。 -/
theorem z075_runs_odd_smallL (n c : ℕ) (hn : n ≤ 5) (hc : 1 ≤ c) (hcn : c ≤ n) :
    z075Count (2 * n) (2 * c + 1) = 2 * Nat.choose (n - 1) c * Nat.choose (n - 1) (c - 1) := by
  interval_cases n <;> interval_cases c <;> decide

/-- **L = 12 抽样（偶数段，`c = 3`）**：`b_6 = 2·C(5, 2)² = 200`，`n = 6`。
（`L = 12` 的枚举要 4096 个词，内核归约需要放开递归深度，见上方文件级选项。） -/
theorem z075_runs_L12_even_c3 :
    z075Count 12 6 = 2 * (Nat.choose (6 - 1) (3 - 1)) ^ 2 := by decide

/-- **L = 12 抽样（奇数段，`c = 3`）**：`b_7 = 2·C(5, 3)·C(5, 2) = 200`，`n = 6`。 -/
theorem z075_runs_L12_odd_c3 :
    z075Count 12 7 = 2 * Nat.choose (6 - 1) 3 * Nat.choose (6 - 1) (3 - 1) := by decide

#print axioms z075_runs_even_smallL
#print axioms z075_runs_odd_smallL
#print axioms z075_runs_L12_even_c3
#print axioms z075_runs_L12_odd_c3

end Zero