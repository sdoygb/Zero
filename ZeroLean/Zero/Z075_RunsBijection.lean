/-
  Z075 §3 段结构闭式 —— 「词 ↦ 段长向量」的**索引搬运**（续作轮）。

  前一轮（`Z075_RunsGeneral.lean`）已备好 stars-and-bars：
    `z075Compositions n c` + `z075_mem_compositions` + `z075_compositions_card`（一般情形）。
  本轮把剩下的两处索引搬运中的**第一处**做掉，并给出它的三条结构内容：

    · `z075RunLens`：**段长向量的提取**（极大常值段长，按出现次序）；
    · `z075_runLens_sum`    ：段长之和 = 词长 `L`；
    · `z075_runLens_length` ：段长向量长度 = 段数 `z075Runs w`。

  实现口径：`z075RunLens w = z075RunLensList (List.ofFn w)`，扫描相邻「转向点」做块分解；
  「长度 = 段数」这一步把 `List` 扫描与 `z075Runs` 的 `Finset` 濾波对接
  （`z075TurnsAux_ofFn`：`Fin m` 上「相邻不同」的个数，靠 `Fin.sum_univ_succ` 与
  `Fin.castSucc_succ` 搬索引）。**无 `sorry`/`admit`/`native_decide`**；`#print axioms` 见文末。

  **本轮做到哪、没做到哪（诚实边界，必读）**：

  · 已证（无 `sorryAx`）：上面三条结构定理。
  · **未证**：`z075_runs_even_general` / `z075_runs_odd_general` / `z075_runs_sum_general`，
    以及题面 `z075_even_run_lengths`（见下）。它们要的是**计数**双射
    `{w : Fin (2n) → Bool // z075Runs w = r ∧ Σw = 0} ≃ (分拆对) × {首字母}`，
    本轮只搬动了「词 → 段长向量」这一半，**还没搬**「段长向量的奇偶下标拆分」
    （`finProdFinEquiv` / `Fin.sum_univ_add` 那一处），故三条计数式与合计式不写入本文件
    （不许 `sorry`、不许削弱陈述冒充）。`z075_even_run_lengths` 亦不写入：见下。
  · **关于 `z075_even_run_lengths`（必须留意）**：题面该条的结论
    `∃ v : Fin c → ℕ, (∀ i, 1 ≤ v i) ∧ ∑ i, v i = n` **量词里不出现 `w`**，
    而 `hc : 1 ≤ c`、`hcn : c ≤ n` 已足以显式构造 `v = (n−c+1, 1, …, 1)`，
    即该陈述**与词无关**（不是假的，是弱于本意）。本文件因此**不**把它当「偶数段结构」来交付：
    真正把词与段长挂钩的命题应当是「偶数下标段长之和 = n」（需要 `∑` 按符号加权 = `wordSum`），
    那一层本轮没做。留待续作的正是这一层与上面的计数双射。
-/

import Mathlib
import Zero.Z075_RunsGeneral

namespace Zero

open Finset

/-! ### 一、段长向量：扫描相邻转向点做块分解 -/

/-- 段长累加器：当前段符号 `b`、当前段已累计长度 `acc`，继续扫描剩余符号 `bs`；
读到不同符号就把 `acc` 记为一个段长、并开启新段（新段从 `1` 起算）。 -/
def z075RunLensAux : Bool → ℕ → List Bool → List ℕ
  | _, acc, [] => [acc]
  | b, acc, b' :: bs =>
      if b' = b then z075RunLensAux b (acc + 1) bs else acc :: z075RunLensAux b' 1 bs

/-- 列表版段长向量：空表无段；非空表从首符号起算。 -/
def z075RunLensList : List Bool → List ℕ
  | [] => []
  | b :: bs => z075RunLensAux b 1 bs

/-- **段长向量的提取**：词的极大常值段长（按出现次序）。 -/
def z075RunLens {L : ℕ} (w : Fin L → Bool) : List ℕ :=
  z075RunLensList (List.ofFn w)

/-! ### 二、段长之和 = 词长；段长向量长度 = 段数 -/

/-- 累加器版：段长之和 = 已累计长度 + 剩余符号数。 -/
theorem z075_runLensAux_sum (b : Bool) (acc : ℕ) (bs : List Bool) :
    (z075RunLensAux b acc bs).sum = acc + bs.length := by
  induction bs generalizing b acc with
  | nil => simp [z075RunLensAux]
  | cons b' bs ih =>
      simp only [z075RunLensAux]
      by_cases h : b' = b
      · rw [h, ite_eq_left rfl, ih b (acc + 1), List.length_cons]
        omega
      · rw [ite_eq_right h, List.sum_cons, ih b' 1, List.length_cons]
        omega

/-- 列表版：段长之和 = 列表长度。 -/
theorem z075_runLensList_sum (l : List Bool) : (z075RunLensList l).sum = l.length := by
  cases l with
  | nil => rfl
  | cons b bs =>
      simp only [z075RunLensList, z075_runLensAux_sum, List.length_cons]
      omega

/-- **Z075 · 段长向量之和 = 词长** -/
theorem z075_runLens_sum {L : ℕ} (w : Fin L → Bool) :
    (z075RunLens w).sum = L := by
  simp only [z075RunLens, z075_runLensList_sum, List.length_ofFn]

/-- 列表版转向计数：`z075TurnsAux b bs` 数出 `bs` 内部相邻符号不同之处，
外加首符号与「前一个符号」`b` 的比较。 -/
def z075TurnsAux : Bool → List Bool → ℕ
  | _, [] => 0
  | b, b' :: bs => (if b' = b then 0 else 1) + z075TurnsAux b' bs

/-- 列表版转向数（空表 0 个转向）。 -/
def z075TurnsList : List Bool → ℕ
  | [] => 0
  | b :: bs => z075TurnsAux b bs

/-- 累加器版：段长向量的长度 = 转向数 + 1。 -/
theorem z075_runLensAux_length (b : Bool) (acc : ℕ) (bs : List Bool) :
    (z075RunLensAux b acc bs).length = z075TurnsAux b bs + 1 := by
  induction bs generalizing b acc with
  | nil => simp [z075RunLensAux, z075TurnsAux]
  | cons b' bs ih =>
      simp only [z075RunLensAux, z075TurnsAux]
      by_cases h : b' = b
      · rw [h, ite_eq_left rfl, ite_eq_left rfl, ih b (acc + 1)]
        omega
      · rw [ite_eq_right h, List.length_cons, ih b' 1, ite_eq_right h]
        omega

/-- 列表版：`b :: bs` 的段长向量长度 = 转向数 + 1。 -/
theorem z075_runLensList_length_cons (b : Bool) (bs : List Bool) :
    (z075RunLensList (b :: bs)).length = z075TurnsList (b :: bs) + 1 := by
  simp only [z075RunLensList, z075TurnsList, z075_runLensAux_length]

/-! ### 三、段长为正（下一轮双射的输入：段长向量是正整数向量） -/

/-- 累加器版：段长向量里每一项都 `≥ 1`（`acc + 1 ≥ 1`，故无需额外假设）。

写法上把累计量取成 `acc + 1`，免得归纳时把一个依赖 `acc` 的假设一起搬来搬去。 -/
theorem z075_runLensAux_succ_pos (b : Bool) (acc : ℕ) (bs : List Bool) :
    ∀ x ∈ z075RunLensAux b (acc + 1) bs, 1 ≤ x := by
  induction bs generalizing b acc with
  | nil =>
      intro x hx
      simp only [z075RunLensAux, List.mem_singleton] at hx
      omega
  | cons b' bs ih =>
      intro x hx
      simp only [z075RunLensAux] at hx
      by_cases h : b' = b
      · rw [h, ite_eq_left rfl] at hx
        exact ih b (acc + 1) x hx
      · rw [ite_eq_right h] at hx
        rcases List.mem_cons.mp hx with hx0 | hx'
        · omega
        · exact ih b' 0 x hx'

/-! ### 四、`List` 扫描 ↔ `Fin` 索引扫描（`z075Runs` 的对接）

`z075Runs w = 1 + #{i : Fin m | w i.castSucc ≠ w i.succ}`（`L = m+1`）；
下面把左边那个 `List` 扫描的转向数搬到同一把尺子上。 -/

/-- `Finset.card_filter` 的 `univ` 版：基数写成指示函数之和（便于对 `Fin` 做归纳拆分）。 -/
private theorem z075_card_filter_eq_sum {m : ℕ} (p : Fin m → Prop) [DecidablePred p] :
    ((Finset.univ : Finset (Fin m)).filter p).card = ∑ i, (if p i then 1 else 0) := by
  rw [Finset.card_filter]

/-- **列表扫描 ↔ 有限索引扫描**：`z075TurnsAux` 在 `List.ofFn w` 上的计数，
等于 `Fin m` 上「相邻不同」的位置数（即 `z075Runs w − 1`）。

证明：对 `m` 归纳；把 `Fin (k+1)` 的滤子按 `Fin.sum_univ_succ` 拆出头一项，
再用 `Fin.castSucc_succ` 把 `(i.succ).castSucc` 搬成 `i.castSucc.succ`，与归纳假设对齐。 -/
theorem z075TurnsAux_ofFn (m : ℕ) (w : Fin (m + 1) → Bool) :
    z075TurnsAux (w 0) (List.ofFn (fun i : Fin m => w i.succ))
      = ((Finset.univ : Finset (Fin m)).filter (fun i => w i.castSucc ≠ w i.succ)).card := by
  revert w
  induction m with
  | zero =>
      intro w
      simp only [List.ofFn_zero, z075TurnsAux, Finset.univ_eq_empty, Finset.filter_empty,
        Finset.card_empty]
  | succ k ih =>
      intro w
      have hsplit : ((Finset.univ : Finset (Fin (k + 1))).filter
            (fun i => w i.castSucc ≠ w i.succ)).card
          = (if w 0 ≠ w 1 then 1 else 0)
            + ((Finset.univ : Finset (Fin k)).filter
                (fun i => w (i.castSucc.succ) ≠ w (i.succ.succ))).card := by
        have h1 : (∑ i : Fin (k + 1), (if w i.castSucc ≠ w i.succ then (1 : ℕ) else 0))
            = (if w (0 : Fin (k + 1)).castSucc ≠ w (0 : Fin (k + 1)).succ then 1 else 0)
              + ∑ i : Fin k, (if w (i.succ).castSucc ≠ w (i.succ).succ then (1 : ℕ) else 0) :=
          Fin.sum_univ_succ (fun i : Fin (k + 1) => (if w i.castSucc ≠ w i.succ then (1 : ℕ) else 0))
        rw [z075_card_filter_eq_sum (fun i : Fin (k + 1) => w i.castSucc ≠ w i.succ), h1]
        rw [z075_card_filter_eq_sum (fun i : Fin k => w (i.castSucc.succ) ≠ w (i.succ.succ))]
        have h2 : (∑ i : Fin k, (if w (i.succ).castSucc ≠ w (i.succ).succ then (1 : ℕ) else 0))
            = ∑ i : Fin k, (if w (i.castSucc.succ) ≠ w (i.succ.succ) then (1 : ℕ) else 0) :=
          Finset.sum_congr rfl fun i _ => by rw [Fin.castSucc_succ]
        rw [h2]
        have h3 : (if w (0 : Fin (k + 1)).castSucc ≠ w (0 : Fin (k + 1)).succ then (1 : ℕ) else 0)
            = (if w 0 ≠ w 1 then (1 : ℕ) else 0) := by
          have hc : (w (0 : Fin (k + 1)).castSucc ≠ w (0 : Fin (k + 1)).succ) ↔ (w 0 ≠ w 1) :=
            Iff.rfl
          by_cases h : w 0 ≠ w 1
          · exact (ite_eq_left (hc.mpr h)).trans (ite_eq_left h).symm
          · exact (ite_eq_right (fun hh => h (hc.mp hh))).trans (ite_eq_right h).symm
        rw [h3]
      have hih := ih (fun i : Fin (k + 1) => w i.succ)
      simp only [Fin.succ_zero_eq_one] at hih
      rw [hsplit]
      simp only [List.ofFn_succ, z075TurnsAux, Fin.succ_zero_eq_one]
      rw [hih]
      have hlead : (if w 1 = w 0 then 0 else 1) = (if w 0 ≠ w 1 then 1 else 0) := by
        by_cases h : w 0 = w 1
        · rw [ite_eq_left h.symm, ite_eq_right (fun hh => hh h)]
        · rw [ite_eq_right (fun hh : w 1 = w 0 => h hh.symm), ite_eq_left h]
      rw [hlead]

/-- **Z075 · 段长向量长度 = 段数** -/
theorem z075_runLens_length {L : ℕ} (w : Fin L → Bool) :
    (z075RunLens w).length = z075Runs w := by
  cases L with
  | zero =>
      simp only [z075RunLens, List.ofFn_zero, z075RunLensList, List.length_nil, z075Runs]
  | succ m =>
      simp only [z075RunLens, List.ofFn_succ, z075RunLensList, z075Runs]
      rw [z075_runLensAux_length, z075TurnsAux_ofFn]
      omega

/-! ### 五、词版段长为正 + 一条审计用的机器核验事实 -/

/-- 词版：段长向量里每一项都 `≥ 1`。 -/
theorem z075_runLens_pos {L : ℕ} (w : Fin L → Bool) :
    ∀ x ∈ z075RunLens w, 1 ≤ x := by
  cases L with
  | zero =>
      intro x hx
      simp only [z075RunLens, List.ofFn_zero, z075RunLensList, List.not_mem_nil] at hx
  | succ m =>
      intro x hx
      simp only [z075RunLens, List.ofFn_succ, z075RunLensList] at hx
      exact z075_runLensAux_succ_pos (w 0) 0 (List.ofFn (fun i : Fin m => w i.succ)) x hx

/-- **把 `n` 拆成 `c` 个正整数总是可能的**（显式构造 `(n−c+1, 1, …, 1)`）。

**审计用途（重要）**：题面 `z075_even_run_lengths` 的**结论**
`∃ v : Fin c → ℕ, (∀ i, 1 ≤ v i) ∧ ∑ i, v i = n` 就是这一条，
而它只需 `hc : 1 ≤ c`、`hcn : c ≤ n` —— **量词里不出现 `w`**。
故那条陈述的三个假设 `w, hsum, hr` 对结论没有贡献：它不是假，是**弱于本意**
（真正与词挂钩的「偶数段结构」应当是「偶数下标段长之和 = n」，
需要先把 `∑` 按符号加权搬成 `wordSum`，本轮没做）。 -/
theorem z075_exists_pos_composition (n c : ℕ) (hc : 1 ≤ c) (hcn : c ≤ n) :
    ∃ v : Fin c → ℕ, (∀ i, 1 ≤ v i) ∧ ∑ i, v i = n := by
  have h0 : Fin c := ⟨0, hc⟩
  refine ⟨fun i : Fin c => if i = h0 then n - c + 1 else 1, ?_, ?_⟩
  · intro i
    by_cases h : i = h0
    · dsimp only
      rw [ite_eq_left h]
      omega
    · dsimp only
      rw [ite_eq_right h]
  · rw [← Finset.add_sum_erase _ _ (Finset.mem_univ h0)]
    rw [ite_eq_left rfl]
    have hz : ∀ i ∈ (Finset.univ : Finset (Fin c)).erase h0,
        (if i = h0 then n - c + 1 else 1) = 1 := by
      intro i hi
      exact ite_eq_right (Finset.mem_erase.mp hi).1
    rw [Finset.sum_congr rfl hz]
    rw [Finset.sum_const, Finset.card_erase_of_mem (Finset.mem_univ h0),
      Finset.card_univ, Fintype.card_fin]
    simp only [smul_eq_mul, mul_one]
    omega

#print axioms z075_runLensAux_sum
#print axioms z075_runLensList_sum
#print axioms z075_runLens_sum
#print axioms z075_runLensAux_length
#print axioms z075_runLensList_length_cons
#print axioms z075_runLensAux_succ_pos
#print axioms z075_runLens_pos
#print axioms z075_exists_pos_composition
#print axioms z075TurnsAux_ofFn
#print axioms z075_runLens_length

end Zero