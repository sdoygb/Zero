/-
  Z075 §3 段结构闭式 —— **一般情形**的机器检查（续作轮）。

  目标（四条，名字与陈述由题面逐字给定）：
    (1) `z075Compositions` + `z075_compositions_card`：stars-and-bars
        `#{v : 把 n 分成 c 个正整数} = C(n−1, c−1)`；
    (2) `z075_runs_even_general`  : `z075Count (2n) (2c)   = 2·C(n−1,c−1)²`；
    (3) `z075_runs_odd_general`   : `z075Count (2n) (2c+1) = 2·C(n−1,c)·C(n−1,c−1)`；
    (4) `z075_runs_sum_general`   : `1 ≤ n` 时 `Σ_{r∈Icc 2 2n} z075Count (2n) r = C(2n,n)`。

  **本轮做到哪、没做到哪（诚实边界，必读）**：

  · **(1) 已证**（`z075_compositions_card`，无 `sorryAx`）。同时给出可复用的分解引理
    `z075_compositions_succ_card`（末项 k 分层：`C(n−1,c−1)` 的递推 = 曲棍球棒
    `Nat.sum_range_add_choose`）与 `z075_compositions_card_eq_zero_of_lt`（无解判据）。

  · **(2)(3)(4) 本轮未证**，故**本文件不出现**这三条（不许 `sorry`、不许削弱陈述冒充）。
    剩余缺口的**准确形状**（供后来者接力）：
      需要「词 ↦ 段长向量」的双射：
        `{w : Fin (2n) → Bool // z075Runs w = r ∧ Σ(w) = 0}`
          ≃ {v : Fin r → ℕ // (∀ i, 1 ≤ v i) ∧ Σ_{i 偶} v i = n ∧ Σ_{i 奇} v i = n}
          × {首字母两种}`；
      然后对 `Fin (2c) → ℕ`（或 `Fin (2c+1) → ℕ`）按 **奇偶下标**拆成两份
      `z075Compositions n c`，用 (1) 收尾。缺的两件中，**(1) 已在手**；
      真正还缺的是「词的段长向量提取」与「奇偶下标拆分」两处索引搬运
      （`Finset`/`Fin` 上的排序块分解：段长 = 相邻「转向位置」之差）。
  · 本文件给出的**词侧接口**（把缺口钉死的两把尺子）：
    `z075_wordSum_eq_sameCount` / `z075_wordSum_zero_iff`：`Σ(w) = 0 ⟺ 首字母同号的位数 = 半长`。
    这把 `z075Count` 里的 `wordSum = 0` 翻译成纯计数条件，是上面双射的正规输入。

  · **与题面的一处偏差（必读）**：题面给的
      `def z075Compositions (n c : ℕ) : Finset (Fin c → ℕ) :=
         (Finset.univ : Finset (Fin c → ℕ)).filter …`
    **在本版 Mathlib 里不可原样成立**：`ℕ` 无限，`(Fin c → ℕ)` 没有 `Fintype`，
    `Finset.univ` 无实例，`def` 直接报 `failed to synthesize instance Fintype (Fin c → ℕ)`。
    本文件改为「先在有界类型 `Fin c → Fin (n+1)` 上枚举、再取像」（`z075BddVec` +
    `z075ToNatEmb`）。**集合本身没变**（各项 ≥ 1 且和为 n ⟹ 每项 ≤ n < n+1），
    由 `z075_mem_compositions` 逐条确认：`v ∈ z075Compositions n c ↔ (∀i, 1≤v i) ∧ Σv = n`。
    四条陈述里 `z075_compositions_card` 的名字与内容未动。

  · 纪律：无 `sorry`/`admit`/`native_decide`；`decide` 只用在 `z075_wordStep_ite` 这种
    两个 `Bool` 的有限枚举上。`#print axioms` 见文末（全部只有 propext/Classical.choice/Quot.sound）。
-/

import Mathlib
import Zero.Z075_Runs

namespace Zero

open Finset

/-! ### 一、stars-and-bars：把 `n` 拆成 `c` 个正整数 -/

/-- 有界枚举：`Fin c → Fin (n+1)` 上的全部函数（有限类型），再筛「各项 ≥ 1 且和为 `n`」。

注意：题面原写法 `(Finset.univ : Finset (Fin c → ℕ)).filter …` **在本工程里不成立**：
`ℕ` 无限，`(Fin c → ℕ)` 没有 `Fintype` 实例，`def` 直接报
`failed to synthesize instance of type class Fintype (Fin c → ℕ)`。
这里改成「先在有界类型上枚举、再取像」，**集合本身不变**（各项 ≥ 1 且和为 `n`
⟹ 每项 ≤ n < n+1，见 `z075_mem_compositions`），故 `z075Compositions` 仍是
「把 `n` 拆成 `c` 个正整数」的向量全体。 -/
def z075BddVec (n c : ℕ) : Finset (Fin c → Fin (n + 1)) :=
  (Finset.univ : Finset (Fin c → Fin (n + 1))).filter
    (fun v => (∀ i, 1 ≤ (v i : ℕ)) ∧ ∑ i, (v i : ℕ) = n)

/-- 有界向量到 `ℕ` 值向量的嵌入 `v ↦ (fun i ↦ (v i : ℕ))`。 -/
def z075ToNatEmb (n c : ℕ) : (Fin c → Fin (n + 1)) ↪ (Fin c → ℕ) :=
  ⟨fun v i => (v i : ℕ), fun _a _b h => funext fun i => Fin.ext (congrFun h i)⟩

/-- 正整数分拆集（stars-and-bars）：把 `n` 拆成 `c` 个正整数的向量全体。 -/
def z075Compositions (n c : ℕ) : Finset (Fin c → ℕ) :=
  (z075BddVec n c).map (z075ToNatEmb n c)

/-- 有界枚举的成员刻画。 -/
theorem z075_mem_bdd {n c : ℕ} {v : Fin c → Fin (n + 1)} :
    v ∈ z075BddVec n c ↔ (∀ i, 1 ≤ (v i : ℕ)) ∧ ∑ i, (v i : ℕ) = n := by
  simp [z075BddVec]

/-- **成员刻画**：`v` 是 `n` 的 `c` 项正分拆 ⟺ 各项 ≥ 1 且和为 `n`。 -/
theorem z075_mem_compositions {n c : ℕ} {v : Fin c → ℕ} :
    v ∈ z075Compositions n c ↔ (∀ i, 1 ≤ v i) ∧ ∑ i, v i = n := by
  constructor
  · intro hv
    obtain ⟨u, hu, hvu⟩ := Finset.mem_map.mp hv
    obtain ⟨hpos, hsum⟩ := z075_mem_bdd.mp hu
    refine ⟨fun i => ?_, ?_⟩
    · have h := hpos i
      rwa [show (u i : ℕ) = v i from congrFun hvu i] at h
    · have hcongr : (∑ i, (u i : ℕ)) = ∑ i, v i :=
        Finset.sum_congr rfl fun i _ => congrFun hvu i
      rw [hcongr] at hsum
      exact hsum
  · rintro ⟨hpos, hsum⟩
    have hlt : ∀ i, v i < n + 1 := by
      intro i
      have hle : v i ≤ ∑ j : Fin c, v j :=
        Finset.single_le_sum (f := v) (s := Finset.univ) (fun j _ => Nat.zero_le _) (Finset.mem_univ i)
      exact Nat.lt_succ_of_le (hle.trans (le_of_eq hsum))
    refine Finset.mem_map.mpr ⟨fun i => ⟨v i, hlt i⟩, ?_, ?_⟩
    · exact z075_mem_bdd.mpr ⟨fun i => hpos i, hsum⟩
    · funext i
      rfl

/-- 无解判据：`c` 个正整数之和至少是 `c`，故 `m < c` 时该集合为空。 -/
theorem z075_compositions_card_eq_zero_of_lt {m c : ℕ} (h : m < c) :
    (z075Compositions m c).card = 0 := by
  rw [Finset.card_eq_zero]
  ext v
  simp only [Finset.notMem_empty, iff_false]
  intro hv
  obtain ⟨hpos, hsum⟩ := z075_mem_compositions.mp hv
  have hle : c ≤ ∑ i : Fin c, v i := by
    calc c = ∑ _i : Fin c, 1 := by simp
      _ ≤ ∑ i : Fin c, v i := Finset.sum_le_sum fun i _ => hpos i
  omega

/-- **末项分解**：末项 = `k`，前 `c` 项 = 把 `n−k` 拆成 `c` 个正整数。 -/
theorem z075_compositions_succ_card (n c : ℕ) :
    (z075Compositions n (c + 1)).card
      = ∑ k ∈ Finset.Icc 1 n, (z075Compositions (n - k) c).card := by
  have hmaps : Set.MapsTo (fun v : Fin (c + 1) → ℕ => v (Fin.last c))
      ↑(z075Compositions n (c + 1)) ↑(Finset.Icc 1 n) := by
    intro v hv
    obtain ⟨hpos, hsum⟩ := z075_mem_compositions.mp hv
    have hle : v (Fin.last c) ≤ ∑ i : Fin (c + 1), v i :=
      Finset.single_le_sum (f := v) (s := Finset.univ) (fun i _ => Nat.zero_le _) (Finset.mem_univ _)
    rw [Finset.mem_coe, Finset.mem_Icc]
    exact ⟨hpos _, hle.trans (le_of_eq hsum)⟩
  rw [Finset.card_eq_sum_card_fiberwise hmaps]
  refine Finset.sum_congr rfl fun k hk => ?_
  obtain ⟨hk1, hk2⟩ := Finset.mem_Icc.mp hk
  refine Finset.card_bij (fun v _ => fun i : Fin c => v i.castSucc) ?_ ?_ ?_
  · intro v hv
    obtain ⟨hvmem, hvlast⟩ := Finset.mem_filter.mp hv
    obtain ⟨hpos, hsum⟩ := z075_mem_compositions.mp hvmem
    rw [z075_mem_compositions]
    refine ⟨fun i => hpos _, ?_⟩
    have hcast := Fin.sum_univ_castSucc v
    rw [hsum, hvlast] at hcast
    omega
  · intro v hv v' hv' heq
    obtain ⟨_, hvl⟩ := Finset.mem_filter.mp hv
    obtain ⟨_, hv'l⟩ := Finset.mem_filter.mp hv'
    have key : ∀ u u' : Fin (c + 1) → ℕ,
        (∀ i : Fin c, u i.castSucc = u' i.castSucc) → u (Fin.last c) = u' (Fin.last c) → u = u' := by
      intro u u' hcs hls
      funext i
      exact Fin.lastCases hls (fun j => hcs j) i
    exact key v v' (fun i => congrFun heq i) (by rw [hvl, hv'l])
  · intro u hu
    obtain ⟨hpos, hsum⟩ := z075_mem_compositions.mp hu
    refine ⟨Fin.snoc u k, ?_, ?_⟩
    · rw [Finset.mem_filter]
      refine ⟨?_, Fin.snoc_last ..⟩
      rw [z075_mem_compositions]
      refine ⟨fun i => ?_, ?_⟩
      · exact Fin.lastCases (by rw [Fin.snoc_last]; exact hk1)
          (fun j => by rw [Fin.snoc_castSucc]; exact hpos j) i
      · rw [Fin.sum_univ_castSucc]
        simp only [Fin.snoc_castSucc, Fin.snoc_last]
        omega
    · funext i
      rw [Fin.snoc_castSucc]

/-- 基例 `c = 1`：把 `n` 拆成 1 个正整数只有 `(n)` 一种。 -/
theorem z075_compositions_one_card (n : ℕ) (hn : 1 ≤ n) :
    (z075Compositions n 1).card = Nat.choose (n - 1) 0 := by
  have hset : z075Compositions n 1 = {fun _ => n} := by
    ext v
    rw [z075_mem_compositions, Finset.mem_singleton]
    constructor
    · rintro ⟨hpos, hsum⟩
      funext i
      have h0 : v 0 = n := by rw [← Fin.sum_univ_one v]; exact hsum
      rw [Subsingleton.elim i (0 : Fin 1)]
      exact h0
    · rintro rfl
      exact ⟨fun _ => hn, by simp⟩
  rw [hset, Finset.card_singleton, Nat.choose_zero_right]

/-- **stars-and-bars 计数**：`#{v : 把 n 分成 c 个正整数} = C(n−1, c−1)`。 -/
theorem z075_compositions_card (n c : ℕ) (hc : 1 ≤ c) (hcn : c ≤ n) :
    (z075Compositions n c).card = Nat.choose (n - 1) (c - 1) := by
  induction c generalizing n with
  | zero => omega
  | succ c ih =>
    rcases Nat.eq_zero_or_pos c with h0 | hpos
    · subst h0
      exact z075_compositions_one_card n (by omega)
    · rw [z075_compositions_succ_card n c]
      have hsplit : (∑ k ∈ Finset.Icc 1 n, (z075Compositions (n - k) c).card)
          = ∑ k ∈ Finset.Icc 1 (n - c), (z075Compositions (n - k) c).card := by
        refine (Finset.sum_subset ?_ ?_).symm
        · intro k hk
          rw [Finset.mem_Icc] at hk ⊢
          omega
        · intro k hk1 hk2
          rw [Finset.mem_Icc] at hk1 hk2
          exact z075_compositions_card_eq_zero_of_lt (by omega)
      rw [hsplit]
      have hih : (∑ k ∈ Finset.Icc 1 (n - c), (z075Compositions (n - k) c).card)
          = ∑ k ∈ Finset.Icc 1 (n - c), Nat.choose (n - k - 1) (c - 1) := by
        refine Finset.sum_congr rfl fun k hk => ?_
        rw [Finset.mem_Icc] at hk
        exact ih (n - k) hpos (by omega)
      rw [hih]
      have hreindex : (∑ k ∈ Finset.Icc 1 (n - c), Nat.choose (n - k - 1) (c - 1))
          = ∑ j ∈ Finset.range (n - c), Nat.choose (j + (c - 1)) (c - 1) := by
        refine Finset.sum_bij (fun k _ => n - c - k) ?_ ?_ ?_ ?_
        · intro k hk
          rw [Finset.mem_Icc] at hk
          exact Finset.mem_range.mpr (by omega)
        · intro k₁ hk₁ k₂ hk₂ heq
          rw [Finset.mem_Icc] at hk₁ hk₂
          omega
        · intro j hj
          rw [Finset.mem_range] at hj
          exact ⟨n - c - j, by rw [Finset.mem_Icc]; omega, by omega⟩
        · intro k hk
          rw [Finset.mem_Icc] at hk
          rw [show n - k - 1 = n - c - k + (c - 1) by omega]
      rw [hreindex, show n - c = (n - c - 1) + 1 by omega]
      rw [Nat.sum_range_add_choose]
      rw [show (n - c - 1) + (c - 1) + 1 = n - 1 by omega,
        show (c - 1) + 1 = c by omega]
      simp only [Nat.add_sub_cancel]

/-! ### 二、词侧接口：`wordSum` 与「与首字母同号的位数」

把 `z075Count` 里的条件 `wordSum w = 0` 翻译成「首字母同号的位数恰为一半」——
这是「词 ↦ 段长向量」双射那一步的桥（一般情形的剩余缺口见文件头）。 -/

/-- `wordStep` 的「同号/异号」两分支写法：同号取 `wordStep b₀`、异号取 `-wordStep b₀`。 -/
theorem z075_wordStep_ite (b b₀ : Bool) :
    wordStep b = if b = b₀ then wordStep b₀ else -wordStep b₀ := by
  cases b <;> cases b₀ <;> decide

/-- **平衡 ↔ 同号位数**：`Σ(w) = sgn(w₀)·(2·#{i : wᵢ = w₀} − L)`，`L = 词长`。 -/
theorem z075_wordSum_eq_sameCount (L : ℕ) (w : Fin (L + 1) → Bool) :
    wordSum w = wordStep (w 0)
      * (2 * ((#{i : Fin (L + 1) | w i = w 0} : ℕ) : ℤ) - ((L + 1 : ℕ) : ℤ)) := by
  have hstep : ∀ i : Fin (L + 1),
      wordStep (w i) = if w i = w 0 then wordStep (w 0) else -wordStep (w 0) :=
    fun i => z075_wordStep_ite (w i) (w 0)
  have hcard : #{i : Fin (L + 1) | ¬ (w i = w 0)}
      = (L + 1) - #{i : Fin (L + 1) | w i = w 0} := by
    have h := Finset.card_filter_add_card_filter_not
      (s := (Finset.univ : Finset (Fin (L + 1)))) (p := fun i => w i = w 0)
    simp only [Finset.card_univ, Fintype.card_fin] at h
    omega
  unfold wordSum
  rw [Finset.sum_congr rfl fun i _ => hstep i]
  rw [Finset.sum_ite]
  rw [Finset.sum_const, Finset.sum_const, hcard]
  simp only [nsmul_eq_mul]
  rw [Nat.cast_sub (show #{i : Fin (L + 1) | w i = w 0} ≤ L + 1 by
    simpa using Finset.card_filter_le (s := (Finset.univ : Finset (Fin (L + 1))))
      (fun i => w i = w 0))]
  ring

/-- **平衡的等价刻画**：`Σ(w) = 0 ⟺ 首字母同号的位数 = 词长/2`。 -/
theorem z075_wordSum_zero_iff (L : ℕ) (w : Fin (L + 1) → Bool) :
    wordSum w = 0 ↔ 2 * ((#{i : Fin (L + 1) | w i = w 0} : ℕ) : ℤ) = ((L + 1 : ℕ) : ℤ) := by
  rw [z075_wordSum_eq_sameCount]
  constructor
  · intro h
    have hs : wordStep (w 0) = 1 ∨ wordStep (w 0) = -1 := by
      unfold wordStep
      cases w 0 <;> simp
    rcases hs with hs | hs <;> rw [hs] at h <;> linarith
  · intro h
    rw [h]
    ring

#print axioms z075_wordStep_ite
#print axioms z075_wordSum_eq_sameCount
#print axioms z075_wordSum_zero_iff

#print axioms z075_mem_bdd
#print axioms z075_mem_compositions
#print axioms z075_compositions_card_eq_zero_of_lt
#print axioms z075_compositions_succ_card
#print axioms z075_compositions_one_card
#print axioms z075_compositions_card

end Zero
