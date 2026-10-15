/-
  Z061 / Z021 §8：Heawood 图（Fano 平面的关联图）的**围长 = 6**。

  Z065_Heawood.lean 已分条核验了「无三角形」「无 C₄」「存在 6-环」「二部 ⟹ 无奇闭道」，
  但没有单条的围长定理（审计所指的缺口）。本篇把这几条合成一条：

    `z065HeawoodGraph.girth = 6`。

  证明路线：
  · 上界（girth ≤ 6）：把 Z065 的显式 6-环 p₀ — L₀ — p₁ — L₁ — p₂ — L₆ — p₀
    （顶点序列 0,7,1,8,2,13,0）写成显式道 `z061SixCycleWalk`，逐字段证它是 `IsCycle`
    （`edges.Nodup` 与 `support.tail.Nodup` 由 `decide` 在具体列表上算掉；`≠ nil` 由长度 6 ≠ 0），
    再用 `Walk.IsCycle.girth_le_length`。
  · 下界（6 ≤ girth）：`SimpleGraph.le_egirth` 把 `6 ≤ egirth` 化归为
    「任何 `IsCycle` 道 w 都有 6 ≤ w.length」。设 L = w.length：
    `IsCircuit.three_le_length` 给 L ≥ 3；把 w 的顶点序列 `w.getVert`（总函数 ℕ → Fin 14）喂给
    Z065 的 `z065_heawood_no_odd_closed`（二部 ⟹ 闭道长度必偶）得 L 为偶。
    若 L < 6，则 L = 4：此时 v₀, v₁, v₂, v₃ 互异 —— 因为 `IsCycle.support_nodup` 说
    `support.tail = [v₁, v₂, v₃, v₄]` 无重复，而 v₄ = v₀（见 `z061_getVert_succ_eq_tail`
    把 v_{m+1} 与 `support.tail` 的第 m 项对接）；于是 v₁ 与 v₃ 是 v₀ 与 v₂ 的两个互异公共邻居，
    与 `z065_heawood_no_four_cycle`（不同顶点至多一个公共邻居）矛盾。
  · 合成：`le_antisymm`。

  注：本篇只用 `decide`（内核归约，不动用外部求值器）、`rw`、`omega`，无 `sorry` / `native_decide`。
-/

import Mathlib
import Zero.Z065_Heawood

namespace Zero

/-! ## 一、显式 6-环见证（围长 ≤ 6） -/

/-- 6-环 p₀ — L₀ — p₁ — L₁ — p₂ — L₆ — p₀ 对应的**道**：顶点序列 0,7,1,8,2,13,0
    （点 0,1,2；线 L₀ = 7、L₁ = 8、L₆ = 13）。 -/
private def z061SixCycleWalk : z065HeawoodGraph.Walk (0 : Fin 14) 0 :=
  SimpleGraph.Walk.cons (by decide : z065HeawoodGraph.Adj (0 : Fin 14) 7)
    (SimpleGraph.Walk.cons (by decide : z065HeawoodGraph.Adj (7 : Fin 14) 1)
      (SimpleGraph.Walk.cons (by decide : z065HeawoodGraph.Adj (1 : Fin 14) 8)
        (SimpleGraph.Walk.cons (by decide : z065HeawoodGraph.Adj (8 : Fin 14) 2)
          (SimpleGraph.Walk.cons (by decide : z065HeawoodGraph.Adj (2 : Fin 14) 13)
            (SimpleGraph.Walk.cons (by decide : z065HeawoodGraph.Adj (13 : Fin 14) 0)
              SimpleGraph.Walk.nil)))))

/-- 该道长度 = 6（6 条边，具体列表内核可算）。 -/
private theorem z061SixCycleWalk_length : z061SixCycleWalk.length = 6 := by decide

/-- 该道是环：`edges` 无重复、`support.tail = [7,1,8,2,13,0]` 无重复、长度 6 ≠ 0
    （三处都是具体列表上的 `decide`；`IsCycle` 本身没有 `Decidable` 实例，故逐字段交证据）。 -/
private theorem z061SixCycleWalk_isCycle : z061SixCycleWalk.IsCycle :=
  ⟨⟨⟨by decide⟩, by
      intro h
      have hlen := congrArg SimpleGraph.Walk.length h
      rw [z061SixCycleWalk_length] at hlen
      exact absurd hlen (by decide)⟩, by decide⟩

/-- **Z061 · 围长 ≤ 6**：显式 6-环见证给上界。 -/
theorem z061_girth_le_six : z065HeawoodGraph.girth ≤ 6 := by
  calc z065HeawoodGraph.girth ≤ z061SixCycleWalk.length :=
        z061SixCycleWalk_isCycle.girth_le_length
    _ = 6 := z061SixCycleWalk_length

/-! ## 二、环上顶点的下标账（下界用） -/

/-- 环上的顶点序 `v i := w.getVert i` 满足 `v_{m+1}` 就是 `support.tail` 的第 m 项：
    因为 `support = [v₀, v₁, …, v_L]`，其 `tail` 即 `[v₁, …, v_L]`（`List.get_tail`）。 -/
private theorem z061_getVert_succ_eq_tail {a : Fin 14} (w : z065HeawoodGraph.Walk a a) (m : ℕ)
    (hm : m + 1 ≤ w.length) :
    w.getVert (m + 1) =
      w.support.tail[m]'(by rw [List.length_tail, SimpleGraph.Walk.length_support]; omega) := by
  have hlen : m < w.support.tail.length := by
    rw [List.length_tail, SimpleGraph.Walk.length_support]
    omega
  rw [SimpleGraph.Walk.getVert_eq_support_getElem w hm]
  exact (List.get_tail w.support m hlen _).symm

/-! ## 三、下界：长度为 3/4/5 的环都不存在 -/

/-- **长度下界**：任何环道的长度 ≥ 6。

    设 L = w.length。L ≥ 3（`three_le_length`）且 L 为偶（二部 ⟹ 闭道长度为偶，
    把顶点序 `w.getVert` 直接喂给 Z065 的 `z065_heawood_no_odd_closed`）。
    若 L < 6 则 L = 4；此时 v₀,v₁,v₂,v₃ 互异（`support.tail` 无重复 + v₄ = v₀），
    而 v₁、v₃ 都在 N(v₀) ∩ N(v₂) 里，与该交的基数 ≤ 1（`z065_heawood_no_four_cycle`）矛盾。 -/
private theorem z061_six_le_length (a : Fin 14) (w : z065HeawoodGraph.Walk a a)
    (hw : w.IsCycle) : 6 ≤ w.length := by
  -- L ≥ 3（电路的长度下界）
  have h3 : 3 ≤ w.length := hw.three_le_length
  -- L 为偶：二部染色的奇偶推进（Z065 的 no_odd_closed），顶点序取 w.getVert
  have hev : Even w.length :=
    z065_heawood_no_odd_closed w.length w.getVert (fun i hi => w.adj_getVert_succ hi)
      (by rw [SimpleGraph.Walk.getVert_length, SimpleGraph.Walk.getVert_zero])
  by_contra h6
  -- 3 ≤ L < 6 且 L 偶 ⟹ L = 4
  obtain ⟨k, hk⟩ := hev
  have h4 : w.length = 4 := by omega
  -- v₄ = v₀（环回到起点）
  have h40 : w.getVert 4 = w.getVert 0 := by
    rw [← h4, SimpleGraph.Walk.getVert_length, SimpleGraph.Walk.getVert_zero]
  -- support.tail 无重复
  have hnodup : w.support.tail.Nodup := hw.support_nodup
  -- 把 v₀…v₃ 都搬到 support.tail 上（下标 3,0,1,2）
  have hv0 : w.getVert 0 = w.support.tail[(3 : ℕ)]'(by
      rw [List.length_tail, SimpleGraph.Walk.length_support]; omega) := by
    rw [show w.getVert 0 = w.getVert (3 + 1) from by rw [h40]]
    exact z061_getVert_succ_eq_tail w 3 (by omega)
  have hv1 : w.getVert 1 = w.support.tail[(0 : ℕ)]'(by
      rw [List.length_tail, SimpleGraph.Walk.length_support]; omega) :=
    z061_getVert_succ_eq_tail w 0 (by omega)
  have hv2 : w.getVert 2 = w.support.tail[(1 : ℕ)]'(by
      rw [List.length_tail, SimpleGraph.Walk.length_support]; omega) :=
    z061_getVert_succ_eq_tail w 1 (by omega)
  have hv3 : w.getVert 3 = w.support.tail[(2 : ℕ)]'(by
      rw [List.length_tail, SimpleGraph.Walk.length_support]; omega) :=
    z061_getVert_succ_eq_tail w 2 (by omega)
  have hinj := hnodup.injective_get
  -- v₀ ≠ v₂（在 tail 里是下标 3 与 1）
  have h02 : w.getVert 0 ≠ w.getVert 2 := by
    intro hcon
    have hget : w.support.tail[(3 : ℕ)]'(by
        rw [List.length_tail, SimpleGraph.Walk.length_support]; omega) =
        w.support.tail[(1 : ℕ)]'(by
        rw [List.length_tail, SimpleGraph.Walk.length_support]; omega) := by
      rw [← hv0, ← hv2]
      exact hcon
    have hfin := hinj hget
    have hval : (3 : ℕ) = 1 := by simpa using congrArg Fin.val hfin
    omega
  -- v₁ ≠ v₃（在 tail 里是下标 0 与 2）
  have h13 : w.getVert 1 ≠ w.getVert 3 := by
    intro hcon
    have hget : w.support.tail[(0 : ℕ)]'(by
        rw [List.length_tail, SimpleGraph.Walk.length_support]; omega) =
        w.support.tail[(2 : ℕ)]'(by
        rw [List.length_tail, SimpleGraph.Walk.length_support]; omega) := by
      rw [← hv1, ← hv3]
      exact hcon
    have hfin := hinj hget
    have hval : (0 : ℕ) = 2 := by simpa using congrArg Fin.val hfin
    omega
  -- 邻接（v₀—v₁、v₂—v₁、v₂—v₃、v₃—v₀）
  have ha01 : z065HeawoodGraph.Adj (w.getVert 0) (w.getVert 1) :=
    w.adj_getVert_succ (i := 0) (by omega)
  have ha21 : z065HeawoodGraph.Adj (w.getVert 2) (w.getVert 1) :=
    (w.adj_getVert_succ (i := 1) (by omega)).symm
  have ha23 : z065HeawoodGraph.Adj (w.getVert 2) (w.getVert 3) :=
    w.adj_getVert_succ (i := 2) (by omega)
  have ha03 : z065HeawoodGraph.Adj (w.getVert 0) (w.getVert 3) := by
    have h := w.adj_getVert_succ (i := 3) (by omega)
    rw [h40] at h
    exact h.symm
  -- v₁、v₃ 是 v₀ 与 v₂ 的两个互异公共邻居
  have hm1 : w.getVert 1 ∈
      (z065HeawoodGraph.neighborFinset (w.getVert 0)) ∩
        (z065HeawoodGraph.neighborFinset (w.getVert 2)) := by
    rw [Finset.mem_inter, SimpleGraph.mem_neighborFinset, SimpleGraph.mem_neighborFinset]
    exact ⟨ha01, ha21⟩
  have hm3 : w.getVert 3 ∈
      (z065HeawoodGraph.neighborFinset (w.getVert 0)) ∩
        (z065HeawoodGraph.neighborFinset (w.getVert 2)) := by
    rw [Finset.mem_inter, SimpleGraph.mem_neighborFinset, SimpleGraph.mem_neighborFinset]
    exact ⟨ha03, ha23⟩
  exact h13 (Finset.card_le_one.mp
    (z065_heawood_no_four_cycle (w.getVert 0) (w.getVert 2) h02) _ hm1 _ hm3)

/-- 上界方向对 `egirth` 也成立（`girth = egirth.toNat`，下界要借 `egirth` 传递）。 -/
private theorem z061_egirth_le_six : z065HeawoodGraph.egirth ≤ 6 := by
  calc z065HeawoodGraph.egirth ≤ (z061SixCycleWalk.length : ℕ∞) :=
        z061SixCycleWalk_isCycle.egirth_le_length
    _ = 6 := by rw [z061SixCycleWalk_length]; norm_num

/-- **Z061 · 6 ≤ 围长**：无 3-环（二部 ⟹ 无奇闭道）、无 5-环（同上）、无 4-环
    ⇒ 一切圈长 ≥ 6。经 `SimpleGraph.le_egirth` 把 `6 ≤ egirth` 化归到环道的长度下界，
    再经 `ENat.toNat_le_toNat` 落回 `girth = egirth.toNat`。 -/
theorem z061_six_le_girth : 6 ≤ z065HeawoodGraph.girth := by
  have heg : (6 : ℕ∞) ≤ z065HeawoodGraph.egirth :=
    SimpleGraph.le_egirth.mpr fun a w hw => by exact_mod_cast z061_six_le_length a w hw
  -- egirth ≠ ⊤：图非无环（6-环见证）
  have hne : z065HeawoodGraph.egirth ≠ ⊤ := by
    intro htop
    rw [SimpleGraph.egirth_eq_top] at htop
    exact htop z061SixCycleWalk z061SixCycleWalk_isCycle
  have h := ENat.toNat_le_toNat heg hne
  simpa [SimpleGraph.girth, ENat.toNat_natCast] using h

/-- **Z061 · 围长 = 6**：Heawood 图（Fano 关联图）的围长恰为 6。 -/
theorem z061_girth_eq_six : z065HeawoodGraph.girth = 6 :=
  le_antisymm z061_girth_le_six z061_six_le_girth

#print axioms z061_girth_le_six
#print axioms z061_six_le_girth
#print axioms z061_girth_eq_six

end Zero