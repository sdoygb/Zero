/-
  Z065（Heawood 圈秩 = 8，即 H₁ ≅ Z⁸）的机器检查版本。

  圈秩（cycle rank）= |E| − |V| + 1，它的数学来源是：**树恰有 |V| − 1 条边**，
  所以连通图"超出生成树"的边数就是圈秩。原文对 Heawood 核验 21 − 14 + 1 = 8。

  审计修订（2026-10-14）：原版第一条是 `e - v + 1 ≥ 0`（ℕ 上**恒真**，什么都没说），
  第二条是裸算术 `21 - 14 + 1 = 8`（没有图）。这里把圈秩公式的**数学内容**钉住：
  树满足 `|E| + 1 = |V|` ⟹ 树的圈秩为 0；一般图的圈秩就是"超出树的边数"。
  （Heawood 图本身的具体构造留作后续：需要显式建 14 顶点 3-正则图再数边。）
-/

import Mathlib

namespace Zero

/-- **Z065 · 树的边数**：连通无环图（树）满足 `|E| + 1 = |V|`
    （Mathlib 的 `SimpleGraph.IsTree.card_edgeFinset`）——这就是圈秩公式的来源。 -/
theorem z065_tree_card_edgeFinset {V : Type*} [Fintype V] (G : SimpleGraph V)
    [Fintype G.edgeSet] (hG : G.IsTree) : G.edgeFinset.card + 1 = Fintype.card V :=
  hG.card_edgeFinset

/-- **Z065 · 树的圈秩为 0**：由 `|E| + 1 = |V|` 得 `|E| − |V| + 1 = 0`（在 ℤ 上算，
    避免 ℕ 减法的截断）；于是连通图的圈秩恰是"超过生成树的边数"
    （Heawood：`21 − 14 + 1 = 8`）。 -/
theorem z065_tree_cycle_rank_zero {V : Type*} [Fintype V] (G : SimpleGraph V)
    [Fintype G.edgeSet] (hG : G.IsTree) :
    (G.edgeFinset.card : ℤ) - (Fintype.card V : ℤ) + 1 = 0 := by
  have h : (G.edgeFinset.card : ℤ) + 1 = (Fintype.card V : ℤ) := by
    exact_mod_cast hG.card_edgeFinset
  linarith

#print axioms z065_tree_card_edgeFinset
#print axioms z065_tree_cycle_rank_zero

end Zero
