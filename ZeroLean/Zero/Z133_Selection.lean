/-
  Z133（选择律：为何三代而非 2cos 族全部 11 个）的机器检查版本 —— 群论那一半。

  Z133 的核心判定：**槽数 3 不是从 2cos 族（11 个频率）里选出来的，而是由群的子群结构固定的**
  —— S₃ 只有一个 3 阶子群（即 Z₃ = A₃）；2cos 族只是频率库。

  这里把"只有一个 Z₃"这条在 S₃ = `Equiv.Perm (Fin 3)` 上钉死。两条命题合起来就是唯一性：
    · 阶 3 的元素恰 **2** 个（两个 3-循环）；
    · 这两个元素互为平方 ⟹ 它们生成**同一个**子群 ⟹ 3 阶子群唯一。

  工程注记：`orderOf` 在 Mathlib 里是 noncomputable 的，`decide` 归约不动它；
  因此这里用**等价的可计算刻画** `g³ = 1 ∧ g ≠ 1`（3 阶群里它就是"阶为 3"，
  因为元素阶整除群阶 3）。这样仍走内核 `decide`，不引入 `native_decide` 的 `Lean.ofReduceBool`。
-/

import Mathlib

namespace Zero

/-- **Z133 · S₃ 里阶 3 的元素恰 2 个**（用可计算刻画 `g³ = 1 ∧ g ≠ 1`；
    在 3 阶群中它等价于 `orderOf g = 3`）。 -/
theorem z133_order_three_elements_count :
    (Finset.univ.filter (fun g : Equiv.Perm (Fin 3) => g ^ 3 = 1 ∧ g ≠ 1)).card = 2 := by
  decide

/-- **Z133 · 唯一 Z₃ 子群**：S₃ 中任意两个不同的阶 3 元素互为平方。
    于是 ⟨g⟩ = {1, g, g²} ∋ h，即它们生成同一个子群 —— 3 阶子群只有 A₃ = Z₃ 一个。
    （配合上一条计数：阶 3 元素只有 g 与 g²，没有第三种 3 阶结构。） -/
theorem z133_three_cycles_are_squares :
    ∀ g h : Equiv.Perm (Fin 3), g ^ 3 = 1 → h ^ 3 = 1 → g ≠ 1 → h ≠ 1 → g ≠ h → h = g ^ 2 := by
  decide

#print axioms z133_order_three_elements_count
#print axioms z133_three_cycles_are_squares

end Zero
