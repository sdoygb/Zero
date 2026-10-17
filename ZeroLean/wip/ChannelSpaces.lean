/-
  D0 的通道空间维数（图/组合版）。

  一体空间：顶点空间，维数 = 槽数（16 槽载体 `Fin 4 × Fin 2 × Fin 2` → 16）。
  两体空间：2 元子集空间，维数 = C(n,2)（8 槽载体 → 28）；
            质量单元是振幅²，故取 (C(8,2))² = 784。

  这些是 Z136 载体上的组合事实；"哪个算子取哪个空间"仍是推导稿的建模读法。
-/

import Mathlib

namespace Zero

/-- **两体通道数**：$n$ 个槽的 2 元子集数 $=\binom n2$。 -/
theorem two_body_channels (n : ℕ) :
    ((Finset.univ : Finset (Fin n)).powersetCard 2).card = Nat.choose n 2 := by
  rw [Finset.card_powersetCard, Finset.card_univ, Fintype.card_fin]

/-- **8 槽的两体通道数**：$\binom82=28$。 -/
theorem two_body_channels_8 : Nat.choose 8 2 = 28 := by decide

/-- **质量单元的通道数**：$\binom82^2=784$。 -/
theorem mass_unit_channels : (Nat.choose 8 2) ^ 2 = 784 := by decide

/-- **一体通道数**（16 槽载体）：$|Fin 4 × Fin 2 × Fin 2| = 16$。 -/
theorem one_body_channels :
    Fintype.card (Fin 4 × Fin 2 × Fin 2) = 16 := by decide

#print axioms two_body_channels
#print axioms two_body_channels_8
#print axioms mass_unit_channels
#print axioms one_body_channels

end Zero
