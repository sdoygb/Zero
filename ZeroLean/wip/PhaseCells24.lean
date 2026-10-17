/-
  24 的上游结构：闭合 Z8 与代际 Z3 互素，生成的相位格是 Z24（CRT）。

  Z121：闭合把 $2\pi$ 分成 8 份（每份 $2\pi/8$）。
  Z090：$K_4$ 的 $C_3$ 在三代上给出 $1\oplus\omega\oplus\omega^2$（旋转 $2\pi/3$）。
  两者都是圆上的旋转，且 $\gcd(8,3)=1$，故生成的群是 $\mathbb Z_{24}$
  （中国剩余定理）：$2\pi$ 分成 24 格，每格 $2\pi/24=\pi/12$——正是 $\delta_0$ 的角量子。
-/

import Mathlib

namespace Zero

/-- $8$ 与 $3$ 互素（闭合 8 分与三代 $3$ 分的互素性）。 -/
theorem coprime_8_3 : Nat.Coprime 8 3 := by decide

/-- **CRT**：相位格 $\mathbb Z_{24}\cong\mathbb Z_8\times\mathbb Z_3$——闭合与代际生成 $24$ 格。 -/
noncomputable def phase_cells_24 : ZMod 24 ≃+* (ZMod 8 × ZMod 3) :=
  ZMod.chineseRemainder (m := 8) (n := 3) (by decide)

/-- 相位格数：$|Z_{24}|=24$。 -/
theorem phase_cells_card : Fintype.card (ZMod 24) = 24 := by
  simp

/-- 每八分圆的格数：$24/8=3$。 -/
theorem cells_per_octant : 24 / 8 = 3 := by norm_num

#print axioms coprime_8_3
#print axioms phase_cells_24
#print axioms phase_cells_card

end Zero
