/-
  Zb092 · 分子手性判据的代数核

  对应篇：Zb092 §2（命题 2.1、定理 2.3、推论 2.4）。
  本文件只核验代数与有限性部分：
    · 自同构集合有限（命题 2.1）；
    · 相位场在镜面自同构下按共轭变换 ⟹ 乘积为实 ⟹ Im ∏z = 0（定理 2.3）；
    · 逆否：Im ∏z ≠ 0 ⟹ 不存在那样的镜面自同构（推论 2.4）。
  几何侧（det O = −1 的等距实现、Kabsch/SVD 的最优解）属数值判定，
  由 `ZeroGPU/src/zb092_chirality_topology.py` 机器核验，不在本文件形式化。
-/

import Mathlib

namespace Zero

open scoped BigOperators

variable {V : Type*} [Fintype V]

/-- 图 `G` 的自同构集合：保持邻接关系的顶点置换（Zb092 定义 1.1、命题 2.1）。 -/
def z092AutSet (G : SimpleGraph V) : Set (V ≃ V) :=
  {π | ∀ u v : V, G.Adj (π u) (π v) ↔ G.Adj u v}

/-- **Zb092 命题 2.1**：自同构集合有限（它是有限置换群的子集）。 -/
theorem zb092_aut_finite (G : SimpleGraph V) : (z092AutSet G).Finite :=
  Set.toFinite _

/-- **Zb092 定理 2.3（第一步）**：若相位场 `z` 在置换 `π` 下按共轭变换，
则整积等于自身的共轭。 -/
theorem zb092_prod_mirror_eq_conj (z : V → ℂ) (π : V ≃ V)
    (h : ∀ v : V, z (π v) = starRingEnd ℂ (z v)) :
    (∏ v, z v) = starRingEnd ℂ (∏ v, z v) := by
  calc
    ∏ v, z v = ∏ v, z (π v) := by rw [Equiv.prod_comp]
    _ = ∏ v, starRingEnd ℂ (z v) := by simp [h]
    _ = starRingEnd ℂ (∏ v, z v) := by simp

/-- **Zb092 定理 2.3**：镜面自同构把整积强制为实数，故手征荷的虚部为零。 -/
theorem zb092_im_prod_eq_zero (z : V → ℂ) (π : V ≃ V)
    (h : ∀ v : V, z (π v) = starRingEnd ℂ (z v)) :
    (∏ v, z v).im = 0 := by
  have hc := zb092_prod_mirror_eq_conj z π h
  exact Complex.conj_eq_iff_im.mp hc.symm

/-- **Zb092 推论 2.4**：若整积的虚部非零，则不存在那样的镜面自同构
（即判据必判手性）。 -/
theorem zb092_no_mirror_of_im_ne_zero (z : V → ℂ)
    (hne : (∏ v, z v).im ≠ 0) :
    ∀ π : V ≃ V, ¬∀ v : V, z (π v) = starRingEnd ℂ (z v) := by
  intro π h
  exact hne (zb092_im_prod_eq_zero z π h)

#print axioms zb092_aut_finite
#print axioms zb092_prod_mirror_eq_conj
#print axioms zb092_im_prod_eq_zero
#print axioms zb092_no_mirror_of_im_ne_zero

end Zero
