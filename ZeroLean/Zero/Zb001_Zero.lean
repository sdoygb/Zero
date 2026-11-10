/-
  Zb001 §2 的机器检查版本：K₄ 上的「零 = 核，场 = 零和超平面」。

  正文定理 2.2 / 推论 2.3 的实质是两条有限维线性代数事实：
    · `K4L v = 0 ↔ v` 是常向量（核空间 = 齐平方向）；
    · `K4L v = 4 v ↔ Σ v = 0`（λ = 4 的特征空间 = 零和超平面）。

  证明复用 `Zero.Basic` 里已经过内核检查的同一组命题（那里用显式四元组
  `V4`，避免依赖 Mathlib 的矩阵库；本文件只把名字对到 Zb001 的篇号上）。
-/

import Zero.Basic

namespace Zero

/-- **Zb001 · 零 = 核**：`K4L v = 0` 当且仅当 `v` 的四个分量全相等。 -/
theorem zb001_kernel_eq_const (v : V4) :
    K4L v = V4.zero ↔ (v.a = v.b ∧ v.b = v.c ∧ v.c = v.d) :=
  K4L_eq_zero_iff v

/-- **Zb001 · 场 = 零和超平面**：`v` 是 `K4L` 的 λ = 4 特征向量，当且仅当 `Σ v = 0`。 -/
theorem zb001_eigen_four_iff (v : V4) :
    K4L v = V4.smul 4 v ↔ v.a + v.b + v.c + v.d = 0 :=
  K4L_eigen_four_iff v

/-- **Zb001 · 推论 2.3 的直接形式**：λ = 4 的模式一定落在零和超平面里。 -/
theorem zb001_eigen_four_zero_sum (v : V4)
    (h : K4L v = V4.smul 4 v) : v.a + v.b + v.c + v.d = 0 :=
  eigen_four_zero_sum v h

#print axioms zb001_kernel_eq_const
#print axioms zb001_eigen_four_iff
#print axioms zb001_eigen_four_zero_sum

end Zero
