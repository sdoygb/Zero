/-
  Zero 推导的可机器检查核心 —— 第一片：K₄ 谱

  为什么先证这个：Z115 / Z116 的波数 k = √λ_max(K₄) = 2 全部建立在这条谱事实之上，
  而 Python 侧只做了「数值验算」（脚本算一遍特征值），不是「证明」。这里把它变成
  kernel 检查过的定理。

  工程取舍：本文件只依赖 Lean 核心（不 import Mathlib）。核心库的记号比 Mathlib 少
  ——没有 `ℤ`（要用 `Int`）、没有 `fin_cases` / `ring` / `linarith` / `Fintype`，
  只有 `omega`（Presburger 线性整数算术）。所以 K₄ 的向量用显式四元组 `V4` 表示而
  不是 `Fin 4 → Int`；换成 Mathlib 版（`Matrix (Fin n) (Fin n) ℤ` + 一般图 Laplacian
  + 重数论证）是后续步骤。

  对 Zero 审计最有用的不是"证了几个定理"，而是文末的 `#print axioms`：
  它列出每条定理实际依赖的**非导出前提**，那就是机械版的「要素来源表」。
-/

namespace Zero

/-- K₄ 上的一条"向量"（四个顶点各一个整数分量）。 -/
structure V4 where
  /-- 顶点 0 的分量 -/ a : Int
  /-- 顶点 1 的分量 -/ b : Int
  /-- 顶点 2 的分量 -/ c : Int
  /-- 顶点 3 的分量 -/ d : Int
deriving DecidableEq, Repr

/-- 零向量。 -/
def V4.zero : V4 := ⟨0, 0, 0, 0⟩

/-- K₄ 的图 Laplacian（每点度 3，完全图）：
    `(L v) i = 3 v i − Σ_{j≠i} v j`。 -/
def K4L (v : V4) : V4 :=
  ⟨3 * v.a - v.b - v.c - v.d,
   3 * v.b - v.a - v.c - v.d,
   3 * v.c - v.a - v.b - v.d,
   3 * v.d - v.a - v.b - v.c⟩

/-- 数乘（用于陈述特征向量）。 -/
def V4.smul (k : Int) (v : V4) : V4 := ⟨k * v.a, k * v.b, k * v.c, k * v.d⟩

/-- 定理 1：**核空间 = 常向量**（λ = 0 的特征空间就是"全同"方向）。 -/
theorem K4L_eq_zero_iff (v : V4) :
    K4L v = V4.zero ↔ (v.a = v.b ∧ v.b = v.c ∧ v.c = v.d) := by
  constructor
  · intro h
    have ha := congrArg V4.a h
    have hb := congrArg V4.b h
    have hc := congrArg V4.c h
    have hd := congrArg V4.d h
    simp only [K4L, V4.zero] at ha hb hc hd
    refine ⟨?_, ?_, ?_⟩ <;> omega
  · rintro ⟨ab, bc, cd⟩
    simp only [K4L, V4.zero, V4.mk.injEq]
    refine ⟨?_, ?_, ?_, ?_⟩ <;> omega

/-- 定理 2：**λ = 4 的特征空间 = 零和超平面**（`Σ v = 0`）。
    定理 1 与定理 2 合起来，就是 Z021/Z115 用到的那句「λ_max(K₄) = 4」的实质内容
    （重数 3 = 零和超平面的维数；维数论证留给 Mathlib 版）。 -/
theorem K4L_eigen_four_iff (v : V4) :
    K4L v = V4.smul 4 v ↔ v.a + v.b + v.c + v.d = 0 := by
  constructor
  · intro h
    have ha := congrArg V4.a h
    simp only [K4L, V4.smul] at ha
    omega
  · intro h
    simp only [K4L, V4.smul, V4.mk.injEq]
    refine ⟨?_, ?_, ?_, ?_⟩ <;> omega

/-- 推论（审计要用的那句话）：λ = 4 的模式**一定**带零和结构，
    所以"拿 K₄ 的谱去给另一张图（方格）上的波定波数"在形式上就缺一步 ——
    上面的定理只说了 K₄，没有说任何别的图。 -/
theorem eigen_four_zero_sum (v : V4)
    (h : K4L v = V4.smul 4 v) : v.a + v.b + v.c + v.d = 0 :=
  (K4L_eigen_four_iff v).mp h

/-- 具体实例：`(1,−1,0,0)` 确实是 λ = 4 的特征向量（Z021 里"3 个非平凡模式"之一）。 -/
example : K4L ⟨1, -1, 0, 0⟩ = V4.smul 4 ⟨1, -1, 0, 0⟩ := rfl

#print axioms K4L_eq_zero_iff
#print axioms K4L_eigen_four_iff
#print axioms eigen_four_zero_sum

end Zero
