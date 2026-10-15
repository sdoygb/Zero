/-
  Z117 §一E7（亚 Nyquist 混叠）的机器检查版本。

  背景：Z115 / Z116 取 k = √λ_max(K₄) = 2，于是强度条纹的空间频率是 2k = 4 rad/格，
  超过格点的 Nyquist 极限 π —— 也就是说：**在整数采样下，频率 4 与频率 4 − 2π 给出
  完全相同的采样序列**，图上无法区分它们（这正是"1.571 间距的条纹在图上不可表示"的
  算术内核）。

  这里把"不可区分"证成一条定理：对任意整数 n，两个复指数的取值逐点相等。
-/

import Mathlib

namespace Zero

open Complex

/-- **Z117 · 亚 Nyquist 混叠**：整数采样下，空间频率 `4` 与 `4 − 2π` 的采样序列恒等。
    （`n` 是格点序号；`exp (ω n I)` 是频率 ω 的平面波在整数点的取值。） -/
theorem z117_aliasing (n : ℤ) :
    Complex.exp (4 * n * Complex.I) = Complex.exp ((4 - 2 * Real.pi) * n * Complex.I) := by
  have hsplit : (4 - 2 * Real.pi) * n * Complex.I
      = 4 * n * Complex.I - 2 * Real.pi * n * Complex.I := by ring
  rw [hsplit, Complex.exp_sub]
  have htwo : 2 * Real.pi * n * Complex.I = (n : ℂ) * (2 * Real.pi * Complex.I) := by ring
  rw [htwo, Complex.exp_int_mul, Complex.exp_two_pi_mul_I]
  simp

/-- 同样的内容用"采样序列相等"表述：把平面波当成 `ℤ → ℂ` 的函数，两者相等。 -/
theorem z117_aliasing_funext :
    (fun n : ℤ => Complex.exp (4 * (n : ℂ) * Complex.I))
      = fun n : ℤ => Complex.exp ((4 - 2 * Real.pi) * (n : ℂ) * Complex.I) := by
  funext n
  exact z117_aliasing n

/-- **Z117 · Nyquist 判据本身**：条纹空间频率 `2k = 4`（Z115/Z116 取 k = 2）超过格点的
    Nyquist 极限 `π`。与上面的混叠定理合起来，就是"1.571 间距的条纹在图上不可表示"的完整理由：
    频率已越过 π，采样后与 `4 − 2π` 无法区分。 -/
theorem z117_beyond_nyquist : Real.pi < 4 := Real.pi_lt_four

#print axioms z117_aliasing
#print axioms z117_aliasing_funext
#print axioms z117_beyond_nyquist

end Zero
