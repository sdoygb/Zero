/-
  Z117 §一E6（度量错配）的机器检查版本 —— 最小反例。

  背景：Z115 §2④ 把"图上路径差"取成**曼哈顿距离差** `|x−6| − |x−14|`。实测（CK4a）
  图上波的相位走的是**欧氏程函**：5/5 条非平凡割线欧氏胜。

  这里把两者的差别证成一条定理，用最小的几何：
    · 两个缝在 (0,0) 与 (2,0)（缝距 2）；
    · 屏点在 (0,h)：`h ≠ 0` 就是"离开缝所在的那一行"（Z115 恰好取 h = 0，两制度量在那
      一行重合，所以错误被掩盖）；
    · 到两缝的曼哈顿距离之差 = |0−0| + |h| − (|0−2| + |h|) = −2 （与 h 无关）；
    · 到两缝的欧氏距离之差 = h − √(4 + h²)，它等于 −2 **当且仅当** h = 0。

  于是：只要屏不在缝所在行，两种"图上距离差"就给出不同的相位差。
-/

import Mathlib

namespace Zero

/-- **Z117 · L1 ≠ 欧氏（离轴即分道）**：`h ≠ 0` 时，欧氏程差 `h − √(4+h²)` 不等于曼哈顿
    路径差 `−2`。 -/
theorem z117_euclid_ne_l1 (h : ℝ) (hh : h ≠ 0) :
    h - Real.sqrt (4 + h ^ 2) ≠ -2 := by
  intro hcontra
  -- 由 h − √(4+h²) = −2 得 √(4+h²) = h + 2
  have hsqrt : Real.sqrt (4 + h ^ 2) = h + 2 := by linarith
  -- 右边非负（左边是 √），于是可平方
  have hsq : 4 + h ^ 2 = (h + 2) ^ 2 := by
    rw [← hsqrt, Real.sq_sqrt (by positivity)]
  -- 展开即得 4h = 0，与 h ≠ 0 矛盾
  have hzero : h = 0 := by nlinarith [hsq]
  exact hh hzero

/-- 曼哈顿距离差与 h 无关（这正是"模型说每一行条纹都一样"的形式表述）。 -/
theorem z117_l1_diff_const (h : ℝ) :
    (|(0 : ℝ) - 0| + |h - 0|) - (|(0 : ℝ) - 2| + |h - 0|) = -2 := by
  simp

/-- 两制度量在共线行（h = 0）重合 —— 这就是 Z115 把屏放在 h = 0 时错误被掩盖的原因。 -/
theorem z117_metrics_agree_on_axis :
    (0 : ℝ) - Real.sqrt (4 + 0 ^ 2) = -2 := by
  norm_num

#print axioms z117_euclid_ne_l1
#print axioms z117_l1_diff_const

end Zero
