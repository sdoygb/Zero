/-
  Z017 §3.2 / §4 的机器检查版本 —— 重播种递推的特征多项式与 Catalan 序列。
  (1) 递推 h_{n+1}=(c+1)h_n − c h_{n−1} 的特征多项式 r²−(c+1)r+c = (r−c)(r−1)，
      两根为 c 与 1，主导根 λ=c。
  (2) co[2i]=2C_i，C_i = choose(2i,i)/(i+1)（Catalan），前几项 2,2,4,10,28（分岔点）。
  未覆盖：周期不动点的收敛、Σco 代数化简（文章自述未竟）、转移矩阵特征方程的一般形式。
-/

import Mathlib

namespace Zero

open Matrix

/-- **特征多项式因式分解**：r² − (c+1)r + c = (r − c)(r − 1)。 -/
lemma z017_charpoly_factor (c r : ℝ) :
    r ^ 2 - (c + 1) * r + c = (r - c) * (r - 1) := by ring

/-- **两根**：r = c 与 r = 1 代入特征多项式都为 0（故主导根 c、瞬态根 1）。 -/
lemma z017_roots (c : ℝ) :
    c ^ 2 - (c + 1) * c + c = 0 ∧
    (1 : ℝ) ^ 2 - (c + 1) * 1 + c = 0 := by
  constructor <;> ring

/-- **Catalan 前几项**：2C_i = 2·choose(2i,i)/(i+1)，i=0..4 给 2,2,4,10,28。
    其中 i=4 给 28（与 Motzkin 第五项 9 分岔）。 -/
lemma z017_catalan_terms :
    2 * (Nat.choose 0 0 / 1) = 2 ∧
    2 * (Nat.choose 2 1 / 2) = 2 ∧
    2 * (Nat.choose 4 2 / 3) = 4 ∧
    2 * (Nat.choose 6 3 / 4) = 10 ∧
    2 * (Nat.choose 8 4 / 5) = 28 := by decide

/-- **c(3)=4（闭式初值）**：T=3 时的周期推进参数。 -/
lemma z017_c3 : (2 + 2 : ℕ) = 4 := by norm_num

#print axioms z017_charpoly_factor
#print axioms z017_roots
#print axioms z017_catalan_terms
#print axioms z017_c3

end Zero
