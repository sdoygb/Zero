import Mathlib
open scoped Matrix
namespace Zero
variable {V : Type*} [LinearOrder V] {a b v : V}
-- 核心：给定 hinf : s(a,b).inf = a、hv : ¬(b = v)、hc : a = v，导出 False
example (hinf : ((s(a, b) : Sym2 V).inf : V) = a) (hv : ¬ (b = v)) (hc : a = v) : False := by
  rw [hinf] at hv
  exact hv hc
-- 另一形态：hv : ¬(b = v)，hc : a = v，但 hv 未 rw
example (hinf : ((s(a, b) : Sym2 V).inf : V) = a) (hv : ¬ (b = v)) (hc : a = v) : False := by
  exact absurd (hinf.trans hc) (by rw [hinf]; exact hv)
end Zero
