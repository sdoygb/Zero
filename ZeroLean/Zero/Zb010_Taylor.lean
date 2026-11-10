/-
  Zb010 §1 的 Taylor 余项：离散 Laplacian 的 `O(a²)` 展开。

  目标（正文定理 1.2）：
    `(f(x+a) − 2f(x) + f(x−a)) / a² = f''(x) + (a²/12) f⁗(x) + O(a⁴)`。

  证明路线：
    · 用 Mathlib 的 `taylor_isLittleO_univ`（`C⁶` 的 Taylor 余项是 `o(a⁶)`）；
    · 把 `T₆(x+a)` 与 `T₆(x−a)` 相加，奇次项相消，偶次项加倍；
    · 余项 `o(a⁶)` 与显式 `a⁶/360·f⁶(x)` 相加，得主差 `O(a⁶)`；
    · 除以 `a²`（在 `𝓝[≠] 0` 上）得 `O(a⁴)`。
-/

import Mathlib

open Finset Filter
open scoped Topology

set_option maxHeartbeats 800000

namespace Zero

/-- `T₆(x+a)+T₆(x−a)−2f(x)` 的显式偶次项：
    `a² f'' + a⁴/12 f⁗ + a⁶/360 f⁶`。 -/
theorem zb010_taylor_even_terms (f : ℝ → ℝ) (x a : ℝ) :
    taylorWithinEval f 6 Set.univ x (x + a)
      + taylorWithinEval f 6 Set.univ x (x - a) - 2 * f x
      = a ^ 2 * iteratedDeriv 2 f x
        + a ^ 4 / 12 * iteratedDeriv 4 f x
        + a ^ 6 / 360 * iteratedDeriv 6 f x := by
  rw [taylor_within_apply, taylor_within_apply]
  simp only [iteratedDerivWithin_univ]
  norm_num [Finset.sum_range_succ]
  ring

/-- **Zb010 · 主差是 `O(a⁶)`**：
    `f(x+a) − 2f(x) + f(x−a) − (a² f'' + a⁴/12 f⁗) = O(a⁶)`。 -/
theorem zb010_taylor_main_bigO (f : ℝ → ℝ) (hf : ContDiff ℝ 6 f) (x : ℝ) :
    (fun a : ℝ => f (x + a) - 2 * f x + f (x - a)
        - (a ^ 2 * iteratedDeriv 2 f x + a ^ 4 / 12 * iteratedDeriv 4 f x))
      =O[𝓝 (0 : ℝ)] (fun a : ℝ => a ^ 6) := by
  have hT := taylor_isLittleO_univ (f := f) (x₀ := x) (n := 6) hf
  have hplus : (fun a : ℝ => f (x + a) - taylorWithinEval f 6 Set.univ x (x + a))
      =o[𝓝 (0 : ℝ)] (fun a : ℝ => a ^ 6) := by
    have h1 : Tendsto (fun _ : ℝ => x) (𝓝 (0 : ℝ)) (𝓝 x) := tendsto_const_nhds
    have h2 : Tendsto (fun a : ℝ => a) (𝓝 (0 : ℝ)) (𝓝 (0 : ℝ)) := tendsto_id
    have hk : Tendsto (fun a : ℝ => x + a) (𝓝 (0 : ℝ)) (𝓝 x) := by simpa using h1.add h2
    have h := hT.comp_tendsto hk
    change (fun a : ℝ => f (x + a) - taylorWithinEval f 6 Set.univ x (x + a))
      =o[𝓝 (0 : ℝ)] (fun a : ℝ => ((x + a) + -x) ^ 6) at h
    have hfun : (fun a : ℝ => ((x + a) + -x) ^ 6) = (fun a : ℝ => a ^ 6) := by
      funext a; ring
    rwa [hfun] at h
  have hminus : (fun a : ℝ => f (x - a) - taylorWithinEval f 6 Set.univ x (x - a))
      =o[𝓝 (0 : ℝ)] (fun a : ℝ => a ^ 6) := by
    have h1 : Tendsto (fun _ : ℝ => x) (𝓝 (0 : ℝ)) (𝓝 x) := tendsto_const_nhds
    have h2 : Tendsto (fun a : ℝ => a) (𝓝 (0 : ℝ)) (𝓝 (0 : ℝ)) := tendsto_id
    have hk : Tendsto (fun a : ℝ => x - a) (𝓝 (0 : ℝ)) (𝓝 x) := by simpa using h1.sub h2
    have h := hT.comp_tendsto hk
    change (fun a : ℝ => f (x - a) - taylorWithinEval f 6 Set.univ x (x - a))
      =o[𝓝 (0 : ℝ)] (fun a : ℝ => ((x - a) + -x) ^ 6) at h
    have hfun : (fun a : ℝ => ((x - a) + -x) ^ 6) = (fun a : ℝ => a ^ 6) := by
      funext a; ring
    rwa [hfun] at h
  have hrem : (fun a : ℝ => (f (x + a) - taylorWithinEval f 6 Set.univ x (x + a))
        + (f (x - a) - taylorWithinEval f 6 Set.univ x (x - a)))
      =o[𝓝 (0 : ℝ)] (fun a : ℝ => a ^ 6) := by
    exact hplus.add hminus
  have hbig : (fun a : ℝ => a ^ 6 / 360 * iteratedDeriv 6 f x)
      =O[𝓝 (0 : ℝ)] (fun a : ℝ => a ^ 6) := by
    have h : (fun a : ℝ => a ^ 6) =O[𝓝 (0 : ℝ)] (fun a : ℝ => a ^ 6) := Asymptotics.isBigO_refl _ _
    simpa [div_eq_mul_inv, mul_comm, mul_left_comm, mul_assoc] using
      h.const_mul_left (iteratedDeriv 6 f x / 360)
  have hmain : (fun a : ℝ => f (x + a) - 2 * f x + f (x - a)
        - (a ^ 2 * iteratedDeriv 2 f x + a ^ 4 / 12 * iteratedDeriv 4 f x))
      = (fun a : ℝ => ((f (x + a) - taylorWithinEval f 6 Set.univ x (x + a))
        + (f (x - a) - taylorWithinEval f 6 Set.univ x (x - a)))
        + a ^ 6 / 360 * iteratedDeriv 6 f x) := by
    funext a
    have h := zb010_taylor_even_terms f x a
    linarith [h]
  rw [hmain]
  exact hrem.isBigO.add hbig

#print axioms zb010_taylor_even_terms

/-- **Zb010 · Taylor 余项**：`C⁶` 函数上
    `(f(x+a) − 2f(x) + f(x−a))/a² = f''(x) + (a²/12) f⁗(x) + O(a⁴)`。 -/
theorem zb010_taylor_remainder (f : ℝ → ℝ) (hf : ContDiff ℝ 6 f) (x : ℝ) :
    (fun a : ℝ => (f (x + a) - 2 * f x + f (x - a)) / a ^ 2
        - iteratedDeriv 2 f x - a ^ 2 / 12 * iteratedDeriv 4 f x)
      =O[𝓝[≠] (0 : ℝ)] (fun a : ℝ => a ^ 4) := by
  have hmain : (fun a : ℝ => f (x + a) - 2 * f x + f (x - a)
        - (a ^ 2 * iteratedDeriv 2 f x + a ^ 4 / 12 * iteratedDeriv 4 f x))
      =O[𝓝[≠] (0 : ℝ)] (fun a : ℝ => a ^ 6) :=
    (zb010_taylor_main_bigO f hf x).mono (nhdsWithin_le_nhds (s := {0}ᶜ) (a := 0))
  have hmul : (fun a : ℝ => (f (x + a) - 2 * f x + f (x - a)
        - (a ^ 2 * iteratedDeriv 2 f x + a ^ 4 / 12 * iteratedDeriv 4 f x)) * (a ^ 2)⁻¹)
      =O[𝓝[≠] (0 : ℝ)] (fun a : ℝ => a ^ 6 * (a ^ 2)⁻¹) :=
    hmain.mul (Asymptotics.isBigO_refl (fun a : ℝ => (a ^ 2)⁻¹) (𝓝[≠] (0 : ℝ)))
  refine hmul.congr' ?_ ?_
  · filter_upwards [self_mem_nhdsWithin] with a ha
    have ha0 : a ≠ 0 := ha
    have ha2 : a ^ 2 ≠ 0 := pow_ne_zero 2 ha0
    field_simp [ha0, ha2]
    ring
  · filter_upwards [self_mem_nhdsWithin] with a ha
    have ha2 : a ^ 2 ≠ 0 := pow_ne_zero 2 ha
    field_simp [ha, ha2]

#print axioms zb010_taylor_remainder

#print axioms zb010_taylor_main_bigO

end Zero
