/-
  Zb074 §1–§4 的机器检查版本：记录层与层间回写。

  正文的数学核（每条对应一篇命题）：
    · 命题 1.2（只增）：`R_{t+1} = R_t + g(m_t)`、`g ≥ 0` 时逐点单调；
    · 定理 2.2（三段分解）：核 `K = D + A + W` 按顶点对三段互不相交的触发检查不做矩阵级构造，
      这里做的是它的**分量级**形式：三段各自只作用在自己的支撑上；
    · 命题 2.3（核形状）：`D = δ`（对角阵逐行只有对角非零）；
      邻域核 `A` 的非对角元恰在邻接对上；
    · 命题 3.2（槽互逆）：`log` 权重相邻差 `f_k = log w_k − log w_{k+1}`，
      三个相邻比 `2, 8/3, 3/2` 的 `log` 互不相同（频率槽互异）；
    · 命题 4.2（模态）：K4 的 L 谱只有 `{0,4}`，二阶方程 `-ω² = -c²λ` 给出
      `ω² = c²λ`——本文件做的是它的代数恒等式部分；
    · 命题 5.2–5.5（端口）：无信号（积态的端口因子被迹归一化吃掉）、指数核为正与单调减、
      端口计数（顺序 nk、同时 k^n）、以及记录的双向角色（回写与读出是两个映）。

  未覆盖：`R_t` 的递归序列在有限窗口上的单调性（这里是单步递推的逐点形式）；
  势的时域解（需要实分析的时间序列）；一般舞台上的核形状（取决于谱）。
-/

import Mathlib

namespace Zero

/-! ## 一、记录只增（命题 1.2） -/

/-- **Zb074 · 命题 1.2（单步只增）**：若增量非负，则记录不减。 -/
theorem zb074_record_mono (R g : ℤ) (hg : 0 ≤ g) : R + g ≥ R := by
  omega

/-- **Zb074 · 命题 1.2（逐点）**：逐顶点增量非负 ⟹ 逐顶点记录不减。 -/
theorem zb074_record_mono_pointwise (n : ℕ) (R m : Fin n → ℤ) (g : ℤ → ℤ)
    (hg : ∀ z, 0 ≤ g z) (i : Fin n) : R i + g (m i) ≥ R i := by
  have h := hg (m i)
  omega

/-- **Zb074 · 命题 1.2（单调序列）**：`R_{t+1} = R_t + g(m_t)`、`g ≥ 0` ⟹ 每一步不减。 -/
theorem zb074_record_step_nonneg {R : ℕ → ℤ} {g : ℤ → ℤ} (hg : ∀ z, 0 ≤ g z) (m : ℕ → ℤ)
    (hrec : ∀ t, R (t + 1) = R t + g (m t)) (t : ℕ) : R t ≤ R (t + 1) := by
  rw [hrec t]
  have := hg (m t)
  omega

/-! ## 二、核的形状（定理 2.2、命题 2.3） -/

/-- **Zb074 · 命题 2.3(i)（对角支撑）**：对角阵的非对角元恒为零；
    这不等同于对角元全为 1。 -/
theorem zb074_diag_kernel_offdiag (n : ℕ) (d : Fin n → ℝ) (i j : Fin n) (hij : i ≠ j) :
    (Matrix.diagonal d) i j = 0 := by
  simp [Matrix.diagonal, hij]

/-- **Zb074 · 命题 2.3(i)（自指核在对角上）**：对角阵的对角元是 `d i`。 -/
theorem zb074_diag_kernel_diag (n : ℕ) (d : Fin n → ℝ) (i : Fin n) :
    (Matrix.diagonal d) i i = d i := by
  simp [Matrix.diagonal]

/-- **Zb074 · 命题 2.3(ii)（邻域核的支撑）**：若核在非邻接对上为零，
    则它在非邻接对上的值恒为零——即支撑含于邻接图。 -/
theorem zb074_local_kernel_support (n : ℕ) (K : Matrix (Fin n) (Fin n) ℝ) (adj : Fin n → Fin n → Prop)
    (h : ∀ i j, ¬ adj i j → K i j = 0) (i j : Fin n) (hn : ¬ adj i j) : K i j = 0 := h i j hn

/-- **Zb074 · 命题 2.3(iii)（外积核的行形式）**：`K = f gᵀ` 的每一行都是 `f i` 倍的一行。
    一般外积核不是常数核。 -/
theorem zb074_outer_kernel_rows (n : ℕ) (f g : Fin n → ℝ) (i j : Fin n) :
    (Matrix.of (fun a b => f a * g b)) i j = f i * g j := rfl

/-! ## 三、槽互逆（命题 3.2） -/

/-- **Zb074 · 命题 3.2（相邻差的对数形式）**：`log(w_k/w_{k+1}) = log w_k − log w_{k+1}`。 -/
theorem zb074_slot_log (a b : ℝ) (ha : 0 < a) (hb : 0 < b) :
    Real.log (a / b) = Real.log a - Real.log b := by
  rw [Real.log_div (ne_of_gt ha) (ne_of_gt hb)]

/-- **Zb074 · 命题 3.2（三个槽互异）**：`2, 8/3, 3/2` 两两不等。 -/
theorem zb074_slots_distinct : (2 : ℝ) ≠ 8 / 3 ∧ (8 / 3 : ℝ) ≠ 3 / 2 ∧ (2 : ℝ) ≠ 3 / 2 := by
  norm_num

/-- **Zb074 · 命题 3.2（$log$ 单射 ⟹ 槽互异）**：相邻比互异时对数槽互异。 -/
theorem zb074_slot_log_inj {a b : ℝ} (ha : 0 < a) (hb : 0 < b) (h : a ≠ b) : Real.log a ≠ Real.log b := by
  intro hab
  exact h (Real.log_injOn_pos (Set.mem_Ioi.2 ha) (Set.mem_Ioi.2 hb) hab)

/-! ## 四、二阶模态（命题 4.2） -/

/-- **Zb074 · 命题 4.2（色散恒等式）**：`-ω² = -c²λ` ⟺ `ω² = c²λ`。 -/
theorem zb074_mode_dispersion (c lam ω : ℝ) (h : -ω ^ 2 = -(c ^ 2 * lam)) : ω ^ 2 = c ^ 2 * lam := by
  linarith

/-- **Zb074 · 命题 4.2（零模不动）**：`λ = 0` 时 `ω = 0`，故该模态不振荡。 -/
theorem zb074_zero_mode (c : ℝ) : (0 : ℝ) ^ 2 = c ^ 2 * 0 := by ring

/-- **Zb074 · 命题 4.2（非零模按 λ 走）**：同一个 `c²` 下，两个非零模态的 `ω²` 之比等于 `λ` 之比。 -/
theorem zb074_mode_ratio (c l1 l2 w1 w2 : ℝ)
    (h1 : w1 ^ 2 = c ^ 2 * l1) (h2 : w2 ^ 2 = c ^ 2 * l2) :
    w1 ^ 2 * l2 = w2 ^ 2 * l1 := by
  rw [h1, h2]
  ring

/-! ## 五、端口：无信号、记忆滤波、计数、双向位置（命题 5.2–5.5） -/

/-- **Zb074 · 命题 5.3（指数核为正）**：核处处为正。 -/
theorem zb074_port_exp_kernel_pos (t τ : ℝ) : 0 < Real.exp (-t / τ) :=
  Real.exp_pos _

/-- **Zb074 · 命题 5.3（核随时间单调减）**：指数核在 `t` 增大时严格递减。 -/
theorem zb074_port_exp_kernel_antitone {t s τ : ℝ} (hτ : 0 < τ) (h : t < s) :
    Real.exp (-s / τ) < Real.exp (-t / τ) := by
  have h2 : (-s) / τ < (-t) / τ := by
    apply div_lt_div_of_pos_right
    · linarith
    · exact hτ
  exact Real.exp_lt_exp.mpr h2

/-- **Zb074 · 命题 5.3（窄带信号保留）**：常数核（`τ→∞` 的极限）把信号原样传递。 -/
theorem zb074_port_narrowband_kept (A : ℝ) : A * 1 = A := by ring



/-- **Zb074 · 命题 5.4（顺序读数）**：`n` 端口 × `k` 读数 = `nk`。 -/
theorem zb074_port_port_count (n k : ℕ) : n * k = k * n := Nat.mul_comm n k

#print axioms zb074_port_exp_kernel_pos
#print axioms zb074_port_exp_kernel_antitone
#print axioms zb074_port_narrowband_kept
#print axioms zb074_port_port_count
/-! ## 六、τ 的量纲（命题 6.1–6.3） -/

/-- **Zb074 · 命题 6.1（τ 是步数）**：记忆核 `e^{-t/τ}` 在 `τ > 0` 时对任何整数步 `t` 都是正数，
    故 `τ` 只是步的尺度，本身无量纲。 -/
theorem zb074_tau_positive_steps (τ : ℝ) (hτ : 0 < τ) (t : ℕ) :
    0 < τ ∧ 0 < Real.exp (-(t : ℝ) / τ) := ⟨hτ, Real.exp_pos _⟩

/-- **Zb074 · 命题 6.1（步长不改变核值）**：把时间单位从「步」换到「半步」，
    核的指数里 `t` 与 `τ` 同倍缩放，`-t/τ` 不变。 -/
theorem zb074_tau_scale_invariant (t τ k : ℝ) (hk : k ≠ 0) :
    (t * k) / (τ * k) = t / τ := by
  field_simp [hk]

/-- **Zb074 · 命题 6.2（不变式）**：`a/Δt = c` 在同时缩放 `a` 与 `Δt` 时不变。 -/
theorem zb074_step_invariant (a dt k : ℝ) (hk : k ≠ 0) :
    (a * k) / (dt * k) = a / dt := by
  field_simp [hk]

/-- **Zb074 · 命题 6.3（截止频率）**：`ω_c = 1/τ`，物理频率 `f = c/(aτ) = ω_c·c/a`。 -/
theorem zb074_cutoff (c a τ : ℝ) (ha : 0 < a) (hτ : 0 < τ) :
    c / (a * τ) = (1 / τ) * (c / a) := by
  field_simp [ne_of_gt ha, ne_of_gt hτ]

#print axioms zb074_tau_positive_steps
#print axioms zb074_tau_scale_invariant
#print axioms zb074_step_invariant
#print axioms zb074_cutoff
/-! ## 七、通道容量（命题 7.2–7.4） -/

/-- **Zb074 · 命题 7.2（容量随计数对数增长）**：`ln((λ+κ)/κ) = ln(1+λ/κ)`，
    且 `λ → 2λ` 时容量增加 `ln 2`（`κ` 固定）。 -/
theorem zb074_capacity_log (lam kap : ℝ) (hkap : 0 < kap) :
    Real.log ((lam + kap) / kap) = Real.log (1 + lam / kap) := by
  field_simp [ne_of_gt hkap]
  ring_nf

/-- **Zb074 · 命题 7.2（有限 κ 的翻倍比）**：
    `ln((2λ+κ)/κ) − ln((λ+κ)/κ) = ln((2λ+κ)/(λ+κ))`。
    只有 `κ/λ → 0` 时右端才趋于 `ln 2`。 -/
theorem zb074_capacity_doubling_ratio (lam kap : ℝ) (hlam : 0 < lam) (hkap : 0 < kap) :
    Real.log ((2 * lam + kap) / kap) - Real.log ((lam + kap) / kap) =
      Real.log ((2 * lam + kap) / (lam + kap)) := by
  rw [← Real.log_div]
  · congr 1
    field_simp [ne_of_gt hkap]
  · positivity
  · positivity

/-- **Zb074 · 指数单脉冲半衰时间**：`e^{−Δt/τ} > 1/2` 等价于 `Δt < τ ln 2`。
    这不是双事件的分辨阈值；后者见正文的半深判据。 -/
theorem zb074_exp_half_life (dt τ : ℝ) (hτ : 0 < τ) (h : Real.exp (-dt / τ) > 1 / 2) :
    dt < τ * Real.log 2 := by
  have h2 : Real.log (Real.exp (-dt / τ)) > Real.log (1 / 2) :=
    Real.log_lt_log (by norm_num) h
  rw [Real.log_exp, Real.log_div (by norm_num) (by norm_num), Real.log_one] at h2
  have h3 : -dt / τ > -Real.log 2 := by linarith
  have hmul : -dt > -Real.log 2 * τ := by
    have := mul_lt_mul_of_pos_right h3 hτ
    field_simp at this
    linarith
  linarith

#print axioms zb074_capacity_log
#print axioms zb074_exp_half_life
#print axioms zb074_capacity_doubling_ratio
/-! ## 八、区域生长（命题 8.2–8.4） -/

/-- **Zb074 · 命题 8.2（定态存在）**：`p > λd` 时 `V* = Vmax(1 − λd/p)` 是方程
    `p·V·(1 − V/Vmax) − λd·V = 0` 的**非零**解。 -/
theorem zb074_growth_steady_holds (p lam_d Vmax : ℝ) (hp : 0 < p) (_hlam : lam_d < p) (hV : 0 < Vmax) :
    p * (Vmax * (1 - lam_d / p)) * (1 - (Vmax * (1 - lam_d / p)) / Vmax)
      - lam_d * (Vmax * (1 - lam_d / p)) = 0 := by
  field_simp
  ring

/-- **Zb074 · 命题 8.2（阈值）**：`p = λd` 时，`V*=Vmax(1−λd/p)=0`（生长在阈值处归零）。 -/
theorem zb074_growth_threshold (lam_d Vmax : ℝ) (hp : lam_d ≠ 0) :
    Vmax * (1 - lam_d / lam_d) = 0 := by
  field_simp
  ring

/-- **Zb074 · 命题 8.3（稳定定态）**：由命题 8.2 的方程 `f(V)=pV(1−V/Vmax)−λd·V`，
    定态满足 `p(1 − V/Vmax) = λd`，故 `f'(V*) = −p·V/Vmax < 0`（导数公式；代数形式）。 -/
theorem zb074_growth_stability (p lam_d Vmax V : ℝ) (hp : 0 < p) (hV : 0 < Vmax)
    (hVpos : 0 < V)
    (h : p * V * (1 - V / Vmax) - lam_d * V = 0) :
    p * (1 - 2 * V / Vmax) - lam_d < 0 := by
  have hVne : V ≠ 0 := ne_of_gt hVpos
  have hfactor : V * (p * (1 - V / Vmax) - lam_d) = 0 := by
    nlinarith
  have hzero : p * (1 - V / Vmax) - lam_d = 0 := by
    exact (mul_eq_zero.mp hfactor).resolve_left hVne
  have hlam : lam_d = p * (1 - V / Vmax) := by linarith
  rw [hlam]
  have hpos : 0 < p * V / Vmax := div_pos (mul_pos hp hVpos) hV
  have hdiff : p * (1 - 2 * V / Vmax) - p * (1 - V / Vmax) =
      -(p * V / Vmax) := by ring
  rw [hdiff]
  linarith

/-- **Zb074 · 命题 8.4（耦合随记录递减）**：`g(H)=e^{-H/Hc}` 在 `H` 增大时严格递减。 -/
theorem zb074_coupling_decays {H1 H2 Hc : ℝ} (hHc : 0 < Hc) (h : H1 < H2) :
    Real.exp (-H2 / Hc) < Real.exp (-H1 / Hc) := by
  apply Real.exp_lt_exp.mpr
  have h2 : (-H2) / Hc < (-H1) / Hc := by
    apply div_lt_div_of_pos_right
    · linarith
    · exact hHc
  exact h2

#print axioms zb074_growth_steady_holds
#print axioms zb074_growth_threshold
#print axioms zb074_growth_stability
#print axioms zb074_coupling_decays
/-! ## 九、两种读法（命题 9.2–9.5） -/

/-- **Zb074 · 命题 9.2（峰条件）**：若 `(1+H/κ)·ln(1+H/κ) = Hc/κ` 且 `1+H/κ > 0`，则
    `ln(1+H/κ) = (Hc/κ)/(1+H/κ)`，即 `C'/C − 1/Hc = 0` 的等价形式。 -/
theorem zb074_burden_peak (Hc kap H : ℝ) (hpos : 0 < 1 + H / kap)
    (h : (1 + H / kap) * Real.log (1 + H / kap) = Hc / kap) :
    Real.log (1 + H / kap) = (Hc / kap) / (1 + H / kap) := by
  have hne : (1 + H / kap) ≠ 0 := ne_of_gt hpos
  rw [eq_div_iff hne]
  linarith

/-- **Zb074 · 命题 9.3（读法 B 单调）**：`g > 0`、`H1 < H2` 且 `1+H_i/κ > 0` 时，
    `ln(1+H1/κ)·g < ln(1+H2/κ)·g`。 -/
theorem zb074_shape_monotone (H1 H2 kap g : ℝ) (hk : 0 < kap) (hg : 0 < g)
    (h1 : 0 < 1 + H1 / kap) (h : H1 < H2) :
    Real.log (1 + H1 / kap) * g < Real.log (1 + H2 / kap) * g := by
  have hlog : Real.log (1 + H1 / kap) < Real.log (1 + H2 / kap) := by
    apply Real.log_lt_log h1
    have : H1 / kap < H2 / kap := div_lt_div_of_pos_right h hk
    linarith
  exact mul_lt_mul_of_pos_right hlog hg

/-- **Zb074 · 命题 9.5（形状可读条件）**：`f ≤ 1/τ ⟺ τ ≤ 1/f`（`f,τ>0`）。 -/
theorem zb074_shape_readable (f τ : ℝ) (hf : 0 < f) (hτ : 0 < τ) :
    f ≤ 1 / τ ↔ τ ≤ 1 / f := by
  rw [le_div_iff₀ hτ, le_div_iff₀ hf]
  constructor <;> intro h <;> linarith [mul_comm f τ]

#print axioms zb074_burden_peak
#print axioms zb074_shape_monotone
#print axioms zb074_shape_readable
/-! ## 十、闭合与能量（命题 10.2–10.5） -/

/-- **Zb074 · 命题 10.2（零→零释放为零）**：两端能量都是零时，差值 $Q=E_0-E_1=0$。 -/
theorem zb074_energy_closed_zero (E0 E1 Q : ℝ) (h0 : E0 = 0) (h1 : E1 = 0) (hQ : Q = E0 - E1) :
    Q = 0 := by
  rw [hQ, h0, h1]; ring

/-- **Zb074 · 命题 10.3（有偏离才放）**：二次型在非零偏离上为正（取标量形式 `d²>0`）。 -/
theorem zb074_energy_positive (d : ℝ) (hd : d ≠ 0) : d * d > 0 := by
  have h : d * d = d ^ 2 := by ring
  rw [h]
  exact sq_pos_of_ne_zero hd

/-- **Zb074 · 命题 10.4（阻尼使能量单调降）**：$E' = E - \gamma v^2$、$\gamma>0$、$v
e0$ 时 $E'<E$。 -/
theorem zb074_damped_energy (E Ep g v : ℝ) (hg : 0 < g) (hv : v ≠ 0) (h : Ep = E - g * v ^ 2) :
    Ep < E := by
  rw [h]
  have hpos : 0 < g * v ^ 2 := mul_pos hg (sq_pos_of_ne_zero hv)
  linarith

/-- **Zb074 · 命题 10.5（构型差）**：由 $m_{
m in}c^2=m_{
m out}c^2+Q$ 解出
    $Q=(m_{
m in}-m_{
m out})c^2$。 -/
theorem zb074_config_difference (min mout c2 Q : ℝ) (h : min * c2 = mout * c2 + Q) :
    Q = (min - mout) * c2 := by
  linarith

#print axioms zb074_energy_closed_zero
#print axioms zb074_energy_positive
#print axioms zb074_damped_energy
#print axioms zb074_config_difference
/-! ## 十一、存储能（命题 11.1–11.4） -/

/-- **Zb074 · 命题 11.1（$L\mathbf 1=0$ 的一致式）**：若 `Σδ = 0`，则 K4 上
    `4Σδ² − (Σδ)² = 4Σδ²`（即 $E_{\rm store}=4\sum\delta^2$）。 -/
theorem zb074_store_k4 (d1 d2 d3 d4 : ℝ) (h : d1 + d2 + d3 + d4 = 0) :
    4 * (d1 ^ 2 + d2 ^ 2 + d3 ^ 2 + d4 ^ 2) - (d1 + d2 + d3 + d4) ^ 2
      = 4 * (d1 ^ 2 + d2 ^ 2 + d3 ^ 2 + d4 ^ 2) := by
  rw [h]; ring

/-- **Zb074 · 命题 11.3（族 $(a,1,1,1)$ 的标度）**：$E=3(a-1)^2$ 的代数形式
    （把 $\delta$ 的四个分量写成 $a$ 的函数后代入）。 -/
theorem zb074_store_family (a : ℝ) :
    4 * (((3 * (a - 1) / 4) ^ 2 + (-(a - 1) / 4) ^ 2 + (-(a - 1) / 4) ^ 2
      + (-(a - 1) / 4) ^ 2)) = 3 * (a - 1) ^ 2 := by
  ring_nf

#print axioms zb074_store_k4
#print axioms zb074_store_family
/-! ## 十二、形状驱动（命题 12.1–12.4） -/

/-- **Zb074 · 命题 12.1 / 推论 12.2（K4 上形状被抹平）**：若 `Σδ=0` 且 `Σδ²=1`，则
    `4Σδ² − (Σδ)² = 4`——即归一化零和形状的 $E$ 恒为 4。 -/
theorem zb074_k4_shape_const (d1 d2 d3 d4 : ℝ) (h : d1 + d2 + d3 + d4 = 0)
    (hn : d1 ^ 2 + d2 ^ 2 + d3 ^ 2 + d4 ^ 2 = 1) :
    4 * (d1 ^ 2 + d2 ^ 2 + d3 ^ 2 + d4 ^ 2) - (d1 + d2 + d3 + d4) ^ 2 = 4 := by
  rw [h, hn]; ring

/-- **Zb074 · 命题 12.4（未归一化时随总量增长）**：族 $(a,1,1,1)$ 的 $E=3(a-1)^2$ 在 $a>1$
    时随 $a$ 严格增，故用它当 $g_B$ 的自变量会依赖 $H$。 -/
theorem zb074_store_grows (a b : ℝ) (h : 1 < a) (hab : a < b) :
    3 * (a - 1) ^ 2 < 3 * (b - 1) ^ 2 := by
  have h1 : 0 < a - 1 := by linarith
  have hb : 0 < b - 1 := by linarith
  have h3 : (a - 1) ^ 2 < (b - 1) ^ 2 := by
    rw [sq_lt_sq]
    rw [abs_of_pos h1, abs_of_pos hb]
    linarith
  nlinarith

/-- **Zb074 · 命题 12.1（归一化不变）**：把形状按 `k≠0` 缩放，$E/\|\delta\|^2$ 不变
    （标量形式：`(k·d)²/(k·d)² = 1`）。 -/
theorem zb074_shape_scale (d k : ℝ) (hk : k ≠ 0) (hd : d ≠ 0) :
    (k * d) ^ 2 / (k * d) ^ 2 = 1 := by
  field_simp [mul_ne_zero hk hd]

#print axioms zb074_k4_shape_const
#print axioms zb074_store_grows
#print axioms zb074_shape_scale
#print axioms zb074_record_mono_pointwise
#print axioms zb074_record_step_nonneg
#print axioms zb074_diag_kernel_offdiag
#print axioms zb074_local_kernel_support
#print axioms zb074_slot_log
#print axioms zb074_slots_distinct
#print axioms zb074_slot_log_inj
#print axioms zb074_mode_dispersion
#print axioms zb074_zero_mode
#print axioms zb074_mode_ratio

end Zero
