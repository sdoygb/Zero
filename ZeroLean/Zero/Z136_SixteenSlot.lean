/-
  Z136 的机器检查版本：**十六槽载体上「ν_R 由轨道完备化逼出」** ——
  对应 Z136《右手中微子的强制性：十六槽载体与锚定完备化》§3–§4（核验 ①–④）。

  载体：`Z136Slot = Fin 4 × Fin 2 × Fin 2`（16 槽）= K₄ 顶点 × 同位旋重数 × 手征。
    · `z136Iso` / `z136Chir`：两条**自由**对合（核验 ①2、①3）；
    · `z136F`：SM 每代 **15 场**（**不含 ν_R**）的填充集（核验 ③1）；
    · `z136Orbit v c`：「顶点 v × 手征扇区 c」的 ℤ₂^{iso} 轨道（恰 2 元素）。

  主结果（核验 ③1–③4）：
    · `z136_filled_card`       : |填充集| = 15；
    · `z136_filled_rule`       : 色顶点四槽全填；色单态顶点填 `(0,↑,L)`、`(0,↓,L)`、`(0,↓,R)`；
    · `z136_leak_orbit_card`   : **恰 1 条轨道半满**（顶点 0、右手扇区）；
    · `z136_orbits_full_of_ne` : 其余 **7 条全满**；
    · `z136_complement_eq_nuR` : 补集恰是单点集 `{(0,↑,R)}`；
    · `z136_unique_completion` : **唯一**补位使 8 条轨道全满 ⟺ 补位 = `(0,↑,R)` = ν_R 槽。
  关键「输入侧无 ν_R」：`z136_nuR_not_in_F`（`(0,↑,R) ∉ z136F`）——ν_R 是本文件的**结论**，不是输入。

  量子数（核验 ④1、④2；一律取「六分之一」为单位，全整数算术）：
    · `z136_nuR_hypercharge`  : Y_L = −6、T₃^L = +3 ⟹ Y_R = Y_L + 2T₃^L = 0；
    · `z136_nuR_charge`       : T₃^R = 0、Y_R = 0 ⟹ Q = 2T₃ + Y = 0（sterile）；
    · `z136_eR_charge`、`z136_uR_hypercharge`、`z136_dR_hypercharge`：
      同一条规则复算 e_R、u_R、d_R ⟹ −12、+8、−4（即 −2、+4/3、−2/3）全对 PDG。

  与被取代的 `wip/_orphan_Z097_QNumbers.lean`（原 `Zero/Z097_QNumbers.lean`）的关系（**本文件取代它**，理由逐条可核）：
    · `z097_gmn (t3 y : ℚ) : t3 + y/2 = t3 + y/2 := rfl` 是 `x = x` 型**重言式**，不陈述任何事；
    · `z097_right_left_hypercharge` 在右手态上**仍带着 T₃**（其 `qR := t3 + yR/2`），
      于是证出 `Q_R = Q_L + T₃`；但右手是 SU(2) 单态，`T₃^R = 0`，正确关系是 `Q_R = Q_L`
      （`z136_qR_eq_qL`）。代入 e_R（T₃ = −1/2）应为 Q = −1，而前者给 −3/2。
      该式作为"符号恒等式"为真、作为"e_R 电荷"为假 —— 机器检查只保证前者
      （见 `z136_qR_ne_qL_without_T3R_zero`：不置 T₃^R = 0 时对 t3 = 1、yL = 0 已经不成立）。
-/

import Mathlib

namespace Zero

/-- Z136 的槽位：K₄ 顶点 × 同位旋 × 手性。 -/
abbrev Z136Slot := Fin 4 × Fin 2 × Fin 2

/-- `Fin 2` 上的翻转（对合）。 -/
def z136Flip (a : Fin 2) : Fin 2 := if a = 0 then 1 else 0

/-- ℤ₂^{iso}：`(v,a,c) ↦ (v,ā,c)`。 -/
def z136Iso (s : Z136Slot) : Z136Slot := (s.1, z136Flip s.2.1, s.2.2)

/-- ℤ₂^{chir}：`(v,a,c) ↦ (v,a,c̄)`。 -/
def z136Chir (s : Z136Slot) : Z136Slot := (s.1, s.2.1, z136Flip s.2.2)

/-- ν_R 的槽位 `(v=0, ↑, R)`。**本文件中它只作结论出现，不进入填充集。** -/
def z136NuRSlot : Z136Slot := (0, 0, 1)

/-- `e_R` 的槽位 `(v=0, ↓, R)`（属于 SM 15 场，是输入的已填元素）。 -/
def z136ERSlot : Z136Slot := (0, 1, 1)

/-- SM 每代 15 场（**不含 ν_R**）的填充判据：
    色顶点（`v ≠ 0`）四槽全填；色单态顶点只填 `(0,↑,L)`、`(0,↓,L)`、`(0,↓,R)`。
    记法：第二分量 `0/1` = 同位旋 ↑/↓，第三分量 `0/1` = 手征 `L/R`。 -/
def z136Filled (s : Z136Slot) : Prop := s.1 ≠ 0 ∨ s.2.2 = 0 ∨ s.2.1 = 1

/-- SM 每代 15 场的填充集（**输入**）。
    注意：判据必须**内联**进 lambda —— 具名 `def` 谓词不会被 instance search 展开，
    `Finset.filter z136Filled` 会因合成不出可计算的 `DecidablePred` 而报
    "failed to synthesize instance"（本文件第一版即栽在这里）。 -/
def z136F : Finset Z136Slot :=
  Finset.univ.filter fun s => s.1 ≠ 0 ∨ s.2.2 = 0 ∨ s.2.1 = 1

/-- 内联判据与具名谓词一致（把两者对上，供后文引用）。 -/
theorem z136_mem_F_iff : ∀ s : Z136Slot, s ∈ z136F ↔ z136Filled s := by
  simp [z136F, z136Filled]

/-- 轨道「顶点 `v` × 手征扇区 `c`」（ℤ₂^{iso} 轨道，恰 2 元素）。 -/
def z136Orbit (v : Fin 4) (c : Fin 2) : Finset Z136Slot :=
  Finset.univ.filter (fun s => s.1 = v ∧ s.2.2 = c)

/-- **Z136 · 载体槽数**：16 槽 = 4 顶点 × 2 同位旋 × 2 手征。 -/
theorem z136_slot_card : Fintype.card Z136Slot = 16 := by decide

/-- **Z136 · ℤ₂^{iso} 自由**：无不动点 ⟹ 每条轨道恰 2 元素。 -/
theorem z136_iso_free : ∀ s : Z136Slot, z136Iso s ≠ s := by decide

/-- **Z136 · ℤ₂^{iso} 对合**。 -/
theorem z136_iso_involutive : ∀ s : Z136Slot, z136Iso (z136Iso s) = s := by decide

/-- **Z136 · ℤ₂^{chir} 自由**。 -/
theorem z136_chir_free : ∀ s : Z136Slot, z136Chir s ≠ s := by decide

/-- **Z136 · 填充集大小**：SM 每代 15 场（不含 ν_R）。 -/
theorem z136_filled_card : z136F.card = 15 := by decide

/-- **Z136 · e_R 在输入里**（半满轨道上已填的那一个元素）。 -/
theorem z136_eR_in_F : z136ERSlot ∈ z136F := by decide

/-- **Z136 · ν_R 不在输入里**（关键：本文件的 ν_R 是结论，不是输入）。 -/
theorem z136_nuR_not_in_F : z136NuRSlot ∉ z136F := by decide

/-- **Z136 · 轨道半满计数**：8 条轨道里，恰 1 条半满——顶点 0、右手扇区。 -/
theorem z136_leak_orbit_card : (z136Orbit 0 1 ∩ z136F).card = 1 := by decide

/-- **Z136 · 其余 7 条轨道全满**（含色顶点左右手 6 条 + 色单态左手 1 条）。 -/
theorem z136_orbits_full_of_ne :
    ∀ v : Fin 4, ∀ c : Fin 2, (v, c) ≠ (0, 1) → (z136Orbit v c ∩ z136F).card = 2 := by
  decide

/-- **Z136 · 补集唯一**：未填的槽位恰是 `{(0,↑,R)}`。 -/
theorem z136_complement_eq_nuR : Finset.univ \ z136F = ({z136NuRSlot} : Finset Z136Slot) := by
  decide

/-- **Z136 · 完备化唯一（主结果）**：补一个槽位 `s` 能使 **8 条轨道全部**满，
    当且仅当 `s = (0,↑,R)` —— 即右手中微子的槽位。
    （`s ∈ z136F` 时补的是重复元素，半满轨道仍是 1，故不含在内。） -/
theorem z136_unique_completion :
    ∀ s : Z136Slot,
      (∀ v : Fin 4, ∀ c : Fin 2, (z136Orbit v c ∩ insert s z136F).card = 2) ↔
        s = z136NuRSlot := by
  decide

/-- **Z136 · 修正后的手征配对关系**：右手是 SU(2) 单态（`T₃^R = 0`），
    故 `Q_R = T₃^R + Y_R/2 = (Y_L + 2T₃^L)/2 = T₃^L + Y_L/2 = Q_L`。 -/
theorem z136_qR_eq_qL (t3 yL : ℚ) :
    let t3R : ℚ := 0
    let yR := yL + 2 * t3
    t3R + yR / 2 = t3 + yL / 2 := by
  simp only
  ring

/-- **Z136 · 反面**：不置 `T₃^R = 0` 时，`Q_R = Q_L` 已经不成立（`t3 = 1, yL = 0` 即反例）。
    这是 `z097_right_left_hypercharge`（得 `Q_R = Q_L + T₃`）不能当"右手电荷"用的机器证据。 -/
theorem z136_qR_ne_qL_without_T3R_zero : ∃ t3 yL : ℚ, t3 + (yL + 2 * t3) / 2 ≠ t3 + yL / 2 :=
  ⟨1, 0, by norm_num⟩

/-- **Z136 · ν_R 的超荷为零**（六分之一为单位：Y_L = −6、T₃^L = +3 ⟹ Y_R = 0）。 -/
theorem z136_nuR_hypercharge : (-6 : ℤ) + 2 * 3 = 0 := by decide

/-- **Z136 · ν_R 的电荷为零**（六分之一为单位：Q = 2T₃ + Y = 0 + 0 = 0 ⟹ sterile）。 -/
theorem z136_nuR_charge : 2 * 0 + ((-6 : ℤ) + 2 * 3) = 0 := by decide

/-- **Z136 · 对偶核对 e_R**：Y_R = −6 + 2(−3) = −12（即 −2）⟹ Q = −12/6 = −1 ✓ PDG。 -/
theorem z136_eR_charge : 2 * 0 + ((-6 : ℤ) + 2 * (-3)) = -12 := by decide

/-- **Z136 · 对偶核对 u_R**：Y_R = +2 + 2(+3) = +8（即 +4/3）✓ PDG。 -/
theorem z136_uR_hypercharge : (2 : ℤ) + 2 * 3 = 8 := by decide

/-- **Z136 · 对偶核对 d_R**：Y_R = +2 + 2(−3) = −4（即 −2/3）✓ PDG。 -/
theorem z136_dR_hypercharge : (2 : ℤ) + 2 * (-3) = -4 := by decide

#print axioms z136_unique_completion

end Zero
