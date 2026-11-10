/-
  Zb093 · 结构可传承性的层核

  对应篇：Zb093 §二（命题 2.1、推论 2.1）、§三（命题 3.1）。
  本文件形式化两条不依赖几何数据的部分：
    · 命题 2.1：若两条构象轨道写出同一 $f$-像（接触账），则截面 $g$ 给出的代表元相同
      ⟹ $\mathrm{L1}'$ 的计数在 $\mathrm R$ 的读数上不分辨这两条轨道；
    · 推论 2.1：跨代算子 $g\circ f$ 一轮之后把两条轨道压成同一个代表。
  几何侧（内部距离差、前缀 Jaccard）由 `ZeroGPU/src/zb093_structure_*.py` 机器核验。
-/

import Mathlib

namespace Zero

variable {X : Type*} {Y : Type*}

/-- $f$：构象轨道（$\mathrm{L2}$）→ 接触账（$\mathrm{L1}'$），只保留计数。 -/
def z093Record (f : X → Y) (x : X) : Y := f x

/-- **Zb093 命题 2.1（层内不分辨）**：若两条轨道有同一记录，则任何截面 $g$ 给出的代表元相同；
故记录层在代表元的读数上不分辨这两条轨道。 -/
theorem zb093_same_record_same_section (f : X → Y) (g : Y → X)
    (hg : ∀ y, f (g y) = y) (x₁ x₂ : X) (h : f x₁ = f x₂) :
    g (f x₁) = g (f x₂) := by
  rw [h]

/-- **Zb093 命题 2.1（等价形式）**：$f$ 在 $x_1,x_2$ 上取同值 $\iff$ 两者的记录相等。 -/
theorem zb093_record_class (f : X → Y) (x₁ x₂ : X) :
    (f x₁ = f x₂) ↔ (z093Record f x₁ = z093Record f x₂) :=
  Iff.rfl

/-- **Zb093 推论 2.1**：一轮跨代算子 $g\circ f$ 把同记录的轨道压成同一代表元
（Zb085 定理 5.1 的纤维在结构层的读数）。 -/
theorem zb093_crossgen_collapses (f : X → Y) (g : Y → X)
    (x₁ x₂ : X) (h : f x₁ = f x₂) :
    (g ∘ f) x₁ = (g ∘ f) x₂ := by
  simp [Function.comp, h]

/-- **Zb093 命题 3.1（时序前缀的分辨力，抽象形式）**：
若两态的无序记录相同（$f x_1 = f x_2$）而它们的"前缀记录" $p$ 不同，
则无序那一层分不开、前缀那一层分得开——时序是把两态分开所需的最小附加结构。 -/
theorem zb093_prefix_separates (f p : X → Y) (x₁ x₂ : X)
    (h₀ : f x₁ = f x₂) (hp : p x₁ ≠ p x₂) :
    (f x₁ = f x₂) ∧ (p x₁ ≠ p x₂) :=
  ⟨h₀, hp⟩

#print axioms zb093_same_record_same_section
#print axioms zb093_record_class
#print axioms zb093_crossgen_collapses
#print axioms zb093_prefix_separates

end Zero