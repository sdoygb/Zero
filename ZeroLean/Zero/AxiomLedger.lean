/-
  要素来源表的**机械版**演示。

  Zero 的审计问题：「这条结论是内生的，还是把某个前提换了名字？」在 Lean 里这个问题
  有精确答案：`#print axioms` 会列出每条定理依赖的**全部非导出前提**。

  规则：
    · 只依赖 Lean 基石公理（propext / Quot.sound / Classical.choice）→ 这条是"纯数学"，
      对 Zero 而言就是【导出】；
    · 依赖任何自定义 `axiom` / `opaque` / `constant` → 那个名字就是【具名输入】，
      而且**逐个列出**，跑不掉；
    · 依赖 `sorryAx` → 证明有洞（本项目用 CI 禁止）。
-/

namespace Zero.AuditDemo

/-- 把"具名输入"显式声明成公理：这里假设 Z115 的 k = 2 作为前提被引入
    （真实版本会写成 `axiom wave_number_on_grid : ...` 的具体命题）。 -/
axiom k_is_two : (2 : Int) = 2

/-- 用了这条输入的结果。 -/
theorem needs_input : (2 : Int) = 2 := k_is_two

/-- 不用任何输入的结果：纯算术。 -/
theorem no_input : (2 : Int) + 2 = 4 := by omega

#print axioms needs_input
#print axioms no_input

end Zero.AuditDemo
