# ZeroLean 形式化待办账（LEAN BACKLOG）

**用途**：与 `FORMALIZATION_COVERAGE.md`（覆盖账：什么**可**形式化）配套 ——
这里管**已经试过但没做完**的项。每项写成可重启的四元组：

> **命题** ｜ **卡在哪** ｜ **需要什么** ｜ **下次从哪开始**

**纪律（2026-10-14 审计教训）**：不在 Lean 报错循环里把命题削到能过。
若失败的性质是"缺构造／缺定理／命题要回原文核对"，就记录到这里后推；
硬啃的必然结局是空心定理（`2*3=6`、`True := trivial`、`P↔P` 那类）。

---

## 一、wip/ 草稿（命题真、实现未完成）

见 `ZeroLean/wip/README.md`，三份重写要点：

| 项 | 命题 | 卡在哪 | 需要什么 | 下次从哪开始 |
|---|---|---|---|---|
| Z022 饱和判据 | `C(T,⌊T/2⌋) ≤ 2(T+1) ⟺ T ≤ 5` | ~~`Nat.choose_le_choose` 参数名~~ | — | **部分完成**（2026-10-14）：`Zero/Z022_Saturation.lean` 3 条 —— 中心二项式**单调**（Pascal 递推，真有内容）、`T ≤ 5` 成立、`6 ≤ T ≤ 20` 不成立（逐值 decide）。**未证**：一般的 `T ≥ 6 ⟹ 不成立`，需 (a) 同 n 下 `C(n,3) ≤ C(n,⌊n/2⌋)`（Mathlib 无 k 方向单调件）、(b) `2(T+1) < C(T,3)`。坑记：`Nat.even_or_odd` 的偶数支给**加法形式** `m+1+(m+1)`、奇数支给 `2*m+1`，且 `2*x ≠ x+x`（非 defeq）—— 两边陈述要分别对齐 |
| Z047 对合 | 有限集上无不动点对合 ⟹ 基数为偶 | 归纳（挖掉 `{x, σx}` 后递归）的目标管理有类型错 | `Finset` 的 `erase` 基数引理，或改用 `Fintype` 上的强归纳 | 先证 `s.card = s'.card + 2`（`s' = s.erase x |>.erase (σ x)`），再 `ih` 递归 |
| Z052 gcd | `b` 偶 ⟹ `gcd(b+1, b−1) = 1` | ~~原写法 `simpa [Nat.gcd_rec] using rfl` 不能证那个等式~~ | — | ✅ **已完成**（2026-10-14）：走"公因数整除两数之差 2，而 `b−1` 奇 ⟹ 只能是 1"，并**手工构造整除**避开 `dvd_sub` 的 API 细节；已接入构建 |

---

## 二、需要新构造／新定理的项（不是打错字，是数学缺口）

| 项 | 命题 | 卡在哪 | 需要什么 | 下次从哪开始 |
|---|---|---|---|---|
| Z065 Heawood 图 | 14 顶点、3-正则、21 条边，圈秩 = 8 | `SimpleGraph` 结构字段的 `symm`/`loopless` 是隐式 binder，`intro i j h` 失败；`edgeFinset` 的 `Fintype`/`DecidableRel` 实例要显式 | 用 `SimpleGraph.fromRel`（自动对称化）或补 `classical` + 显式实例 | 改用 `fromRel (fun i j => j = i+1 ∨ j = i+5)`，loopless 用 `(1 : ZMod 14) ≠ 0` 与 `(5 : ZMod 14) ≠ 0`；边数用 `decide` |
| Z096 A₂ 根系统 | `R = {±1, ±ω, ±ω²}` 的 Cartan 矩阵 `[[2,−1],[−1,2]]` | ~~`Complex.normSq` 不是 `(z*conj z).re`；`Complex.sin` 干扰 simp；`ω` 的 `re/im` 展开要手写~~ | — | ✅ **已完成**（2026-10-14）：`Zero/Z096_Cartan.lean` 9 条定理。关键坑：①`Complex.conj` 这个**函数名**在当前 Mathlib 里不存在（共轭要用 `star`），但引理名仍是 `Complex.conj_re/im`；②`(star z).re = z.re`、`(star z).im = -z.im` 直接 `rfl` 给（别等 simp）；③`Complex.normSq` 的定义不是 `(z*conj z).re`，别走那条路 |
| Z080 一般 ballot | `m_j = C(L,L/2−j) − C(L,L/2−j−1)` 对**任意 L** | Mathlib 无 ballot 定理；反射双射需要"首次达到 j+1 的位�置"（`Nat.find`）与翻转的互逆证明 | 代数核心**已证**（`z080_flip_prefix_sum`）；缺"首次达到位置"的良定义 + 双射互逆 | 定义 `firstHit w j := Nat.find (...)`，证翻转后 `M⁺` 首次达到位置不变，再建 `Finset` 双射 |
| Z132 2cos 族判定 | `3/2` 与 `8/3` **不是** `2cos(2πm/n)` 的值 | 需要代数数论（Niven 型：`2cos(2πm/n)` 的有理值只有 0, ±1, ±2） | Mathlib 的 `Niven` 或代数整数工具 | 先查 `Mathlib/NumberTheory/Niven.lean`（已存在！）能否直接推有理值分类 |
| Z021 Γ=K₄ 唯一性 | 4 顶点 3-正则连通图 ≅ K₄（`SimpleGraph.Iso` 版） | 与 Z089 的 `Finset` 版本重复度高；`Iso` 的构造成本 | `SimpleGraph.Iso` + 度数与完全图的刻画 | 先看 Z089 的 `z089_four_vertex_three_regular_is_complete` 能否直接升级到 `SimpleGraph` 层 |

---

## 三、Mathlib 缺口（先立陈述，别硬做）

覆盖账判定 4 已登记：Tomita–Takesaki／有限维因子分类（Z005/Z034/Z068/Z085 定理 34）、
有限群表示论（Z090/Z011/Z093/Z099）、Stirling 型渐近（Z026/Z029）、Lorentz 表示论（Z028 墙 IV）、
Clebsch–Gordan（Z009）。**这些项的正确做法是写"降级版陈述"**（例如 Z068 的模流平凡性已经退化为
交换性引理并证掉），而不是在 Lean 里空等 Mathlib。

---

## 四、来自 ZeroGPU/doc 的待形式化裁定（新目标源，2026-10-14 由用户提示）

`ZeroGPU/README.md` 的"一页结论"表里有一批**已裁决但从未形式化**的结论，每一条都带 `doc/*.md`
推导记录 —— 这是比"从覆盖账里挑"更省力的目标源（裁定已完成，只差机器检查）：

| 裁定（ZeroGPU/doc） | 数学内容 | 形式化起点 |
|---|---|---|
| 饱和判据／`T` 的唯一解（`saturation_law_rediscovered.md`、`algebraic_capacity.md`） | `C(T,⌊T/2⌋) ≤ 2(T+1) ⟺ T ≤ 5`；再用**词数** `\|W_T\| ≤ 2(T+1)` 收到 `{2,4}`、四维窗再收到 `{4}` | 先做 `wip/Z022` 的饱和判据（命题已在账上），再补"词数版"与四维窗两步 |
| `Γ = K₄` 的来源（`gamma_origin_verdict.md`） | 4 顶点、度 3、6 边、3 对边对，**由链条推出而非申报** | 已部分形式化（`z089_four_vertex_three_regular_is_complete`）→ 可升级为 `SimpleGraph.Iso` 版 |
| Heawood 图（`scale_alignment_verdict.md`） | 14 顶点、3-正则、围长 6、二部；3-边着色给 3 个轴 | 就是 backlog §二里的 Z065 显式构造 |
| 耦合构造与映射（`coupling_verdict.md`、`coupling_map.md`） | 正则 `D`-单纯形；`ρ = −1/D` 在等角族内唯一；夹角 = 四面体角；空间维数 = 3 | 内积/夹角计算 —— **Z096 的 Cartan 那条路已通**（`innerC` + `exp_ofReal_mul_I`） |
| 维数字典（`dict_verdict.md`） | `D = m − 1` 是代数唯一解；四维窗 `(1/2, 3/5)` | 与 Z021 的 `C(D,2) = D` 相邻，可直接扩 |
| 长时演化（`longterm_eval.md`） | **非周期演化被公理排除**（定理） | 与 Z020 的鸽巢原理同源，可直接形式化 |

**用法**：下一轮起，队列优先从这张表取（裁定已定、证据在 doc、只差 kernel 检查）；
每完成一条，回来把该行标 ✅ 并记下卡点（供后来者复用）。

---

## 五、已完成的（审计清理，2026-10-14）

- 16 条重言式删除、9 条纯算术改写为真命题、`sorry` 文件删除；
- Z020/Z024/Z076 重写并接入构建；**Z052 当场修好并接入**（判据验证：局部实现错 → 当场修）；
- 终态：39 个文件、**153 条定理**、守则全过。
