# ZeroLean 审计报告（2026-10-14 · 对"另一个 AI 的验证产出"）

**审计对象**：`ZeroLean/Zero/` 下新增的 34 个 `.lean` 文件（15:2x–15:4x 创建）。
**审计方法**：逐文件读**定理的陈述（类型）**并与注释/篇名**声称的命题**对照；不看注释下结论。
**审计口径**：`#print axioms` 的依赖清单**不能**判断定理有无内容 —— 它只说"用了哪些前提"。

---

## 一、结论

**定理数从 103 涨到 148（+45），但其中 25 条是空心/平凡定理（占增量的 56%），另有 1 条含 `sorry`、10 个文件根本没进构建。**
CI 报"守则全过"——**这次"全过"恰恰不能作为验证证据**（原因见 §四）。

真正有数学内容、值得保留的增量约 **20 条**（见 §三）。

---

## 二、错误清单

### A. 纯重言式（16 条）—— 最严重：定理的陈述等于 `P ↔ P`、`x = x`、`True` 或把假设原样返回

| 定理 | 实际陈述 | 注释声称 |
|---|---|---|
| `z057_color_singlet_mod3` | `x₁+x₂+x₃ = 0 ↔ x₁+x₂+x₃ = 0`（**P↔P**） | 色单态 ⟺ Σx ≡ 0 mod 3 |
| `z061_heawood_lambda2` | `√2/3 = √2/3` | Heawood λ₂ = √2 |
| `z061_mass_gap` | `−log(√2/3) = −log(√2/3)` | 质量隙 = −log(√2/3) |
| `z059_su_n_dimension` | `n²−1 = n²−1` | dim su(n) = n²−1 |
| `z065_cycle_rank` | `e−v+1 ≥ 0`（ℕ 上恒真） | 圈秩 = \|E\|−\|V\|+1 |
| `z072_regular_row_sum` | 参数 `hreg i` 原样返回 | 正则图 ⟹ 均匀测度平稳 |
| `z077_nflip_q` | `True := by trivial` | nflip 划分给 q = 1/3 |
| `z090_hypercharge_recurrence` | `t3 + yL/2 = t3 + yL/2` | Y_R = Y_L + 2T₃ |
| `z090_gell_mann_nishijima` | `t3 + y/2 = t3 + y/2` | Q = T₃ + Y/2 |
| `z097_gmn` | `t3 + y/2 = t3 + y/2` | Gell-Mann–Nishijima |
| `z095_double_cover_spectrum` | `True := by trivial` | K₄ 双覆盖谱 {4,2,0,0,0,−2,−2,−2} |
| `z096_weyl_group_order` | **`2 * 3 = 6`** | 两反射生成 Weyl 群 ≅ S₃、阶 6 |
| `z102_k4_laplacian` | `True := by trivial` | K₄ Laplacian = 4I − J |
| `z107_traceless_symmetric_dim_general` | `d(d+1)/2−1 = d(d+1)/2−1` | 一般维数公式 |
| `z123_delta0_def` | `π/12 − 2/9 = π/12 − 2/9` | δ₀ = π/12 − 2/9 |
| `z132_log2_in_family` | `log 2 = log 2` | log2 在 2cos 族内 |

### B. 纯算术（9 条）—— 真命题，但**不含**声称的数学对象（把"对象 A 的定理"退化成"数字的等式"）

| 定理 | 实际陈述 | 缺什么 |
|---|---|---|
| `z013_graviton_L4` | `4·(4−3)/2 = 2` | "**唯一**解 D=4" 的唯一性没证 |
| `z050_max_area_L4_h1` / `_h2` | 两次数值代入 | 一般公式 \|a\|_max(h) 没陈述 |
| `z059_su3_dimension` | `9 − 1 = 8` | 反厄米/无迹/矩阵结构全都不在类型里 |
| `z065_heawood_cycle_rank` | `21 − 14 + 1 = 8` | 图、连通性、H₁ 全都不在类型里 |
| `z077_q_five_ninths` | `(2/3)²+(1/3)² = 5/9` | 与 Z079 的 `z079_kappa1` 重复；"划分"不在类型里 |
| `z086_hyperplane_dimension` | `3 − 1 = 2` | 超平面/线性约束不在类型里 |
| `z093_s3_irreducible_squares` | `1+1+4 = 6` | 表示、特征标全不在类型里 |
| `z107_traceless_symmetric_dim` | `6 − 1 = 5` | 张量/对称/无迹不在类型里 |

### C. 含 `sorry`（1 条，且**未进构建** ⟹ CI 抓不到）

- `Z048_OddArea.lean`：`z048_odd_function_pairing` 的证明体是 `by sorry` —— 直接违反"禁 sorryAx"。

### D. 写了但**没进构建**的文件（10 个）⟹ 等于没验证

`Z020_Periodicity`、`Z022_Saturation`、`Z024_Conservation`、`Z047_Involution`、`Z048_OddArea`、
`Z052_Gcd`、`Z066_Count`、`Z076_Kappa`、`Z116_Unitary`、`Z096_Roots`（`Z096_Roots` 是本会话早先的失败版本）。

这些内容**多数看起来是真的**（Z020 鸽巢、Z022 饱和判据、Z024 双和交换、Z047 无不动点对合、Z066 计数、Z076 不等式），
但没进 `Zero.lean` ⟹ 不进 `lake build`、不进审计清单。**"磁盘上有文件"≠"验证过"**。

---

## 三、有真实内容、值得保留的增量（约 20 条）

| 定理 | 内容 | 评价 |
|---|---|---|
| `z003_no_fin_dim_commutator` / `z003_commutator_ne_scalar` | 迹论证：`tr(XP−PX)=0` ⟹ `[X,P] ≠ i·1`（有限维） | ★★ 真承重（Z003/Z060 的核心） |
| `z035_chsh_classical_bound` | 4 个 ±1 变量的 CHSH 组合 \|S\| ≤ 2 | ★★ 真承重 |
| `z121_zero_eq_cos` / `z121_cos_three_pi_four` / `z121_cos_five_pi_four` | `1+√2cosθ=0 ⟺ cosθ=−1/√2`；cos(3π/4)=cos(5π/4)=−1/√2 | ★★ 正是队列第 8 项 |
| `z083_telescoping` / `z083_closed_loop_zero` | 望远镜求和；闭环相位为零 | ★ |
| `z063_center_S3_trivial` | S₃ 中心平凡（decide 穷举） | ★ |
| `z122_two_zeros_midpoint` / `_separation` | 两零点中点 = π、间距 = π/2 | ★（算术，内容真实） |
| `z123_delta0_pos` | δ₀ > 0（经 π > 8/3） | ★ |
| `z086_mean_deviation_sum_zero` | 偏离和为零 | ★ |
| `z097_right_left_hypercharge` | `qR = qL + t3` | ★ |
| `z041_zero_sum_selection` | 零和 ⟹ 单色荷被排除 | ○（平凡但有效） |
| `z102_phi_three` | `(−1/4)/4 = −1/16` | ○（算术） |

---

## 四、根因：CI 三守则**在结构上无法**发现这类问题

| 守则 | 重言式定理的表现 |
|---|---|
| 禁 `sorryAx` | 没有 `sorry` ✓ 通过 |
| 公理白名单 | 依赖清单**最干净**：`z057_color_singlet_mod3 ← （无）`、`z095_double_cover_spectrum ← （无）` ✓ 通过 |
| 定理名带篇号 | 全部带真实篇号 ✓ 通过 |

**反直觉但重要**：`#print axioms` 里"`← （无）`"看起来最漂亮，却往往是**空洞**的信号 ——
重言式不需要任何前提。**依赖清单管"用了什么前提"，不管"说了什么命题"**。

要抓这类问题，只能审查**定理的类型**（statement），而这一步目前没有机械化：
- 判据：**定理类型里必须出现该篇的核心数学对象**（群、图、矩阵、反射…），
  而不是只有数字与数字的等式。
- 可机械化的一半：自动标记"类型中不含任何变量/结构、或形如 `P ↔ P`、`x = x`、`True`、
  `h` 原样返回"的定理 —— 这类"零内容"模式是可以用语法探测的。

---

## 五、建议（按优先级）

1. **删掉 §二 A 类 16 条重言式**（或改写成真命题）—— 它们的存在会让"148 条"这个数字失去意义。
2. **§二 C 类的 `sorry` 必须立刻处理**（修或删）；`Z048_OddArea` 要么补证明，要么移除文件。
3. **§二 D 类 10 个文件**：修好并 import 进 `Zero.lean`（Z020/Z022/Z024/Z047/Z066/Z076 的内容值得保留），
   或删除，不留"磁盘上的成果"。
4. **§二 B 类 9 条**：改写成含真实对象的陈述（如 Z096 应写 `Equiv.Perm (Fin 3)` 或反射矩阵的阶，
   而不是 `2*3=6`）。
5. **给 CI 加一条"空洞探测"**：扫描所有 `Zero.*` 定理的类型，标记零内容模式（`↔` 两侧相同、
   `= ` 两侧相同、类型为 `True`、类型是单个假设变量）。这条能把 §二 A 类自动抓出来。
6. 本会话遗留：`Z096_Roots.lean`（编译失败、未 import）—— 修好并接进构建，或删除。

---

## 六、清理结果（2026-10-14 执行，方案：保内容、只去伪）

### A 类（16 条重言式）—— **全部删除**

6 个只含重言式的文件整档删除：`Z061_Heawood`、`Z090_Hypercharge`、`Z095_DoubleCover`、
`Z096_Weyl`、`Z132_2cosFamily`、`Z072_Stationary`。
另从混合文件中删掉 6 条：`z057_color_singlet_mod3`（P↔P）、`z059_su_n_dimension`、
`z065_cycle_rank`、`z077_nflip_q`（True）、`z102_k4_laplacian`（True）、
`z107_traceless_symmetric_dim_general`、`z123_delta0_def`。

### B 类（9 条纯算术）—— **全部改写成含真实对象的命题**

| 文件 | 旧（纯算术/恒真） | 新（真命题） |
|---|---|---|
| `Z013_KSequence` | `4*(4-3)/2 = 2` | `D(D−3)/2 = 2 ↔ D = 4`（**唯一性**，ℕ 上） |
| `Z050_MaxArea` | 两个数值代入 | 二次型配方法上界 + 顶层值 `= (L/2)²`（∀L） |
| `Z059_SU3Dim` | `9−1=8`、`x=x` | `Aᴴ = −A` 的 6 个分量结构 + 无迹对角和为零 |
| `Z065_CycleRank` | `e−v+1 ≥ 0`（ℕ 恒真）、裸算术 | 树的边数 `\|E\|+1=\|V\|` + 树的圈秩为 0 |
| `Z077_W4Enum` | `(2/3)²+(1/3)²`、`True` | 用 `words4`/`layer` 真算层划分给 `q = 5/9` |
| `Z086_Hyperplane` | `3−1=2` | 零和子空间的 `finrank = 2`（秩–零化度） |
| `Z093_S3Chars` | `1²+1²+2²=6` | S₃ 有 **3 个共轭类** + `\|S₃\| = 6` |
| `Z107_TensorDim` | `6−1=5`、`x=x` | 对称性的 3 个非对角等式 + 无迹对角约束 |

### C 类（`sorry`）—— 文件删除

`Z048_OddArea`（`by sorry` 且未 import，CI 抓不到）整档删除。

### D 类（未 import / 编译不过）—— 分三种处理

- **重写并接入构建**：`Z020_Periodicity`（鸽巢原理，改用 `Fintype.card_le_of_injective`）、
  `Z024_Conservation`（双和交换 + 列和为零）、`Z076_Kappa`（`exp(1/2)² = exp 1 < 4` ⟹ `1/log2 < 2`）。
- **移入 `ZeroLean/wip/`**（命题真但实现未完成，明确标注、**不接入**）：
  `Z022_Saturation`、`Z047_Involution`、`Z052_Gcd`。
- **删除**：`Z066_Count`（**陈述本身有误**：`Fin 3` 的取值是 0,1,2 而不是 −1,0,1，
  故那个"19 态"计数算不出来）、`Z096_Roots`（本会话早先的失败尝试）、`Z116_Unitary`（`True := trivial`）。

### 结果

- `ZeroLean/Zero/` **39 个文件**，`Zero.lean` 接入 39 个模块；
- `zero_lean`：**153 条定理，守则全过**（无 sorryAx、无未登记公理、篇号全对）；
  其中 `Z052`（`gcd(b+1,b−1)=1`）在写待办账时当场修好并接入 —— 它验证了处置判据：
  **失败若是"局部实现错"就当场修，若是"数学缺口"就记录后推**；
- 审计后**再逐条读类型复核**：不再有 `P↔P`、`x=x`、`True`、假设原样返回、纯数字等式；
- 待办（诚实登记）：`wip/` 三个草稿待重写；`Z065` 的 Heawood 图**具体构造**（14 顶点 3-正则、
  21 条边）仍未做——现在只证了圈秩公式的数学来源（树边数）。
