# Z002 §6.7 维数账：`dim ker B = dim im B ⟺ m = 4` 的钉死与机器核验

**日期**：2026-10-16　**性质**：核查 + 形式化（承接 Lean 验证线）
**触发**：Z002 当天大改（现 987 行，新增 §6.7「闭合＝运动读法」、§6.8 三条等价形式），
而 ZeroLean 里对应它的 `Z002_Quadform.lean` 只覆盖旧版的三条二次型定理 —— 新写的承重线没有机器核验。

---

## 一、结论（先行）

1. **§6.7 的核验句本身是对的**：$K_m$ 上"圈空间维数 = 运动空间维数"**仅当** $m=4$。
2. **但公式里的对象必须先钉死**：文档写的 `dim ker B`／`dim im B` 只有在
   **`B` 的上边界读法**下才给出这句话。原样的"关联矩阵"读法给不出：
   | 读法 | $\dim\ker B$ | $\dim\operatorname{im}B$ | $m=4$ 时相等？ |
   |:--|:--|:--|:--|
   | $B$ = 关联映射（边↦顶点，流量取边界） | $\binom m2-(m-1)$（**圈空间**） | $m-1$（零和态） | ✅ 相等 |
   | $B$ = 上边界（顶点↦边，$\mathrm df$） | $1$（常数函数） | $\binom m2-1$ | ❌ |
   ⟹ 文档那句要么补上"读法"说明，要么把 $\ker B$ 改写成 **`ker` 上边界**（等价地：圈空间 $= \ker(B^\top\cdot)$）。
3. **数值（本机算，非引用）**：

   | $m$ | 2 | 3 | **4** | 5 | 6 | 7 | 8 |
   |:--|--:|--:|--:|--:|--:|--:|--:|
   | 圈空间 $\binom m2-(m-1)$ | 0 | 1 | **3** | 6 | 10 | 15 | 21 |
   | 运动空间 $m-1$ | 1 | 2 | **3** | 4 | 5 | 6 | 7 |
   | 关系 | < | < | **=** | > | > | > | > |

   与 Z002 §6.7「$m=2$ 给 $0\ne1$、$m=3$ 给 $1\ne2$、$m\ge5$ 一律 $\dim\ker B>\dim\operatorname{im}B$」**逐字对上**。

---

## 二、机器核验（Lean）

**文件**：`ZeroLean/wip/Z002_Dimension.lean`（设计稿，待并入 `Zero/Z002_Dimension.lean`）
**内容**：关联矩阵 $B$、关联映射 $T$、圈空间、零和超平面，以及
`z002_dim_eq_iff : finrank (圈空间) = finrank (range T) ↔ m = 4`。

**已跑通（无 sorry，`lake env lean` 退出码 0）**：`wip/Z002_Dimension.lean`
- `z002B`／`z002T`／`z002Closure`／`z002ZeroSum` 四个定义（`T` 用**显式 λ**，绕开矩阵取向记法）；
- `z002B_sum_col`：**列和为零**（每条边两端各 $\mp1$）；
- `z002T_range_le_zeroSum`：**像 ⊆ 运动空间**（列和为零 + `Finset.sum_comm` 交换求和）。

**尚欠（草稿在 `wip/Z002_Dimension_probe.lean`，是技术性簿记、不是数学缺口）**：
1. `z002ZeroSum_finrank`：`finrank H_Q = card V − 1`
   （路子已通：`Module.Dual.finrank_ker_add_one_of_ne_zero` + `Module.finrank_fintype_fun_eq_card`）；
2. `z002T_ker_eq_zero`：连通图上 `T y = 0 ⟹ y = 0`
   （路子：以 `v` 为超端点的边上流量和为零 ⟹ 每条边为零；欠 `Sym2` 端点见证的簿记）；
3. 由 1+2 经秩-零化度得 `z002_finrank_range`／`z002_finrank_closure`；
4. `z002_card_edge`（$=\binom m2$）与主定理 `z002_dim_eq_iff`（算术核已写好，等 3 落地后接上）。

**另外两条可独立核验的算术（不依赖线性代数）**：
$\binom m2-(m-1)=m-1\iff m=4$ 与 $\binom{m-2}2=1\iff m=4$ —— 上面表里的数字就是它们。

**取向上踩到的坑（值得记）**：本库 `Matrix.mulVecLin` 的取向由探针定死 ——
`M : Matrix m n R ⟹ M.mulVecLin : (n → R) →ₗ (m → R)`（与 `mulVec_apply_eq_sum` 一致）；
中间我被 `#check` 的显示形式误导过一次，最后改用**显式 λ 定义** `T` 绕开记法歧义。

---

## 三、这条线对 §6.8 的意义

§6.8 的三条等价形式里，**（一）配对唯一性** $\binom{m-2}2=1\iff m=4$ 是**纯算术**，
可以独立机器核验（不需线性代数）：

| $m$ | 2 | 3 | **4** | 5 | 6 | 7 |
|:--|--:|--:|--:|--:|--:|--:|
| 每边的不相邻边数 $\binom{m-2}2$ | 0 | 0 | **1** | 3 | 6 | 10 |

⟹ 与 §6.7 的维数账是**同一条算式的两种包装**（$\binom m2-(m-1)=m-1$ 与 $\binom{m-2}2=1$ 同解）。

**待办**：把 `ZeroLean/wip/Z002_Dimension.lean` 收成正式文件（补 4 处欠账、挂 `Zero.lean`、
`zero_lean` 跑守则），并在 Z002 §6.7 就地补一句"读法"限定（不改结论）。

**记忆**：本次结论本想写进项目记忆（`~/.dsh/memory/projects/--Users-oygb-Downloads-Zero--.md`），
但写入被文件沙箱挡住（`Operation not permitted`，工作区外/策略），故只落在本报告；
**下次可在允许的会话里补记一条**（要点：读法钉死 + 圈空间 vs 常数核 + `mulVecLin` 取向坑）。
