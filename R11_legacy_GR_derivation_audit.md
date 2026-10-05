# R11 · 旧理论与历史 GR 推导审计

**日期**：2026-10-02  
**范围**：`cosmos-construct`、`modular-equilibrium`，以及 `lh` 中早期 `G*`／`Z*`／`D2xx` 的 GR 推导记录。  
**性质**：历史与外部路线审计。只回收可审查的证明部件，不把旧结论计入当前项目进度。  
**当前状态入口**：[`STATUS.md`](STATUS.md)。  
**核验**：[`R11_check.py`](R11_check.py)。

$$

\begin{aligned}
&\text{旧理论没有无条件导出四维 GR。}\\
&\text{最强遗留物是有限维模代数、K32 类时锥代数引理、面积系数预算和条件化几何组装。}\\
&\text{这些材料只能补 R8／R10 的条件桥，不能关闭 J1／J5，也不计入 O1–O5。}
\end{aligned}
$$

---

## §0 审计结论

这次审计回答三个问题：

1. 旧理论有没有已经完成的 4D GR 推导？
2. 哪些证明部件可以迁移到当前 `R8`／`R10` 条件桥？
3. 哪些旧“已完成”说法必须降级为历史？

第一问的答案是**没有**。旧 `cosmos-construct` 的最强结果是 2D CFT 条件推导和 4D 条件结构闭合；旧 `modular-equilibrium` 把条件恢复层拆得更严格，但其恢复账中的连续网、几何 boost、四维 Lorentz、面积／`G`、全域平衡和状态到几何仍全部未建。

第二问的答案是有若干独立部件可迁移，但没有一件能单独补上当前决定性的 L1 或 L5。

第三问的答案集中在 §8。最重要的三条是：

- `K17` 的“BW boost 已确认”已被 `K19` 撤回；
- `GRCOMPLETE` 的“完整非线性 EFE”是条件结构推导，不是从旧公理无条件导出；
- 旧谱作用量／内空间路线是独立唯象路线，不能并入 Zero 的 `Z/G` 进度。

---

## §1 三条旧路线必须分开记账

| 路线 | 旧目录 | 真正得到什么 | 当前判决 | 能否关闭 R8／R10 新缺口 |
|:--|:--|:--|:--|:--|
| 谱作用量／内空间几何 | `cosmos-construct` | `a₂ ⇒ EH`、局部 PPN、条件面积熵 | **历史／外部基准**：内空间是唯象 ansatz，KK 约化未做 | 不能；没有连续区域网、模流桥或物理站点映射 |
| 纠缠平衡 Jacobson／FAZ | `cosmos-construct` | 2D CFT 条件推导；4D 条件结构闭合；K22 负结果；K32 代数引理 | **条件证成＋已撤回项**：全族与面积系数仍是输入 | 只能迁移代数与失败边界；不能关闭 J1／J5 |
| `U1–U4` 条件恢复 | `modular-equilibrium` | 有限维模代数、精确第一定律前身、条件面积系数和条件几何组装 | **严格前身**：恢复链条仍有多项未建 | 提供检查表和缺口分解；不能关闭 J1／J5 |

这三条路线不能相加：

$$
\text{谱作用量路线}
\not\equiv
\text{Jacobson／FAZ 条件路线}
\not\equiv
\text{Zero 主 }Z/G\text{ 路线}.
$$

---

## §2 `cosmos-construct`：谱作用量／内空间几何

### 2.1 旧结果

该路线给出一个可复现的唯象框架：

1. 选定内空间陪集几何；
2. 计算谱作用量；
3. 由热核系数 $a\_2$ 得到 Einstein–Hilbert 项；
4. 检查局部牛顿势与 PPN。

来源：

- [`GR_LIMIT.md`](../cosmos-construct/GR_LIMIT.md)；
- [`GR_FROM_GEOMETRY.md`](../cosmos-construct/GR_FROM_GEOMETRY.md)；
- [`submission/GR_FROM_GEOMETRY.md`](../cosmos-construct/submission/GR_FROM_GEOMETRY.md)。

### 2.2 必登记的输入与缺口

| 项 | 旧文件承认的状态 | 当前判决 |
|:--|:--|:--|
| 内空间陪集 | **phenomenological ansatz**，不是从零和原语导出 | 历史／外部输入 |
| `a₂ ⇒ EH` | 标准热核结果；前提是先接受谱作用量 | 条件证成，不算 Zero 新定理 |
| 10D → 4D KK 约化 | 从未计算 | **开放**；不能声称零和宇宙完成该约化 |
| KK 塔与观测 | 有严重冲突 | **排除将旧谱路线直接并入当前主链** |
| Lorentzian 底空间 | 是输入 | 不能补 R8 的连续 Lorentzian 网 |
| 体积模与零模 | 体积模无有界极小；零模靠公设排除 | 不能当背景稳定性定理 |

### 2.3 判决

这条路线可以作为独立的“谱作用量唯象交叉检查”，但不能回答：

1. 类为何对应物理站点；
2. 四维为何被选中；
3. 连续局域代数网如何出现；
4. 小球模流如何等于几何 boost；
5. 面积密度如何从 Zero 状态普适化。

因此它**不并入** `O1–O5`，也不改变 `STATUS.md` 的主判定。

---

## §3 `cosmos-construct`：Jacobson／FAZ 条件路线

### 3.1 真正有用的负结果

`K22` 证明单个球的第一定律只给球平均的 `l=0` 投影，不能单独给出 4D 张量方程。它防止把球面平均误读为“完整张量源”。

来源：[`NEW_UNIVERSE_DESIGN_GR10.md`](../cosmos-construct/NEW_UNIVERSE_DESIGN_GR10.md)。

### 3.2 BW boost 的历史张力

早期 `K17` 声称格点上确认了 Bisognano–Wichmann boost 结构。随后 `K19` 在同一项目内订正：

- 实际测到的是有隙链的局域均匀模 Hamiltonian；
- 它不是 BW boost；
- BW boost 是无隙真空的性质；
- 必须走解析方法或更换正则化，不能在原模型上继续重复。

来源：

- [`NEW_UNIVERSE_DESIGN_GRFINAL.md`](../cosmos-construct/NEW_UNIVERSE_DESIGN_GRFINAL.md)；
- [`NEW_UNIVERSE_DESIGN_GR7.md`](../cosmos-construct/NEW_UNIVERSE_DESIGN_GR7.md)。

当前处置：**`K17` 的强声明标为已撤回；L1 仍开放。**

### 3.3 全张量条件推导

旧终局稿的逻辑是：

1. `K30` 要求所有球心、半径和类时参考系的平衡；
2. `K31` 引入 Jacobson 精确面积展开和外部系数；
3. `K32` 给出代数引理：若对称张量 $M$ 对所有类时单位向量 $u$ 满足

$$
M_{ab}u^a u^b=0,
$$

则 $M=0$；

4. 因而把局部分量关系升为全张量方程。

来源：

- [`NEW_UNIVERSE_DESIGN_GRSTUCK.md`](../cosmos-construct/NEW_UNIVERSE_DESIGN_GRSTUCK.md)；
- [`NEW_UNIVERSE_DESIGN_GRDERIVE.md`](../cosmos-construct/NEW_UNIVERSE_DESIGN_GRDERIVE.md)；
- [`NEW_UNIVERSE_DESIGN_GRCOMPLETE.md`](../cosmos-construct/NEW_UNIVERSE_DESIGN_GRCOMPLETE.md)。

### 3.4 为什么不能称为“4D GR 已导出”

旧终局稿自己保留的关键条件仍在：

| 条件 | 状态 | 对当前结论的影响 |
|:--|:--|:--|
| 所有球 × 所有类时参考系 | `K30` 只定位族的要求，没有构造全族 | 仍在 R8 的 J3／J4 区域 |
| Jacobson 面积系数 | 外部输入或面积律一致性识别 | 不能算 Zero 已导出 $G$ |
| 连续 Lorentzian 区域代数 | 旧模型没有一般连续网 | 连续局域化仍开放 |
| $K\_B\to2\pi B\_B$ | 没有算子级极限定理 | **L1 仍开放** |
| 低维相关算符控制 | 未做 | **L5 仍开放** |
| $\Lambda$ | 仍为输入 | 不改变 E4 的处置 |

所以正确措辞是：

> 在 Jacobson 系数、固定体积平衡、全球球族和所有观察者均作前提时，旧路线给出 4D EFE 的条件结构闭合。

不能写成：

> 旧理论已经从自身公理无条件导出完整非线性 4D GR。

---

## §4 `modular-equilibrium`：最严格的条件恢复前身

### 4.1 公理与恢复层

`U1–U4` 只规定上游零和结构；几何、连续网、模流几何极限、面积熵和 GR 都在恢复层。其全库恢复账中的 `G1–G8` 仍全部记“未建”，包括：

1. 连续局域代数网；
2. Type III 结构；
3. Hadamard 态；
4. 几何 boost；
5. 四维 Lorentz 结构；
6. 面积密度与 `G`；
7. 全域平衡；
8. 状态到几何映射。

来源：[`AXIOMS.md`](../modular-equilibrium/AXIOMS.md)。

### 4.2 已有严格部件的层级

| 部件 | 旧文件 | 当前可回收形式 |
|:--|:--|:--|
| 有限维模第一定律 | [`D23_modular_first_law.md`](../modular-equilibrium/derivations/D23_modular_first_law.md) | `R8.1` 的独立历史确认 |
| 条件 EFE 桥 | [`D24_equilibrium_to_einstein.md`](../modular-equilibrium/derivations/D24_equilibrium_to_einstein.md) | 条件桥失败树；明确 C1–C4 是输入 |
| 正能平移不等于 boost | `D42–D44` | 可直接补强 L1 的排除边界 |
| 面积系数预算 | [`D152_entanglement_area_coefficient_and_G.md`](../modular-equilibrium/derivations/D152_entanglement_area_coefficient_and_G.md) | $G=1/(4\kappa c I\_1)$ 的缺口分解 |
| 张量 vs 直接和 | [`D225_tensor_vs_direct_sum_factorization_gap.md`](../modular-equilibrium/derivations/D225_tensor_vs_direct_sum_factorization_gap.md) | 排除把代数因子直接当空间站点 |
| 模密度剖面缺口 | `D231`、`D233`、`D234` | 候选核与“支持不能唯一选 boost”的边界 |
| ADM／Dirichlet 条件几何 | `D253–D256` | 条件恢复空间几何；仍缺时钟、lapse、shift 与细化 |
| 电阻度量非局域性 | [`D257_resistance_metric_fixed_point_and_locality_gap.md`](../modular-equilibrium/derivations/D257_resistance_metric_fixed_point_and_locality_gap.md) | 排除把它直接当局部度规极限 |
| D259 条件汇流 | [`D259_z_layer_condition_dictionary_and_four_dimensional_confluence.md`](../modular-equilibrium/derivations/D259_z_layer_condition_dictionary_and_four_dimensional_confluence.md) | 五个条件通道仍是输入 |

### 4.3 关键开放项

`modular-equilibrium` 不能补上当前以下缺口：

1. **连续四维 Type III 局域代数网**：只有抽象因果网和有向归纳网，没有连续 QFT 拓扑；
2. **几何模流极限**：正能平移不能推出 Lorentz boost；平移与代数交换；
3. **普适面积密度**：支持结构本身不给面积律；加入模边界条件后仍依赖谱隙、状态类和几何定标；
4. **固定体积平衡**：`D24` 直接把它作为外部几何引理；
5. **低维相关算符污染**：旧库没有 $R^{2\Delta}$ 或 $\Delta\le d/2$ 首阶污染分析；
6. **Cao–Carroll 的 Radon 反演**：条件 ADM 组装不是从互信息切割数据反演空间度规；
7. **四维Lorentzian 动力学**：旧库只条件组装几何，不给连续动力学或 $D=4$ 的物理选择。

因此它不是“已经被更严格版本取代的历史错误”，而是**当前条件桥的有限维前身和失败边界总表**。

---

## §5 `lh` 早期 `G*`／`Z*` 的迁移结果

早期 `lh` 推导并不是全部作废。已有两项已进入当前主线：

1. `R1`／`R6`／`R7` 把 E1 的连续极限压成 H1–H7；
2. `R8.1` 把有限维第一定律升级为任意有限维忠实态的精确熵差恒等式。

但旧文档中仍有多处写作时点冲突，必须在引用时按当前状态重读：

| 旧说法 | 旧出处 | 当前处置 |
|:--|:--|:--|
| 层结构对 GR 帮助很小 | [`G21`](G21_do_the_layers_help_derive_GR.md) | `G22` 更正：层结构给局域性与时间定向，不给度规与号差 |
| “唯一堵点是 Z0③（历史标号 A3）” | [`G55`](G55_dynamics_line_degeneration_to_GR.md) | `G56`／`G59` 撤回；I7 不承重，精确锥仍开放 |
| $m\_{\rm stag}=\Theta/(2L)$ | [`G78`](G78_area_law_in_3d.md) | `G77` 的口径更正为 $m\_{\rm stag}=|\Theta|/L$；数字须用 $0.23534171$ |
| $\xi\simeq1.078L$ 与缺口 $8.498$ | `G77`／`G78` | `R8-L2` 修正：漏了 $v\_F$，物理关联长度约 $2.157L$；$8.498$ 是有限窗口算术差距 |
| $1.14$ 与 $m=1.75$ | `G80`／`G83` | 已撤回 |
| I2a 是唯一真瓶颈 | [`Z2`](Z2_zero_to_gr_direct_route.md) | 当前视为 E1 条件链；无条件 Γ-收敛仍未证 |
| 数学缺口 0、承重 no-go 0 | [`Z6`](Z6_stall_autopsy_and_released_ledger.md) | 只在“输入可有”的政策下成立，不是绝对无条件状态 |
| E1 显式字典即完成 | [`Z7`](Z7_embedding_input_explicit_dictionary.md) | 只完成条件显式化；I5b、GDL、站点嵌入仍缺 |

这些不是要删除的旧文本，而是历史记录。当前引用必须以 [`STATUS.md`](STATUS.md) §8 为准。

---

## §6 对 R8 的逐项映射

| R8 缺口 | 旧理论可提供 | 审计判决 |
|:--|:--|:--|
| C1 连续局域完成 | `D35`／`D37` 抽象网；旧谱路线假设 Lorentz 底空间 | 输入或开放；没有连续 Type III 网 |
| C2 几何模流极限 | `D42–D44` 排除正能平移生成 boost；`K17→K19` 撤回 BW 确认 | 开放；旧材料只给负结果 |
| C3 普适面积密度 | `D149–D153`、`K15/K18` | 条件预算；面积系数三因子仍不导出 |
| C4 固定体积平衡 | `D24`、`K30`、`K31` | 外部几何引理或条件输入 |
| J1 几何 boost | `D227`／`D228` 有限生成元族；`D231`／`D234` 候选剖面 | **开放，不能补** |
| L2 面积密度 | `D152`；`K15/K18` | 可迁移为缺口预算；强 gap 无关性不成立 |
| L3 四维 Lorentzian | `D193–D196`、`D253–D259` | 条件组装；不足以为连续四维 Lorentzian 定理 |
| L4 固定体积变分 | `D24`、`K31` | 外部系数与几何引理；不能由旧有限模型推出 |
| L5 低维相关算符污染 | 无直接对应材料 | **开放，旧库不能补** |

结论：

$$

\text{旧理论只把 J1–J5 的失败边界写详细了，没有关掉其中任何决定性缺口。}

$$

---

## §7 对 R10 的逐项映射

| Cao–Carroll 缺口 | 旧理论材料 | 判决 |
|:--|:--|:--|
| CC1 首选局域张量完成 | `D224`／`D225`；`G62` 的有限维量子侧 | 抽象代数因子不等于空间 Hilbert 因子；开放 |
| CC2 近似 RC 与割函数 | `D126`／`D128` 双分划精确分解 | 面积或互信息恒等式不等于任意切割 RC；开放 |
| CC3 跨切割面积–互信息比例 | `D129`／`D135`／`D152` | 规则区域数值与条件容量不足；统一 $\alpha$ 未证 |
| CC4 背景度规与 Radon 反演 | `D257`、`D253–D256` | 条件几何组装不是 Radon 反演；无材料 |
| CC5 MEEC、EFT 与 Rindler 第一定律 | `D23`、`D116`、`R8.1` | 有限维第一定律可用；连续 generic Rindler EFT 未建 |
| CC6 Lorentzian 组装 | `D195`／`D196`、`D253`／`D259` | 仅条件空间几何；连续动力学与 lapse／shift 未解 |
| CC7 局部 Lorentz 完成 | `D44`、`D194` | 仍归结为 L1；$D=4$ 仍为输入 |

---

## §8 可迁移证明部件

以下部件值得保留，但不能扩大它们的结论：

| 部件 | 可迁移内容 | 目标 | 不能声称 |
|:--|:--|:--|:--|
| K32 类时锥代数引理 | 对称张量对所有类时单位向量收缩为零则 $M=0$ | R8 的全张量收缩步骤 | 不能构造全球球族 |
| K30 小球探测 | 球心、半径与偏心探测的参数化 | 局部场重构的辅助恒等式 | 不是“所有球”族的构造 |
| K31 面积展开 | Jacobson 面积展开的内部一致性检查 | R8-L4 的独立复核 | 不导出 Jacobson 系数 |
| K22 负结果 | 单球第一定律只给 $l=0$ | 防止把迹投影当全张量 | 不代替全族条件 |
| D152 面积系数 | $G=1/(4\kappa c I\_1)$ 的缺口预算 | R8-L2 | 三个因子仍未导出 |
| D23 模第一定律 | 有限维精确熵差恒等式的前身 | R8.1 交叉核验 | 不升级到连续 QFT |
| D42–D44 | 正能平移不生成 Lorentz boost | L1 的排除边界 | 不证明 boost 存在 |
| D225 | 抽象代数因子不选择空间张量分解 | CC1 | 不排除特殊动力学选择 |
| D231／D233／D234 | 支持不决定剖面；候选球核 | J1／J5 的目标函数 | 不把候选核写成 Zero 推导 |
| D253–D256 | Dirichlet 张量到条件空间度规 | J3／R7 的交叉检查 | 不关闭时钟、lapse、shift 与细化 |

---

## §9 历史冲突总表

| 编号 | 冲突 | 历史正确读法 |
|--:|:--|:--|
| C1 | `GRFINAL` 声称 BW boost 确认，`GR7` 随后撤回 | `K17` 是历史；`K19` 的负结果有效；L1 开放 |
| C2 | `GRCOMPLETE` 写“完整非线性 EFE”，但把全族与 Jacobson 系数当前提 | 只能记作条件结构闭合 |
| C3 | 谱作用量路线写“恢复 GR”，但内空间是 ansatz，KK 约化未做 | 外部唯象交叉检查，不并入 Zero |
| C4 | `D213` 说当前骨架不能直接推 GR；后文条件链仍存在 | 直接路线 no-go 成立；条件路线不因此自动关闭 |
| C5 | `G78` 仍保留旧 $m\_{\rm stag}=\Theta/(2L)$ 附注 | 当前数字以 `G77`／`R8-L2` 为准 |
| C6 | `G55` 称唯一堵点是 Z0③（历史标号 A3）；`G56`／`G59` 撤回 | Z0③ 唯一性作废；I7 不承重 |
| C7 | `G21` 称层结构帮助很小；`G22` 更正 | 采纳 `G22`：有限局域性与时间定向贡献，但仍不给度规 |
| C8 | `G80`／`G83` 的 `1.14`、`m=1.75` 后被撤回 | 只作历史；不得进入当前数字 |
| C9 | `Z6` 写“数学缺口 0” | 仅是“输入可有”政策下的相对账本 |
| C10 | `lh` 的 `D2xx` 副本与 `modular-equilibrium` 原文不完全相同 | 引用时区分“旧理论原文”和“当前归一化副本” |

---

## §10 这次审计对主路线的影响

### 10.1 已经得到

1. 一份把旧理论分成三条互不相加路线的总表；
2. 一组可迁移的严格或条件部件；
3. 一份把旧“已解决”与“未解决”冲突统一到当前状态的历史账；
4. 对 R8 和 R10 的旧材料映射。

### 10.2 没有改变

`R11` **不改变**以下当前判定：

- 当前仍是条件恢复，不是无条件导出；
- `L1` 仍开放；
- `L5` 仍开放，且应优先于几何 boost；
- `D=4` 仍没有物理独立的选维原则；
- 主 `Z/G` 路线与 `D259` 路线仍没有完整等价；
- `G,\Lambda` 的绝对单位仍是 E4；
- Cao–Carroll 只到弱场， Jacobson 条件桥没有完成。

### 10.3 下一步顺序

1. 把 `R11` 的迁移部件接到 `R8` 的 J1–J5 检查表；
2. 先做 L5 的低维相关算符污染判定；
3. 再做 L1 的算子级 $K\_B-2\pi B\_B$ 收敛或反例；
4. 同时用 `R11` 的历史冲突表阻止旧“已完成”文本回流；
5. 保持 `R11` 为外部／历史审计，不计入 `R0` 的 O1–O5 进度。

---

## §11 核验边界

`R11_check.py` 只核验：

1. 三条旧路线被分开登记；
2. 已撤回与未建项没有被写成当前完成；
3. 可迁移部件与“不可声称”边界成对出现；
4. K32 代数引理在四维类时锥上的独立代数证书；
5. 面积系数分解的算术关系；
6. 本文不进入 O1–O5、不关闭 J1／J5。

它**不**核验旧目录中每个数值脚本，也**不**把旧程序的退出码当作连续极限证明。
