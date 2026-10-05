# INDEX · 零和宇宙文档总清单

**生成方式**：由文件系统实际内容生成（标题取自各文件首行），非手写。  
**核验**：[`INDEX_check.py`](INDEX_check.py) —— 无孤儿文档 + 计数一致 + 范围声明。

---

## 0 当前状态源（先读）

**唯一当前进度与解决状态**：[`STATUS.md`](STATUS.md)。
本清单只负责文档归属、编号、范围与核验入口；任何“已解决／未解决／输入／no-go／缺口 0”都必须以 `STATUS.md` 为当前状态。
`G*`、`Z*`、`D2xx` 内的状态文字是写作时点的推导记录，除非被 `STATUS.md` 收录，否则不得当作当前结论。

| 当前项 | 状态入口 |
|:--|:--|
| 主 `Z/G` 路线 | 条件恢复；R1-A／R1-B＋R6／R7 条件定理；R8 给出 Jacobson 条件桥但 L1／L5 开放；R9／R10 把次选定为 Cao-Carroll 且只到弱场；R11 统一旧理论历史、R12 收紧 L1／L5／RC 边界、R13 证伪 L1 原样强预解并修正 G79 截断读数、R14 把 Z-CRIT-DER 重组装为三块已证＋Z-CAR 候选＋单费米点残留识别、R15 把 Z-CAR 升级为双覆盖并给出 Z-STRESS 常数 ≈π²/3、R16 证明开放具名簇停在 8 且选择型稳定为 5、R17 判定 Z14–Z16 不是 L1 临界路径并把主线切换为 L5 的 go/no-go、R18 把 L5 收成低维通道维数与归一化约束数的双因子判据并证明 Z-STRESS 约束数不足时必败、R19 重排 L1 上游并证明旋转双覆盖不能提供洛伦兹 boost、R20 判定 A1 形状成立，R21 将 R20 的常数额校正为 2π/v_F（本模型 v_F=2，故连续系数为 π），R22 分离主符号层与 R20 估计器并排除 0.9255 的有限尺寸解释，R23 给出后代优势选维的乘积 no-go 与条件模型，并证明若 L 自由则后代优势不选四维（单周期无有限峰；每步增长率选偶数 L=8）；演化层条件存活率另由精确组合计数 S_D^evo=(binom(L,L/2)/2^L)^D 证明，在 GR 的 D>=4 筛选下四维在 GR 兼容演化层扇区中条件存活率唯一最大，且不声明整个宇宙或其它层，但绝对后代数／长期占比仍需 EVO-NORM，DIM-SECTOR 仍未导出；R24 纠正路线顺序：GR-LB 只是 TEMP-GR，最终要撤掉；乘积存活率在全局候选上唯一峰在 D=1，SURV4-GLOBAL 未证，正路线是从 Zero 补 DIM-COST-Q（3/4<q<4/5）或 DIM-INTERACT；R25 把 DIM-INTERACT 收窄为 PAIR-CARRIER：原生旋转类账本给 q_L<3/4，单方向路线不选四维；成对模型 C(D,2)q^D 的四维窗口是 1/2<q<3/5，L=4 的 q=5/9 条件证成全局唯一 D=4，PAIR-CARRIER-DER 仍开放；Z14 把闭合循环序与双覆盖群论部分登记为无偏好基础扩展，Z15 证明 Z-CAR 不能由双覆盖唯一选择，并在具名 Z-READ 下由 Jordan-Wigner 条件构造 CAR 与费米宇称，Z16 证明平衡正则模块原则 Z-UNIF 不能由 Z0/Z-E* 自动推出，但在该具名输入下可由正则表示与 Jordan-Wigner 条件构造 CAR，自旋结构与循环切口仍开放，均不计入 O1–O5；E1 余 I5b／GDL，E2／E3／E4 仍为条件或具名输入 |
| 借入 `D259` 路线 | 条件图册；五通道、时间线、体积与汇流仍为路线本地未解输入 |
| 路线关系 | **未证明等价；进度不可相加** |
| 历史状态 | 见 `STATUS.md` §8；旧“唯一账本”或“未建立”清单只作历史引用 |

---

## 0.5 五项推进、外部基准与可发表主定理接口（R0–R48）

### Zero 演进与专题审计

- [`LAYER_LEDGER.md`](LAYER_LEDGER.md) — 分层账本：断言逐条带层指标（L0／L1／L1′／L2／读出）
- [`SYNTHESIS_zero_to_standard_model.md`](SYNTHESIS_zero_to_standard_model.md) — 到标准物理模型的路线图 ＋ 外部文献清单
- [`LIT_SURVEY.md`](LIT_SURVEY.md) — 外部文献调查
- [`E1_NN_verdict.md`](E1_NN_verdict.md) — Nielsen–Ninomiya 检验：平移不变性不成立；真实对称为反射
- [`README.md`](README.md) — 仓库入口

**L2 演化层专题（毁灭–重播种周期）**

| 文件 | 内容 |
|:--|:--|
| [`L2_C_recursion_verdict.md`](L2_C_recursion_verdict.md) | L2 · $C(T)$ 递推：**已对齐闭合**（含闭式） |
| [`L2_T_absolute_verdict.md`](L2_T_absolute_verdict.md) | L2 · $T$ 的绝对值：**容量封顶给出 $T=5$** |
| [`L2_anchor_verdict.md`](L2_anchor_verdict.md) | L2 · 锚定：毁灭周期（年）的**结构闭合 ＋ 一个自由锚** |
| [`L2_cat_vs_motzkin_verdict.md`](L2_cat_vs_motzkin_verdict.md) | L2 · Catalan 与 Motzkin：术语更正 ＋ 参考文献定位 |
| [`L2_catalan_verdict.md`](L2_catalan_verdict.md) | L2 · $\text{co}[s]$ 的闭式：**Catalan 数**（已定），部分和公式未竟 |
| [`L2_decoupling_verdict.md`](L2_decoupling_verdict.md) | L2 · 「相位解耦 ⇒ 无预兆、瞬间毁灭」的核查 |
| [`L2_imprint_verdict.md`](L2_imprint_verdict.md) | L2 · 路线 A：毁灭事件的可观测印记 —— **印记存在，且 T 可反解** |
| [`L2_layer_retention_corrected.md`](L2_layer_retention_corrected.md) | L2 · 层保留猜想（修正版）：**与 `D222` 完全一致** |
| [`L2_layer_retention_verdict.md`](L2_layer_retention_verdict.md) | L2 · 层保留猜想核查：**L0 完整保留 ✓，L1 保留一层 vs 两层是可区分的** |
| [`L2_period_abs_verdict.md`](L2_period_abs_verdict.md) | L2 · 毁灭周期 $T$ 的**绝对取值**试算：容量封顶给出 $T\in\{5,6\}$ |
| [`L2_period_verdict.md`](L2_period_verdict.md) | L2 · 毁灭周期 $T$：可辨识的不动点（**不是自由参数，但也不是一个数**） |
| [`L2_rho_closed_verdict.md`](L2_rho_closed_verdict.md) | L2 · $\rho(T)$ 的闭式与「反解唯一整数 $T$」的判定：**闭式存在，但反解不唯一** |
| [`L2_transfer_verdict.md`](L2_transfer_verdict.md) | L2 · 补全状态后的转移矩阵：$\lambda(T)$ 的精确代数方程 |
| [`L2_years_verdict.md`](L2_years_verdict.md) | L2 · 毁灭周期「多少年」：**可算部分已闭合，年数需要一个锚** |

**主线接口**：[`R0_publication_theorem.md`](R0_publication_theorem.md) 把“条件恢复四维 GR”写成可逐项证明、可被审稿人否定的主定理，固定五条证明义务 O1–O5，并列出验收门槛 M0–M4 与反过度主张句。以下为该分支的文档与核验脚本（状态以 `STATUS.md` 为准）。

`R8` 登记用 Zero 基础补 Jacobson 2016 前提的条件桥；`R9` 对外部强路线排序，`R10` 登记 Cao-Carroll 2018 的弱场条件桥；`R11` 统一 `cosmos-construct`／`modular-equilibrium` 与早期 `G*`／`Z*` 的历史 GR 推导；`R12` 对 L5／L1／Cao-Carroll RC 三个缺口做定向判定，`R13` 集中攻坚该 L1 正路由并证伪其原样强预解形式、修正 G79 截断读数，`R14` 把 `Z-CRIT-DER` 重组装为三块已证＋Z-CAR 候选＋单费米点残留识别，`R15` 把 Z-CAR 更正为旋转群双覆盖并给出 `Z-STRESS` 归一化常数 ≈π²/3，`R16` 审计归约树并判定开放具名簇停在 8、选择型输入稳定为 5，`R17` 把 Z14–Z16 移出 L1 临界路径并切回 L5 前置门，`R18` 把 L5 收成低维通道维数与归一化约束数的双因子判据，并证明 `Z-STRESS` 的约束数不足时 L5 必败，`R19` 重排 L1 上游并证明旋转双覆盖不能提供洛伦兹 boost，`R20` 在自由费米子支线上判定 A1 形状成立，`R21` 校正其常数额为 2π/v_F（本模型 v_F=2，故为 π），`R22` 分离主符号层与 R20 估计器，排除 0.9255 的有限尺寸解释，`R23` 把后代优势选维写成乘积 no-go 与条件模型，补上自由寿命归一化 no-go，并由演化层精确组合计数与 GR 的 D≥4 筛选证明 GR 兼容演化层扇区中四维条件存活率唯一最大；该结论只关于演化层，绝对后代数／长期占比另需 EVO-NORM，DIM-SECTOR 仍未导出，O3 未关闭。`R24` 纠正路线优先级：GR-LB 只是 TEMP-GR，最终要撤掉；SURV4-GLOBAL 未证，下一步是从 Zero 补 DIM-COST-Q 或 DIM-INTERACT。`R25` 把后一条收窄为 `PAIR-CARRIER`：单方向原生账本 no-go，成对窗口为 `1/2<q<3/5`，`L=4` 的 `q=5/9` 条件证成全局唯一 `D=4`，但 `PAIR-CARRIER-DER` 仍开放。它们只作外部基准、条件候选与历史／方向审计，不替代 O1–O5，也不把工程进度与外部论文或旧理论进度相加。

`R26` 对 R25 的成对载体作对抗审计并把 `PAIR-CARRIER-DER` 拆成六项：`DIR-DICT-BETA`、`PAIR-GRAPH-KD`、`PAIR-ID-EDGE`、`PAIR-COUNT-1`、`PAIR-COST-FACTORIZATION`、`PAIR-NO-EXTRA-MULT`。已证 no-go 为：Z1 定理 1 连通性只给 `D-1≤|E(Γ_D)|≤C(D,2)`，合法环图 `C_D` 不给 `C(D,2)`；D194/D259 的 `D=m-1` 字典把自然标签改成 `C(D+1,2)`，`q=5/9` 的唯一峰移到 `D=3`；成对重数本身不推出 `q^D` 代价，只付两维支撑时无有限峰、每条边独立付代价时峰在 `D=2`。它们只作 O3 条件候选审计，不关闭 O3。

`R27` 把 R26 的六项接口重排为相位 1-上链商路线：采用 `D=m-1` 时，原始标签仍有 `C(D+1,2)`，但边相位模顶点相位后的规范商维数为 `C(m,2)-(m-1)=C(D,2)`。若 `PHASE-1-COCHAIN`、`PAIR-ID-QUOTIENT`、`FULL-SUPPORT-LEDGER`、`WIPE-RESET-LEDGER` 与 `L=4` 同时成立，则每代峰与共同毁灭代际下的代际增长率峰都可条件落到 `D=4`；但这些桥均未从 Zero 导出，D211 的共同 `T`／Seed 仍是恢复层输入，D222 局部异步与绝对层占比另需共同代际账本及 `EVO-NORM`。旧 `T_D=τD` 反例不适用于 D211 演化层。R27 不关闭 O3。

`R28` 对 R27 的身份簇作归约：`dim(C^1/dC^0)=C(D,2)` 只给目标空间维数，不能保证实际记录已经给出这些身份。若所有边相位都是顶点势的梯度 `a_r=dθ_r`，则它们全在 `dC^0` 中，商身份子空间为零；因此 `PHASE-1-COCHAIN+PAIR-ID-QUOTIENT` 单独不足，必须补成 `PHASE-IDENTITY-DER = EDGE-CONNECTION + HOLONOMY-FULL-SPAN + INHERITANCE-IDENTITY`。`HOLONOMY-FULL-SPAN` 管身份秩，`FULL-SUPPORT-LEDGER` 管 `q^D` 代价，两者正交。R28 完成 no-go 与归约，不关闭 O3。

`R29` 对 R27 的代价簇作归约：全支撑不推出乘积记录，联合记录可以只有二维张量支撑而重叠保持 `q`，也可以把所有方向合并成一笔共同记录。代价侧最小输入改写成 `LEDGER-FACTORIZATION = DIR-SUPPORT-D + RECORD-FAMILY-D + PRODUCT-LEDGER + SAME-Q`；四项均未从 Zero 原生导出。身份簇与代价簇独立，条件四维峰必须同时采用 `PHASE-IDENTITY-DER + LEDGER-FACTORIZATION + WIPE-RESET-LEDGER + L=4`。R29 不关闭 O3。

`R30` 对「用原生循环相位充当无质量极化小群」这条选维路线作对抗审计并**封死**它：字面读法「全小群 `ISO(D−2)` 交换」在 `D≥4` 上给空解（`ISO(2)` 已非交换，数值核验 `max|[R,T]|=1.000000`）；可修复读法「旋转部分 `SO(D−2)` 交换」与 `G89` §4 第 1／2 条逐点同真（`dim SO(D−2)=(D−2)(D−3)/2=dim P_D^+`，且「`SO` 交换」⟺「`dim P_D^-≤1`」⟺`D≤4`），故只是同义改写；「循环加强」两读法全灭（有限 `Z_p⊂SO(n)` 对每个 `n≥2` 成立因而不选维；`U(1)` 可除故抽象循环读法连 `D=4` 也排除）。R30 同时撤回自己的前一稿主张：依赖没有降低（进口**不更少**，与 `(Z₂)′` 共用无质量／极化／Lorentz 三件外部结构，以 Wigner 分类替换同一轨道与 `S_m` 提升；且过筛前提 `P_grav` 已含「无质量自旋 2」），`Z0` §4.3 明写「Zero 层没有相位」，Zero 的原生非交换是**有限**的（`D_L`、`M_2(C)`）。唯一出口写成可否证的函子 `LG-FUNCTOR`（须自行给出满射到 `SO(D−2)`，不许默认）；另新登记 `Z4-EVIDENCE-SCOPE`（`Z0` §2.6 终端款依据是位形空间图 `G_T`，对 `Γ` 的效力未证）。R30 不关闭 O3。

`R31` 把 Z14 的循环序（相位）接到 R25 的成对账本：旋转类的**轨道大小分布**给出 `q_L=\sum_c\omega_c^2`，`\omega_c=o_c/N_L`，`N_L=C(L,L/2)`。**定理 R31.1**：在 `LEDGER-ROT` ＋ 成对账本 `F_D=B·C(D,2)·q_L^D` 下，全部偶寿命中**恰有 `L=4`** 使有限峰落在引力子允许域 `D\ge4`，且此时峰 `={4}`（`L=2` 给 `q_2=1`、账本无有限峰；`L\ge6` 由 `q_L\le L/N_L\le r_6=3/10<1/3` 知唯一峰为 `D=2`；`L=4` 给 `q_4=5/9\in(1/2,3/5)`）。因此 `L=4` 不再依赖 `G61` 的**最小性**（R3 §1.2 已判最小性不是导出），而由具名生存要求 `PEAK-IN-GRAVITON-DOMAIN` 给出；`q=5/9` 由轨道分布 `{4,2}` 固定。R30 与 R31 不冲突：相位**做不了小群**，但**能做账本**。未解决的张力：`R23.6` 在另一比较量（每步增长率）下给 `L=8`。R31 不关闭 O3、`SURV4-GLOBAL`、`DIM-SECTOR` 或 `EVO-NORM`。\n
`R32` 把 R31 的生存要求 `PEAK-IN-GRAVITON-DOMAIN` 施加到 `G72` §4 列出的**三条账本路线**上，证明只有**闭类旋转轨道**（路线 A）能给出落在引力子域的主导维数，且只在 `L=4`：路线 B/C（时间残类，`q=1/L`）在 `L\ge4` 上 `q\le1/4<1/3` 故峰恒为 `D=2`；路线 D12（`q=e^{-1}`）峰为 `D=3`（另有账本纯度必为有理数、`e^{-1}` 超越的 no-go）；路线 A 由 R31.1 只在 `L=4` 给峰 `{4}`。故**定理 R32.1**：唯一存活组合是 `(A, L=4)`，`argmax={4}`——同一个生存要求**同时**选出账本寄存器与寿命，`LEDGER-ROT` 从【具名输入】升为【条件选择】。**定理 R32.2** 同时消解 `R23.6` 的 `L=8` 对照：其 `q_L=e^{-1/L}` 不是账本纯度（超越 vs 有理 `M_2/M^2`），其 `L=8` 依赖的「每步增长率」正是 R27 §6 `WIPE-RESET-LEDGER` 已撤回的比较量，而在册每代比较下 R23.6 自陈无有限极大点。边界：`L\ge4` 由 `G61` (A)∧(B) 给出（不用最小性）；`L=2` 时 B/C 会给并列峰 `{3,4}`，故唯一性依赖排除 `L=2`；路线表穷尽性未证。**命题 R32.3（有条件弱化）**：在**身份计数已固定**为 `C(D,2)` 时，判据可弱化为 `PEAK-NOT-GRAVITY-FREE`（只要求峰 `\ne2`），结论不变（`{D\ge3}` 下只多留 `D12`，由 `G72` §4 有理纯度 no-go 独立排除）。**命题 R32.4（联合唯一性）**：在身份计数 `{D, C(D,2), C(D+1,2)}` × 路线 `{A, B/C, D12}` 的叉积中，`S={D\ge4}` 下**恰有一个**组合存活——`(C(D,2), A, L=4)`，峰 `{4}`（单方向峰 `D=2`；字典 `C(D+1,2)` 峰 `D=3`）。故 `PAIR-CARRIER` 也由【具名结构输入】升为【条件选择】，且弱判据不能用于联合选择。R32 不关闭 O3／`SURV4-GLOBAL`。\n
`R33` 把「**振幅的相位是否就是几何作用量的相位**」立为独立目标 `ACTION-PHASE-MATCH`（**立项，未证**），三种等价形式：(T1) 模流＝几何流、(T2) 振幅相位 `=e^{iS_geo}`、(T3) `K_ω=c·G_geo` 且 `c=2π`（区域特例即 L1 的 `K_B→2πB_B`）。材料侧：`K=−log ω` 由整数计数唯一确定、无自由参数（`G72`），复振幅／干涉／Born／模流均已导出（`G62`／`G68`）。缺口侧四个可否证子目标：S1 `2π` 归一化（现为识别）、S2 非恒定剖面＋局域性（`R12` 的 blocker）、S3 相位可加性、S4 经典极限复现 `V=κ₁^N`；三种失败形态 F1 非几何／F2 非局域／F3 归一化自由。**第一击**：boost 必须从**可逆／不可逆分裂**（`G1` 引理 5 的 `R_τ⊕H_Q`）来，**不能**从旋转双覆盖来（`R19` 已排除；代数内容 `[J,J]⊂so(3)` 生不出 boost，需混合生成元）。R33 不改动任何既有判定，也不关闭 L1。\n
`R34` 报告 R33 第一击（P2）的探针结果：**混合生成元在原生 `M_2(C)` 里存在**（`so(1,3)` 三组关系 `[J,J]=iJ`、`[J,K]=iK`、`[K,K]=-iJ` 全部实现），**但有限维不可能承载 boost**——`K_i=i\sigma_i/2` 必反 Hermitian（不存在实系数 Hermitian 解），故 `e^{i\theta K}` 不酉且范数无界（`1→1.65→12.2→148.4`）；而有限 `T` 下原生模流由 `K_\omega=-\log\omega` 生成（Hermitian、谱有限离散、跨度 `=log 50`），故为**内**自同构、闭包**紧**。**定理 R34.1**：非紧半单 Lie 群无非平凡有限维酉表示 ⟹ `(T1)`／`(T3)` 在任何有限 `T` 为假。失败原因是**维数**而非 Zero（任何有限维代数都缺 boost）。目标因此被唯一化：**P2′** 细化极限须含 **type III** 因子（BW 的几何模流需 type III₁）；**P2″** 给出二值数值指纹——`spec(log Δ_T)` 有限 vs type III₁ 需 `R`。R34 让 `R19`（旋转太紧→维数不够）、`R12`（常数剖面）、`R13`（谱半径发散）三条 no-go 合流，不关闭 L1。\n
`R35` 执行 R34 的 P2″ 数值纲领，得到**判据性**结论：极限因子的类型由 $G=\langle\log(w_i/w_j)\rangle$ 的形状决定——`{0}`／`cZ`／稠密，对应平凡／III$_\lambda$（$\lambda=e^{-c}$）／III$_1$。数值：均匀轮廓给 `{0}`（R12 的机制）；**原生满分支 `2^a` 给 III$_{1/2}$**（残差 `2.9e-12`）——**不是** BW 所需的 III$_1$；而幂律 `(a+1)^2`、阶乘 `(a+1)!`、甚至 `2^a(a+1)` **全部给 III$_1$**，即**指数因子不足以强制 III$_\lambda$**，只要轮廓不是**恰好等差**。文档里的原生例（G29 的块 `2,5,20,100`，比值 `2.5,4,5` 递增）非等差，指向 III$_1$。**净效果**：`ACTION-PHASE-MATCH` 的 (T1)/(T3) 反过来给缺失输入 `π`（`E5`）一个可否证约束——**渐近块轮廓须非等差**；并给出二值否证判据：若极限**恰为** III$_{1/2}$，则其模论非 QFT 那一支（局域代数 III$_1$）⇒ 无 BW。原生轮廓本身仍未算出。R35 不关闭 L1，也不解决 `2π`。\n
`R36` 对 Zero 现有结构做 **CHSH 检验**：机制校验通过（Bell 态 `2.8284=2√2`、直积态 `2`、Werner 门槛 `0.7075≈1/√2`）；但 `G68` 的多路径设定**不是双体**（单系统多路径，干涉 ≠ Bell 违反），而 `G82` 的两个独立 `Z₂`（`B=2×2`）因**独立动机**而联合态为**直积**，故 `S=2.0000`、`T` 秩 1——**单方向关联再强也不违反**。结论：**现有可检验的每一个二分割上 CHSH≤2（Bell 局域）**。违反门槛被量化：`S>2 ⟺ u₁²+u₂²>1 ⟺ T 秩 ≥2`（两比特时等价于纠缠；Werner 参考线 `p>1/√2`）。缺口因此被准确命名：需要的不是更强干涉，而是**空间二分割＋跨它的纠缠**——正是 `D31` 早在 `G88` §0 登记为未恢复的「**复合系统与张量积**」，与 R35 的 `III₁` 问题共享同一缺失结构。R36 不否定量子性（是「未能检验」，非「已排除」）。\n
`R37` 在**单体**层面做 **KCBS/语境性检验**（不需要 `E1`、空间分离或精确类空）。关键是把「五角星取向自由」这个陷阱用 **von Neumann 迹不等式**正确处理：对全部取向取最大有闭式 `S_max(rho)=sum_k lambda_k(rho) mu_k(A)`，`mu(A)=(sqrt5, 1.3819660, 1.3819660)`，非语境界 `2`。对照通过（相邻正交误差 `2.8e-16`；完全混合 `5/3<2`；最优纯态 `sqrt5`；3000 随机态抽查无一起界）。**结论：原生态违反 KCBS** —— 以 `G72`/`G29` 推前权重 `2,5,20,100` 的 3 维归约，谱 `(0.7407,0.1852,0.0741)` 给 `S_max=2.0146>2`（去最小块给 `2.0652`），阈值 `lambda*=0.7236`、原生余量 `+0.0171`。这是量子栏**第一个正面硬证据**：Zero 的单体统计**是语境的（量子）**，不是经典概率。与 `R36` 合起来：**单体量子性成立；多体（Bell）层面只是尚无场地**。边界：检验是 state-dependent（测量集按态选）；3 维归约＝`pi`/`E5`；谱用文档例，**若真实谱 `lambda_1<0.7236` 则结论翻转**（二值可证伪入口）。\n
`R38` 检验纠缠的涌现机制 **H：共享的闭合起因（历史层记录）＋ 一个历史层未记录的自由度 ⇒ 跨位点纠缠**，用 `Z15` 的 Jordan-Wigner 构造与 `R36` 的 CHSH 判据。结果：单位点 `c_j†|0>` 给 `S=2.000`（不纠缠）；**跨位点相干叠加 `(c_0†+c_1†)|0>/sqrt2` 给 `S=2sqrt2=2.8284`（最大纠缠）**；位点被记录（`50/50` 经典混合）退回 `S=2.000`。机制成立的原因是 **JW 字符串把「位点」非局域化**，故「同一个费米子」这一共享事实与「哪个位点」这一未记录自由度可并存。判据 **(R38-1)**：`纠缠 = 共同起因 ∧ 历史层未记录的自由度`；**(R38-2)** 它有两条实现路径：A′ 给关系加相位（`EDGE-CONNECTION`）或 **B′ 多一个未记录自由度（格点，不需新相位）**。这把 `R36` 的「没有纠缠」与 `R37` 的「单体量子性成立」缝起来。决定性下游问题 **(R38-3)**：**Zero 的历史层对「位点」究竟是失明还是记录？**（二值）。边界：`Z-READ` 未导出、链是 `1+1` 维、共同起因在探针中是给定的。\n
`R39` 裁决 `R38` 的下游问题 (R38-3)：闭合记录 `(精确词 w, 闭合类 [w])` 里有没有位点信息。**结论：完全失明** —— 环图 `C_m` 上闭合性与起始位点无关，故每个记录与**全部 m 个位点**一致，`I(位点;记录)=0`（10 组 `(m,L)` 全部 `|I|<1e-15`）；而记录仍**非平凡**（区分词形：`m=6,L=6` 给 6 类）。**结构性理由**：位点标签在 Zero 里**非原生**（`E1`/`I5b` 是第一号承重项），故记录**不可能**含它 —— `E1` 缺失这件原第一号卡点，在此**第一次变成资源**。**裁决**：走 **B′**（多一个未记录自由度），不需要 `EDGE-CONNECTION`。**唯一翻盘入口**：播种点 `P_i`（D 系列未解选择器）——若播种把位点写进记录，失明被打破，须改走 A′。\n
`R40` 裁决 `R39` 的唯一翻盘入口：播种点 `P_i` 是否把位点写进记录。**结论：不写。** 实现层：`zero_sum_cycle_evolution.py` 的播种输入是 `seed_word`（词），经 `canonical_cycle` 映到**旋转类**（L257-265）——形状层面而非位置层面；脚本自述 seed 是显式输入（L566）。结构层：等变性探针显示位点旋转下记录**相同**（6/6 组），`I(位点;记录)=0`，故位点信息在任何环节都注入不进去。于是这条线从机制到落地闭合：`R37`（单体量子性）→ `R36`（现有二分割无纠缠）→ `R38`（纠缠＝共同起因＋未记录自由度，`S=2sqrt2`）→ `R39`（记录对位点失明）→ `R40`（播种不注入位点）。**R38 的路线 B′ 保持有效，不需要改走 A′**（`EDGE-CONNECTION`）。残余：播种『用哪个词』仍是具名输入，但那是形状选择。\n
`R41` 更正 `R38`/`R40` 的一处不精确注脚并立成命题：**纠缠机制（共同起因 ＋ 未记录标签 ⇒ 纠缠）不含任何空间维度量，故维数无关**。数值：1D(`L=4`)、2D(`2x2`,`2x4`)、3D(`2x2x2`) 格子上沿哈密顿路径做 Jordan-Wigner，单粒子跨位点相干叠加**一律给 `S=2sqrt2=2.828427`**。原文「仍为 1+1 维、接到 3 维空间仍需 E1」**不精确**：`1+1` 只是所选**实现**（JW 需路径排序）的属性。真正受几何限制的是：① 把「两方」识别为「两处空间」（`E1`）；② Bell 的**类空**前提（精确锥）；③ 面积律标度（含 `D-1`）——与选维（`D=4`）无关。**反向洞见**：机制不受限 ⇒ **纠缠是默认的，退相干（完整记录）才是需要解释的那个**，与真实物理一致。\n
`R42` 把 `pi` **显式取成旋转类**（`R31` 路线 A：块＝平衡词的 `Z_L` 轨道，权重 `omega_c=o_c/N_L`），于是 `R35`/`R37` 的两项二值判据**同时变成可算**。**交叉校验**：`q_L=sum omega_c^2` 与 `R31` 的表**四个值全吻合**（`1, 5/9, 7/25, 19/175`）。**裁决一（语境性，负）**：`lambda_1` 从 `0.667` 掉到 `0.001`，`L>=4` 全部 `<0.7236` ⇒ **非语境**。**裁决二（模论类型，正）**：轨道大小 `o_c=L/d (d|L)` ⇒ 比值 `log p`；**一般 `L` 秩>=2 ⇒ `III_1`**（如 `L=12` 秩 4），素数幂子列 ⇒ `III_lambda` ⇒ 类型还取决于**沿哪条细化序列**取极限。**收窄 `R37`**：其正面结论条件于文档里的**粗分块**（`2,5,20,100`，`lambda_1≈0.79`），不是 `pi` 无关。**本轮核心**：`pi` 的第一条**双侧约束**——① 主导块 `lambda_1>0.7236`（要粗）；② 尾部对数比秩>=2（要细）；三个自然候选**没有一个同时满足**。\n
`R43` 换路——不枚举自然分割，而是**扫描分割族并构造**。**7 个自然族没有一个同时满足** `R42` 的双侧约束（最接近者：闭合长度 `L=20`，`lambda_1=0.7362` 过阈值但 `S_max=1.9849`，**差 0.8%**，且 `L->oo` 渐近失败）；而**两层族**（主导类 `p` ＋ 幂律尾 `a^{-alpha}`）在 60 组网格中 **46 组同时满足**：`p≳0.8`、`alpha≳1`、任意 `K` ⇒ `S>2`（语境）且秩≥2（`III_1`）。**最漂亮的一点**：仓库的文档例 `(2,5,20,100)` **本身就是两层形状**（`p=0.7874≈0.8`），只差**无限尾**。于是 `pi` 的形状要求被两条尺子夹出来：**主导类（要粗）＋ 幂律尾（要细）**。未导出：该形状仍需从 Zero 层结构产生，且未检查与 `R32` 生存要求（另一个泛函）的相容性。\n
`R44` 执行 `R43` 的第二笔硬账，得到**解析 no-go**：**选维（生存）与单体语境性互斥**。两者都是同一个泛函 `q=sum omega^2` 的函数：生存窗口 `F_D=C(D,2)q^D` 峰在 `D=4 <=> q in (1/2,3/5)`；语境性上界 `S_max(q)=mu2+(mu1-mu2)(1+sqrt(2q-1))/2`，**`S_max(3/5)=2.000000` 恰好**（`lambda1*=0.723607`）。于是 `q<3/5` 给 `D=4` 峰但 `S<2`，`q>3/5` 给 `S>2` 但峰 `>=D=5`——**两者在 `q=3/5` 相接不重叠**。两条既有正面结论正落在互斥两侧：路线 A `L=4`（`q=5/9`）生存✅/语境✗；文档例 `(2,5,20,100)`（`q=0.6466`）语境✅/生存✗（峰 `D=5`）。no-go 只在 `R31`/`R32` 账本框架内成立，故**逃逸出口＝改账本形式 `F_D`**（新目标）。\n
`R45` 执行 `R44` 的新目标：扫描账本族 `F_D=M(D)q^{E(D)}`（8x8=64 个形式），找使「`D=4` 峰窗口」覆盖 `q>3/5` 者。**结果：39 个形式都能——`R44` 的 no-go 可逃。** 最干净的一行：把「对」的计数从**洛伦兹对** `C(D,2)=dim so(D-1,1)` 换成**欧氏对** `C(D+1,2)=dim so(D+1)`，窗口即从 `[0.502,0.600]` 移到 `[0.602,0.666]`，**整个落进语境性区**；此时文档例（`q=0.6466, S=2.0328` ✅）与两层族（`q=0.6585, S=2.0150` ✅）**同时**给出 `D=4` 峰与语境性。**代价照实说**：39 个形式都行 ⇒ 「`D=4`」不再唯一钉住账本形式 ⇒ `R32` 的「唯一存活者」收窄为「**洛伦兹对字典内**唯一」。**新目标**：从零和输运的对结构**导出** `C(D+1,2)`（或排除它）。\n
`R46` 执行 `R45` 的目标：**从零和输运的对结构导出 `C(D+1,2)`**。关键观察：`C(k,2)` 里的 `k` 是**对象个数**，而 Zero 的原语对象是**通道**（补偿移动 `T_ex: x -> x+e_j-e_i` 由通道对 `(i,j)` 指标化），故多重度 `=C(|C|,2)`；若通道是 `D`-单纯形的顶点（`|C|=D+1`），多重度即 `C(D+1,2)`（单纯形边数）。**仓库自身的证据**：`G29` 核验三的 `single_cut: M=1+r` **逐值等于** `r`-单纯形顶点数（`r=1..7` 全对），而 `all_cuts: 2^r` 是超立方、非单纯形。**价格**：单纯形是**欧氏**的、没有签名，故选它就必须让**因果结构**供出洛伦兹签名（`G59` 锥／`R33` T1）。**二值判据**：供得出 ⇒ `D=4` 与单体量子性同时到手；供不出 ⇒ `R44` 的互斥重新生效。\n
`R47` 执行 `R46` 的判据：**因果结构能否供出洛伦兹签名？裁决：能。** 论证链：(1) 光滑锥场 ⟺ 共形洛伦兹结构（锥＝二次型零集，符号差 `(1,D-1)`；数值 `D=2..6` 全对）；(2) **签名是离散不变量** ⇒ `G59` 的「有效锥」（锥外指数小而非零）无害——指数尾巴只把边界抹糊（锥外占比随阈值 `0.395 -> 0.052`），改不了锥的拓扑（pointed/convex）。于是 `R44` 的互斥被**完全绕过**：账本多重度取单纯形边 `C(D+1,2)`（`R46`）、签名由因果锥供出（`R47`）、峰位窗口 `q in (0.6,2/3)`（`R45`）、语境性 `q>3/5`（`R44`）⇒ **`D=4` 与单体量子性可以同时到手**（文档例 `q=0.6466, S=2.0328`、两层族 `q=0.6585, S=2.0150` 均满足）。残留：共形因子/尺度不导出（与 `G57` 一致，非新缺口）、精确光锥仍缺、动力学仍需 `G1`+L1。\n
`R48` 关闭 `G59` 的开放项「有效锥 -> 精确锥」，并更正 `R47` §3 残留 2 与 `G59` §3.3 的措辞：**支持锥（精确锥）从来就是精确的**——Z1 定理 1 的最近邻传输每步最多外扩一格，故 `supp rho_n subset [-n,n]`（`B` 无关；`B=2,3,3.9,4,5` 核验零违反，`max(tip-t)=0`）。缺的不是锥，而是**锥边的填充／可见性** `eps(B) = lim_t (1/t) log rho_t(t)`：基准实测给出 `eps = 0.346579, 0.235004, 0.143841, 0.066766, 0.012659`（`B=2,2.5,3,3.5,3.9`），与闭式 `(1/2)log(4/B)` 一致到 `1e-5`（独立无饱和积分器复算同值）。于是三区制：`B<4` 锥边**指数不可见**（观察到的是被选锥 `c_*`，深阈锥边速度与 `c_*` 相对差 `<0.35%`）；`B=4` 锥边**可见**（`rho_t(t)*t -> 7.99`，幂律指数 `0.9961`，尖端速度精确 `1.000000`）；`B>4` 锥边被**饱和填满**（`O(1)`，尖端速度仍精确 `1`，与初值振幅 `1 -> 1e-10` 无关）。而 `B<=4` 是闭式（`f(+inf)=log 2`）⇒ **可见锥边 <=> B=4**，与 `G56` 的 `c_*=1` 是同一条件，故 `B=4` 从「模型参数」升级为**自洽条件**（不是新公理），`G56` §4 的「`B` 无原生约束」被收窄为「`B<=4` ＋ 绝对来源仍未导出」。

诚实边界：`eps(B)` 的**常数 `1/2` 目前是数值律（`1e-5` 一致），无解析推导**；精确的是 `eps(4)=0`、`B<=4`、以及 `lambda(mu*)=mu*c*` 与 KPP 定义。残留：`B=4` 的绝对来源、尺度（`G57`）、动力学（`G1`）、高维常数因子，以及 `R45` 的账本唯一性。核验：`R48_check.py`。

### 0.5.1 术语迁移档案（非理论文档）

[`A2Z_MIGRATION_RECORD.md`](A2Z_MIGRATION_RECORD.md)：记录 2026-10-03 的 A0–A5 → Z0 条款术语迁移（用户裁决、迁移红线、事故与纪律）。**家谱映射的权威位置是 [`G0`](G0_bottom_layer_and_derivation_route.md) §0.1**；回归闸门 [`Z0_axiom_hygiene_check.py`](Z0_axiom_hygiene_check.py)。

`R49` 执行 `R48` §5 的量子侧目标（让 $\pi$ 非等差），把 `R35` 的类型指纹与 `R42`/`R43` 的双侧约束合并成**一个有限可判定的条件**：**判据**：块权重的对数比生成子群 $G$ 稠密 $\iff$ 相邻比的对数在 $\mathbb Q$ 上线性无关 $\iff$ **素数指数差向量秩 $=k-1$**（有限、可判定）。**不可能**：$k\le3$ ⇒ 秩 $\le2$ ⇒ 任何 3 块轮廓（含 3 维归约）永远是 $III_\lambda$——这也解释了 `R42` 表里「$\lambda_1$ 越大秩越小」不是巧合。**可行**：4 块上语境性（$\lambda_1>0.723607$）与稠密性（秩 3）**解耦**，显式解 $(10^4,2,3,5)$ 给 $S_{\max}=2.2349$、$(10^3,1,3,15)$ 给 $2.2188$、$(2,246,1,5)$ 给 $2.2037$。**新障碍**：旋转类的轨道权重全是 2 的幂 ⇒ 秩恒为 1 ⇒ **旋转类恒 $III_\lambda$**（除非 $L$ 含非 2 素因子，而 $L=12$ 已被 `R31.1` 排除）。未做：从 Zero 原生生成该形状的 $\pi$。核验：`R49_check.py`。\n
`R50` 立一条**推导纪律**并据此会诊全库矛盾：**每个量、定理、常数都带层指标 $\ell$**，断言写成 $P_\ell(v)$ 才完整——默认解释不是矛盾而是**层不同**，只有**同层相反**才是真矛盾。层：L0 底层（`Z0`/`Z1`–`Z5`）／L1 历史层（$\mathcal P$）／L1′ 全局闭合类层（$\mathcal Z_\ast$）／L2 演化层（活动层）／$\mathcal R$ 读出面（$\pi,\omega,K$）。**会诊 10 条**：真矛盾 1（`L=8`，已由 `R32.2` 撤回）、符号碰撞 3（$K$／$N$／局部编号）、其余 6 条全是**层坍塌**——例如「不设概率」(L0) vs「必须靠概率」(L2/L3)、「$B$ 不可导出」(L0) vs「$B=4$ 钉住」($\mathcal R$)。并把 `R48` 的 $\varepsilon(B)$ 从单侧公式**重算为双侧定律**：$B<4$ 指数衰减、$B=4$ 幂律、$B>4$ 指数增长被 **L2 容量饱和**截断（实测斜率 $\approx0$）。三条操作规则：①断言带层指标；②比较前对齐层；③**层坍塌优先于改结论**。**纪律 ① 已落地**：6 条层坍塌＋3 条符号碰撞已逐条补层指标（`G28`／`G73`／`G59`／`R48`／`Z1`／`Z17`／`R32`／`G56`／`G61`／`G32`／`G72`／`G33`／`Z0` §0.5），只加指标不改结论，落地后 170/170 通过。核验：`R50_check.py`。\n
| 分支 | 文档 | 标题 |
|:--|:--|:--|
| `R0` | [`R0_publication_theorem.md`](R0_publication_theorem.md) | 可发表主定理接口与五项验收门槛 |
| `R1` | [`R1_gamma_convergence_theorem.md`](R1_gamma_convergence_theorem.md) | Γ-收敛：从数值二阶到可审查的条件定理 |
| `R2` | [`R2_site_identification.md`](R2_site_identification.md) | 旋转类到物理独立站点：**单值识别不可行；平衡对应可构造** |
| `R3` | [`R3_dimension_selection.md`](R3_dimension_selection.md) | 四维选维：不循环原则的审查与当前 no-go |
| `R4` | [`R4_falsifiable_prediction.md`](R4_falsifiable_prediction.md) | 条件恢复之外的候选可检验预言 |
| `R5` | [`R5_route_equivalence.md`](R5_route_equivalence.md) | 主 `Z/G` 路线与 `D210–D259` 恢复路线的等价性审计 |
| `R6` | [`R6_refutation_attempt.md`](R6_refutation_attempt.md) | H5 字典误差命题的对抗性证伪审计 |
| `R6` | [`R6_h5_dictionary_error_bound.md`](R6_h5_dictionary_error_bound.md) | H5 闭合：实际闭环权到字典多项式的一致误差界 |
| `R7` | [`R7_refutation_attempt.md`](R7_refutation_attempt.md) | H3/H7 对抗性审计：离散系数、尺度限制与正则性迁移 |
| `R7` | [`R7_h3_h7_regularity.md`](R7_h3_h7_regularity.md) | H3/H7 条件闭合：类测度到连续系数场的 mollification 构造 |
| `R8` | [`R8_L1_refutation_attempt.md`](R8_L1_refutation_attempt.md) | R8 L1 对抗审计：C2/R8-GMA 能否由现有 Zero 基础补上 |
| `R8` | [`R8_L2_area_density_audit.md`](R8_L2_area_density_audit.md) | J2 独立审计（外部阶梯 J2；与 Zero 的层 L2 无关）：普适面积密度 |
| `R8` | [`R8_jacobson_entanglement_equilibrium_completion.md`](R8_jacobson_entanglement_equilibrium_completion.md) | 用 Zero 层补 Jacobson 2016：纠缠平衡路线的条件闭合尝试 |
| `R9` | [`R9_external_GR_derivations_landscape.md`](R9_external_GR_derivations_landscape.md) | 外部“恢复 Einstein 方程或 GR”路线再排序 |
| `R10` | [`R10_cao_carroll_bulk_entanglement_completion.md`](R10_cao_carroll_bulk_entanglement_completion.md) | Cao-Carroll 2018 的 Zero 条件桥审计 |
| `R11` | [`R11_legacy_GR_derivation_audit.md`](R11_legacy_GR_derivation_audit.md) | 旧理论与历史 GR 推导审计 |
| `R12` | [`R12_zero_native_gap_filling.md`](R12_zero_native_gap_filling.md) | 用 Zero 基础补外部论文缺口：L5、L1 与 Cao-Carroll 三个定向尝试 |
| `R13` | [`R13_L1_strong_resolvent_attempt.md`](R13_L1_strong_resolvent_attempt.md) | L1 正路由：从 Zero 临界链到强预解收敛的诚实状态 |
| `R13` | [`R13_refutation_attempt.md`](R13_refutation_attempt.md) | L1 正路由对抗审计：临界自由费米链能否强预解收敛到 BW |
| `R13` | [`R13_external_limit_lemmas.md`](R13_external_limit_lemmas.md) | L1 临界自由费米链到 Bisognano–Wichmann 模 Hamiltonian：外部严格结果与缺口 |
| `R14` | [`R14_L1_from_zero_assembly.md`](R14_L1_from_zero_assembly.md) | 从 Zero 组装 L1：把 `Z-CRIT-DER` 从“识别”拆成五个部件 |
| `R15` | [`R15_zcar_double_cover_and_zstress_scale.md`](R15_zcar_double_cover_and_zstress_scale.md) | L1 推进：`Z-CAR` 的双覆盖升级，与 `Z-STRESS` 的归一化常数 |
| `R16` | [`R16_direction_audit_reduction_tree.md`](R16_direction_audit_reduction_tree.md) | 方向检查：L1 归约树是否收敛 |
| `R17` | [`R17_L1_critical_path_and_L5_gate.md`](R17_L1_critical_path_and_L5_gate.md) | L1 必要性审计与 L5 前置门：停止继续下钻 `Z-CAR` |
| `R18` | [`R18_L5_cert_dimension_normalization_gate.md`](R18_L5_cert_dimension_normalization_gate.md) | L5-CERT 第一关：低维通道的目录与归一化判据 |
| `R19` | [`R19_L1_upstream_probability_phase_and_missing_boost.md`](R19_L1_upstream_probability_phase_and_missing_boost.md) | L1 上游重排：概率与相位在位，缺的是洛伦兹 boost |
| `R20` | [`R20_A1_verdict_shape_holds_constant_fails.md`](R20_A1_verdict_shape_holds_constant_fails.md) | A1 判定：形状成立，常数不成立 |
| `R21` | [`R21_vf_normalization_resolves_the_r20_factor.md`](R21_vf_normalization_resolves_the_r20_factor.md) | R20 归一化缺口的判定：因子是 $2\pi/v\_F$，不是新增的 $1.85$ |
| `R22` | [`R22_principal_symbol_vs_r20_estimator.md`](R22_principal_symbol_vs_r20_estimator.md) | A1 主符号与 R20 比值口径的分离 |
| `R23` | [`R23_dimension_descendant_selection.md`](R23_dimension_descendant_selection.md) | 从后代优势选维（仅限演化层）：乘积 no-go 与 `DIM-DESC` 条件模型 |
| `R24` | [`R24_global_four_survival_gate.md`](R24_global_four_survival_gate.md) | 全局四维生存峰门槛：先去掉 GR 假设 |
| `R25` | [`R25_native_pair_cost_and_four_dim_peak.md`](R25_native_pair_cost_and_four_dim_peak.md) | 成对维度增益与旋转类读出账本成本：四维全局峰的首个非 GR 正候选 |
| `R26` | [`R26_pair_carrier_reduction_no_go.md`](R26_pair_carrier_reduction_no_go.md) | 成对载体约化：`C(D,2)` 不能由连通性直接得到 |
| `R27` | [`R27_phase_cochain_quotient_and_dimension_dictionary.md`](R27_phase_cochain_quotient_and_dimension_dictionary.md) | 相位上链商与维数字典汇流：从 `D=m-1` 得到 `C(D,2)` |
| `R28` | [`R28_phase_identity_cluster_reduction.md`](R28_phase_identity_cluster_reduction.md) | 相位身份簇归约：纯规范 no-go 与三重最小条件 |
| `R29` | [`R29_full_support_ledger_factorization_no_go.md`](R29_full_support_ledger_factorization_no_go.md) | 全支撑账本归约：全支撑不推出乘积记录 |
| `R30` | [`R30_little_group_phase_route_audit.md`](R30_little_group_phase_route_audit.md) | 小群—相位路线的对抗审计：一条空解 no-go、一条同义改写，与唯一幸存的桥 |
| `R31` | [`R31_phase_ledger_and_lifetime_selection.md`](R31_phase_ledger_and_lifetime_selection.md) | 相位接通账本：旋转类轨道分布唯一选出 `L=4`，并让 `D=4` 成为结论 |
| `R32` | [`R32_ledger_readout_selection_and_L8_resolution.md`](R32_ledger_readout_selection_and_L8_resolution.md) | 账本读出的选择与 `L=8` 张力的消解：一个生存要求同时选出寄存器与寿命 |
| `R33` | [`R33_action_phase_match_project.md`](R33_action_phase_match_project.md) | 作用量相位立项：`ACTION-PHASE-MATCH` 的目标、材料、可否证子目标与第一击 |
| `R34` | [`R34_finite_dimensional_boost_obstruction.md`](R34_finite_dimensional_boost_obstruction.md) | P2 探针结果：混合生成元**存在**，但有限维**不可能**承载 boost |
| `R35` | [`R35_type_iii_classification.md`](R35_type_iii_classification.md) | P2″ 结果：极限因子的**类型**由粗粒化轮廓决定，且几何模流反过来**约束 π** |
| `R36` | [`R36_chsh_bell_locality.md`](R36_chsh_bell_locality.md) | CHSH 探针：现有二分割上 **Bell 局域**，且违反所需的门槛已量化 |
| `R37` | [`R37_kcbs_contextuality.md`](R37_kcbs_contextuality.md) | KCBS 探针：Zero 的单体统计**是语境的**（量子性第一个硬证据） |
| `R38` | [`R38_entanglement_from_shared_closure_origin.md`](R38_entanglement_from_shared_closure_origin.md) | 纠缠的涌现机制：**共同起因 ＋ 未记录的自由度** |
| `R39` | [`R39_history_layer_site_blindness.md`](R39_history_layer_site_blindness.md) | 历史层对"位点"是否失明？（`R38-3` 的裁决） |
| `R40` | [`R40_seeding_does_not_inject_site.md`](R40_seeding_does_not_inject_site.md) | 播种是否把位点写进记录？（`R39` §3 残余条件的裁决） |
| `R41` | [`R41_dimension_independence.md`](R41_dimension_independence.md) | 维数无关性：纠缠机制不含空间维度量 |
| `R42` | [`R42_explicit_pi_rotation_class.md`](R42_explicit_pi_rotation_class.md) | 显式 `π`（旋转类）：两项二值判据的实际裁决 |
| `R43` | [`R43_pi_two_layer_construction.md`](R43_pi_two_layer_construction.md) | `π` 的构造：**两层族**同时满足双侧约束 |
| `R44` | [`R44_survival_vs_contextuality_no_go.md`](R44_survival_vs_contextuality_no_go.md) | no-go：**选维（生存）与单体语境性互斥** |
| `R45` | [`R45_ledger_form_scan.md`](R45_ledger_form_scan.md) | 账本形式扫描：`R44` 的 no-go **可逃**，代价是 `R32` 的唯一性被削弱 |
| `R46` | [`R46_pair_ledger_objects_and_simplex.md`](R46_pair_ledger_objects_and_simplex.md) | 对账本的"对象"是什么？`C(D+1,2)` 的**导出**与它的价格 |
| `R47` | [`R47_signature_from_causal_cone.md`](R47_signature_from_causal_cone.md) | 因果结构供出洛伦兹签名：`R46` 的价格**可以付** |
| `R48` | [`R48_exact_cone_and_effective_cone.md`](R48_exact_cone_and_effective_cone.md) | 精确锥的结算：它**一直是精确的**，缺的是**锥边的可见性** |
| `R49` | [`R49_pi_nonarithmetic_criterion.md`](R49_pi_nonarithmetic_criterion.md) | 让 $\pi$ 非等差：**精确判据、3 块不可能、4 块可构造** |
| `R50` | [`R50_layer_discipline.md`](R50_layer_discipline.md) | 分层纪律：为什么同一个量在不同层"对、反、全局、局部"都对 |
| `R54` | [`R54_L3_removed_quantum_emergence.md`](R54_L3_removed_quantum_emergence.md) | 撤掉 $\mathcal R$ 之后量子力学从哪来：**$L=16$ 上的原生推导** |
| `R57` | [`R57_evolution_layer_quantum_chain.md`](R57_evolution_layer_quantum_chain.md) | 演化层逼出量子力学：**猜想串起来的完整链** |
| `R58` | [`R58_evolution_layer_quantum_emergence.md`](R58_evolution_layer_quantum_emergence.md) | Zero 量纲内量子力学的涌现 |
| `R59` | [`R59_conjecture_inventory.md`](R59_conjecture_inventory.md) | 猜想总盘点：证毕、否证与修正 |
| `R64` | [`R64_lorentzian_stable_phase.md`](R64_lorentzian_stable_phase.md) | 3+1 维洛伦兹度量为何是稳定相 |
| `R68` | [`R68_causal_metric_emergence.md`](R68_causal_metric_emergence.md) | 方向 2：闭环演化的因果度规 —— 洛伦兹号差是**涌现**的 |
| `R71` | [`R71_continuum_and_dimension.md`](R71_continuum_and_dimension.md) | 连续极限与维数：因果度规的两个缺口 |
| `R72` | [`R72_layer_attribution.md`](R72_layer_attribution.md) | 维数与 GR 的层归属修正 |
| `R73` | [`R73_pair_carrier_der.md`](R73_pair_carrier_der.md) | PAIR-CARRIER-DER：通道对数、身份秩、与和乐秩的三方辨析 |
| `R74` | [`R74_euclidean_simplex_settlement.md`](R74_euclidean_simplex_settlement.md) | 欧氏单纯形账本结算：四维峰与单体语境性同时闭环 |
| `R77` | [`R77_length_grading_necessity.md`](R77_length_grading_necessity.md) | 长度分级的必然性：根稳定子轨道 |
| `R80` | [`R80_beta_lock_verdict.md`](R80_beta_lock_verdict.md) | 最后一块骨牌：$T^{ab}$ 能否把 $\beta\varepsilon$ 钉死在 $4/3$？ |
| `R82` | [`R82_conformal_gap_memo.md`](R82_conformal_gap_memo.md) | 卡点备忘录：不是"度规构造失败"，而是"共形因子无来源" |
| `R85` | [`R85_local_isotropy_verdict.md`](R85_local_isotropy_verdict.md) | 方向 3 的裁决：局部各向同性在单纯形上**必然失败** |
| `R86` | [`R86_reseed_class_verdict.md`](R86_reseed_class_verdict.md) | 重播种是否保类号？——**判据落空，但缺口 1 的真身被抓到** |
| `R87` | [`R87_class_weight_consistency.md`](R87_class_weight_consistency.md) | 类级账本权重的一致性（(a) 项）——`R25` 的权重**是动力学的推前**，F4 是口径错配 |
| `R88` | [`R88_path_yi_falsification.md`](R88_path_yi_falsification.md) | 路径乙的证伪（(b) 项）——符号对称幸存，但代价是整条 `R31`/`R32` 链 |
| `R89` | [`R89_layer_audit.md`](R89_layer_audit.md) | 层审计：`R86`–`R88` 的对象各自在哪一层？ |
| `R90` | [`R90_layer_distribution.md`](R90_layer_distribution.md) | 量子力学是**分布在各层**的（用户猜想的判定） |
| `R91` | [`R91_L3_status_contradiction.md`](R91_L3_status_contradiction.md) | `$\mathcal R$` 的地位：语料内的一处**真矛盾**（并结算它对量子栏的影响） |
| `R92` | [`R92_layer_inventory_and_dictionary.md`](R92_layer_inventory_and_dictionary.md) | 层清单与物理对应 |
| `R93` | [`R93_L3_sublayer.md`](R93_L3_sublayer.md) | $\mathcal R$ 是不是 L2 的亚层？ |
| `R94` | [`R94_GR_layer_compatibility.md`](R94_GR_layer_compatibility.md) | 把 GR 套上分层方法：兼容性检验 |
| `R95` | [`R95_layer_table_and_discipline.md`](R95_layer_table_and_discipline.md) | 统一层表与层纪律（单一来源） |
| `R96` | [`R96_quantum_inventory.md`](R96_quantum_inventory.md) | Zero 量子栏成果与缺口盘点 |
| `R97` | [`R97_sublayer_structure.md`](R97_sublayer_structure.md) | L0／L1／L2 的亚层细分 |
| `R98` | [`R98_phase_group.md`](R98_phase_group.md) | L1 的相位缺口：原生相位群是**L 次单位根** |
| `R99` | [`R99_L3_status_and_phase_grading.md`](R99_L3_status_and_phase_grading.md) | $\mathcal R$ 还存在吗？＋ 相位分级靶子的否定结果 |
| `R100` | [`R100_symbol_resolution.md`](R100_symbol_resolution.md) | 符号决议：$L$ 只留给层，参数改用 $n$／$\tau$ |
| `R101` | [`R101_L3_disposition_and_L_cleanup.md`](R101_L3_disposition_and_L_cleanup.md) | $\mathcal R$ 的处置 ＋ $L$-编号系统的清场 |
| `R102` | [`R102_quantum_gap_recount.md`](R102_quantum_gap_recount.md) | 量子栏缺口再盘点：最大缺口是**作用量**，不是多体结构 |

核验：[`R0_check.py`](R0_check.py)、[`R10_check.py`](R10_check.py)、[`R11_check.py`](R11_check.py)、[`R12_check.py`](R12_check.py)、[`R13_check.py`](R13_check.py)、[`R14_check.py`](R14_check.py)、[`R15_check.py`](R15_check.py)、[`R16_check.py`](R16_check.py)、[`R17_check.py`](R17_check.py)、[`R18_check.py`](R18_check.py)、[`R19_check.py`](R19_check.py)、[`R1_check.py`](R1_check.py)、[`R20_check.py`](R20_check.py)、[`R21_check.py`](R21_check.py)、[`R22_check.py`](R22_check.py)、[`R23_check.py`](R23_check.py)、[`R24_check.py`](R24_check.py)、[`R25_check.py`](R25_check.py)、[`R26_check.py`](R26_check.py)、[`R27_check.py`](R27_check.py)、[`R28_check.py`](R28_check.py)、[`R29_check.py`](R29_check.py)、[`R2_check.py`](R2_check.py)、[`R30_check.py`](R30_check.py)、[`R31_check.py`](R31_check.py)、[`R32_check.py`](R32_check.py)、[`R33_check.py`](R33_check.py)、[`R34_check.py`](R34_check.py)、[`R35_check.py`](R35_check.py)、[`R36_check.py`](R36_check.py)、[`R37_check.py`](R37_check.py)、[`R38_check.py`](R38_check.py)、[`R39_check.py`](R39_check.py)、[`R3_check.py`](R3_check.py)、[`R40_check.py`](R40_check.py)、[`R41_check.py`](R41_check.py)、[`R42_check.py`](R42_check.py)、[`R43_check.py`](R43_check.py)、[`R44_check.py`](R44_check.py)、[`R45_check.py`](R45_check.py)、[`R46_check.py`](R46_check.py)、[`R47_check.py`](R47_check.py)、[`R48_check.py`](R48_check.py)、[`R49_check.py`](R49_check.py)、[`R4_check.py`](R4_check.py)、[`R50_check.py`](R50_check.py)、[`R58_check.py`](R58_check.py)、[`R5_check.py`](R5_check.py)、[`R6_check.py`](R6_check.py)、[`R6_independent_check.py`](R6_independent_check.py)、[`R74_check.py`](R74_check.py)、[`R7_check.py`](R7_check.py)、[`R7_independent_check.py`](R7_independent_check.py)、[`R8_L1_check.py`](R8_L1_check.py)、[`R8_L2_check.py`](R8_L2_check.py)、[`R8_check.py`](R8_check.py)、[`R95_check.py`](R95_check.py)、[`R9_check.py`](R9_check.py)。

---

## 1 三个来源与归属

| 前缀 | 来源 | 归属 | `lh/` 中的数量 |
|:--|:--|:--|--:|
| **`G*`** | 本会话独立推导（自底层条款（Z0 条款 ＋ Z1–Z5 定理）到 GR） | **零和宇宙** | 90 篇 + 88 个核验脚本 |
| **`D2xx`** | 零层弧（D210–D259） | **零和宇宙**（其中 8 篇含非原生桥接） | 50 篇 |
| **`zero_sum_*`** | 零和宇宙的仿真验证笔记与程序 | **零和宇宙** | 11 篇 + 10 个脚本 |
| **`Z*`** | **基础层**：单一公理「零不断乱动」＋ 取代 A1–A5（历史命名）＋ Zero 结构扩展 | **零和宇宙** | 18 篇 + 18 个核验脚本 |
| `D1–D209` | 混合语料 | **不属零和宇宙** | **0**（未拷入） |

$$
\boxed{\text{zero-sum universe}\ =\ \texttt{zero\_sum\_*}\ \to\ \texttt{Z*}\ \to\ \texttt{G*}\ +\ \texttt{D2xx}}
$$

---

## 2 编号规则（**不改字头**）

**结论：保留 `zero_sum` / `Z` / `G` / `D` 四个字头，不统一。** 字母序即**读序**：`zero_sum`（仿真）→ `Z`（基础层）→ `G`（推导）→ `D`（对照）。理由：

1. **D 编号是稳定 ID**——被正文、母项目、D 文档间交叉引用、以及核验脚本共同引用；改字头全断；
2. **三个来源确实不同**（作者／纪律／核验基础设施），合成一个字头会丢失来源信息；
3. **改动量巨大、收益为零**：G 系列内部有数百条链接，D 之间有数百条交叉引用；
4. **真正的歧义只有 D**（D1–D209 vs D210–D259），而在 `lh/` 里**只有 D210–D259**。

**消歧方式＝本清单（文档）而不是改名。** 如需物理分目录，见 §8。

---

## 3 `zero_sum_*` 仿真验证（11 篇笔记 / 10 个程序 ＋ 4 个核验脚本）——**Zero 原始层**

这是体系的**最底层**：先有仿真，后有公理。分诊见 [`G90`](G90_zero_series_reference_triage.md)。

**本层由两部分构成：文章（笔记）与程序（可执行沙盒）**。程序的定理已写成文章：[`旋转类代数`](zero_sum_rotation_class_algebra.md)、[`繁殖转移定理`](zero_sum_reproduction_transition_theorems.md)、[`闭合图定理`](zero_sum_closure_graph_theorems.md)、[`持续性定理`](zero_sum_persistence_theorems.md)。

### （甲）原始笔记（4 篇）——仿真记录

| 笔记 | 标题 |
|:--|:--|
| [`zero_sum_closure_exit.md`](zero_sum_closure_exit.md) | 闭合路径退出演化层：程序验证 |
| [`zero_sum_global_R_local_P.md`](zero_sum_global_R_local_P.md) | 全局 $R$、局部 $P$ 验证记录 |
| [`zero_sum_open_reservoir.md`](zero_sum_open_reservoir.md) | 持续开放储层与少量闭合：程序验证 |
| [`zero_sum_periodic_destruction.md`](zero_sum_periodic_destruction.md) | 演化层定期毁灭与再生：程序验证 |

### （乙）由程序写成的定理文章（7 篇）

程序里装着定义与定理；这些文章把它们写成可读形式（公式＋数值核验＋与体系的接口）。

| 文章 | 内容 |
|:--|:--|
| [`zero_sum_closure_graph_theorems.md`](zero_sum_closure_graph_theorems.md) | 闭合图定理：度无界与"无稳定谱维数" |
| [`zero_sum_closure_time_selection.md`](zero_sum_closure_time_selection.md) | 闭合时间的选择 |
| [`zero_sum_persistence_theorems.md`](zero_sum_persistence_theorems.md) | 持续性定理：熄灭判据与重播种映射 |
| [`zero_sum_reproduction_audit.md`](zero_sum_reproduction_audit.md) | 繁殖公式审计 |
| [`zero_sum_reproduction_transition_theorems.md`](zero_sum_reproduction_transition_theorems.md) | 繁殖转移定理：五条规则、谱半径与幂零灭绝 |
| [`zero_sum_rotation_class_algebra.md`](zero_sum_rotation_class_algebra.md) | 旋转类代数：项链计数、内部归零与三套闭合重数 |
| [`zero_sum_tri_layer_universe.md`](zero_sum_tri_layer_universe.md) | 零和三层宇宙 |

### （丙）本层的核验脚本（4 个）

核验与程序分离：程序给数据，脚本复核文章里的公式与数字。

| 脚本 | 核验的文章 |
|:--|:--|
| [`zero_sum_closure_graph_theorems_check.py`](zero_sum_closure_graph_theorems_check.py) | `zero_sum_closure_graph_theorems.md` |
| [`zero_sum_persistence_theorems_check.py`](zero_sum_persistence_theorems_check.py) | `zero_sum_persistence_theorems.md` |
| [`zero_sum_reproduction_transition_theorems_check.py`](zero_sum_reproduction_transition_theorems_check.py) | `zero_sum_reproduction_transition_theorems.md` |
| [`zero_sum_rotation_class_algebra_check.py`](zero_sum_rotation_class_algebra_check.py) | `zero_sum_rotation_class_algebra.md` |
另有 **4 个动力学沙盒**（`zero_sum_living_universe.py`、`zero_sum_geometry_probe.py`、`zero_sum_cycle_evolution.py`、`zero_generative_selection.py`）——已在 [`G28`](G28_dynamics_audit.md) 审计。

---

## 4 `Z*` 基础层：单一公理与取代 A1–A5（历史命名）（18 篇 ＋ 18 个核验脚本）

**读序**：本层在 Zero 原始层**之上**、G 推导系列**之前**。它把 Zero 的架构确立为 Z0 条款（公理只有 Z0），并把 A1–A5 降为定理。

| 文档 | 标题 |
|:--|:--|
| [`Z0_zero_never_rests_single_axiom.md`](Z0_zero_never_rests_single_axiom.md) | 零不断乱动：**单一公理**与 Z1–Z5 的导出 |
| [`Z1_zero_layer_as_the_foundation.md`](Z1_zero_layer_as_the_foundation.md) | Zero 层作为基础：A1–A5的**取代**与 Zero 结构的**扩充** |
| [`Z2_zero_to_gr_direct_route.md`](Z2_zero_to_gr_direct_route.md) | 从 Zero **直连** GR：免去 A0–A5桥接的路线与两处残余输入 |
| [`Z3_i2a_dimension_drift_verdict.md`](Z3_i2a_dimension_drift_verdict.md) | I2a 的维数漂移判定：**固有**还是**截断假象** |
| [`Z4_p1_weight_scaling_and_shape.md`](Z4_p1_weight_scaling_and_shape.md) | P1：G40 权重在细化下的标度与形状 —— I2a 到底缺什么 |
| [`Z5_finite_k_locality_escape.md`](Z5_finite_k_locality_escape.md) | 有限 $k$：逃出三角困境，与三笔代价（其中一笔是新账） |
| [`Z6_stall_autopsy_and_released_ledger.md`](Z6_stall_autopsy_and_released_ledger.md) | 卡点解剖：为什么"退化出 GR"老停在同一处，以及**放开输入**后的终态账本 |
| [`Z7_embedding_input_explicit_dictionary.md`](Z7_embedding_input_explicit_dictionary.md) | E1 的显式化：嵌入输入＝**一个局部标度场**；度规的**闭式字典**（1D／2D／3D 二阶验证） |
| [`Z8_native_scale_field_candidate.md`](Z8_native_scale_field_candidate.md) | $c$ 的**原生候选**：$\pi$ 的推前重数（附一次**被证伪**的候选） |
| [`Z9_pi_filter_and_lifetime_fork.md`](Z9_pi_filter_and_lifetime_fork.md) | $\pi$ 的筛子：**年龄筛**与**生存筛**，以及"局部寿命"的**分叉** |
| [`Z10_position_field_and_amplitude_criterion.md`](Z10_position_field_and_amplitude_criterion.md) | 位置型原生场：D214 的环与**归零概率**——逃出引理 79，但**死在 $k=L$ 的放大**上 |
| [`Z11_correction_k_is_fixed_and_age_is_not_site.md`](Z11_correction_k_is_fixed_and_age_is_not_site.md) | **更正 [`Z10`](Z10_position_field_and_amplitude_criterion.md) §3**：$k=L$ 是**固定步数**（不是格子大小）——并由此排除一条识别（**年龄 ≠ 站点**） |
| [`Z12_geometric_input_closed_binary_labeling.md`](Z12_geometric_input_closed_binary_labeling.md) | E1 的收口：几何输入＝**一个原生分类标号**（$L=4$ 时**二值**）——类规模的**因子律**与均匀化 |
| [`Z13_zero_foundation_missing_principle.md`](Z13_zero_foundation_missing_principle.md) | Zero 基础缺什么：不是公理，而是选择／读出原理 |
| [`Z14_closure_cyclic_order_base_theorem.md`](Z14_closure_cyclic_order_base_theorem.md) | 闭合的循环序与双覆盖：基础扩展 Z-E* |
| [`Z15_zcar_no_go_and_jordan_wigner_readout.md`](Z15_zcar_no_go_and_jordan_wigner_readout.md) | `Z-CAR` 的物理读出：无唯一性定理与 Jordan–Wigner 条件构造 |
| [`Z16_zunif_balanced_regular_module.md`](Z16_zunif_balanced_regular_module.md) | 读出无偏性与平衡正则模块：`Z-UNIF` 条件构造 |
| [`Z17_A0_A5_retirement_vacancy_ledger.md`](Z17_A0_A5_retirement_vacancy_ledger.md) | A0–A5 退场后的空缺账本：no-go 的 Z0 化与维数缺口的定位 |

核验：[`Z0_check.py`](Z0_check.py)、[`Z1_check.py`](Z1_check.py)、[`Z2_check.py`](Z2_check.py)、[`Z3_check.py`](Z3_check.py)、[`Z4_check.py`](Z4_check.py)、[`Z5_check.py`](Z5_check.py)、[`Z6_check.py`](Z6_check.py)、[`Z7_check.py`](Z7_check.py)、[`Z8_check.py`](Z8_check.py)、[`Z9_check.py`](Z9_check.py)、[`Z10_check.py`](Z10_check.py)、[`Z11_check.py`](Z11_check.py)、[`Z12_check.py`](Z12_check.py)、[`Z13_check.py`](Z13_check.py)、[`Z14_check.py`](Z14_check.py)、[`Z15_check.py`](Z15_check.py)、[`Z16_check.py`](Z16_check.py)、[`Z17_check.py`](Z17_check.py)。表头 `Z0` = **唯一公理**；`Z1` = **取代 A1–A5 与结构扩展**。


**四条结构扩展**（均在 G 系列之前，各带价签）：

| 条款 | 买回 | 处置 |
|:--|:--|:--|
| **Z-E1** 全局读出耦合 | Z1 定理 1 的**连通性**；G23 (a) 的全局读出 $R_Z$／层间核 $K_{ij}$ | **识别 U**（不是公理）；与 `tri_layer` 笔记的局部重建相反，属改选 |
| **Z-E2** 局部寿命与年龄 | G23 **(b) 年龄簇（22 篇）**主体 | **导出**：年龄＝词长；局部 $\tau_i$ 为参数异质性 |
| **Z-E3** 周期与相位 | G23 **(c) 周期／相位簇（9 篇）**主体 | 周期＝**参数**；相位＝**定义**（需先有 $T$） |
| **Z-E4** 符号与配对 | G23 **(d) 符号／配对簇（9 篇）**主体 | **导出**：符号＝双向步；配对＝补偿步对 $(+,-)$ |

（G23 的 (e) 模簇 4 篇**明确不并入**——属 D1–D209 上游侧，登记【开放】。）

---

## 5 `G*` 推导系列（90 篇）

### 底层与路线

- [`G0_bottom_layer_and_derivation_route.md`](G0_bottom_layer_and_derivation_route.md) — 底层**条款表**与到 GR 的推导路线（**公理只有 [`Z0`](Z0_zero_never_rests_single_axiom.md) 一条**）

### 几何骨架（引理 1–8）

- [`G1_derivations_from_the_bottom_layer.md`](G1_derivations_from_the_bottom_layer.md) — 从零和底层导出几何骨架
- [`G2_local_continuum_limit.md`](G2_local_continuum_limit.md) — 局域连续极限：输入 I2 的三分解与两个障碍
- [`G3_admittance_fixed_point.md`](G3_admittance_fixed_point.md) — 导纳自洽不动点：显式解、尺度含义与一个负结果
- [`G4_assembly_route_obstruction.md`](G4_assembly_route_obstruction.md) — 代数装配路线的排除：等边化、刚性与度量不可导出
- [`G5_stress_lift_and_conservation.md`](G5_stress_lift_and_conservation.md) — 应力提升：守恒流的来源、正则尘埃提升与 (C) 的不足
- [`G6_geodesy_of_the_coarse_grained_flow.md`](G6_geodesy_of_the_coarse_grained_flow.md) — 粗粒化流的测地性：I3c 判定与物质类别的选择
- [`G7_one_operator_and_the_dissipation_obstruction.md`](G7_one_operator_and_the_dissipation_obstruction.md) — 同一算子的两个后果：结构性耦合，与它制造的耗散障碍
- [`G8_dimension_selection.md`](G8_dimension_selection.md) — 维数筛选：Z0 条款不选维数，但结构把 $D\le3$ 排除

### D 系列对照与收官

- [`G9_d_series_reference_triage.md`](G9_d_series_reference_triage.md) — D 系列参考对照：哪些可承接、哪些要打问号
- [`G10_final_derivation_and_input_ledger.md`](G10_final_derivation_and_input_ledger.md) — 收官：从零和底层到四维 GR 的完整推导与输入总账

### 维数（引理 35–39）

- [`G11_dimension_as_consistency.md`](G11_dimension_as_consistency.md) — 维数：定位更正、no-go 与「极化两标签无偏好」条件

### 规范扇区评估

- [`G12_gauge_sector_minimal_extension.md`](G12_gauge_sector_minimal_extension.md) — 评估：规范扇区的最小扩展——需要什么、代价是什么、会不会破坏已导出的结论

### 叶层与洛伦兹不变性

- [`G13_foliation_and_lorentz_invariance_gap.md`](G13_foliation_and_lorentz_invariance_gap.md) — 叶层与洛伦兹不变性缺口（账本第 6 条 I6）
- [`G14_causal_closure_and_lorentz_emergence.md`](G14_causal_closure_and_lorentz_emergence.md) — I6 的物质层：因果闭合把抛物与双曲统一成一个参数族
- [`G15_bare_ax3_has_no_characteristic_speed.md`](G15_bare_ax3_has_no_characteristic_speed.md) — 裸 Z0③（无偏好）下不存在有限特征速度：I6 物质层的 no-go

### 不增扩充条款清算

- [`G16_repair_audit_without_new_axioms.md`](G16_repair_audit_without_new_axioms.md) — 「不增扩充条款」的清算：哪些缺口能在 Z0 条款内修复

### 优先帧的正面物理

- [`G17_positive_physics_of_the_preferred_frame.md`](G17_positive_physics_of_the_preferred_frame.md) — 优先帧的正面物理：不增扩充条款下能算出什么

### 连续极限的可攻性

- [`G18_attackability_of_the_continuum_limit.md`](G18_attackability_of_the_continuum_limit.md) — I2a 的可攻性评估：化归为 Γ-收敛命题并数值刻画

### 条款精简


### 语料、动力学与概率

- [`G20_axiom_audit_extended_to_zero_and_D.md`](G20_axiom_audit_extended_to_zero_and_D.md) — 把 Z0 条款审计扩展到 zero/D 语料（原题：公理审计）
- [`G21_do_the_layers_help_derive_GR.md`](G21_do_the_layers_help_derive_GR.md) — 层级结构对推 GR 有帮助吗？它能被严格导出吗？
- [`G22_correction_sublayers_and_local_layers.md`](G22_correction_sublayers_and_local_layers.md) — 修正 G21：补上亚层（stratification）与局部／全局层结构
- [`G23_zero_layer_structure_inventory.md`](G23_zero_layer_structure_inventory.md) — 零层基础结构清点：Z0 条款漏了什么
- [`G24_age_structure_and_its_conflicts.md`](G24_age_structure_and_its_conflicts.md) — 年龄结构：零层弧最大簇的完整结构，以及它与 I2a / I8 的冲突
- [`G25_age_to_geometry_channel_is_obstructed.md`](G25_age_to_geometry_channel_is_obstructed.md) — 年龄 → 几何通道被堵住：D234/D235/D236 与四条独立的 no-go
- [`G26_scope_and_non_native_structures.md`](G26_scope_and_non_native_structures.md) — 范围与「非原生结构」清点
- [`G27_purification_attempt.md`](G27_purification_attempt.md) — 净化尝试：7 类非原生结构里 6 类可从 Z0 条款重新导出
- [`G28_dynamics_audit.md`](G28_dynamics_audit.md) — 动力学审计：零和宇宙的动力学建起来没有？
- [`G29_probability_as_derived_not_postulated.md`](G29_probability_as_derived_not_postulated.md) — 概率作为导出量：计数测度的粗粒化推前（不增扩充条款）
- [`G30_memory_kernel_test.md`](G30_memory_kernel_test.md) — 记忆核检验：推前动力学是**非马尔可夫**的
- [`G31_characteristic_speed_and_saturation.md`](G31_characteristic_speed_and_saturation.md) — 特征速度：因果锥 vs 被选速度，及**饱和**的作用
- [`G32_native_origin_of_saturation.md`](G32_native_origin_of_saturation.md) — 饱和的原生来源：局部代数的有限维给出有限容量
- [`G33_macro_master_equation_and_mz_kernel.md`](G33_macro_master_equation_and_mz_kernel.md) — 宏观主方程的导出：精确 Mori–Zwanzig 核与可集块判据
- [`G34_exit_rule_absorbed_into_pi.md`](G34_exit_rule_absorbed_into_pi.md) — 退出规则并入 pi：R-Z-CLOSURE-EXIT-RULE 不是独立的动力学规则
- [`G35_reseeding_and_chirality.md`](G35_reseeding_and_chirality.md) — 再播种机制（I8）：D233 的障碍是**小周期**现象，$L\ge8$ 被原生结构打破
- [`G36_objective_completion_audit.md`](G36_objective_completion_audit.md) — 目标完成审计：零和宇宙动力学的推导
- [`G37_reseeding_law_and_sign_symmetry_theorem.md`](G37_reseeding_law_and_sign_symmetry_theorem.md) — 定量重播种律与符号对称定理（更正 G35）
- [`G38_d232_verdict_and_g35_withdrawal.md`](G38_d232_verdict_and_g35_withdrawal.md) — `D232`/`D233` 读法定案：$f(a)$ 就是相位比 —— **G37 确认、G35 撤回**
- [`G39_is_Z2_designed_for_the_theorems.md`](G39_is_Z2_designed_for_the_theorems.md) — Z2 的「无偏好」（即 Z0③）是为后面那些定理而设的吗？
- [`G40_metric_from_closed_walk_counting.md`](G40_metric_from_closed_walk_counting.md) — 从其余条款导出度规：**闭环计数度规**
- [`G41_lovelock_premises_under_nonuniform_weight.md`](G41_lovelock_premises_under_nonuniform_weight.md) — 闭环计数度规是否满足 Lovelock 前提？—— **$(L')$（$g$-locality）不成立；$(L)$ 自动成立**
- [`G42_third_route_local_weights.md`](G42_third_route_local_weights.md) — 第三条路：接受度规为输入，但要求权重**局域** —— **Einstein 方程到手**
- [`G43_dynamics_is_metric_free.md`](G43_dynamics_is_metric_free.md) — 动力学是 **metric-free** 的：度规问题可以搁置
- [`G44_metric_needs_a_scale_not_an_origin.md`](G44_metric_needs_a_scale_not_an_origin.md) — 度规缺的不是「起源」，而是「一个**尺度** $k$
- [`G45_units_vs_scales_is_the_ruler_human.md`](G45_units_vs_scales_is_the_ruler_human.md) — 「这把尺子跟人有关吗？」—— **单位 vs 尺度**
- [`G46_k_is_the_lifetime.md`](G46_k_is_the_lifetime.md) — $k$ 由 Z4 的**寿命**唯一确定：$k=L$
- [`G47_refinement_limit_of_the_effective_metric.md`](G47_refinement_limit_of_the_effective_metric.md) — 有效度规在格距细化下收敛吗？（**I2a 的新形态**）
- [`G48_per_layer_scales_and_metrics.md`](G48_per_layer_scales_and_metrics.md) — 每层有自己的标度与度规 —— **分层是解药**
- [`G49_four_boundaries_advanced.md`](G49_four_boundaries_advanced.md) — 推进 [`G48`](G48_per_layer_scales_and_metrics.md) 的四条诚实边界
- [`G50_coarsening_axis_vs_length_axis.md`](G50_coarsening_axis_vs_length_axis.md) — 粗粒化轴 vs 长度轴：**两轴独立**
- [`G51_three_exponents_and_the_two_faces_verdict.md`](G51_three_exponents_and_the_two_faces_verdict.md) — 三个指数与"两个脸孔"的**判定：未实现**（资格限定）
- [`G52_spectrum_preserving_rg_nontrivial_fixed_point.md`](G52_spectrum_preserving_rg_nontrivial_fixed_point.md) — 保谱（保测度）重整化：**非平凡不动点出现了**（限定 [`G51`](G51_three_exponents_and_the_two_faces_verdict.md) 的 (d)）
- [`G53_connecting_the_rg_axis_to_lattice_refinement.md`](G53_connecting_the_rg_axis_to_lattice_refinement.md) — 把 [`G52`](G52_spectrum_preserving_rg_nontrivial_fixed_point.md) 的粗粒化接到 [`G47`](G47_refinement_limit_of_the_effective_metric.md) 的细化上：**I2a 判定**

### 定量剖面、退化与归一化

- [`G54_quantitative_profile_age_measure.md`](G54_quantitative_profile_age_measure.md) — 定量剖面 $f(a)$：非恒定性只能来自「参考测度」，不能来自「配对」
- [`G55_dynamics_line_degeneration_to_GR.md`](G55_dynamics_line_degeneration_to_GR.md) — 动力学线能否退化出 GR？—— **形式能、因果不能**，以及到 GR 的路线图
- [`G56_degeneration_attempt2_six_slots.md`](G56_degeneration_attempt2_six_slots.md) — 退化尝试 #2：以 `D213` 的六槽位预算重估
- [`G57_unreachability_of_absolute_normalization.md`](G57_unreachability_of_absolute_normalization.md) — 不可达定理：绝对归一化（$C\_{\rm norm}$ 收官）

### Γ-收敛与 I2a 判定

- [`G58_I2a_resolved_as_embedding_input.md`](G58_I2a_resolved_as_embedding_input.md) — I2a 判定：**Γ-收敛，不是 UV 不动点** —— I2a 归并为嵌入输入

### 因果锥与 I7 结算

- [`G59_I7_settled_native_cone_and_its_residue.md`](G59_I7_settled_native_cone_and_its_residue.md) — I7 结算：因果锥是**原生**的，残余只在 $B=4$ 处准精确实现

### 无量纲账本与自由单位

- [`G60_dimensionless_ledger_and_one_free_unit.md`](G60_dimensionless_ledger_and_one_free_unit.md) — 无量纲账本 ＋ 一个自由单位

### 参数锁定、量子扇区与对照清单

- [`G61_locking_the_five_integers.md`](G61_locking_the_five_integers.md) — 锁定五个整数参数（$L$, SPAWN, $N$, $K$, $B$）
- [`G62_quantum_sector_from_GNS_modular_flow_gleason.md`](G62_quantum_sector_from_GNS_modular_flow_gleason.md) — 量子扇区的推导：GNS ＋ 模流 ＋ Gleason（并**限定 G11**）
- [`G63_target_list_and_audit.md`](G63_target_list_and_audit.md) — 对照清单：零和宇宙**本该**预言的量

### 自旋 1/2 的核验

- [`G64_is_spin_half_native.md`](G64_is_spin_half_native.md) — 自旋 1/2 是否原生？—— **载体已到，双值性差一个识别**
- [`G65_spin_half_needs_3d_rotations.md`](G65_spin_half_needs_3d_rotations.md) — 自旋 1/2 的第二条路：为什么 **2 维旋转不够、3 维旋转才够**

### 几何扇区的 SU(2) 双覆盖

- [`G66_SU2_double_cover_from_geometry.md`](G66_SU2_double_cover_from_geometry.md) — 几何扇区的 $SU(2)$ 双覆盖：把 [`G11`](G11_dimension_as_consistency.md) 的反射 $\mathbb Z\_2$ 接上自旋 1/2
- [`G67_reflection_generates_spin_Z2.md`](G67_reflection_generates_spin_Z2.md) — 反射 $\mathbb Z\_2$ **生成** 自旋 $\mathbb Z\_2$：$(b)=(a)^2$

### 干涉与 Born 的第二条路

- [`G68_interference_from_coarse_graining.md`](G68_interference_from_coarse_graining.md) — 干涉的导出（附 Born 的第二条路：Schur）

### 洛伦兹自旋、B/tau 结算、退相干

- [`G69_lorentz_spin_from_SL2C.md`](G69_lorentz_spin_from_SL2C.md) — 洛伦兹自旋：$SU(2)$ 升级到 $SL(2,\mathbb{C})$
- [`G70_B_and_tau_closing.md`](G70_B_and_tau_closing.md) — 清单里剩下两项的结算：$B$ 的识别 与 $\tau$ 的初条件依赖
- [`G71_decoherence_from_the_terminal_ledger.md`](G71_decoherence_from_the_terminal_ledger.md) — 退相干：干涉**何时**消失（把"$\pi$ 合并哪些路径"变成可判据的动力学）

### kappa1、B 入册、I2a 与量子

- [`G72_kappa1_from_the_ledger.md`](G72_kappa1_from_the_ledger.md) — $\kappa\_1$：D12 的"率"是**约定**，Z3 能给的是**闭式＋界**（**降级声明**）
- [`G73_B_is_an_input.md`](G73_B_is_an_input.md) — $B$ 的识别：**推不出来**——并把 `D17` 的教训用在我自己身上
- [`G74_I2a_meets_quantum_structure.md`](G74_I2a_meets_quantum_structure.md) — I2a 与量子/自旋结构的接口：**代数半免费，尺度才是硬输入**

### 量子几何：模读出与面积律

- [`G75_quantum_geometry_modular_readout.md`](G75_quantum_geometry_modular_readout.md) — 量子几何：**几何是态的模读出**（面积律由宇称打开）
- [`G76_area_law_in_2d.md`](G76_area_law_in_2d.md) — 高维面积律：2D 的 $S\sim L\log L$（临界）vs $S\propto L$（有 gap）

### 交错耦合、3D 面积律、视界热力学

- [`G77_staggered_coupling_from_A5.md`](G77_staggered_coupling_from_A5.md) — 交错耦合**从 Z3 导出**（把 [`G75`](G75_quantum_geometry_modular_readout.md)／[`G76`](G76_area_law_in_2d.md) 的"识别"变成"推导"）
- [`G78_area_law_in_3d.md`](G78_area_law_in_3d.md) — 3D 面积律：$S\propto L^2$（判据是**拟合残差 ＋ 系数稳定性**）
- [`G79_horizon_thermodynamics.md`](G79_horizon_thermodynamics.md) — 视界热力学：**视界 = 纠缠面**（模流的温度零点）

### 全对全耦合、中央荷、B 的识别

- [`G80_all_to_all_age_coupling.md`](G80_all_to_all_age_coupling.md) — 补那个因子 8：**全对全年龄耦合**（闭合词全历史补偿）
- [`G81_central_charge_from_area_density.md`](G81_central_charge_from_area_density.md) — 中央荷：从**面积密度**接 [`D87`](../modular-equilibrium/derivations/D87_brown_henneaux_central_charge.md) 的 $c=6k=\frac{3\ell}{2G}$
- [`G82_B_from_two_independent_Z2.md`](G82_B_from_two_independent_Z2.md) — $B$ 的识别：两条**有独立动机**的约束（年龄奇偶 × 词的取向）

### v_F 记账、整数 level、全对全的导出

- [`G83_the_missing_1_14.md`](G83_the_missing_1_14.md) — 补 G80 剩下的 1.14：**判据人为 ＋ $v\_F$ 记账**
- [`G84_integer_level_resolved.md`](G84_integer_level_resolved.md) — 整数 level 张力的**重述**：$k=\rho v\_F/m^2$，$k=1$ 需 $v\_F=m^2/\rho$
- [`G85_all_to_all_from_closed_walks.md`](G85_all_to_all_from_closed_walks.md) — 从 Z3 显式导出**全对全耦合**（填上 `R-Z-LONG-RANGE-MEMORY-GAP`）

### 配对核直接定义、xi 判据、测量即典型性

- [`G86_pairing_kernel_from_A5_directly.md`](G86_pairing_kernel_from_A5_directly.md) — 配对核的**直接定义**（Z3 的账本），绕过 `D244` 的 no-go
- [`G87_xi_criterion_in_2d_large.md`](G87_xi_criterion_in_2d_large.md) — 大尺寸 2D 面积律：判据是**比值** $\xi/L$（更正 [`G83`](G83_the_missing_1_14.md) 的绝对形式）
- [`G88_measurement_as_typicality.md`](G88_measurement_as_typicality.md) — 测量诠释的重新审视：**典型性（导出）＋ 第一人称（输入）**

### 维数 no-go 与平衡条件

- [`G89_dimension_no_go_and_the_balance_condition.md`](G89_dimension_no_go_and_the_balance_condition.md) — 维数的分类级 no-go 与「极化两标签无偏好」条件

### Zero 系列参考分诊

- [`G90_zero_series_reference_triage.md`](G90_zero_series_reference_triage.md) — Zero 系列参考分诊：哪些已承接、哪些要打问号、哪些是纯历史

### 核验脚本（88 个）

[`ledger_sync.py`](ledger_sync.py) 是**同步入口**：按 mtime 缓存逐个跑 `G*_check.py` 与 `Z*_check.py`（未变者不重跑），写回账本合计，再重生本清单。

```
python3 ledger_sync.py            # 增量同步（用缓存）
python3 ledger_sync.py --all      # 忽略缓存，全部重跑
```

账本合计：**独立实断言 3027 / 不符 0**（依赖上文的结论行以 `[i]` 单列，不计入合计；见 [`G10_final_derivation_and_input_ledger.md`](G10_final_derivation_and_input_ledger.md) §6）。

---

## 6 `D2xx` 零层弧（50 篇）

编号**保持不变**（= 母项目中的同一编号）。★ 标记的 8 篇含**非原生桥接结构**（见 [`G26`](G26_scope_and_non_native_structures.md)／[`G27`](G27_purification_attempt.md)）。

**范围**：本层是**借入的旧理论零层弧**。其中的外部公理体系（U 系）字样一律为**否证性声明**（"这不是从它推出的"）或**接口审计标的**（[`D212`](D212_zero_universe_to_u_interface.md)／[`D215`](D215_cycle_local_semantic_u_interface.md)，二者已自行放弃"无条件导出"的主张）。**本体系的推理链不使用它们**；本层的公理只有 [`Z0`](Z0_zero_never_rests_single_axiom.md) 一条。

| 文档 | 标题 |
|:--|:--|
| [`D_arc/D210_zero_sum_closure_filter.md`](D_arc/D210_zero_sum_closure_filter.md) | 零和闭合筛选：无有限预算与无通道权重的谱系增长 |
| [`D_arc/D211_global_static_closure_zero_layer.md`](D_arc/D211_global_static_closure_zero_layer.md) | 全局静态闭合零层：多重零记录、局部历史与全局读回边界 |
| [`D_arc/D212_zero_universe_to_u_interface.md`](D_arc/D212_zero_universe_to_u_interface.md) | 周期零宇宙的载体–态–支持三层接口 |
| [`D_arc/D213_periodic_skeleton_to_gr_direct_audit.md`](D_arc/D213_periodic_skeleton_to_gr_direct_audit.md) | 周期零宇宙骨架到 GR 的直接桥审计 |
| [`D_arc/D214_local_zero_sum_transport.md`](D_arc/D214_local_zero_sum_transport.md) | 局域零和传输：闭环邻接上的最小补偿 |
| [`D_arc/D215_cycle_local_semantic_u_interface.md`](D_arc/D215_cycle_local_semantic_u_interface.md) | 周期局域环的语义接口 |
| [`D_arc/D216_history_phase_state.md`](D_arc/D216_history_phase_state.md) | 精确历史相位计数：忠实态的最小选择器 |
| [`D_arc/D217_history_phase_asymptotic_no_go.md`](D_arc/D217_history_phase_asymptotic_no_go.md) | 历史相位态的渐近障碍：暂态非中心态回到中心 |
| [`D_arc/D218_active_phase_occupancy_and_d_boundary_state.md`](D_arc/D218_active_phase_occupancy_and_d_boundary_state.md) | 活动相位占用与 $D$ 层边界：稳定非中心态候选 |
| [`D_arc/D219_active_phase_seed_independence.md`](D_arc/D219_active_phase_seed_independence.md) | 活动相位结构的种子无关性：从任意闭合历史到 $(1,1,2)$ |
| [`D_arc/D220_general_period_active_phase_spectrum.md`](D_arc/D220_general_period_active_phase_spectrum.md) | 一般周期的活动相位谱：中心二项系数与忠实态族 |
| [`D_arc/D221_cyclic_naturality_and_primitive_phase_gap.md`](D_arc/D221_cyclic_naturality_and_primitive_phase_gap.md) | 载体代数的循环自然性：交叉积表示与原始相位障碍 |
| [`D_arc/D222_stratified_destruction_and_local_memory.md`](D_arc/D222_stratified_destruction_and_local_memory.md) | 分层毁灭与局部记忆：活动亚层全清、历史层保留最高两层 |
| [`D_arc/D223_age_carrier_tensor_interface.md`](D_arc/D223_age_carrier_tensor_interface.md) | 年龄载体与矩阵载体的张量接口：毁灭周期负责支持协变，矩阵因子负责载体代数／忠实态 |
| [`D_arc/D224_minimal_m2_age_bridge_and_modular_support_independence.md`](D_arc/D224_minimal_m2_age_bridge_and_modular_support_independence.md) | 最小年龄矩阵桥：$M\_2$ 因子、D220 年龄权重与模流支持寿命分离 |
| [`D_arc/D225_tensor_vs_direct_sum_factorization_gap.md`](D_arc/D225_tensor_vs_direct_sum_factorization_gap.md) | 张量因子与直接和替代：三层条款不选择局部耦合方式 |
| [`D_arc/D226_uniform_matrix_coexistence_selector.md`](D_arc/D226_uniform_matrix_coexistence_selector.md) | 矩阵均匀共存原则：在张量与直接和之间选出张量载体 |
| [`D_arc/D227_age_local_modular_witness_selector.md`](D_arc/D227_age_local_modular_witness_selector.md) | 年龄扇区局部模见证：把均匀共存改写为 GR 侧可检验选择器 |
| [`D_arc/D228_age_local_modular_family_support_covariance.md`](D_arc/D228_age_local_modular_family_support_covariance.md) | 年龄局部模生成元族：从逐扇区见证到支持协变 |
| [`D_arc/D229_continuous_age_support_obstruction.md`](D_arc/D229_continuous_age_support_obstruction.md) | 连续年龄支持的原子障碍：为什么不能把有限扇区直接取连续极限 |
| [`D_arc/D230_linfinity_age_support_and_projection_complement.md`](D_arc/D230_linfinity_age_support_and_projection_complement.md) | 连续年龄支持的 $L^\infty$ 构造：投影压缩加补单位 |
| [`D_arc/D231_modular_density_profile_gap.md`](D_arc/D231_modular_density_profile_gap.md) | 局部模密度剖面缺口：区间支持不等于几何 boost |
| [`D_arc/D232_profile_as_matrix_age_correlation.md`](D_arc/D232_profile_as_matrix_age_correlation.md) | 模密度剖面就是矩阵相位的年龄相关性 |
| [`D_arc/D233_sign_age_symmetry_no_go_for_profile.md`](D_arc/D233_sign_age_symmetry_no_go_for_profile.md) | 符号年龄对称性无解：现有重播种不能产生非恒定剖面 |
| [`D_arc/D234_geometric_ball_profile_candidate.md`](D_arc/D234_geometric_ball_profile_candidate.md) | 几何球模核参照：抛物型剖面候选与未闭合的年龄映射 |
| [`D_arc/D235_age_radial_reparametrization_no_go.md`](D_arc/D235_age_radial_reparametrization_no_go.md) | 年龄到径向映射的重参数化障碍：二次年龄剖面不是上游结论 |
| [`D_arc/D236_source_operator_identification_embedding.md`](D_arc/D236_source_operator_identification_embedding.md) | 源算子识别的条件嵌入与不唯一性 |
| [`D_arc/D237_contact_channel_single_profile_reduction.md`](D_arc/D237_contact_channel_single_profile_reduction.md) | 接触项与单剖面归约：中央通道、平行通道与独立通道 |
| [`D_arc/D238_linear_sign_hazard_mechanism.md`](D_arc/D238_linear_sign_hazard_mechanism.md) | 线性符号危险率机制：抛物型年龄比例的最小条件实现 |
| [`D_arc/D239_quadratic_sign_potential_audit.md`](D_arc/D239_quadratic_sign_potential_audit.md) | 线性危险率的二次势审计：中心临界、常曲率与年龄配对 |
| [`D_arc/D240_additive_record_no_go.md`](D_arc/D240_additive_record_no_go.md) | 可加记录无解：静态零层为什么不能生成二次符号势 |
| [`D_arc/D241_sign_blind_pairing_no_go.md`](D_arc/D241_sign_blind_pairing_no_go.md) | 符号盲配对无解：时间方向不等于正负符号耦合 |
| [`D_arc/D242_finite_memory_pairing_no_go.md`](D_arc/D242_finite_memory_pairing_no_go.md) | 有限记忆配对无解：两层历史不能生成二次符号势 |
| [`D_arc/D243_global_readback_not_interaction.md`](D_arc/D243_global_readback_not_interaction.md) | 全局读回不是相互作用：零层可见性不产生全对全耦合 |
| [`D_arc/D244_zero_sum_matching_linearity.md`](D_arc/D244_zero_sum_matching_linearity.md) | 零和匹配线性性：闭合词配对本身不给全对全 |
| [`D_arc/D245_fixed_pair_kernel_uniqueness.md`](D_arc/D245_fixed_pair_kernel_uniqueness.md) | 固定二体核唯一性：二次势把全对全核逼成常数 |
| [`D_arc/D246_collective_zero_mode_rank_one.md`](D_arc/D246_collective_zero_mode_rank_one.md) | 全对全核的集体零模表示：从 $k^2$ 条边到单个 $Q^2$ |
| [`D_arc/D247_zero_defect_stiffness_scale_audit.md`](D_arc/D247_zero_defect_stiffness_scale_audit.md) | 零缺陷刚度尺度审计：形状与数值必须分开 |
| [`D_arc/D248_interlayer_readout_dynamics.md`](D_arc/D248_interlayer_readout_dynamics.md) | 层间读出动力学：从历史与零层到守恒源 |
| [`D_arc/D249_minimal_stress_lift_from_interlayer_scalar.md`](D_arc/D249_minimal_stress_lift_from_interlayer_scalar.md) | 最小应力提升：层间标量流、刚性物态与应力类别缺口 |
| [`D_arc/D250_pure_exchange_selects_stiff_scalar_source.md`](D_arc/D250_pure_exchange_selects_stiff_scalar_source.md) | 纯层间交换选择刚性标量源：无势 Dirichlet、常量零模与 onsite 缺口 |
| [`D_arc/D251_layer_type_complex_and_local_readout_presheaf.md`](D_arc/D251_layer_type_complex_and_local_readout_presheaf.md) | 层类型复合与局域读出预层：物理局域接口的最小底座 |
| [`D_arc/D252_layer_time_atlas_and_global_time_potential.md`](D_arc/D252_layer_time_atlas_and_global_time_potential.md) | 分层时间图册与全局时间势：定向可行条件与 Lorentz 锥缺口 |
| [`D_arc/D253_adm_metric_assembly_and_lapse_shift_gap.md`](D_arc/D253_adm_metric_assembly_and_lapse_shift_gap.md) | 叶状 ADM 度规组装与 lapse/shift 缺口 |
| [`D_arc/D254_dirichlet_tensor_leaf_metric_and_clock_gauge.md`](D_arc/D254_dirichlet_tensor_leaf_metric_and_clock_gauge.md) | Dirichlet 张量到叶层度规与读回时钟规范 |
| [`D_arc/D255_tetrahedral_dirichlet_assembly_and_gluing.md`](D_arc/D255_tetrahedral_dirichlet_assembly_and_gluing.md) | 交换权重的四面体 Dirichlet 组装与面胶合 |
| [`D_arc/D256_clock_gauge_intrinsic_simplex_and_refinement_convergence.md`](D_arc/D256_clock_gauge_intrinsic_simplex_and_refinement_convergence.md) | 时钟规范拆分、内禀单纯几何与细化收敛条件 |
| [`D_arc/D257_resistance_metric_fixed_point_and_locality_gap.md`](D_arc/D257_resistance_metric_fixed_point_and_locality_gap.md) | 有效电阻度量、交换权重固定点与局部性缺口 |
| [`D_arc/D258_full_exchange_projection_and_local_isotropy_obstruction.md`](D_arc/D258_full_exchange_projection_and_local_isotropy_obstruction.md) | 全边交换投影、$K\_4$ 普遍固定点与局部各向异性障碍 |
| [`D_arc/D259_z_layer_condition_dictionary_and_four_dimensional_confluence.md`](D_arc/D259_z_layer_condition_dictionary_and_four_dimensional_confluence.md) | Z 层条件字典与四维时空的汇流构造 |

---

## 7 工具

| 文件 | 用途 |
|:--|:--|
| [`latex_lint.py`](latex_lint.py) | 扫描 13 类 LaTeX 问题（C1–C12） |
| [`latex_fix.py`](latex_fix.py) | 修复器（默认 dry-run，`--apply` 才写盘） |
| [`LATEX_AUDIT.md`](LATEX_AUDIT.md) | LaTeX 审计记录 |
| [`INDEX_check.py`](INDEX_check.py) | 核验本清单：无孤儿文档 / 计数一致 / 范围声明 |
| [`STATUS_check.py`](STATUS_check.py) | 核验唯一状态源、关键页指向、旧状态措辞清零 |
| [`ledger_sync.py`](ledger_sync.py) | 同步入口：增量重跑 `G*_check.py`、写回账本合计、重生本清单 |
| [`gen_index.py`](gen_index.py) | 重生成本清单（标题与合计数自动取自文件系统与账本） |

---

## 8 建议的阅读顺序（**自底向上**）

| 顺序 | 读什么 | 为什么 |
|--:|:--|:--|
| 0 | [`STATUS.md`](STATUS.md) | **唯一当前状态源**：先分清主路线、D259 路线与历史替代关系 |
| 1 | `zero_sum_*`（见 §3 的 7 篇笔记） | **Zero 原始层**：先有仿真，后有公理 |
| 2 | [`Z0`](Z0_zero_never_rests_single_axiom.md) | **唯一公理**：零不断乱动（三款各买一件东西） |
| 3 | [`Z1`](Z1_zero_layer_as_the_foundation.md) | 取代 A1–A5 ＋ Zero 结构的扩充（含散度字典） |
| 4 | [`Z2`](Z2_zero_to_gr_direct_route.md) | **直连路线**：从 Zero 到 GR（免去 A0–A5（历史命名）桥接；标出两处残余输入） |
| 5 | [`G0`](G0_bottom_layer_and_derivation_route.md) | **底层条款表**（Z0 条款 ＋ Z1–Z5 定理；A0–A5 历史命名）＋ 路线 ＋ 账本 I1–I12 |
| 6 | [`G19`](G19_axiom_reduction.md) | 条款精简：6 条 A 条款 → 3 条独立定理（公理只有 Z0） |
| 7 | [`G1`](G1_derivations_from_the_bottom_layer.md) → [`G8`](G8_dimension_selection.md) | 几何骨架 → 维数 |
| 8 | [`G13`](G13_foliation_and_lorentz_invariance_gap.md) → [`G15`](G15_bare_ax3_has_no_characteristic_speed.md) | 叶层与洛伦兹不变性（含 no-go） |
| 9 | [`G16`](G16_repair_audit_without_new_axioms.md) | 不增扩充条款清算：**主定理存活** |
| 10 | [`G23`](G23_zero_layer_structure_inventory.md) → [`G25`](G25_age_to_geometry_channel_is_obstructed.md) | 语料审计、范围 |
| 11 | [`G27`](G27_purification_attempt.md) → [`G29`](G29_probability_as_derived_not_postulated.md) | 净化与概率的导出 |
| 12 | [`G28`](G28_dynamics_audit.md) | 动力学完成度（含闭环图平均度无界） |
| 13 | [`G90`](G90_zero_series_reference_triage.md) | Zero 系列分诊：哪些已承接、哪些待造 |

---

## 9 如需物理分目录（可选方案）

本清单已解决**归属**问题。若还要物理分目录，最小改动方案是：

```
lh/
  INDEX.md            ← 本文件
  G/                  ← 全部 G*.md + G*_check.py
  D_zero_arc/         ← 全部 D2xx_*.md
  zero_sum_sims/      ← 全部 zero_sum_* / zero_*
  tools/              ← latex_lint.py / latex_fix.py / LATEX_AUDIT.md
```
代价：需改 2 个核验脚本里的路径（`G20_check.py` 读 `zero_sum_*`、`G22_check.py` 读 `D2xx`），其余脚本用 `MOD` 绝对路径不受影响。
