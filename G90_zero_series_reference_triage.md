# G90 · Zero 系列参考分诊：哪些已承接、哪些要打问号、哪些是纯历史

**日期**：本轮 · **范围**：`zero_sum_*` **7 篇笔记** ＋ **10 个实验**（脚本与其结果文件）。
**判据**：**【承接】**＝它的结论（或其对照实验／其提出的缺口）已在 G 系列里作为**推导依据或登记项**出现；
**【打问号】**＝它提出了 G 系列尚未承接的对象或问题，且该对象至今不存在；
**【纯历史】**＝探索性玩具模型，无笔记、且未被任何 G/D 文档引用。
**等级标签**：【承接】/【打问号】/【纯历史】/【数值证据】/【no-go】。
**核验**：[`G90_check.py`](G90_check.py) —— **独立实断言 47 / 结论行 0 / 不符 0**，退出码 `0`（复算引用数 ＋ 真跑实验）

$$
\ \text{净承接 5 项（含 1 条推导、1 条输入）、打问号 2 项、纯历史 3 项；10 个实验迁入后全部可复跑。}\ 
$$

---

## §0 为什么补这份分诊

D 系列有 [`G9`](G9_d_series_reference_triage.md) 的**逐篇分诊**；Zero 系列**没有对应物**。
[`G20`](G20_axiom_audit_extended_to_zero_and_D.md) 把它当作**语料整体**做公理审计（并据此登记了 **I8**），
但从未逐条结算"7 篇笔记 ＋ 10 个实验各自贡献了什么"。

本文补上这一层——**这是 Zero 系列作为"未结算资产"的最后一块**。

## §1 事实基线（全部可复算）

### §1.1 来源与可复跑性

| 项 | 事实 |
|:--|:--|
| 来源 | 逐字节复制自 `../modular-equilibrium/simulations/`（md5 相同） |
| 原址为何能跑 | 脚本用 `ROOT = Path(__file__).resolve().parents[1]`；在原址 `simulations/` 下解析为 `modular-equilibrium/` |
| 复制到 `lh/` **根目录**后 | 解析为 `/Users/oygb/Downloads/`（不存在）⇒ **10 个实验全部 `rc=1`**（`FileNotFoundError`） |
| 现状 | 已置入 [`simulations/`](simulations/README.md)，`ROOT` 恰好解析为 `lh/`、输出写回 `simulations/` ⇒ **10/10 恢复可跑**（沙箱实测 1–25 秒，退出码全 0） |

### §1.2 被引次数

口径：字符串 `zero_sum_<名>` 在 `G*.md` ＋ `D*.md` 中的出现次数，**不含本文自身**
（本文一旦提及某实验，就会改变被统计量——这一条由 `G90_check.py` F2 强制）。
所列数字是**本版实算值、作为下界**：后续文档（如 [`Z1`](Z1_zero_layer_as_the_foundation.md)）只会使其增加，
`G90_check.py` 按 **≥** 校验（因此虚报会被抓、增长不会误报）。
**另排除"关于 Zero 层本身的文档"**（[`Z1`](Z1_zero_layer_as_the_foundation.md)、[`Z0`](Z0_zero_never_rests_single_axiom.md)）：
它们的职责就是逐条点名 Zero 组件（Zero 成为基础后必然如此），引用它们**不表示理论脊柱在用**。

| 实验 | 有独立笔记 | 被引 | 判定 |
|:--|:--:|--:|:--|
| `closure_exit` | ✅ | **7** | 【承接】 |
| `periodic_destruction` | ✅ | **6** | 【承接】 |
| `global_R_local_P` | ✅ | **6** | 【承接】 |
| `open_reservoir` | ✅ | **4** | 【承接】 |
| `tri_layer_universe` | ✅ | **4** | 【承接】 |
| `closure_time_selection` | ✅ | 0 | 【打问号】 |
| `reproduction_audit` | ✅ | 0 | 【打问号】 |
| `cycle_evolution` | ❌ | 0 | 【纯历史】 |
| `geometry_probe` | ❌ | 0（**但 G28 §4 用了它的数据，原未引其文件**） | **【承接】**（本版更正） |
| `living_universe` | ❌ | 0 | 【纯历史】 |

**提及 `zero_sum` 的 G 文档共 6 篇**：`G0`(6)、`G20`(4)、`G22`(5)、`G21`(3)、`G28`(1)、`G34`(1)。
（另有 `G57`、`G86` 各 1 处，指的是 **D 系列文件名**里的 `zero_sum_transport`／`zero_sum_matching_linearity`，不计入本表。）

**一处结构性事实（本版修正）**：3 个没有笔记的实验中，**2 个**（`cycle_evolution`、`living_universe`）确为零引用；
**第 3 个 `geometry_probe` 是例外**——它零引用是因为 **`G28` 漏引**，不是因为它没被用（本版已补引）。写笔记这个动作本身，就是它们被承接的前提。

## §2 【承接】5 项——它们具体进了哪里

### 2.1 `closure_exit` → `G28` → **`G34`（一条推导）**

笔记 §1 设了两条规则做对照：

| 规则 | 行为 |
|:--|:--|
| `continue_after_closure` | 闭合路径**仍留在**活动层，可继续生成后继 |
| `exit_after_closure` | 闭合路径**退出**活动层，写入 $P\_i$ 与全局 $\mathcal Z\_\ast$ |

**这个对照实验是 [`G34`](G34_exit_rule_absorbed_into_pi.md) 的承重依据**：它显示"继续规则"下闭合路径会**再次闭合**、产生**更多**闭合事件，
由此证明「**闭合即退出**」不是独立的动力学规则，而是 **$\pi$ 的一条定义**。
[`G28`](G28_dynamics_audit.md) 先把这个自述查出并登记，`G34` 再把它归并进 $\pi$。

### 2.2 `periodic_destruction` ＋ `open_reservoir` → `G20` §3 → **`I8`（一条输入）**

- `periodic_destruction` 的三规则对照（`continuous` / `wipe_no_reseed` / `wipe_reseed`）给出**确定性重播种**
  $w\longmapsto\{w+,w-\}$，周期 $T=3$；
- `open_reservoir` 反向指出：`closure_exit` 的活动层之所以会清空，是**两个外加条件**造成的
  （① 固定寿命 $L=4$；② 闭合后不再播种）——**它定位了缺的那个子句**。

两者合起来 ⇒ [`G20`](G20_axiom_audit_extended_to_zero_and_D.md) §3 引理 74：**Z3（历史标号 A5）缺「再播种」子句**（重播种现由 **Z5** 导出）。
按"不增公理"纪律**没有**写进条款表，而是登记为账本 **`I8`**（输入／模型选择）。

### 2.3 `global_R_local_P` → `D211` 的术语更新 → G 系列的 $P\_i$／$\mathcal Z\_\ast$

笔记的猜想是"$R$ 全局、$P$ 局部"，`D211` 把 $R$ 重新解释为**全局静态闭合零层** $\mathcal Z\_\ast$ 的结果索引。
G 系列**确实在用这套词汇**（$P\_i$、$\mathcal Z\_\ast$ 出现在 `G0`、`G20`、`G22`、`G23`、`G25`、`G27`、`G34`、`G49` 等）。

> **但要说清承接路径**：[`G23`](G23_zero_layer_structure_inventory.md) 的基础结构清点是**从 D 系列**抽取的
> （抽取 D 系列每篇的"预先结构"字段），**不是**从 Zero 笔记抽取的。所以这条承接是**经 D211 中转**的。

### 2.4 `tri_layer_universe` → `G22`／`G21` 的 **locality** 问题

笔记结尾自述：*"The next missing structure is **locality**. A single global pool still mixes everything."*
[`G22`](G22_correction_sublayers_and_local_layers.md) 与 [`G21`](G21_do_the_layers_help_derive_GR.md) 均引用这句——
即 **"单全局池"这一诊断被 G 系列接住**（层结构之所以需要局域子层，源头在这里）。

## §3 【打问号】2 项——提出了对象，但至今没有

### 3.1 `reproduction_audit`：**一个被查出的隐藏假设，后继对象始终没造出来**

笔记记录了一次**自我审计**：先前沙盒把

$$
r(w)=\#\{\text{proper returns of }w\text{ to zero}\},\qquad M=1+r\ \text{或}\ M=2^{r}
$$

混用，而 $M=1+r$ 等于**默认每条分支都继承父代的全部繁殖能力**——笔记称之为
*"a hidden copy assumption"*。它给出的结论是：

> *"The next required object is not another dimension estimate. It is a **zero-preserving reproduction or autocatalytic rule** that actually creates new reproductive capacity."*

**这个对象至今不存在**（G 系列没有繁殖／自催化规则这一节）。所以判定【打问号】，不是【承接】。

### 3.2 `closure_time_selection`：**"为什么不同闭合时间能共存"没有后继**

笔记从"全模式在场 ＋ 更快闭合 ⇒ 更多完成周期 ⇒ 更多后代"重建模型，结果是：

> *"mode dominates almost completely."*

它自己指出出路不是再加一个 fitness 参数，而是 *"locality, compatibility, or an autocatalytic network"*。
**locality 这一支被 `G22`／`G21` 接住（见 §2.4）**，但它真正问的那个问题
（**不同闭合时间为何能共存**）在 G 系列里**没有对应结论**。故判【打问号】。

## §4 【纯历史】3 项——无笔记、零引用、自述为独立玩具

| 实验 | 脚本自述（逐字） |
|:--|:--|
| `cycle_evolution` | *"This is an exploratory model. It does not use any D\* derivation."* |
| ~~`geometry_probe`~~ | ~~*"This sandbox asks a narrow question"*~~ → **本版更正：它是【承接】，见下 |
| `living_universe` | 「这是一个独立玩具模型，不依赖 `derivations/D*` 的结论。」 |

**两篇**（`cycle_evolution`、`living_universe`）没有独立笔记、被引次数为 **0**，自述即为探索性模型。

> **本版更正**：`geometry_probe` **不属于【纯历史】**——[`G28`](G28_dynamics_audit.md) §4 的"闭环图平均度无界"（第四条独立障碍）**用的就是它的数据**，只是**原稿没有引用它的文件**。
> 本轮已在 `G28` §4 补上引用，并据此把 [`Z0`](Z0_zero_never_rests_single_axiom.md) §2.6 的**Z4 终端款**从"残余假设"改为【导出】。
> **这个例外恰好证明分诊的价值**：它抓出了一处**未登记依赖**。
其中 `living_universe_results.json` 是 **3.4 MB**，为全语料最大文件，**从未被任何文档消费**。

> **一处必须区分的细节**：`cycle_evolution` 对**理论**零引用，但它同时是**本层的公共基础设施**——
> 另外 **7 个实验**都 `from zero_sum_cycle_evolution import canonical_cycle`（`canonical_cycle(word)`
> 取词的最小旋转，即"闭包类"的规范代表）。这也是 10 个实验**必须成套搬运**的原因：
> 只搬其中几个会在 `import` 处直接失败。所以它的准确定位是：
> **【纯历史】（对理论）＋【层内基础设施】（对 Zero 层）**。

**处置建议**：保留为历史档案（它们仍可复跑，见 §1.1），但**不再计入"理论资产"**；
若日后要用，须先补笔记并给出它与 Z0 条款（A0–A5 历史命名）的关系。

## §5 一条时间线观察（值得记录，但不构成承接）

`closure_exit` 笔记（Sep 30）在模型里**设定**了「有限寿命 $L=4$」：

> 「未闭合路径继续分支。超过有限寿命 $L=4$ 的开放路径进入 $D\_i$。」

而 G 系列后来一度把 $L=4$ 写成锁定候选（[`G46`](G46_k_is_the_lifetime.md) 逼出 $k=L$，[`G61`](G61_locking_the_five_integers.md) 给出最小性候选），
且 `G46`／`G61` **不引用** Zero 系列。当前更正见 [`G73`](G73_B_is_an_input.md)：去掉“取最小”后可行寿命为 `{4,6,8,...}`，故 $L=4$ 是**选择原则／输入**，不是 G 系列的无条件导出。

**所以**：值 $4$ 在 Zero 系列里是**脚本假设**，在 G 系列早期是**条件候选**，最终由 `G73` 明确降为**选择原则／输入**；
但**没有证据**表明后者由前者启发。诚实写法是：**同一个数在两代里出现，但当前都不提升为无条件导出**，
这条时间线只作记录，**不登记为承接**。

## §6 净结论

| 判定 | 项数 | 内容 |
|:--|--:|:--|
| **【承接】** | **6** | `closure_exit`（→ G34 一条推导）、`periodic_destruction` ＋ `open_reservoir`（→ I8 一条输入）、`global_R_local_P`（→ 词汇，经 D211 中转）、`tri_layer_universe`（→ locality 问题）、**`geometry_probe`（→ G28 §4 第四条障碍 ＋ Z0 §2.6 的 Z4 终端款导出）** |
| **【打问号】** | **2** | `reproduction_audit`（零和保持的繁殖／自催化规则，**至今缺**）、`closure_time_selection`（不同闭合时间共存，**无后继**） |
| **【纯历史】** | **2** | `cycle_evolution`、`living_universe` |

$$
\ \text{Zero 系列的真实产出 = 2 条推导（G34＋Z4 终端款）＋ 1 条输入（I8）＋ 1 个问题（locality）＋ 2 个待造对象。}\ 
$$

**开放项**

1. **零和保持的繁殖／自催化规则**（§3.1）——`reproduction_audit` 点名要的那个对象。
2. **不同闭合时间的共存机制**（§3.2）。
3. **3 个纯历史实验的正式归档决定**（§4）——建议标为【纯历史】并停止计入资产。
4. **$L=4$ 的两代关系**（§5）——若要主张启发关系，需要独立证据；目前**不主张**。

**诚实边界**

- 本文的分诊依据是：**笔记自述文本**（7 篇全部通读）＋ **引用计数**（可复算）＋ **脚本自述与实跑**（10 个全部实跑）。
  **我没有**逐行审读这 10 个脚本（合计约 15 万字符）；**也没有**核验它们各自主张的**数值结论**是否成立
  —— 那属于"复算 Zero 系列"，是另一件事（`G20` 已把它们的等级定为**【数值证据】**，不是【导出】）。
- 【承接】的判定标准是"**它的结论或对照实验出现在 G 系列的推导／登记里**"，而不是"它启发了作者"。
  启发关系通常不可考，本文一律不主张。
- 引用计数只统计 `G*.md` 与 `D*.md`，且**排除本文自身**；`simulations/README.md`、`verify/README.md`
  与 `G10` 的说明文字**不计入**。这条排除是必需的：本文正文提到 `zero_sum_cycle_evolution` 一次，
  若不排除，它就会从 0 变成 1，分诊表将无法自洽（`G90_check.py` 首轮就是这样失败的）。
