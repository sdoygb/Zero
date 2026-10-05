# R9 · 外部“恢复 Einstein 方程或 GR”路线再排序

**日期**：2026-10-02  
**范围**：2015-2026 年主流可核验工作，另保留 Jacobson 1995 作为历史控制。  
**性质**：外部路线审计，不修改 `STATUS.md`，不把外部论文的已有证明算作 Zero 层新增定理。  
**核验**：[`R9_check.py`](R9_check.py)。

$$

\begin{aligned}
&\text{R8 把 Jacobson 2016 作为第一接口的选择仍然成立；}\\
&\text{但次选应从 Jacobson 1995 改为 Cao-Carroll 2018；}\\
&\text{Oh-Park-Sin 和 Faulkner 等全息强结果不能当 Zero 的首选接口；}\\
&\text{Zero 目前仍只支持条件恢复，不支持无条件导出四维 GR。}
\end{aligned}
$$

---

## §0 直接结论

### 0.1 R8 的主目标是否需要改变

**不需要。** Jacobson 2016 仍是 Zero 基础最直接、最省跨领域的接口：

| Zero 已有对象 | Jacobson 2016 所需对象 | 对接程度 |
|:--|:--|:--|
| 有限维忠实态、GNS 模流、模 Hamiltonian | 小球约化态与 $K\_B=-\log\rho\_B$ | 直接对应 |
| R8.1 精确熵差恒等式 | 一阶纠缠第一定律 | 已补到有限维任意忠实态 |
| 面积律候选与边界退化 | 小球面积的 UV 项 | 条件对应 |
| 固定体积边界权重 | 固定体积熵平衡 | 只有有限模型接口 |
| 离散细化与局部极限 | 连续 Lorentzian 区域代数网 | 未闭合 |
| 模流候选剖面 | 几何 boost $K\_B\to 2\pi B\_B$ | 决定性缺口 |

Jacobson 2016 的结构把问题集中到一条链上：

$$
\text{有限维态}
\longrightarrow
\text{模流}
\longrightarrow
K_B
\longrightarrow
2\pi B_B
\longrightarrow
\text{小球纠缠平衡}
\longrightarrow
\text{场方程}.
$$

其他路线要么事先加入 AdS/CFT、Ryu-Takayanagi 或连续 CFT，要么只能给出弱场、线性化、离散约束或修改重力。它们不是不能研究，而是会把 Zero 最缺的前提改名成更大的外部输入。

### 0.2 次选应当改选谁

**次选改为 Cao-Carroll 2018**：

> Cao, Carroll, *Bulk Entanglement Gravity without a Boundary: Towards Finding Einstein's Equation in Hilbert Space*, Phys. Rev. D 97, 086003 (2018), [DOI](https://doi.org/10.1103/PhysRevD.97.086003), [arXiv:1712.02803](https://arxiv.org/abs/1712.02803)。

原因是它从抽象 Hilbert 空间的张量分解、互信息图和面积数据出发，而不是从已给定的 Lorentzian 局部量子场论出发。它只推出弱场 Einstein 方程，不冒充完整非线性 GR，但它与 Zero 的离散因子、计数推前态、信息图和度规候选有实质接口。

原先的次选 Jacobson 1995 仍重要，但更适合作为“经典局部 Rindler 极限的历史控制”，不适合作为 Zero 的下一条主接口。因为它一开始就假设局部 Rindler 视界、Unruh 温度、面积熵和 Clausius 平衡，几乎正好绕过本项目最想自己构造的对象。

### 0.3 结论排序

| 排名 | 外部路线 | 判决 |
|--:|:--|:--|
| 1 | Jacobson 2016 | Zero 的最佳主接口；保留 R8 选择 |
| 2 | Cao-Carroll 2018 | Zero 的最佳次接口；从有限维 Hilbert 空间出发，但只到弱场 |
| 3 | Faulkner et al. 2017 | 外部最强非线性纠缠推导之一；连续 CFT 与 HRRT 输入过重 |
| 4 | Oh-Park-Sin 2017/2018 | 能推出完整非线性 EFE，但核心输入是全阶 GDERE，不适合当 Zero 首目标 |
| 5 | Gorard 2020 + Wolfram 2020 | 离散骨架很有价值，但固定维数与连续极限仍是关键输入 |
| 6 | Bianconi 2025 | 最近的最强相对熵作用量路线之一；起点是 Lorentzian 度规，不是 Zero 原语 |

**一句话**：Jacobson 2016 解决“模流怎样变成局部 boost”；Cao-Carroll 2018 解决“有限维纠缠数据怎样变成空间几何”。这两篇正好补 Zero 的两条最薄弱但可攻击的链。

---

## §1 外部来源核查方法

本页只把论文正文、期刊页、Crossref 元数据和 arXiv 元数据当作证据，不把评论文章、博客或二手摘要当作结论来源。

核查步骤：

1. 对每篇目标论文查询 Crossref DOI 记录，确认标题、作者和出版时间。
2. 对每篇有 arXiv 版本的论文查询 arXiv abstract 页及其 citation metadata。
3. 在能取得 PDF 时，直接核对摘要、假设表和结论段。
4. 对“完整非线性 EFE”“弱场”“线性化”“微扰到二阶”“修改重力”等强弱判断，只采用论文自己的表述。
5. 尝试访问 `https://export.arxiv.org/api/query` 时，执行环境返回了 `Rate exceeded`；因此最终逐篇核对改用 arXiv 页面 metadata 与 Crossref DOI 记录，并保留全部可核验链接。

外部论文的关系图如下：

```text
Jacobson 1995
  局部 Rindler + Clausius
        |
        v
Jacobson 2016
  小球模块平衡
        |
        +--> Casini-Galante-Myers 2016 / Speranza 2016
        |        低维相关算符与固定体积假设的批评
        |
        +--> Bueno-Min-Speranza-Visser 2017
        |        高阶引力只到线性化
        |
        +--> Alonso-Serrano-Liska 2020
                 热力学路线更自然得到 unimodular gravity

CFT / 全息支线
  Faulkner et al. 2017
  Oh-Park-Sin 2017/2018
  Dong-Lewkowycz 2018
  Leichenauer et al. 2018

抽象 Hilbert 空间 / 离散支线
  Cao-Carroll 2018
  Gorard 2020 + Wolfram 2020
  Carrasco-Pedraza-Svesko-Weller-Davies 2023
  Bianconi 2025
```

---

## §2 六维比较表

表中的 `显式` 表示该论文在自身前提下直接给出；`输入` 表示它是论文的前提；`条件` 表示论文自己保留条件或猜想；`微扰`、`弱场`、`线性化` 按其实际承诺；`不涉及` 表示论文不处理该项。

| 路线 | 连续流形与区域代数 | modular / boost 桥 | 面积密度 | 固定体积平衡 | 完整非线性 EFE | 可检验预测 | 对 Zero 的主要代价 |
|:--|:--|:--|:--|:--|:--|:--|:--|
| Jacobson 1995 | 局部 Rindler 几何为输入 | Rindler boost 为输入 | $S=\eta A$ 为输入 | 用局部 Clausius 平衡替代 | 在热力学前提下完整 | 无直接新增预言 | 绕过 Zero 要构造的局部几何与温度 |
| Jacobson 2016 | 连续 QFT 与测地球为输入 | CFT 模 Hamiltonian 显式；几何 boost 在球极限 | UV 面积项为输入 | 显式，固定体积 | 共形场完整；非共形留猜想 | 无直接新增预言 | 正对模流链，但要求在算子级完成 boost 极限 |
| Casini-Galante-Myers 2016 | 全息 CFT 计算 | 接受 Jacobson 的模 Hamiltonian 框架 | 用全息熵计算检验 | 正是被检验项 | 不提出完整推导 | 无 | 说明低维相关算符可污染固定体积假设 |
| Speranza 2016 | 共形微扰与球区域 | 接受模 Hamiltonian 框架 | 球熵展开显式 | 检验其适用性 | 只给修正，不给完整 EFE | 无 | 给 $R^{2\Delta}$ 项，Zero 必须解释其抑制 |
| Faulkner et al. 2017 | 连续 CFT 为输入 | 相对熵显式 | HRRT 面积表示 | 不是主轴 | 微扰到二阶 | 无直接定量预言 | 缺少连续 CFT 与全球 HRRT 几何表示 |
| Bueno-Min-Speranza-Visser 2017 | 连续因果菱形为输入 | Wald 熵替代普通面积 | Wald 面积熵显式 | 广义体积为输入 | 只到线性化；论文明确说明不能由小球线性化推出完整非线性 | 无 | 需要高阶曲率作用量与 Wald 形式，Zero 目前没有 |
| Oh-Park-Sin 2017/2018 | AdS/CFT 与 RT 面为输入 | 模 Hamiltonian 与 Holland-Wald 形式 | RT 面积与 Iyer-Wald | 广义第一定律 | 在“全阶 GDERE”假设下完整非线性 | 无直接定量预言 | 只是把真正缺口移到全阶引力相对熵 |
| Cao-Carroll 2018 | 从抽象 Hilbert 因子出发，连续几何为输出 | 不走 boost 桥，走全局切割 | 面积与互信息成比例为输入 | 用跨切割平衡替代局部固定体积 | 弱场 | 无直接定量预言 | 需补首选分解、RC 条件、面积比例与 Lorentzian 组装 |
| Leichenauer-Levine-Shahbazi-Moghaddam 2018 | 全息场论为输入 | 用熵的二阶形状变分 | QNEC 型面积/熵关系 | 通过 Null 变形实现 | 仅到 leading order | QNEC 与量子聚焦方向可检验 | 强依赖全息与 QNEC 饱和，不适合 Zero 首步 |
| Dong-Lewkowycz 2018 | Euclidean 引力路径积分为输入 | 量子极值面与体模块 Hamiltonian | 广义熵显式 | 极值条件 | 给出积分形式 EOM，线性化到任意态 | 无直接定量预言 | 需要引力路径积分，Zero 目前只有离散/有限模型 |
| Alonso-Serrano-Liska 2020 | 局部因果菱形为输入 | Clausius 熵与纠缠熵对照 | 面积项为输入 | 比较两类平衡表述 | 热力学自然得到 unimodular gravity，不是完整 GR | 无直接定量预言 | 是重要反例：热力学第一定律不自动给完整 GR |
| Gorard 2020 + Wolfram 2020 | 离散超图与因果图显式；连续流形为极限 | 不处理 modular/boost 桥 | 离散曲率与体积锥显式 | 不是主轴 | 离散 Ricci 约束；固定维数为输入 | 离散 Lorentz 与几何修正 | 需要固定维数、弱遍历性和连续超图极限 |
| Carrasco-Pedraza-Svesko-Weller-Davies 2023 | 全息复杂体积理论为输入 | 不处理 modular/boost 桥 | CV 字典为输入 | 不是主轴 | 线性化；明确承诺限制在二维膨胀引力 | 无直接定量预言 | 复杂度假设比模流离 Zero 更远 |
| Kumar 2024/2025 | 大因果菱形与 de Sitter 静态片为输入 | 近 Rindler 的纠缠第一定律为输入 | $\delta S\_{\rm out}=\delta A/(4G)$ 为输入 | 用量子膨胀与非平衡熵流 | 宣称完整半经典 EFE，但保留多项假设与因子约定 | 无直接定量预言 | 不能算真正“无假设”；适合作为批评案例 |
| Bianconi 2025 | Lorentzian 度规为输入 | 用相对熵作用量，不走模 boost | 不以边界面积为主轴 | 不是固定体积方案 | 修改 EFE；低耦合极限回到 EFE | 有 $G$-场与暗物质角色，但量化不足 | R8.1 可接相对熵，但缺少度规算符与双度规结构 |

---

## §3 Jacobson 2016 的关键批评

### 3.1 固定体积假设并非免费

Casini, Galante, Myers 计算了被相关算符扰动的全息 CFT。对低维相关算符，他们发现熵变分中出现

$$
R^{2\Delta}\delta\langle O_\Delta\rangle^2
$$

型项。当

$$
\Delta\le \frac d2
$$

时，这一项在小球极限中压过通常的 $R^d\delta\langle T\_{00}\rangle$ 项，于是 Jacobson 所需的 $\delta S\_{\rm IR}$ 形式不再自动成立。

来源：Casini, Galante, Myers, JHEP 03 (2016) 194, [DOI](https://doi.org/10.1007/JHEP03(2016)194), [arXiv:1601.00528](https://arxiv.org/abs/1601.00528)。

### 3.2 Speranza 的独立共形微扰检验

Speranza 对相关算符微扰的球区域熵计算得到同一类标度。论文直接指出，当 $\Delta\le d/2$ 时，新项主导并需要修改 Jacobson 的猜想。

来源：Speranza, JHEP 04 (2016) 105, [DOI](https://doi.org/10.1007/JHEP04(2016)105), [arXiv:1602.01380](https://arxiv.org/abs/1602.01380)。

### 3.3 对 R8 的直接含义

R8 已经写过“无低维相关算符污染该项”属于 C4 的一部分，但 R9 的结论更强：

1. 低维污染不是技术小项，而是可以改变小球熵的首阶标度。
2. Zero 的 gap $m\_{\rm stag}=0.23534171$ 与相关维数控制目前没有建立映射。
3. 若 Zero 的局部相关算符有效维数接近二，固定体积平衡可能失去 Jacobson 所需的首阶形式。
4. 所以 Jacobson 2016 仍是最强接口，但不是“只差 boost”一条；还差相关算符层级与固定体积 UV/IR 分离。

这并不推翻 R8，而是把 R8 的 `L1/L2／$\mathcal R$/L4` 缺口补上一条外部审计要求：

$$

\text{在补几何 boost 之前，先证明低维相关算符不会污染固定体积小球首阶熵变。}

$$

---

## §4 各条强路线逐项判断

### 4.1 Jacobson 1995

Jacobson, PRL 75, 1260 (1995), [DOI](https://doi.org/10.1103/PhysRevLett.75.1260), [arXiv:gr-qc/9504004](https://arxiv.org/abs/gr-qc/9504004)。

**外部强度**：在局部 Rindler 视界、Unruh 温度和面积熵输入下，热力学第一定律给出完整场方程。它是最经典、最简洁的路线。

**Zero 接口**：低。它把局部几何、局部 Killing 结构、温度和面积熵全部预先放入。Zero 若从这条路线开始，实际上必须先把所有这些对象构造出来，再回到与 Jacobson 2016 相同的连续极限问题。

**判决**：历史控制与教学基准，不再作为 Zero 的次选目标。

### 4.2 Jacobson 2016

Jacobson, PRL 116, 201101 (2016), [DOI](https://doi.org/10.1103/PhysRevLett.116.201101), [arXiv:1505.04753](https://arxiv.org/abs/1505.04753)。

**外部强度**：直接连接小球纠缠平衡、模 Hamiltonian 和完整非线性场方程。对共形场达到精确论证；对非共形场保留熵变分猜想。

**Zero 接口**：最高。R8.1 补了有限维任意忠实态的第一定律。GNS 模流和球面积候选正好是该证明的第一段。

**剩余硬缺口**：

1. 有限维模流到连续局部代数网。
2. 算子级 $K\_B\to 2\pi B\_B$。
3. 普适面积密度 $S=\eta A+o(A)$。
4. 固定体积平衡的连续实现。
5. 低维相关算符污染的控制。

**判决**：保留第一目标。

### 4.3 Faulkner et al. 2017

Faulkner, Haehl, Hijano, Parrikar, Rabideau, Van Raamsdonk, JHEP 08 (2017) 057, [DOI](https://doi.org/10.1007/JHEP08(2017)057), [arXiv:1705.03026](https://arxiv.org/abs/1705.03026)。

**外部强度**：非常高。它从 CFT 源微扰和所有球区域熵出发，把二阶源扰动后的熵几何化为 asymptotically AdS 度规，并证明该度规二阶满足 Einstein 方程。论文明确不假设 AdS/CFT 对偶，但仍使用 HRRT 面积表示。

**Zero 接口**：低到中。Zero 没有连续 CFT、没有所有球区域的 HRRT 面积泛函，也没有规则要求熵数据来自共形场论。

**判决**：作为“完整非线性恢复的数学上限”阅读，不作为 Zero 下一条主链。

### 4.4 Oh-Park-Sin 2017/2018

Oh, Park, Sin, PRD 98, 026020 (2018), preprint 2017, [DOI](https://doi.org/10.1103/PhysRevD.98.026020), [arXiv:1709.05752](https://arxiv.org/abs/1709.05752)。

**外部强度**：论文证明，在“全阶引力相对熵表达式”成立时，广义纠缠第一定律等价于完整非线性 Einstein 方程。逻辑链本身很强。

**真正代价**：论文自己把“从 CFT 证明全阶 GDERE”留作未完成。因此，拿它作 Zero 目标只是把“怎样从量子信息得到完整 EFE”改名成“怎样从底层得到全阶 GDERE”。

**判决**：可用于检查 R8 的有限维相对熵恒等式，但不能当 Zero 的第一或第二目标。

### 4.5 Bueno-Min-Speranza-Visser 2017

Bueno, Min, Speranza, Visser, PRD 95, 046003 (2017), [DOI](https://doi.org/10.1103/PhysRevD.95.046003), [arXiv:1612.04374](https://arxiv.org/abs/1612.04374)。

**外部强度**：把 Jacobson 方案扩展到高阶引力，并给出因果菱形上的 Wald 熵变分恒等式。

**关键警告**：论文明确说明，对高阶曲率引力，小球线性化方程不能像 Einstein 情形那样推出完整非线性场方程。

**Zero 接口**：需要高阶曲率作用量与 Wald 形式。当前 Zero 主链先给 Einstein 结构，把它接高阶引力会扩大输入面。

**判决**：不是首选，但必须保留为“不要把线性化误称为完整非线性”的边界条件。

### 4.6 Cao-Carroll 2018

Cao, Carroll, PRD 97, 086003 (2018), [DOI](https://doi.org/10.1103/PhysRevD.97.086003), [arXiv:1712.02803](https://arxiv.org/abs/1712.02803)。

**外部强度**：它从抽象 Hilbert 空间分解、互信息图和冗余约束态出发，用 Radon 变换从切割面积构造空间度规，并在弱场极限推出 Einstein 方程。

**论文自己列出的假设**：

1. 首选张量分解。
2. 冗余约束态。
3. 面积与互信息成比例。
4. 修改版纠缠平衡。
5. 存在 emergent EFT 与其熵变分。
6. 存在生成空间几何序列的动力学。
7. Lorentz 不变性在适当极限出现。

**Zero 可补部分**：

| Cao-Carroll 假设 | Zero 候选来源 | 当前缺口 |
|:--|:--|:--|
| 离散因子与信息图 | 零和图、类、计数推前态 | 首选局部分解与非唯一站点识别 |
| 冗余约束/面积型熵 | Z2、G40、G76、G78 | 跨切割面积比例未证 |
| 度规从面积数据重建 | G40、G47、R7 | 尚需 Radon 型反演与规范固定 |
| 修改版纠缠平衡 | R8.1、G79 | 跨全局切割平衡未构造 |
| 弱场 EFE | G55、G56 的场方程结构 | 四维 Lorentzian 组装仍缺 |

**判决**：次选。它不要求先搬入 AdS/CFT，且明确只承诺弱场。未来若成功，得到的是 Zero 原生时空几何与弱场引力，不会抢称完整非线性 GR。

### 4.7 Leichenauer et al. 2018

Leichenauer, Levine, Shahbazi-Moghaddam, PRD 98, 086013 (2018), [DOI](https://doi.org/10.1103/PhysRevD.98.086013)。

**外部强度**：从全息对偶计算 Von Neumann 熵的二阶形状变分，把 QNEC 饱和、量子聚焦与 leading-order 引力方程联系起来。

**Zero 接口**：低。它需要全息对偶、量子聚焦和 Null 变形，当前 Zero 没有这些对象。

**判决**：作为“熵的形状变分也能控制能量”的对照，不是补前提目标。

### 4.8 Dong-Lewkowycz 2018

Dong, Lewkowycz, JHEP 01 (2018) 081, [DOI](https://doi.org/10.1007/JHEP01(2018)081), [arXiv:1705.08453](https://arxiv.org/abs/1705.08453)。

**外部强度**：Euclidean 路径积分和量子极值面给出任意引力理论的熵极值，并推广到 $G\_N$ 全阶。对积分类场方程给出强结果。

**Zero 接口**：低。它依赖 Euclidean 引力路径积分、极值面和广义熵，而 Zero 目前只有离散/有限模型接口。

**判决**：理论形式漂亮，但会扩大 Zero 的输入面。

### 4.9 Alonso-Serrano-Liska 2020

Alonso-Serrano, Liska, PRD 102, 104056 (2020), [DOI](https://doi.org/10.1103/PhysRevD.102.104056), [arXiv:2008.04805](https://arxiv.org/abs/2008.04805)。

**外部强度**：比较因果菱形的 Clausius 熵与纠缠平衡，并指出热力学路线自然给出 unimodular gravity，而不是完整 GR。

**对 Jacobson 路线的意义**：即使完整非线性场方程出现，也可能只固定无迹部分；宇宙学常数和迹部分需要额外条件。

**Zero 接口**：中。它可以提醒 R8 不要把“得到场方程结构”直接写成“唯一得到 GR”。

**判决**：重要批评，不是首选接口。

### 4.10 Gorard 2020 与 Wolfram 2020

Gorard, Complex Systems 29(2), 599-674 (2020), [DOI](https://doi.org/10.25088/ComplexSystems.29.2.599), [arXiv:2004.14810](https://arxiv.org/abs/2004.14810)。  
Wolfram, Complex Systems 29(2), 107-536 (2020), [DOI](https://doi.org/10.25088/ComplexSystems.29.2.107), [arXiv:2004.08210](https://arxiv.org/abs/2004.08210)。

**外部强度**：离散超图、因果图、离散 Lorentz 和离散 Ricci 约束非常完整，是离散 GR 路线的代表。

**关键限制**：Gorard 论文的摘要明确使用“更新规则在极限中保持因果图维数”这一假设，并只给出离散形式的 EFE 约束。它不是从零和原语自动得到连续四维 Lorentzian 流形。

**Zero 接口**：中。R6/R7 可补一部分连续极限，但与 Wolfram 超图重写、固定维数和弱遍历性没有现成等价定理。

**判决**：离散对照，不作为第一或第二接口。

### 4.11 Carrasco et al. 2023

Carrasco, Pedraza, Svesko, Weller-Davies, JHEP 09 (2023) 167, [DOI](https://doi.org/10.1007/JHEP09(2023)167), [arXiv:2306.08503](https://arxiv.org/abs/2306.08503)。

**外部强度**：从“时空复杂度”原理推出 Einstein 及高阶引力、半经典修正，但论文明确大部分证明在二维膨胀引力中完成，其余维度是提案。

**Zero 接口**：低。复杂度体积字典比模 Hamiltonian 更远，当前 Zero 没有复杂度或全息计算结构。

**判决**：记录为更新路线，不入首选。

### 4.12 Kumar 2024/2025

Kumar, Gen. Rel. Grav. (2024/2025), [DOI](https://doi.org/10.1007/s10714-023-03172-x), [arXiv:2404.16912](https://arxiv.org/abs/2404.16912)。

**论文主张**：用广义熵与非平衡热力学“不需要 Jacobson 式假设”而恢复半经典 EFE。

**直接核对后的保留意见**：

1. 起点是大因果菱形和 de Sitter 静态片。
2. 使用近 Rindler 的纠缠第一定律。
3. 输入 $\delta S\_{\rm out}=\delta A/(4G)$。
4. 假设 $T\_{ab}$ 涨落可忽略。
5. 为恢复方程采用了特殊因子约定。

**判决**：它是值得登记的新提案，但现阶段不能据此宣称 Jacobson 路线已被无假设替代。

### 4.13 Bianconi 2025

Bianconi, PRD 111, 066001 (2025), [DOI](https://doi.org/10.1103/PhysRevD.111.066001), [arXiv:2408.14391](https://arxiv.org/abs/2408.14391)。

**外部强度**：把量子相对熵作为度规与物质诱导度规之间的作用量，得到二阶修改重力；低耦合时回到零宇宙学常数的 Einstein 方程，并可用 $G$-场把作用量改写为带小正宇宙学常数的 dressed Einstein-Hilbert 形式。

**Zero 接口**：中。R8.1 已经给出有限维相对熵恒等式，但 Bianconi 路线的核心不是态的模第一定律，而是两个 Lorentzian 度规之间的相对熵作用量。

**缺少的 Zero 对象**：

1. 连续 Lorentzian 度规算符。
2. 物质场诱导度规。
3. 作用量选择原则。
4. $G$-场的微观来源。

**判决**：最新强路线之一，值得后续跟踪，但当前不适合作为 Zero 的直接补前提目标。

---

## §5 Zero 的相对强项与弱点

### 5.1 相对这些论文的强项

1. **原语更少**：外部强论文多数以连续流形、CFT、全息对偶、路径积分或超图重写为起点。Zero 至少把底层约束压到零和和不变量。
2. **账本更细**：某一步是定理、条件、输入还是排除，在项目中分开登记。外部论文常把连续极限、对偶和正则性写在同一套语言里。
3. **有限维第一定律已闭合**：R8.1 对任意有限维忠实态成立，不依赖 Gaussian 或 CFT。
4. **有可证伪残余**：R4 的面积律离散修正和 Cattaneo 型偏离是外部恢复路线普遍缺少的定量预测。
5. **对抗审计习惯**：R2、R3、R5、R6、R7 已把不可行性、路线不等价和连续性边界写清。

### 5.2 相对这些论文的弱点

1. **连续区域代数网未构造**：这是 Jacobson 2016、Faulkner 2017 和 Dong-Lewkowycz 2018 的共同入口。
2. **几何 boost 未证**：有限维模 Hamiltonian 不能自动等于局部 boost。
3. **固定体积平衡未连续化**：目前只有有限模型和局部微扰。
4. **面积密度尚未普适化**：G76/G78 的面积律候选还带 gap 与尺度障碍。
5. **低维相关算符控制缺失**：R9 新增的外部批评直接命中。
6. **Lorentzian 四维组装不完整**：Cao-Carroll 路线也要求这一项。
7. **绝对单位映射仍是 E4**：外部论文中的 $G$、$\Lambda$、Planck 尺度不是 Zero 原生输出。

### 5.3 “拿得出手吗”的当前回答

可以拿出手的是：

1. 一份把外部强路线按前提拆开、并明确 Zero 能补哪一段的审计；
2. R8.1 这类有限维严格引理；
3. R6/R7 条件连续极限；
4. R4 的可证伪预测框架。

还不能拿出手的是：

1. “Zero 已无条件导出四维 GR”；
2. “Zero 已补完 Jacobson 2016”；
3. “R8 与 Oh-Park-Sin 或 Faulkner 路线等价”；
4. “主 Z/G 路线与 D259 图册路线已经汇流”。

---

## §6 下一步工作顺序

### 6.1 首选线：Jacobson 2016

1. 先攻低维相关算符污染：把 $m\_{\rm stag}$、相关指数和固定体积首阶熵变的关系写成可证或可排除的问题。
2. 再攻算子级几何 boost：在有限维/准局部模型中，构造

$$
K_B-2\pi B_B
$$

的收敛估计，而不是只比较谱剖面。

3. 同时补面积密度：给出与 gap、半径和细化尺度无关的首阶系数。
4. 最后才做固定体积连续变分。

### 6.2 次选线：Cao-Carroll 2018

1. 把 Zero 的零和图与互信息图对应起来。
2. 证明或否证冗余约束条件。
3. 从切割熵恢复局部面积数据，再检验 Radon 型反演。
4. 只在弱场范围内声称 Einstein 方程，不提前抬到完整非线性。

### 6.3 不应先做的事

1. 不先搬入 AdS/CFT 或 Ryu-Takayanagi。
2. 不把 Oh-Park-Sin 的全阶 GDERE 当作已给输入。
3. 不用 Gorard/Wolfram 的固定维数连续极限替代 Zero 的选维问题。
4. 不把 Bianconi 的修改重力直接改名成 Zero 的完整 EFE。
5. 不把外部论文的“条件恢复”写成项目自己的“无条件导出”。

---

## §7 可核验来源

1. Jacobson, *Thermodynamics of Spacetime: The Einstein Equation of State*, PRL 75, 1260 (1995), [DOI](https://doi.org/10.1103/PhysRevLett.75.1260), [arXiv:gr-qc/9504004](https://arxiv.org/abs/gr-qc/9504004)。
2. Jacobson, *Entanglement Equilibrium and the Einstein Equation*, PRL 116, 201101 (2016), [DOI](https://doi.org/10.1103/PhysRevLett.116.201101), [arXiv:1505.04753](https://arxiv.org/abs/1505.04753)。
3. Casini, Galante, Myers, *Comments on Jacobson's "Entanglement Equilibrium and the Einstein Equation"*, JHEP 03 (2016) 194, [DOI](https://doi.org/10.1007/JHEP03(2016)194), [arXiv:1601.00528](https://arxiv.org/abs/1601.00528)。
4. Speranza, *Entanglement Entropy of Excited States in Conformal Perturbation Theory and the Einstein Equation*, JHEP 04 (2016) 105, [DOI](https://doi.org/10.1007/JHEP04(2016)105), [arXiv:1602.01380](https://arxiv.org/abs/1602.01380)。
5. Faulkner, Haehl, Hijano, Parrikar, Rabideau, Van Raamsdonk, *Nonlinear Gravity from Entanglement in Conformal Field Theories*, JHEP 08 (2017) 057, [DOI](https://doi.org/10.1007/JHEP08(2017)057), [arXiv:1705.03026](https://arxiv.org/abs/1705.03026)。
6. Bueno, Min, Speranza, Visser, *Entanglement Equilibrium for Higher Order Gravity*, PRD 95, 046003 (2017), [DOI](https://doi.org/10.1103/PhysRevD.95.046003), [arXiv:1612.04374](https://arxiv.org/abs/1612.04374)。
7. Oh, Park, Sin, *Complete Einstein Equations from the Generalized First Law of Entanglement*, PRD 98, 026020 (2018), [DOI](https://doi.org/10.1103/PhysRevD.98.026020), [arXiv:1709.05752](https://arxiv.org/abs/1709.05752)。
8. Cao, Carroll, *Bulk Entanglement Gravity without a Boundary: Towards Finding Einstein's Equation in Hilbert Space*, PRD 97, 086003 (2018), [DOI](https://doi.org/10.1103/PhysRevD.97.086003), [arXiv:1712.02803](https://arxiv.org/abs/1712.02803)。
9. Leichenauer, Levine, Shahbazi-Moghaddam, *Energy Density from Second Shape Variations of the von Neumann Entropy*, PRD 98, 086013 (2018), [DOI](https://doi.org/10.1103/PhysRevD.98.086013)。
10. Dong, Lewkowycz, *Entropy, Extremality, Euclidean Variations, and the Equations of Motion*, JHEP 01 (2018) 081, [DOI](https://doi.org/10.1007/JHEP01(2018)081), [arXiv:1705.08453](https://arxiv.org/abs/1705.08453)。
11. Alonso-Serrano, Liska, *New Perspective on Thermodynamics of Spacetime: The Emergence of Unimodular Gravity and the Equivalence of Entropies*, PRD 102, 104056 (2020), [DOI](https://doi.org/10.1103/PhysRevD.102.104056), [arXiv:2008.04805](https://arxiv.org/abs/2008.04805)。
12. Gorard, *Some Relativistic and Gravitational Properties of the Wolfram Model*, Complex Systems 29(2), 599-674 (2020), [DOI](https://doi.org/10.25088/ComplexSystems.29.2.599), [arXiv:2004.14810](https://arxiv.org/abs/2004.14810)。
13. Wolfram, *A Class of Models with the Potential to Represent Fundamental Physics*, Complex Systems 29(2), 107-536 (2020), [DOI](https://doi.org/10.25088/ComplexSystems.29.2.107), [arXiv:2004.08210](https://arxiv.org/abs/2004.08210)。
14. Carrasco, Pedraza, Svesko, Weller-Davies, *Gravitation from Optimized Computation: Einstein and Beyond*, JHEP 09 (2023) 167, [DOI](https://doi.org/10.1007/JHEP09(2023)167), [arXiv:2306.08503](https://arxiv.org/abs/2306.08503)。
15. Kumar, *Recovering Semiclassical Einstein's Equation Using Generalized Entropy*, Gen. Rel. Grav. (2024/2025), [DOI](https://doi.org/10.1007/s10714-023-03172-x), [arXiv:2404.16912](https://arxiv.org/abs/2404.16912)。
16. Bianconi, *Gravity from Entropy*, PRD 111, 066001 (2025), [DOI](https://doi.org/10.1103/PhysRevD.111.066001), [arXiv:2408.14391](https://arxiv.org/abs/2408.14391)。

**最终判定**：R8 选择 Jacobson 2016 仍是最强接口。当前项目最适合补的一到两篇是 Jacobson 2016 与 Cao-Carroll 2018。相对外部论文，Zero 的强项是底层对象少、账本细、第一定律和可证伪残余清楚；弱项是连续区域代数、几何 boost、普适面积密度、固定体积连续平衡和四维 Lorentzian 组装都还没有关闭。
