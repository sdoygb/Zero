# ν_R 可测性登记：为什么现在测不出来

**目的**：把「Zero 里的右手中微子为什么还没被测出来」这件事登记清楚——**不是"测不出来"，而是可测性由两个数决定，这两个数目前谁也不知道，而且被 seesaw 锁在一起**。
**问题来源**：用户提问「现在为什么还不能测出来呢」（2026-10-15，Z136 落盘之后）。
**关联**：[Z136](../Z136_右手中微子的强制性_十六槽载体与锚定完备化.md)（ν_R 的存在由结构逼出：色单态、$T_3=0$、$Y=0$、$Q=0$，sterile 为推论）、[`NUR_VERIFICATION.md`](NUR_VERIFICATION.md)（ν_R 独立核查）、[Z098](../Z098_右手中微子深挖_sterileMajorana与跷跷板结构.md)／[Z112](../Z112_右手中微子详细特征_三代质量谱与跷跷板MR与sterile双角色.md)（Zero 侧边界：无 Higgs／汤川／电弱破缺 ⟹ $y_\nu,v,M_R$ 外部输入）。
**日期**：2026-10-15 ｜ **性质**：登记／讨论（不改任何 Z 篇；外部资料按下面的档位标注）
**档位词**（每处数字/结论都标）：【核实】＝本次检索读到原文摘要或其片段；【量级】＝标准公式算出的量级估计；【未核】＝未逐条核实，只作指向。

---

## §0 结论（先行）

1. ν_R 是**完全 sterile**：色单态、$T_3=0$、$Y=0$、$Q=0$ ⟹ 三个规范耦合全零。**这一点 Zero 是导出的**（Z136 定理 2），主流文献的定义与之逐字对应：「as well as being colour and electroweak singlets, they also have zero hypercharge and electric charge」【核实，[King 综述](https://arxiv.org/abs/2502.07877)】。
2. 因此它只有两条可见通道：**(i) 与活性中微子混合**（需汤川＋Higgs）；**(ii) 引力**（量级上无望）。**所有实验实质上都在测同一个量：混合角 $\theta^2$**。
3. 而 seesaw 把 $\theta^2$ 与质量 $M_R$ **锁在一起**：$\theta^2\simeq m_\nu/M_R$【量级，标准单代 seesaw】。⟹ **$M_R$ 越大越"自然"，耦合越小、越不可测**；要可测就得 $M_R$ 小，那又是另一类模型。
4. **Zero 的贡献是"存在性 + 量子数 + 计数"，不是"可测性"**：它给每代 1 个、共 3 个，恰好是 seesaw／SO(10) 16 维旋量需要的最小内容；但 Zero **算不出 $\theta^2$，也算不出 $M_R$**（Z098 边界；Z112 的 $M_R$ 是用 Type-I 公式反推的外部输入）。
5. **当前的零结果并不排除 ν_R**，只排除参数空间的一部分。Zero 的 ν_R 没有被放到任何具体参数点上，所以"没测到"与"Zero 有 ν_R"**不矛盾**。
6. 要让 Zero 变可测，必须补一个 **$M_R$ 的生成机制**——这是下一个真课题（见 §5）。

---

## §1 可测性的结构

**Type-I seesaw（单代）**：$m_\nu\simeq m_D^2/M_R$，重的那个本征值 $\simeq M_R$，轻重混合 $\theta\simeq m_D/M_R$，于是

$$\boxed{\ \theta^2\ \simeq\ \frac{m_\nu}{M_R}\ }\qquad\text{【量级，本文只用它做量级估算】}$$

取 $m_\nu\approx0.05$ eV【量级】：

| $M_R$ | $\theta^2\simeq m_\nu/M_R$ | 可及性 |
|--:|--:|:--|
| $10^{14}$ GeV | $\sim5\times10^{-25}$ | **无任何手段**（这是"自然 seesaw"标度） |
| $1$ GeV | $\sim5\times10^{-11}$ | 仍在 LHC／束流收集器**直接搜寻之外**（那些只够到 $\theta^2\sim10^{-6}$–$10^{-9}$【量级】） |
| $1$ keV | $\sim5\times10^{-5}$ | 才落进 X 射线／结构形成的可及范围（"keV sterile 暗物质"那一支） |

**要点**：可观测性**不是 ν_R 的属性，是模型的属性**——同一个 ν_R，在 sequential dominance（两个右手中微子的最小情形）下预言的 HNL 是「practically unobservable」，而在 inverse seesaw 里「allows observable heavy neutral leptons」【核实，同上综述摘要】。

---

## §2 实验格局（2025–2026）与各自的窗口

| 窗口 | 实验手段 | 测的量 | 状态（含引用） |
|:--|:--|:--|:--|
| **eV 级**（振荡异常那一支） | 氚 β 谱 kink（KATRIN）、反应堆／镓源、加速器消失实验 | $|\Delta m^2|\sim0.1$–$10^2$ eV²、$\theta^2\gtrsim10^{-2}$ | **KATRIN 259 天**：3600 万个氚 β 电子、端点下 40 eV 谱形，探测 $|\Delta m^2|$ 从零点几 eV² 到几百 eV²，**排除混合角大于百分之几**的轻 sterile 态，排除了 **gallium anomaly** 相当一部分空间、并**挑战 Neutrino-4**【核实，[OSTI 摘要](https://www.osti.gov/biblio/3015522)／Nature **648**, 8092 (2025)】；MicroBooNE 双束流搜寻（[Nature 2025](https://www.nature.com/articles/s41586-025-09757-7)）【未核具体数字】 |
| **keV 级**（暖暗物质） | X 射线未认证谱线（含争了十年的 3.5 keV 线）、结构形成（Lyman-α）、TRISTAN | keV 质量窗的 $\theta^2$ 上限 | **XRISM**：10 个星系团堆叠谱，在 **2.5–15 keV** 以 **2.5 eV** 步长搜未认证谱线（含 3.5 keV 线）——**给出的是上限**（标题即 "Constraints"）【核实片段，[ApJL 2025](https://iopscience.iop.org/article/10.3847/2041-8213/ae17ad)；**上限数值未核**】；KATRIN 的 TRISTAN 计划扫全 18.6 keV 谱【未核】 |
| **GeV–TeV 级**（唯一能"直接产生"的窗口） | LHC 位移顶点／长寿命粒子（ATLAS/CMS）、束流收集器、B 工厂 | $\theta^2$（产生率 ∝ $\theta^2$） | 只能覆盖**偏离 seesaw 线**的那部分参数区【核实（结论方向），[King 综述](https://arxiv.org/abs/2502.07877)】；权威汇总见 [PDG 2025 HNL 清单](https://pdg.lbl.gov/2025/listings/rpp2025-list-heavy-neutral-leptons.pdf)、[CMS 13 TeV 汇总图](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SummaryPlotsEXO13TeV)、[PRD 111, 075031](https://journals.aps.org/prd/pdf/10.1103/PhysRevD.111.075031)【三者均未逐条读】 |
| **GUT 标度** $10^{13}$–$10^{15}$ GeV | 无直接手段 | 只能间接：轻子味破坏、$0\nu\beta\beta$、宇宙学 | leptogenesis 需要 $M_R\gtrsim1$ GeV【核实，同上综述】⟹ **"想解释物质-反物质不对称"这个动机本身就把 $M_R$ 推到很难测的地方** |
| **宇宙学** | CMB／BBN 的 $N_{\rm eff}$、结构形成 | "曾经相对论性"的额外自由度 | $N_{\rm eff}$ 的余量很小（Planck 值 $\approx2.99$、误差约 $0.17$【**未核**，只作量级；本次检索未读到原文数值】） |

**一条旁证（只登记标题，未细读）**：ORNL 的新闻稿标题即「New results from physics experiment at ORNL show no sign of sterile neutrinos」【核实标题，[ORNL](https://www.ornl.gov/news/new-results-physics-experiment-ornl-show-no-sign-sterile-neutrinos/)】。

---

## §3 Zero 侧能给什么、不能给什么

| 项 | Zero 状态 | 指针 |
|:--|:--|:--|
| ν_R 的**存在** | **导出**（16 槽载体上 8 条 ℤ₂ 轨道中唯一半满者的补位） | Z136 §4 定理 2；`ZeroLean/Zero/Z136_SixteenSlot.lean` `z136_unique_completion` |
| 量子数：色单态、$T_3=0$、$Y=0$、$Q=0$ | **导出**（由槽位坐标算出） | Z136 式 (4.2)(4.3) |
| **sterile**（三规范全零 ⟹ 只能靠混合可见） | **推论** | Z136 式 (4.3) |
| Majorana 可行 | **条件**（反号不变对任何中性色单态都成立；真条件是实表示＋轻子数破缺） | [`NUR_VERIFICATION.md`](NUR_VERIFICATION.md) §1 第 10 行 |
| 每代 1 个、共 **3 个** | **导出**（×3 代来自 Z089 三频率槽） | Z136 §4、§附录 A |
| 与 SO(10) 16 维旋量的一致性 | **引用性交叉核对**（不作推导） | Z136 §5；[King 综述](https://arxiv.org/abs/2502.07877) 同结构 |
| **$\theta^2$（混合角）** | **给不出**——无 Higgs／汤川／电弱破缺 | Z098 §边界 |
| **$M_R$（质量标度）** | **给不出**——Z112 用 Type-I 公式**反推**，明写外部输入 | Z112 §边界 |
| 绝对质量／跷跷板数值 | 外部输入（$y_\nu,\ v,\ M_R$） | Z112:49 |

**⟹ 一句话**：Zero 现在**既不能说 ν_R 可测、也不能说不可测**；它给出的是「可测性必须经过的那条通道」（混合）与一组结构性预言（sterile、Majorana 可行、3 代、每代 1 个）。

---

## §4 这张登记能立即支持的两件事

1. **回答"为什么还不能测出来"的完整形式**：不是"它太弱所以测不到"这一句话，而是三段式——
   ① sterile ⟹ 所有实验都在测 $\theta^2$；② seesaw 把 $\theta^2$ 与 $M_R$ 反锁（$\theta^2\simeq m_\nu/M_R$）⟹ "自然的" $M_R$ 不可达；③ 于是零结果**只覆盖参数空间的一部分**，与 ν_R 存在不矛盾。
2. **给 Zero 侧划出一条真正的研究课题**（见 §5）。

---

## §5 下一步（两条，可选）

1. **给 $M_R$ 找一个 Zero 内生的标度**：Zero 侧现成的标度只有内禀频率 $\log2$（Z087/Z131）与轻子质量线（Z121–Z123）。把 $M_R$ 挂到哪个内禀量上（例如"$M_R$ 与 $\log2$ 的某个整幂关系"），是让 ν_R 从"存在"走向"可瞄准"的唯一路径。**注意**：这需要新增一笔输入或新机制，按库内纪律应登记价签。
2. **画一张 $\theta^2$–$M_R$ 图**：把 §2 的四个窗口边界与 §1 的 seesaw 线画在同一张图上，标出「Zero 目前能落点的区域」（即整条 seesaw 线，因为 $M_R$ 未定）。图比文字更能回答"为什么还不能测"。

---

## §6 诚实边界（本登记自己的）

1. **§0–§2 讲的是主流物理对 RHN 可测性的结构**，有文献支撑（标【核实】处）；**§3–§5 是 Zero 侧的对接与推断**，其中 $M_R$ 机制仍是**开放项**。
2. **数字分档**：凡标【量级】的是标准公式/常识算出来的量级（$\theta^2\simeq m_\nu/M_R$ 与三档数值），**不是**从检索结果里抄的；凡标【未核】的只作指向，未逐条读原文（含 Planck 的 $N_{\rm eff}$ 数值、XRISM 的上限数值、LHC 各通道的具体灵敏度）。
3. **外部链接是外部内容**（未受信任），本登记只把它们当资料出处引用。
4. 本登记**未改动任何 Z 篇**，也未改 `_staging` 里的其它文件；新增本文件。

---

## 附录 来源清单

| 来源 | 用途 | 档位 |
|:--|:--|:--|
| [King, *Right-handed neutrinos: seesaw models and signatures*, arXiv:2502.07877v3](https://arxiv.org/abs/2502.07877) | sterile 定义；seesaw 参数化；"practically unobservable" vs "observable HNLs"；leptogenesis $M_R\gtrsim1$ GeV；SO(10) 三个 16 | 【核实】 |
| [KATRIN, *Sterile-neutrino search based on 259 days of KATRIN data*, Nature 648, 8092 (2025)](https://www.osti.gov/biblio/3015522) | 259 天／3600 万电子／端点下 40 eV；$|\Delta m^2|$ 覆盖零点几–几百 eV²；排除混合角 > 几 %；排除 gallium anomaly 相当部分、挑战 Neutrino-4 | 【核实】 |
| [MicroBooNE, *Search for light sterile neutrinos with two neutrino beams*, Nature (2025)](https://www.nature.com/articles/s41586-025-09757-7) | eV 级双束流搜寻 | 【未核】 |
| [XRISM, *Constraints on Unidentified X-Ray Emission Lines, Including the 3.5 keV Line…*, ApJL (2025)](https://iopscience.iop.org/article/10.3847/2041-8213/ae17ad) | keV 段 X 射线上限（2.5–15 keV、2.5 eV 步长、10 星系团堆叠） | 【核实片段】 |
| [PDG 2025 · Heavy Neutral Leptons](https://pdg.lbl.gov/2025/listings/rpp2025-list-heavy-neutral-leptons.pdf) | HNL 限值权威汇总 | 【未核】 |
| [CMS Exotica 13 TeV summary plots](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SummaryPlotsEXO13TeV) | LHC 直接搜寻通道 | 【未核】 |
| [*Heavy neutral leptons without prejudice*, PRD 111, 075031](https://journals.aps.org/prd/pdf/10.1103/PhysRevD.111.075031) | LHC 限值的无偏重解释 | 【未核】 |
| [ORNL 新闻稿](https://www.ornl.gov/news/new-results-physics-experiment-ornl-show-no-sign-sterile-neutrinos/) | "no sign of sterile neutrinos" 的气候（只读标题） | 【核实标题】 |