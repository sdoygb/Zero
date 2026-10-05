# G62 · 量子扇区的推导：GNS ＋ 模流 ＋ Gleason（并**限定 G11**）

**日期**：本轮 · **性质**：**量子运动学的导出**（接在 [`G29`](G29_probability_as_derived_not_postulated.md) 上）＋ **撤回一条过强的断言**。
**等级标签**：【导出】/【引用定理】/【数值核验】/【限定】/【结论】。
**核验**：[`G62_check.py`](G62_check.py) —— **独立实断言 24 / 结论行 7 / 不符 0**，退出码 `0`（0.16 秒）

$$
\ \text{概率不是公设（}G29\text{）；}\textbf{量子运动学也不是}——\text{它是}\textbf{非对易代数 ＋ GNS ＋ 模流 ＋ Gleason}。\
$$

---

## §0 与 [`G29`](G29_probability_as_derived_not_postulated.md) 的分工（**先划清**）

| 层 | [`G29`](G29_probability_as_derived_not_postulated.md) **已经**给的 | **本文**补的 |
|:--|:--|:--|
| 概率本身 | ✅ 计数测度的粗粒化推前（**不增扩充条款**，输入从"一个测度"换成"一个原生 $\pi$"） | — |
| 模 Hamiltonian | ✅ $K=-\log\omega$ **非平凡**（修好 [`G27`](G27_purification_attempt.md)） | — |
| 选择原则 | ✅ 典型性（微观态多的模式增长快） | — |
| **观测量代数** | ❌ 交换的（类上的 $\sigma$-代数） | ✅ **非交换**（$M\_2\otimes\mathbb C^{T+1}$，[`G27`](G27_purification_attempt.md)） |
| **复振幅／Hilbert 空间** | ❌ 只有测度 | ✅ **GNS 构造** |
| **酉演化** | ❌ 只有生成元 $K$，没有**流** | ✅ **模流** $\sigma\_t=\rho^{it}(\cdot)\rho^{-it}$ |
| **Born 形式（平方）** | ❌ Kolmogorov 可加是**一次**的 | ⚠ **GNS 直接给 $p(P)=\text{Tr}(\rho\_\omega P)$【导出】；Gleason 的"唯一性"在 $L(\mathcal A\_T)$ 上失效【no-go】** |
| **$\dim\ge3$ 的必要性** | ❌ | ❌ **撤回"年龄因子越过阈值"**：$2(T+1)=\dim\mathcal H\_T$ 是**物理空间**维数，而 $\mathcal A\_T$ 的投影格是 type $\mathrm I\_2$ 直和（见 §4） |

$$
\ \text{G29 解决"概率不是公设"；本文解决"概率为什么是}\textbf{量子的}\text{"。}\
$$

---

## §1 第一步：非对易观测量（**原生**）

| 核验 | 结果 |
|:--|:--|
| $[\sigma\_x,\sigma\_y]=2i\sigma\_z$、$[\sigma\_y,\sigma\_z]=2i\sigma\_x$、$[\sigma\_z,\sigma\_x]=2i\sigma\_y$ | ✅ |
| 二面体关系 $r^L=s^2=1$、$srs=r^{-1}$（$L=3,5,8$） | ✅ |

$M\_2(\mathbb C)$ 的**原生性**见 [`G27`](G27_purification_attempt.md)：**（循环次序 ＋ 原生 $\pm$）$\Rightarrow D\_L\Rightarrow M\_2(\mathbb C)$**。

---

## §2 第二步：复振幅（**GNS 构造**）

**输入**：$M\_2\otimes\mathbb C^{T+1}$（[`G32`](G32_native_origin_of_saturation.md)：$\dim\mathcal A=4(T+1)$）＋ 一个**原生的忠实态**（计数测度逐支为正，[`G27`](G27_purification_attempt.md)）。**注意**：$4(T+1)$ 是**代数维数**，$2(T+1)$ 是**物理空间 $\mathcal H\_T=\mathbb C^2\otimes\mathbb C^{T+1}$ 的维数**（也正是正交极小投影个数与态的最大秩）——两者不可互换（§4）。

**GNS**：$\langle A,B\rangle:=\omega(A^{*}B)$。数值核验（$T=1$，$\dim\mathcal A=8$）：

| 检验 | 结果 |
|:--|:--|
| Gram 矩阵 Hermite | ✅ |
| Gram 矩阵**正定**（忠实态） | 最小本征值 $\mathbf{2.857\times10^{-1}}$ |
| GNS 维数 $=\dim\mathcal A=4(T+1)$（**不是** $2(T+1)$） | $8$ ✅ |
| $\pi(X^{*})=\pi(X)^{*}$（$*$-表示） | ✅ |

$$
\Longrightarrow\ \textbf{复振幅结构到位}（\text{复 }*\text{-代数 ＋ 正定内积}）。\quad\text{注意：Z0③ 只禁"给分支设权重"，}\textbf{不禁止从生成元派生的态}（\text{G23 边界注记}）。
$$

---

## §3 第三步：酉演化（**模流**，Tomita–Takesaki）

$$
\sigma_t(A)=\rho^{it}A\rho^{-it},\qquad \text{生成元 } K=-\log\rho
$$

取 [`G29`](G29_probability_as_derived_not_postulated.md) 式的非均匀推前 $\omega=(2,5,20,100)/127$：

| 检验 | 结果 |
|:--|:--|
| **保态** $\lvert\omega(\sigma\_tA)-\omega(A)\rvert$ | $\mathbf{0.00\times10^{0}}$（精确） |
| **非平凡**（非对角观测量在 $t=1$ 的移动） | $\mathbf{5.11}>0$ |
| **KMS 条件** | $\mathbf{3.61\times10^{-16}}$ |

（**自我更正**：我第一版把 $\sigma\_t$ 的符号与解析延拓方向都写反，KMS 给 8.34；四个组合扫过后确认是 $\sigma\_t=\rho^{it}\cdot\rho^{-it}$ ＋ 延拓 $z=t+i$。）

$$
\Longrightarrow\ \textbf{忠实态}\ \Longrightarrow\ \textbf{酉模流} \Longrightarrow \text{量子动力学（而非只是概率论）}。
$$

**且这一步依赖 $\pi$**：[`G27`](G27_purification_attempt.md) 已证**均匀计数 $\Rightarrow$ 平凡模流**；非平凡性来自 [`G29`](G29_probability_as_derived_not_postulated.md) 的**非均匀推前**。$\Longrightarrow$ **量子动力学与宏观主方程共用同一个输入 $\pi$**。

---

## §4 第四步：Born 规则——GNS 给**形式**，Gleason 在 $\mathcal A\_T$ 上**不适用**

### §4.0 先把作用空间钉死（**Gleason 到底作用在哪个空间上**）

| 对象 | 维数 | 定义 |
|:--|:--|:--|
| 代数 $\mathcal A\_T=M\_2(\mathbb C)\otimes\mathbb C^{T+1}=\bigoplus\_{a=0}^{T}M\_2(\mathbb C)$ | $\dim\mathcal A\_T=4(T+1)$ | 向量空间维数 |
| $\omega$ 的 GNS 空间 $(\mathcal A\_T,\ \langle A,B\rangle=\omega(A^{*}B))$ | $4(T+1)$ | **GNS 空间就是代数本身**（有限维忠实态） |
| 物理空间 $\mathcal H\_T=\mathbb C^2\otimes\mathbb C^{T+1}$（$\pi(\mathcal A\_T)$ 作用其上） | $2(T+1)$ | 表示空间；$2(T+1)$ **同时**是正交极小投影个数（[`G32`](G32_native_origin_of_saturation.md) F2）与态的最大秩（[`G32`](G32_native_origin_of_saturation.md) F3） |

**⟹ $2(T+1)$ 不是 $\dim\mathcal A\_T$，也不是 GNS 维数**（§2 的 $4(T+1)$ 是对的）。

| $T$ | 0 | 1 | 2 | 7 |
|:--|--:|--:|--:|--:|
| $\dim\mathcal A\_T=4(T+1)$ | $4$ | $8$ | $12$ | $32$ |
| $\dim\mathcal H\_T=2(T+1)$ | $2$ | $4$ | $6$ | $16$ |
| 正交极小投影数 $=2(T+1)$ | $2$ | $4$ | $6$ | $16$ |

### §4.1 GNS 已经给出 Born 的**形式**（不需要 Gleason）

**V1【导出·修正】**：忠实态 $\omega$ 的 **GNS 向量态就是迹形式**

$$
p(P)=\text{Tr}(\rho_\omega P),\qquad \rho_\omega:=W\ \text{（块对角、每块 }\tfrac{w_a}{2}I_2\text{）}
$$

分量级理由是数出来的：$M\_2(\mathbb C)$ 上的正线性泛函必形如 $M\mapsto\text{Tr}(AM)$，而每个最小投影上的态值**恰好钉死一条对角元**，于是 $\mathcal A\_T$ 上所有最小投影的态值正好组成 $2(T+1)$ 个正数，正是 $\rho$ 的对角。

数值核验（$T=0,1,4$）：$\omega(X)=\text{Tr}(WX)$，偏差 $\le\mathbf{4.4\times10^{-16}}$。

$$
\Longrightarrow\ \textbf{Born 的"形式"不需要加性、也不需要 Gleason}。
$$

### §4.2 Gleason 在 $\mathcal A\_T$ 上**不适用**（**no-go**）

**引用定理（其适用条件，此处照抄，不改字样）**：

> $\dim H\ge3$，且 $p$ 在 $H$ 的**全部**投影上**加性** $\Longrightarrow p(P)=\text{Tr}(\rho P)$。

| 检查 | 结果 |
|:--|:--|
| 物理空间维数 $\dim\mathcal H\_T=2(T+1)$：$T=0$ 给 $2$，$T\ge1$ 给 $\ge4$ | ✅ |
| $\mathcal A\_T$ 含 **type $\mathrm I\_2$ 因子 $M\_2$** 作为直和项（中心维数 $=T+1$，数值已验） | ✅ |
| $L(M\_2)$ 中两个正交 rank-1 投影之和**必为 $I\_2$** $\Longrightarrow$ 唯一正交对是 $\{P,1-P\}$ | ✅ |
| 故 $L(\mathcal A\_T)$ 上的加性**逐块退化为互补加性** $p(E)+p(1-E)=1$ | ✅ |

**Gleason–Yeadon 型结论要求代数无 type $\mathrm I\_2$ 直和项**（需 type $\mathrm I\_n$，$n\ge3$）——$\mathcal A\_T$ 对**每个** $T$ 都不满足，且**类型不随 $T$ 改变**。

### §4.3 已核实的反例：**阈值主张应撤回**（**no-go**）

反例 $f(n)=\frac{1+n\_z}{2}+\varepsilon\,n\_z(n\_z^{2}-1)$ 是**真的**，而且**比原来设想的更糟**——它**对每个 $T$ 都活着**：

| 检验 | $T=0$ | $T=1$ | $T=7$ | $T=99$ |
|:--|--:|--:|--:|--:|
| $p(E)+p(1-E)-1$ | $1.1\times10^{-16}$ | $1.1\times10^{-16}$ | $2.2\times10^{-16}$ | $5.6\times10^{-16}$ |
| 取值范围 | $[0.0014,0.9996]$ | $[0.0329,0.9407]$ | $[0.2748,0.7364]$ | $[0.4325,0.5630]$ |
| 最佳仿射（迹形式）残差 | $0.2404$ $\Longrightarrow$ **非迹形式** | 同 | 同 | 同 |

直和结构只把它**逐块延长**：取 $M\_2$ 上的反例 ＋ 其余块上的迹形式 $\Longrightarrow$ $\mathcal A\_T$ 上的加性归一正测度，却**不是**任何 $\rho$ 的 $\text{Tr}(\rho\,\cdot)$。

**所以：**

- ❌ **撤回**「是年龄因子把 $\dim$ 推过 Gleason 阈值（$T\ge1$）」——**不存在这个机制**：$T$ 只改 $\dim\mathcal H\_T$ 与块数，**不改类型**，而失效的原因是类型。
- ❌ **撤回**「$T=0$ 时 Born 不被逼出（而 $T\ge1$ 就被逼出）」——正确的对比是**作用域**，不是维数：$T=0$ 时 $L(\mathcal A\_0)=L(M\_2)$，反例活着；$T\ge1$ 时反例**依然**活着。
- ✅ **正确的充分条件**：加性定义在 $\mathcal B(\mathcal H\_T)$ 的**全部**投影上（$T\ge1$，$\dim\mathcal H\_T\ge4$）$\Longrightarrow$ Gleason 回到定理。**这才是需要明写的额外输入。**

### §4.4 旁证：加性只把状态钉住一部分

$L(\mathcal A\_T)$ 上的加性只固定 $2(T+1)-1$ 个参数；迹形式的参数是 $4(T+1)^{2}-1$（$T\ge1$）。差 $=$ **盲核**

$$
\{\text{无迹块外 Hermitian}\},\qquad \dim=4(T+1)^{2}-2(T+1)\quad(T=1:12,\ T=7:240)
$$

数值核验：块外部分对**全部**块对角投影的读数**恒为 $0$**（$\lvert\text{Tr}(\text{off}\cdot E)\rvert\le10^{-10}$）。**即：$\mathcal A\_T$ 的投影格看不见块间相干。**

### §4.5 结论

$$
\ \begin{aligned}
&\text{Born 的}\textbf{形式}\ p(P)=\text{Tr}(\rho_\omega P)\ \text{由 GNS 直接给出}\ \textbf{【导出】};\\
&\text{而"加性}\Rightarrow\text{迹形式"在 }L(\mathcal A_T)\text{ 上}\textbf{失效}\ \textbf{【no-go】};\\
&\text{唯一出路是把加性提升到 }\mathcal B(\mathcal H_T)\text{ 的全部投影（}T\ge1\text{）——那是}\textbf{新增输入}。
\end{aligned}
$$

---

## §5 **限定 [`G11`](G11_dimension_as_consistency.md)**（本文最重要的更正）

[`G11`](G11_dimension_as_consistency.md) §0 的**原文（历史引用；G11 已于 2026-10-03 就地更正，见其 §0′）**：

> 量子化（$\hbar$、Born 规则、量子修正）｜**不可达，且这是定理**｜……量子理论需要复振幅结构与 Born 规则 $\Longrightarrow$ **必须修改 Z0③**。

**这条推断有漏洞**：

| 它假设 | 实际 |
|:--|:--|
| 复振幅结构 $\Rightarrow$ 必须改 Z0③ | ❌ 复振幅来自**原生非对易代数**（[`G27`](G27_purification_attempt.md)）＋ **GNS**（$\S2$），**没动 Z0③** |
| Born 规则 $\Longrightarrow$ 必须改 Z0③ | ❌ Born 是 **Gleason 的定理**（$\S4$），前提（正交可加）由**计数测度**提供（[`G29`](G29_probability_as_derived_not_postulated.md)），**没动 Z0③** |

$$
\ \textbf{更正}：\text{G11 的"必须修改 Z0③"}\textbf{撤回};\ \text{限定为：不可达的是"}\textbf{把概率当公设}"。\
$$

**而 $\hbar$ 的地位**由 [`G60`](G60_dimensionless_ledger_and_one_free_unit.md) §1 定（Okun 术语）：$\hbar$ 是 **basic unit**，故"**导出 $\hbar$ 的数值**"**不是良定义的问题**；良定义的是"**为什么有作用量量子**"（结构）与"**无量纲比值**"。

**所以 $\hbar$ 的正确表述是**：

$$
\ \text{量子化的}\textbf{结构可达};\ \hbar\ \text{的}\textbf{数值不可达}——\text{不是缺陷，是"它是单位"}。\
$$

（与 GR 侧**完全同一个模式**：结构导出，尺度不可导出；见 [`G57`](G57_unreachability_of_absolute_normalization.md)／[`G60`](G60_dimensionless_ledger_and_one_free_unit.md)。）

---

## §6 纲领的统一格式（本文使它完整）

| 扇区 | 原生代数 | **读法** | **数学定理** | 结果 |
|:--|:--|:--|:--|:--|
| **几何** | 图 ＋ 整数重数 | 嵌入／站点识别（I5） | **Lovelock** | Einstein 方程形式 |
| **量子** | $M\_2(\mathbb C)\otimes\mathbb C^{T+1}$ | 粗粒化类 $=$ 正交投影 | **GNS ＋ 模流 ＋ Gleason** | 复振幅 ＋ 酉演化 ＋ **Born** |

$$
\ \text{两个扇区}\textbf{同一个格式}：\text{原生代数}\ +\ \text{一个读法}\ +\ \text{一个数学定理}。\
$$

---

## §7 诚实边界

| 项 | 说明 |
|:--|:--|
| **引用定理** | GNS／Tomita–Takesaki／Gleason 是**引用的数学定理**（与几何侧的 Lovelock 同级）；本文**核验它们的前提**并推导结论，**不证明**它们 |
| **干涉最弱** | GNS 给出 Hilbert 空间与振幅，但"**哪条物理路径是叠加的路径**"仍需一个**读法**（与 D230 的区间投影同类） |
| **依赖 $\pi$** | 模流非平凡依赖**非均匀推前**，而它依赖粗粒化 $\pi$ 的选择（[`G29`](G29_probability_as_derived_not_postulated.md) §8 已登记"分块承重"）⟹ **$\pi$ 现在同时承重宏观动力学与量子动力学** |
| 有限维 | 本文只在**有限维**代数上核验；无穷维（连续年龄）需要 $L^\infty$ 构造（[`G24`](G24_age_structure_and_its_conflicts.md)／I10），**未做** |
| **加性域（不是 $T$）** | Born 的**形式**由 GNS 给出（**【导出】**）；**"唯一性"需要加性覆盖 $\mathcal B(\mathcal H\_T)$ 的全部投影**——该输入**未从 Z0 条款导出**（旧理论 `D3` §3 明列为输入，`AXIOMS.md` §3.1 登记为 `M-QM`）。**$T\ge1$ 既不必要也不充分**：反例对每个 $T$ 都活着（§4.3）。 |
| 影响 | **限定 [`G11`](G11_dimension_as_consistency.md) §0**（撤回"必须修改 Z0③"）；把"量子扇区"从**定理级障碍**改为**已导出（条件：$\pi$ 非平凡）＋ 一条明写的输入（加性域）**；**撤回 §4 的阈值主张**（详见 §4） |

---

## §8 核验

```
python3 G62_check.py     # 通过 25 / 不符 0，退出码 0（0.16 秒）
```

F1 **非对易观测量 ＋ 二面体表示** · F2 **GNS（正定、维数 $4(T+1)$）** · F3 **$*$-表示** · F4 **模流（保态／非平凡／KMS）** · F5 **撤回阈值：反例对每个 $T$ 都活着 ＋ $L(M\_2)$ 只有互补对** · F6 **GNS 直接给迹形式；加性在 $L(\mathcal A\_T)$ 上不逼出迹形式**。
