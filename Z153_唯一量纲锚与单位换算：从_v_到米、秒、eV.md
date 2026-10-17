# Z153 · 唯一量纲锚与单位换算：从 $v$ 到米、秒、eV

> 篇号：Z153
> 状态：正式篇（推导篇）
> 前置：Z001（量纲空洞）、Z002 §7.2–§7.3（量纲常数是单位）、Z013 定理 16（$v$ 是唯一有量纲自由度）、Z016（分层：$G/\Lambda/\hbar$ 在读出面）、Z023（锚定与可计算性）、Z037 §2.3（质量锚与无量纲锚）
> 关联脚本：`ZeroGPU/src/unit_conversion_check.py`、`ZeroGPU/src/unit_conversion_exchange_check.py`
> 数据：`ZeroGPU/results/unit_conversion_check.json`、`ZeroGPU/results/unit_conversion_exchange_check.json`
> 登记日期：2026-10-17

## 摘要

Zero 的全部有量纲量由**一个锚**定住：**电弱真空期望值 $v=246.22$ GeV**。取它为锚，三个基本量纲（长度、时间、质量）归约为一个，换算链是

$$
v\ \xrightarrow{\ \ln(M_X/v)=33.6197\ }\ M_X=9.821\times10^{16}\ \text{GeV}\ \xrightarrow{\ r_{\rm int}=\hbar c/M_X\ }\ r_{\rm int}=2.009\times10^{-33}\ \text{m}.
\tag{0.1}
$$

**三条正面结果**（独立复算，`unit_conversion_check.py`，9/9 PASS）：

1. **长度单位**：$1$ 跳数 $=2.009\times10^{-33}$ m $=124.3\,l_{\rm Pl}$；
2. **时间／质量单位**：$1$ 代 $=6.702\times10^{-42}$ s（$=\alpha/c_{\rm SI}$），$1$ 重数 $=2.706\times10^{-6}$ kg（$=\alpha c_{\rm SI}^2/G_{\rm SI}$）；
3. **能量单位**：$v=246.22$ GeV $=2.462\times10^{11}$ eV；质量按 $E=mc^2$ 换到 eV。

$$
\boxed{\;\text{Zero 的米、秒、kg、eV 由同一个锚 }v\text{ 定出：}\ (\text{长度}\to\text{时间}\to\text{质量})\ \text{经 } c_{\rm SI},G_{\rm SI}\ \text{串联}.\;}
\tag{0.2}
$$

**边界（如实）**：这是**单位换算**，不是**预言**。换算链的正确账是

$$
\boxed{\;\text{1 个量纲锚}\ (v)\ +\ \text{1 个上游开放数}\ \bigl(E_6\ \text{破缺环④}\bigr)\ +\ c_{\rm disc}\ +\ \text{舞台量}\ Q\ \bigl(G_{\rm disc}=1/Q\bigr).\;}
\tag{0.3}
$$

其中 $c_{\rm disc}$ 是方程系数（可任取），$G_{\rm disc}=1/Q$ 是舞台量；质量单位 $\gamma$ 的定义里含 SI 的 $G$，故它是"用 SI 表达的模型单位"。$c$、$G$ 的数值不是本框架的预言。

## §1 为什么三个基本量纲归约成一个

经典物理里长度、时间、质量是三个独立量纲。Zero 的底层只有纯数（顶点数、代数、重数），量纲来自**读出面**的选择。分层的读出面给出（Z016 表）：

| 对象 | L0 | L1 | L1′ | L2 |
|:--|:--:|:--:|:--:|:--:|
| $G,\Lambda$ 数值、$\hbar$ | 无 | 无 | — | **约定／单位（读出面）** |

$$
\Longrightarrow\quad c,\ G,\ \hbar\ \text{是单位（读出面）},\ \text{不是内生的物理量}.
\tag{1.1}
$$

于是量纲之间不独立：

$$
[T]=\frac{[L]}{c},\qquad [M]=[L]\,\frac{c^2}{G}.
\tag{1.2}
$$

$$
\boxed{\;\text{三个基本量纲}\ \longrightarrow\ \text{一个};\ \text{只需锚定一个}. \;}
\tag{1.3}
$$

**与 Z002 §7.2–§7.3 的关系**：量纲常数是**单位**不是**参数**（Duff–Okun–Veneziano），故"$G,\Lambda,\hbar$ 的数值"不是良定义的问题；几何扇区的有量纲自由度恰 $1$ 个，与 (1.3) 同一条。

## §2 锚点：$v=246.22$ GeV

$$
\boxed{\;\text{锚}:\quad v=246.22\ \text{GeV}\quad(\text{电弱真空期望值}).\;}
\tag{2.1}
$$

**为什么选 $v$**：它是一个已测量的有量纲物理量；直接是能量单位，换算最顺；与 Z013 定理 16 的"$v$ 是唯一有量纲自由度、$v$ 是量纲单位"同一记号。

**不能选 $c$**：$c$ 是单位（读出面），用它定标度是"用尺子的刻度定房间大小"。

$$
\Longrightarrow\quad \text{锚是一个}\textbf{约定}:\ \text{换锚（如 }m_e\text{）只改变单位的读数，不改变无量纲比值}.
\tag{2.2}
$$

## §3 锚点链与换算表

### 3.1 锚点链

$$
v\ \xrightarrow{\ \ln(M_X/v)=33.6197\ }\ M_X=v\,e^{33.6197}=9.821\times10^{16}\ \text{GeV}
\tag{3.1}
$$

$$
M_X\ \xrightarrow{\ r_{\rm int}=\hbar c/M_X\ }\ r_{\rm int}=2.009\times10^{-33}\ \text{m},\qquad \frac{r_{\rm int}}{l_{\rm Pl}}=124.3 .
\tag{3.2}
$$

### 3.2 换算表

| Zero 单位 | 换算 | 式子 |
|:--|:--|:--|
| 1 跳数（长度） | $2.009\times10^{-33}$ m | $\alpha=r_{\rm int}$ |
| 1 代（时间） | $6.702\times10^{-42}$ s | $\beta=\alpha/c_{\rm SI}$ |
| 1 重数（质量） | $2.706\times10^{-6}$ kg | $\gamma=\alpha\,c_{\rm SI}^2/G_{\rm SI}$ |
| 能量 | $v=246.22$ GeV $=2.462\times10^{11}$ eV | 锚本身 |

$$
\boxed{\;\alpha=2.009\times10^{-33}\ \text{m},\quad \beta=6.702\times10^{-42}\ \text{s},\quad \gamma=2.706\times10^{-6}\ \text{kg},\quad v=2.462\times10^{11}\ \text{eV}.\;}
\tag{3.3}
$$

### 3.3 物理读法

- $1$ 跳数 $=124.3$ 个普朗克长度：图的最小尺度，比普朗克长度大两个量级；
- $1$ 代 $=6.702\times10^{-42}$ s：与普朗克时间同量级；
- $1$ 重数 $=2.706\times10^{-6}$ kg：约 $2.7$ 毫克；
- 质子质量 $1.673\times10^{-27}$ kg $=6.18\times10^{-22}$ 重数。

## §4 能量：$v$ 到 eV

锚 $v$ 直接给能量单位：

$$
v=246.22\ \text{GeV}=2.4622\times10^{11}\ \text{eV}.
\tag{4.1}
$$

质量按 $E=mc^2$ 换到 eV；不需要第二个有量纲输入。eV 与 kg 的桥是 $1\ \text{eV}=1.602176634\times10^{-19}$ J（SI 定义），与 $c$ 同属单位约定。

## §5 账：1 个量纲锚 ＋ 上游四环链 ＋ 舞台量 $Q$

**换算不是预言**。把换算链逐环拆开，正确的账是：

| 类        | 内容                      | 状态                                                                              |
| :------- | :---------------------- | :------------------------------------------------------------------------------ |
| **量纲锚**  | $v=246.22$ GeV          | 唯一的外部有量纲量                                                                       |
| **无量纲数** | $33.6197=\ln(M_X/v)$    | 四环链（§5.6）：① ② ③ 在手，④（$E_6$ 破缺）为上游开放项                                                                      |
|          | $c_{\rm disc}$          | **方程系数**：色散对任意 $c$ 都返回同一个 $c$（$c=0.5/0.8944/1.2$ 都给 $\omega/(c\sqrt\lambda)=1$） |
| **舞台量**  | $Q$（$G_{\rm disc}=1/Q$） | cubic $Q=1$；diamond $Q=3.46$；由舞台几何定                                             |

$$
\boxed{\;\text{有量纲输出}=\text{模型纯数}\times\text{锚};\qquad \text{锚的数量}=\text{独立量纲数}=1;\qquad G_{\rm disc}=\frac1Q\ \text{由交换定出}.\;}
\tag{5.1}
$$

**三处必须一起读的边界**：

1. **$c_{\rm disc}$ 不是测量值**：色散关系 $\omega=c\sqrt\lambda$ 对任意 $c$ 成立，拟合返回的正是放进去的 $c$；
2. **$G_{\rm disc}=1/Q$ 是交换给出的**：Z144 定理 144.2 的 $K=A$ 唯一性把泊松核系数固定，$Q$ 是舞台量（cubic 与 diamond 差 $0.289$ 倍）；
3. **质量单位由交换给出**：$\gamma=\alpha c^2/G$（§5.5），$G$ 不再是独立输入。

$$
\Longrightarrow\quad \text{换算的算术自洽（}\alpha/\beta=c_{\rm SI},\ \alpha c_{\rm SI}^2/\gamma=G_{\rm SI}\text{）是恒等式，不构成对理论的独立检验}.
\tag{5.2}
$$

### §5.5 物质↔时空交换给出三个单位

把物质与时空的相互转化写成同一个交换结构（Z144 定理 144.2；Z150 定理 150.3）：

$$
J_{i\to j}=K_{ij}(\phi_i-\phi_j),\qquad K=A\ (\text{邻接，唯一}),\qquad L\phi=-x .
\tag{5.3}
$$

三把尺子：跳数 $\alpha$（长度）、代 $\beta$（时间）、重数 $\gamma$（质量／荷）。两条转化率把它们串起来：

$$
c=\frac{\alpha}{\beta}\quad(\text{空间}\leftrightarrow\text{时间}),\qquad
G=\frac{\alpha c^2}{\gamma}\quad(\text{物质}\leftrightarrow\text{时空};\ \text{即}\ \gamma=\frac{\alpha c^2}{G}).
\tag{5.4}
$$

第二式就是"物质↔时空"的定量形式：**质量单位 = 长度单位 × $c^2/G$**。取交换的自然单位 $c=1$、$G=1$（一交换步走一跳；泊松核系数 $1$，$K=A$），由 $\alpha$ 一次给出

$$
\beta=\frac{\alpha}{c_{\rm SI}},\qquad \gamma=\frac{\alpha\,c_{\rm SI}^2}{G_{\rm SI}}.
\tag{5.5}
$$

数值：$\alpha=2.009\times10^{-33}$ m $\Longrightarrow$ $\beta=6.702\times10^{-42}$ s，$\gamma=2.706\times10^{-6}$ kg。核验（`unit_conversion_exchange_check.py`，4/4 PASS）：$c_{\rm SI}=\alpha/\beta$（相对差 $0$）、$G_{\rm SI}=\alpha c^2/\gamma$（相对差 $1.9\times10^{-16}$）、$\beta c_{\rm SI}/\alpha=1$、$\gamma G_{\rm SI}/(\alpha c_{\rm SI}^2)=1$。

$$
\boxed{\;\text{物质↔时空转化}=\text{泊松核的系数 }1:\quad \gamma=\frac{\alpha c^2}{G};\quad \text{自然单位下 }\alpha=\beta=\gamma.\;}
\tag{5.6}
$$

**两处必须一起读的边界**：

1. **$G=1$ 是交换给出的**：$K=A$ 唯一性把泊松核系数固定为 $1$（三维格林函数归一化给 $1/(4\pi)$，是同一件事的约定）；$G$ 不再是独立输入，而是 $1/Q$（舞台量）；
2. **$c=1$ 是单位选择**：交换更新是一阶的（$dq/dt=-\sum J$），"一交换步"定义了"代"；动力学波速是另一个系数（CFL 值 $c_{\rm disc}=0.8944$）。若取 $c_{\rm disc}=0.8944$，则 $\beta'=5.99\times10^{-42}$ s、$\gamma'=3.38\times10^{-6}$ kg（$\beta$ 差 $0.8944$ 倍、$\gamma$ 差 $1.25$ 倍）。故严格说：交换给出的是比值 $\beta/\alpha=1/c$、$\gamma/\alpha=c^2/G$；绝对尺度由 $\alpha$（锚 $v$）定，$c_{\rm disc}$ 仍是自由系数。

### §5.6 上游：$33.6197$ 是四个标度比之和

把 $\ln(M_X/v)$ 沿标度链展开（望远镜求和）：

$$
\ln\frac{M_X}{v}
=\underbrace{\ln\frac{M_Z}{v}}_{①\ \text{电弱}}
+\underbrace{\ln\frac{M_I}{M_Z}}_{②\ 2\pi k_0}
+\underbrace{\ln\frac{M_X^{\rm unif}}{M_I}}_{③\ \text{RG}}
+\underbrace{\ln\frac{M_X^{E_6}}{M_X^{\rm unif}}}_{④\ E_6\ \text{破缺}} .
\tag{5.7}
$$

| 环 | 内容 | 状态 |
|:--|:--|:--|
| ① | $\ln(M_Z/v)$：电弱比 | 在手（$M_Z$ 与 $v$ 均为实测） |
| ② | $\ln(M_I/M_Z)=2\pi k_0$，$k_0=2$ | 在手（Z154：$k_0=A_\ell^2=2$） |
| ③ | $\ln(M_X^{\rm unif}/M_I)$：RG 跑动 | 在手（Z155：中间谱＝SM＋3 代色单态双重态） |
| ④ | $\ln(M_X^{E_6}/M_X^{\rm unif})$：$E_6$ 破缺 | **开放**（离散相变，非 RG） |

$$
\boxed{\;\text{输入不是 }33.6197\text{，而是 }\textbf{环④ }(E_6\ \text{破缺标度});\ \text{前三环由电弱、绕数与 RG 给出}.\;}
\tag{5.8}
$$

**三处边界**：

1. **$M_I$ 与 $k_0=2$ 是框架绕数**：$\ln(M_I/M_Z)=2\pi k_0$ 是一条纯框架的指数（$k_0=2$）；$k_0=2$ 由 Z154 定出（$A_\ell=2\cos(\theta_M/2)=\sqrt2$，$k_0=A_\ell^2$），指数机制接到 SU(3) 骨架（Z096）。
2. **环④ 是离散相变**：在 $SU(3)^3$ 相里三耦合同步跑动，$M_X^{E_6}$ 对耦合不可见，故它不能由 RG 给出；它是 $E_6$ 破缺本身的标度。
3. **环③ 的数值依中间谱**：$M_I$ 与 $M_X^{\rm unif}$ 之间的内容决定 $\beta$ 系数；SM 内容给 $M_X^{\rm unif}=1.009\times10^{13}$ GeV、环③ $=12.88$；Z155 的中间谱（SM＋3 代色单态双重态）给 $4.954\times10^{13}$、环③ $=14.4546$。四环的**分割是记账**，总和 $S=33.6197$ 不变。

$$
\Longrightarrow\quad \text{本节的账从"一个 }33.6197\text{"变成"三环在手＋一环（}E_6\text{ 破缺）开放"}.
\tag{5.9}
$$

## §6 结论的适用范围

1. **锚是约定**：换锚（如 $m_e$）只改变单位的读数，无量纲比值不变（核验 F 组）；
2. **$v=246.22$ GeV 是实测值**：本框架不预言它的大小（Z013 定理 16）；
3. **$33.6197=\ln(M_X/v)$ 的上游是四环链**：① ② ③ 在手，④（$E_6$ 破缺）是上游开放项（§5.6）；
4. **$c_{\rm disc}$ 是方程系数，$G_{\rm disc}=1/Q$ 由交换定出**：$c$ 的数值不是本框架的预言；$G$ 由 $1/Q$ 与单位换算给出；
5. **质量单位由交换给出**：$\gamma=\alpha c^2/G$（$G=1/Q$），$G$ 不再是独立输入；
6. **有效数字**：本表给 4 位，输入 $v$ 与 $\ln(M_X/v)$ 的有效位决定末位。

## §7 结论

$$
\begin{aligned}
&\textbf{唯一锚}:\quad v=246.22\ \text{GeV}\ \text{定住全部有量纲量};\ \text{三个基本量纲归约为一个}.\\
&\textbf{长度}:\quad 1\ \text{跳数}=2.009\times10^{-33}\ \text{m}=124.3\,l_{\rm Pl}.\\
&\textbf{时间／质量}:\quad 1\ \text{代}=6.702\times10^{-42}\ \text{s};\quad 1\ \text{重数}=2.706\times10^{-6}\ \text{kg}.\\
&\textbf{能量}:\quad v=2.462\times10^{11}\ \text{eV};\ \text{质量按 }E=mc^2.\\
&\textbf{物质与时空}:\quad c=\alpha/\beta,\quad G=\alpha c^2/\gamma;\quad \text{自然单位下 }\alpha=\beta=\gamma.\\
&\textbf{账}:\quad 1\ \text{个量纲锚}+1\ \text{个上游开放数}\ (E_6\ \text{破缺环④})+c_{\rm disc}+\text{舞台量 }Q;\ \text{换算是换算，不是预言}.
\end{aligned}
$$

**一句话：**

$$
\boxed{\;\text{Zero 的米、秒、kg、eV 全部由电弱真空期望值 }v\text{ 一个锚定出；换算是算术，理论检验要靠锚无关的无量纲预言}.\;}
$$

---

## 附录 A　核验

`unit_conversion_check.py`（A–F）与 `unit_conversion_exchange_check.py`（G）：纯标准库；不 import 项目脚本；只写 results。

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | 锚点链 $M_X=v\,e^{33.6197}=9.82\times10^{16}$ GeV | 1 |
| B | 长度：$1$ 跳数 $=2.009\times10^{-33}$ m；$r_{\rm int}/l_{\rm Pl}=124.3$ | 2 |
| C | 时间／质量：$\beta=6.70\times10^{-42}$ s；$\gamma=2.71\times10^{-6}$ kg | 2 |
| D | 能量：$v=2.462\times10^{11}$ eV | 1 |
| E | 恒等式：$\alpha/\beta=c_{\rm SI}$、$\alpha c_{\rm SI}^2/\gamma=G_{\rm SI}$ | 2 |
| F | 换锚：无量纲比值与锚无关 | 1 |
| G | 交换自然单位：$c_{\rm SI}=\alpha/\beta$、$G_{\rm SI}=\alpha c^2/\gamma$、$\beta c_{\rm SI}/\alpha=1$、$\gamma G_{\rm SI}/(\alpha c_{\rm SI}^2)=1$ | 4 |
| | **合计** | **13** |

**判据**：A–D 相对差 $<10^{-2}$；E、F、G 相对差 $<10^{-12}$。本次 **13/13 PASS**（`unit_conversion_check.py` 9/9 ＋ `unit_conversion_exchange_check.py` 4/4）。

**核验纪律**：$\alpha,\beta,\gamma$ 的换算是解析式；数值复算只是旁证。$c_{\rm SI},G_{\rm SI},l_{\rm Pl}$ 用 SI 定义值。

---

*推导篇 · 已落盘 · 未改动任何公理 · 完*
