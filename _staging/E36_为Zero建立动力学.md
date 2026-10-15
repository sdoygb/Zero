# E36：为 Zero 建立动力学 —— **公理③ 的唯一实现**

**篇性**：**构造篇（暂存讨论稿）**。用户指示："那我们建立 Zero 自己的动力学吧"。
**结论**：**(一)** 动力学**不需要新公设** —— 它是公理③ 在 Z024 已有的**组态图**上的**唯一实现**；**(二)** 完整**补偿移动**给出的图是**正则的**（度 $=(L/2)^2$），故**均匀计数测度就是平稳分布**，与公理③ **自洽**；**(三)** 得到一个**真实的时标**（谱隙 $\to0$ 为 $1/L$）；**(四)** 但它**不产生净漂移**（时间可逆），且**再生的是度测度，不是账本测度**。
**纪律**：**未改动任何公理**。
**依据**：正式篇 Z001（公理③）、Z004（$\mathcal A_L$）、Z005（定理 7／11）、Z024（定理 2 的补偿移动与连通性）；本系列 E32（为何缺动力学）、E33／E34（L2 的转移）。
**核验**：[`E36_dynamics_check.py`](E36_dynamics_check.py)。

---

## 摘要

$$
\boxed{\;\textbf{（一）}\ \text{动力学}=\text{组态图上的}\ \textbf{均匀随机游走};\ \text{组态图}=\text{Z024 的补偿移动图} .\;}
\tag{0.1}
$$

$$
\boxed{\;\textbf{（二）}\ \text{该图}\ \textbf{正则} \text{（度}=(L/2)^2\text{）} \Longrightarrow \textbf{均匀计数测度是平稳分布} .\;}
\tag{0.2}
$$

$$
\boxed{\;\textbf{（三）}\ \text{谱隙}\ \to0\ \text{为}\ 1/L \Longrightarrow \text{有真实的时标} .\;}
\tag{0.3}
$$

$$
\boxed{\;\textbf{（四）}\ \textbf{平稳测度}\ \textbf{恰好} \text{是账本测度}\ w_h=b_h/\lvert\mathcal W_L\rvert\ \text{——}\textbf{不需任何选择} .\;}
\tag{0.4}
$$

---

## §1 出发点：E32 的缺口定位

$$
\textbf{E32 的结论}:\quad \text{以给定平稳分布}\ p\ \text{为平稳的流有}\ (n-1)^2\ \text{个自由度} \Longrightarrow \text{测度不能定出流} .
\tag{1.1}
$$

$$
\textbf{公理③ 逐字}:\quad \text{"运动不带权重——体系中不存在选择规则"};\quad \text{"每个允许的运动都被实例化（全分支），重数为计数"} .
\tag{1.2}
$$

$$
\Longrightarrow\quad \textbf{本文的主张}:\ \text{这两句话}\ \textbf{合起来} \text{唯一地定出一个动力学} .
\tag{1.3}
$$

---

## §2 §2 动力学不需要新建：图**已经**在 Z024 里

### 2.1 Z024 定理 2 逐字

> **定理 2（能量不守恒）.** 设 $\mathcal R$ 是满足容量约束的**连通**状态集。…
>
> 证明中用到：$\mathcal R$ 在**补偿移动**下**不可约**。

$$
\boxed{\;\textbf{故 Z024 早已}\ \textbf{隐含地} \text{拥有一个组态图}:\ \text{顶点}=\mathcal R,\ \text{边}=\text{补偿移动} .\;}
\tag{2.1}
$$

### 2.2 补偿移动的定义

$$
\textbf{补偿移动}:\quad w\longmapsto w'\quad\text{其中}\quad w_i=+1,\ w_j=-1\ (i\ne j);\quad
w'_i=-1,\ w'_j=+1 .
\tag{2.2}
$$

$$
\text{即：在}\ \textbf{任意} \text{一对}\ (+1,-1)\ \text{之间搬运一个单位} .
\tag{2.3}
$$

$$
\textbf{这}\ \textbf{正是} \text{公理① 的静态表述所允许的}\ \textbf{原生操作} \text{（保长、保零和）} .
\tag{2.4}
$$

---

## §3 关键判据：**均匀性 $\iff$ 正则性**

$$
\textbf{公理③ 的要求}:\quad \text{计数测度}\ \mu_{\rm unif}\ \text{必须是平稳分布} .
\tag{3.1}
$$

$$
\textbf{引理}:\quad \mu_{\rm unif}\ \text{在无向图上的随机游走下平稳}\quad\Longleftrightarrow\quad \text{图}\ \textbf{正则} .
\tag{3.2}
$$

*证明.* 平稳条件 $\sum_u\frac1{\deg u}\cdot\frac{A_{uv}}{\deg u}=\frac1{\deg v}$ 化简为 $\deg v$ 为常数。$\blacksquare$

$$
\boxed{\;\textbf{故判据简化为}:\ \text{补偿移动图是否}\ \textbf{正则}？\;}
\tag{3.3}
$$

### 3.1 一处**必须报告**的排除

$$
\textbf{相邻交换} \text{（只交换相邻的相反对）}:\quad \deg(w)=\#\{k:w_k\ne w_{k+1}\}\in[1,\ L-2]\ \textbf{可变} \Longrightarrow \textbf{不正则} .
\tag{3.4}
$$

$$
\textbf{核验 §A}:\ L{=}4:\deg\in\{1,2,3\};\quad L{=}6:\{1,\dots,5\};\quad L{=}8:\{1,\dots,7\} .
\tag{3.5}
$$

$$
\Longrightarrow\quad \text{相邻交换}\ \textbf{不合格} \text{（其平稳分布}\propto\deg\ne\text{均匀）} .
\tag{3.6}
$$

### 3.2 完整补偿移动**合格**（**核验 §B**）

$$
\textbf{核验 §B}:\quad \deg(w)=(\#\{+1\})\times(\#\{-1\})=\frac L2\cdot\frac L2=\boxed{\Bigl(\frac L2\Bigr)^2}\quad\textbf{对一切}\ w .
\tag{3.7}
$$

| $L$ | $4$ | $6$ | $8$ | $10$ |
|:--|--:|--:|--:|--:|
| $\lvert\mathcal W_L\rvert$ | $6$ | $20$ | $70$ | $252$ |
| $\deg$（实测，全同） | $4$ | $9$ | $16$ | $25$ |
| $(L/2)^2$ | $4$ | $9$ | $16$ | $25$ |

$$
\boxed{\;\Longrightarrow\quad \text{图}\ \textbf{正则} \Longrightarrow \mu_{\rm unif}\ \textbf{是平稳分布} \Longrightarrow \text{与公理③}\ \textbf{自洽} .\;}
\tag{3.8}
$$

$$
\textbf{核验 §B}:\ \text{实测平稳分布 max}\lvert p-1/n\rvert<5\times10^{-17}\ \checkmark .
\tag{3.9}
$$

---

## §4 连通性（**核验 §C**）

$$
\textbf{核验 §C}:\quad \text{补偿图在}\ L\in\{4,6,8,10\}\ \text{上}\ \textbf{连通}（\text{分量数}=1） .
\tag{4.1}
$$

$$
\Longrightarrow\quad \text{随机游走}\ \textbf{遍历} \Longrightarrow \text{平稳分布唯一} .
\tag{4.2}
$$

$$
\textbf{与 Z024 一致}:\ \text{Z024 已证}\ \mathcal R\ \text{在补偿移动下不可约} .
\tag{4.3}
$$

---

## §5 时标：谱隙（**核验 §D**）

| $L$ | $\lambda_1$ | $\lambda_2$ | 谱隙 $1-\lambda_2$ | 混合时间 $\sim1/\text{gap}$ |
|--:|--:|--:|--:|--:|
| $4$ | $1.000000$ | $0.000000$ | $1.000000$ | $1.00$ |
| $6$ | $1.000000$ | $0.333333$ | $0.666667$ | $1.50$ |
| $8$ | $1.000000$ | $0.500000$ | $0.500000$ | $2.00$ |

$$
\Longrightarrow\quad \text{谱隙}\ \textbf{递减} \Longrightarrow \text{混合时间}\ \textbf{随}\ L\ \text{增长} \Longrightarrow \textbf{有真实的时标} .
\tag{5.1}
$$

$$
\textbf{注意}:\ \text{这是}\ \textbf{时间} \text{（步数）};\ \text{要把它换成物理时间，仍需}\ \textbf{一个钟} \text{（Z005 定理 8 的}\ \mathcal C_{\rm clk}\text{）} .
\tag{5.2}
$$

---

## §6 这个动力学**是什么**、**不是什么**

| 特征 | 本文的动力学 | 标准动力学 |
|:--|:--|:--|
| 来源 | **公理③**（唯一实现） | 作用量／Hamiltonian |
| 是否需新公设 | **不需要** | — |
| 平稳分布 | **均匀计数测度**（与公理③自洽） | 由动力学定 |
| 时间可逆 | **是**（无向图） | 视情况 |
| **净漂移** | **无** | 可有 |
| 时标 | 有（谱隙） | 有 |
| **禁闭** | **不能** | QCD 能 |

$$
\boxed{\;\textbf{（核验 §E）}:\ \text{沿}\ \textbf{一切} \text{闭合回路，}\sum_{\rm cycle}\Delta=0 \Longrightarrow \textbf{无净漂移} .\;}
\tag{6.1}
$$

### 6.1 **关键**：平稳测度**恰好**是账本测度（本文的更强结果）

$$
\deg\ \textbf{恒定} \Longrightarrow \text{度测度}=\text{均匀测度}=\text{计数测度}
\Longrightarrow p(w)=\frac1{\lvert\mathcal W_L\rvert}\ \text{逐词均匀} .
\tag{6.2}
$$

$$
\Longrightarrow\quad \text{按层}:\quad p_h=\frac{b_h}{\lvert\mathcal W_L\rvert}=w_h\qquad\textbf{完全相同} .
\tag{6.3}
$$

**核验 §F**：

| $L$ | $p_h=w_h$ |
|--:|:--|
| $4$ | $(0.666667,\ 0.333333)$ |
| $6$ | $(0.4,\ 0.5,\ 0.1)$ |
| $8$ | $(0.228571,\ 0.542857,\ 0.2,\ 0.028571)$ |

$$
\boxed{\;\textbf{故}\ \text{动力学测度}\ \textbf{就是} \text{账本测度} \Longrightarrow \textbf{没有实质选择} .\;}
\tag{6.4}
$$

$$
\textbf{这一点比预期更强}:\quad
\begin{cases}
\text{公理③}\ \textbf{唯一地} \text{定出动力学};\\
\text{该动力学的平稳测度}\ \textbf{恰好} \text{是账本测度} .
\end{cases}
\tag{6.5}
$$

$$
\Longrightarrow\quad \text{账本测度}\ \textbf{不是} \text{"另一个输入"};\ \text{它是}\ \textbf{动力学自己的} \text{平稳测度} .
\tag{6.6}
$$

$$
\textbf{一处更正}:\ \text{本文初稿曾写"度测度}\ne\text{账本测度，需一个实质选择"}\ \textbf{—— 错} .
\tag{6.7}
$$

$$
\text{错因}:\ \text{把}\ \textbf{相邻交换} \text{（度可变）的度测度当成了}\ \textbf{完整补偿} \text{（度恒定）的} .
\tag{6.8}
$$

---

## §7 与前面各篇的关系

| 篇 | 结论 | 本文的关系 |
|:--|:--|:--|
| E32 | 缺动力学（$(n-1)^2$ 自由度） | **本文用公理③ 唯一地填上** |
| E33 | L2 的转移有漂移无转移 | 不同层；本文在 L0 |
| E34 | L2 的转移无局域性 | 不同层 |
| E31 | 约束给方向场 | 本文给**流**；方向仍来自谱序（Z005 定理 7） |
| E35 | 约束给计数力 | 本文给**动力学**；力仍在 |
| Z024 | 补偿移动 ＋ 连通性 | **本文的图就是它** |

$$
\boxed{\;\text{本文补上的正是 E32 所指的那个缺口}:\ \text{一个}\ \textbf{流} .\;}
\tag{7.1}
$$

---

## §8 还差什么（**诚实**）

$$
\textbf{（甲）}\ \text{无净漂移} \Longrightarrow \text{不产生}\ \textbf{定向运动} \text{（需要一个破对称的输入）} .
\tag{8.1}
$$

$$
\textbf{（乙）}\ \text{平稳测度}\ \textbf{已是} \text{账本测度} \Longrightarrow \text{壳权重}\ \textbf{内建} \text{于动力学（不再是独立输入）} .
\tag{8.2}
$$

$$
\textbf{（丙）}\ \text{无标度生成}:\ \text{谱隙}\sim1/L\ \text{是}\ \textbf{幂律} \Longrightarrow \text{仍无}\ \Lambda_{\rm QCD} \text{（E20 不变）} .
\tag{8.3}
$$

$$
\textbf{（丁）}\ \text{无禁闭}:\ \text{时间可逆、无净漂移} \Longrightarrow \text{无弦张力} .
\tag{8.4}
$$

$$
\boxed{\;\text{故本文}\ \textbf{建立了动力学} \text{，但}\ \textbf{不} \text{自动解决禁闭与标度} .\;}
\tag{8.5}
$$

---

## §9 不许主张

1. **不主张**本文"改"了公理——**没有**；动力学是公理③ 的实现（§3.2）；
2. **不主张**这是"唯一可能"的动力学——本文只证：在**补偿移动图**上、**均匀权重**下，它是唯一的（§3）；
3. **不主张** $L=4$ 时平稳分布 $=$ 账本——**不相等**（§6.1）；
4. **不主张**相邻交换"也行"——**不正则**，平稳分布 $\propto\deg$（§3.1）；
5. **不主张**本文解决了禁闭——**没有**（§8.4）；
6. **不主张**本文解决了标度——**没有**（§8.3）；
7. **不主张**谱隙 $=$ 物理时标——它需要**一个钟**（Z005 定理 8）；
8. **不主张**"无净漂移"是缺陷——它是**时间可逆**的等价说法。

---

## §10 结论

$$
\boxed{\;
\begin{aligned}
&\textbf{动力学存在且}\ \textbf{免费}:\ \text{组态图上的均匀随机游走（公理③ 的唯一实现）}.\\[1mm]
&\textbf{图}:\ \text{Z024 的补偿移动图};\quad \textbf{正则} \text{（度}=(L/2)^2\text{）};\quad \textbf{连通}.\\[1mm]
&\textbf{自洽}:\ \mu_{\rm unif}\ \text{是平稳分布} \Longrightarrow \text{与公理③ 一致（实测}<5\times10^{-17}\text{）}.\\[1mm]
&\textbf{时标}:\ \text{谱隙}\to0\ \text{为}\ 1/L;\ \text{混合时间}\sim L .\\[1mm]
&\textbf{更强}:\ \text{平稳测度}\ \textbf{恰好} \text{是账本测度}\ w_h\ \Longrightarrow\ \textbf{不需任何选择}.\\[1mm]
&\textbf{代价}:\ \textbf{无净漂移};\ \textbf{不解决} \text{禁闭与标度}.
\end{aligned}
}
$$

**一句话：**

$$
\boxed{\;\textbf{动力学一直在 Zero 里 —— 它就是公理③ 在 Z024 早已画好的那张组态图上的均匀随机游走。}\;}
$$

---

## 附录 A　核验

`E36_dynamics_check.py`：纯标准库（＋`numpy`）；不 import 项目脚本、不写文件。

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | 相邻交换的度可变（不正则） | 6 |
| B | 完整补偿图的正则性（度 $=(L/2)^2$）与平稳分布均匀 | 9 |
| C | 连通性（$L\le10$） | 5 |
| D | 谱隙与混合时间（$L\le8$） | 6 |
| E | 无净漂移（时间可逆） | 5 |
| F | 平稳测度 **$=$** 账本测度（逐层相同） | 6 |
| | **合计** | **49** |

**核验纪律**：数值只算旁证；**本文更正了初稿一处错误**（§6.7：度测度其实**等于**账本测度）。

---

*构造稿 · 未落盘 · 未改动任何公理 · 完*
