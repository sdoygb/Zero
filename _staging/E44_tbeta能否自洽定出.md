# E44：$t\beta$ 能否由自洽条件定出 —— **不能**（同齐性障碍）

**篇性**：**检验 ＋ 否证篇（暂存讨论稿）**。用户指示："(甲) 检验 $t\beta$ 能否由自洽条件定出"。
**结论**：**(一)** 相位**只依赖**组合 $t\beta$（$=$ 标度不变），故**任何只涉及相位的自洽条件都无法定出 $t$ 或 $\beta$ 的单独值**；**(二)** 权重**依赖 $\beta$ 单独**，故 $\beta$ 原则上可由【分布侧】的条件定；**(三)** 但**所有自然的分布侧条件**（归一化／固定 $\langle r\rangle$／最大熵）都**只给 $\beta=0$**（退回账本态）；**(四)** 故 **(甲) 是否定的**：$t\beta$ **不可由 Zero 的结构定出**；**(五)** 而 $t$ 本身是**具名输入**（Z005 定理 8）。
**纪律**：**未改动任何公理**；**暂存稿**，未落盘到根目录。根目录 Z001–Z079 未改动。
**依据**：正式篇 Z001（§6.1 参数账；§6.2 逐项登记）；Z005（定理 3／定理 8；KMS 条件）；Z004（定理 5）；本系列 E42（最小破缺）、E43（$\theta=n\,t\beta$）。
**核验**：[`E44_tbeta_check.py`](E44_tbeta_check.py)。

---

## 摘要

$$
\boxed{\;\textbf{（一）}\ \text{相位}=t\,\beta\,n \Longrightarrow \textbf{只依赖组合}\ t\beta \Longrightarrow \textbf{标度不变} .\;}
\tag{0.1}
$$

$$
\boxed{\;\textbf{（二）}\ \textbf{权重} \text{依赖}\ \beta\ \textbf{单独} \Longrightarrow \beta\ \text{原则上可由分布侧条件定} .\;}
\tag{0.2}
$$

$$
\boxed{\;\textbf{（三）}\ \text{但所有自然的分布侧条件}\ \textbf{只给}\ \beta=0 .\;}
\tag{0.3}
$$

$$
\boxed{\;\textbf{（四）}\ \text{故}\ t\beta\ \textbf{不可} \text{由结构定};\quad t\ \text{是}\ \textbf{具名输入} .\;}
\tag{0.4}
$$

---

## §1 设置

$$
\textbf{E42 的态}:\quad p_\beta(w)\propto e^{-\beta\,r(w)};\qquad \rho_\beta=\operatorname{diag}(p_\beta) .
\tag{1.1}
$$

$$
\textbf{模生成元}:\quad K=-\log\rho_\beta;\qquad
\textbf{E43 的结果}:\quad \theta=n\,t\,\beta .
\tag{1.2}
$$

$$
\textbf{本文的问题}:\quad \text{是否存在自洽条件定出}\ t\beta\ ?
\tag{1.3}
$$

---

## §2 **同齐性障碍**（**核验 §A**）

$$
\textbf{相位}:\quad \varphi=t\,(K_j-K_i)=t\,\beta\,n .
\tag{2.1}
$$

$$
\textbf{标度变换}:\quad t\to\lambda t,\qquad \beta\to\frac{\beta}{\lambda}
\Longrightarrow t\beta\ \textbf{不变} \Longrightarrow \text{一切相位}\ \textbf{不变} .
\tag{2.2}
$$

$$
\textbf{核验 §A}:\quad (t,\beta)=(1,1)\ \text{与}\ (2,0.5)\ \text{给}\ \textbf{相同} \text{的相位}\ (0,\ 1.0,\ 2.0) .
\tag{2.3}
$$

$$
\textbf{核验 §A}:\quad \text{但两者给}\ \textbf{不同} \text{的权重}:\quad
\begin{cases}(1,1) & p\propto(7.389,\ 1,\ 2.718,\ 2.718,\ 1,\ 7.389)\\
(2,0.5) & p\propto(2.718,\ 1,\ 1.649,\ 1.649,\ 1,\ 2.718)\end{cases}
\tag{2.4}
$$

$$
\boxed{\;\Longrightarrow\quad \textbf{任何只涉及相位的自洽条件}\ \textbf{都无法} \text{定出}\ t\ \text{或}\ \beta\ \text{的单独值} .\;}
\tag{2.5}
$$

$$
\textbf{但}:\quad \textbf{权重} \text{依赖}\ \beta\ \textbf{单独} \Longrightarrow \beta\ \text{可由【分布侧】条件定（若存在）} .
\tag{2.6}
$$

---

## §3 分布侧的自然条件**逐个检验**（**核验 §B**）

$$
L=4:\quad r\ \text{的计数}=(2{:}2,\ 3{:}2,\ 4{:}2);\qquad \langle r\rangle_{\rm unif}=3 .
\tag{3.1}
$$

| 条件 | 结果 |
|:--|:--|
| 归一化 $\sum_wp=1$ | **恒成立**（对一切 $\beta$）$\Longrightarrow$ **不定** $\beta$ |
| 固定 $\langle r\rangle=\langle r\rangle_{\rm unif}=3$ | **解出 $\beta=0$**（退回账本态） |
| 最大熵 | **$\beta=0$**（无约束的极大熵即均匀） |
| KMS 条件 $\omega\circ\sigma_t=\omega$ | **恒成立**（$\rho$ 对角 $\Rightarrow\sigma_t$ 对角）$\Longrightarrow$ **不定** |

$$
\textbf{核验 §B}:\quad \sum_wp_\beta(w)=1\ \text{对全部}\ \beta\ \text{成立} \Longrightarrow \text{归一化不给}\ \beta .
\tag{3.2}
$$

$$
\textbf{核验 §B}:\quad \langle r\rangle_\beta=3\iff\beta=0 .
\tag{3.3}
$$

$$
\boxed{\;\textbf{故}:\ \text{所有自然的分布侧条件}\ \textbf{只给}\ \beta=0 \text{—— 即退回中心态} .\;}
\tag{3.4}
$$

---

## §4 严格陈述：**能定什么、不能定什么**

$$
\textbf{能定的}:\quad \text{给定}\ \beta,\ \text{一切}\ \textbf{分布侧量} \text{（权重、}\langle r\rangle\text{、熵、}\omega([A,B])\text{）都可算} .
\tag{4.1}
$$

$$
\textbf{不能定的}:\quad
\begin{cases}
\beta & \text{无自然条件把它定到}\ \ne0\\
t & \textbf{只能} \text{与}\ \beta\ \text{一起进入物理} \Longrightarrow \textbf{本质上} \text{不可由相位条件定}\\
t\beta & \text{相位侧无条件};\ \text{分布侧条件只给}\ 0
\end{cases}
\tag{4.2}
$$

$$
\boxed{\;\Longrightarrow\quad \textbf{故 (甲) 的检验}\ \textbf{否定}:\ t\beta\ \text{不可由 Zero 的结构定出} .\;}
\tag{4.3}
$$

$$
\textbf{而}\ t\ \text{本身}:\quad \textbf{Z005 定理 8}:\ \text{钟}\ \mathcal C_{\rm clk}\ \text{是}\ \textbf{具名输入} \text{，不由公理推出} .
\tag{4.4}
$$

---

## §5 与 Z001 的**参数账**对照

$$
\textbf{Z001 §6.1 逐字}:\quad Z_0+\underbrace{n}_{\text{自由参数}}+\underbrace{\pi}_{\text{结构输入}}+\underbrace{\kappa}_{\text{尺度输入}}=\text{完整框架} .
\tag{5.1}
$$

$$
\textbf{Z001 §6.2 逐项登记}:\quad n,\ T\ (\text{自由参数});\quad \pi,\ \kappa\ (\text{具名输入});\quad \beta\varepsilon\ (\text{参数}) .
\tag{5.2}
$$

$$
\textbf{本文的发现}:\quad \beta\ \text{应作为}\ \textbf{第五项} \text{登记} \text{（E42 的最小破缺引入它）} .
\tag{5.3}
$$

$$
\Longrightarrow\quad \textbf{参数账更新}:\quad Z_0+n+T+\pi+\kappa+\beta\varepsilon+\mathcal C_{\rm clk}+\beta .
\tag{5.4}
$$

$$
\textbf{但}:\ \beta\ \text{与}\ t\ \textbf{只以乘积进入物理相位} \Longrightarrow \text{两者}\ \textbf{合并为一条具名输入} .
\tag{5.5}
$$

$$
\boxed{\;\textbf{故净增一条具名输入}:\quad \Theta:=t\beta\ (\textbf{无量纲}).\;}
\tag{5.6}
$$

---

## §6 一处**必须写的对照**（证据纪律）

$$
\textbf{本文也对}\ \theta\ \text{做了组合搜索}:\quad \text{从}\ \{L,\ H,\ \tau_{\rm rel},\ \log2,\ \log3,\ \log6,\ 5/9,\ \pi,\ 2\pi\}\ \text{的乘除组合里找}\ \theta/n .
\tag{6.1}
$$

$$
\textbf{核验 §C}:\quad 7\ \text{个组合落在}\ 2\times10^{-3}\ \text{内（最好者差}\ 6.4\times10^{-4}） .
\tag{6.2}
$$

$$
\textbf{但}:\quad \textbf{Q275 逐字}:\ \text{"搜索无效（低分母有理数必然被命中，}\textbf{无判别力}\text{）"} .
\tag{6.3}
$$

$$
\boxed{\;\textbf{故本文}\ \textbf{不主张} \text{任何组合命中} \text{—— 它们}\ \textbf{不构成证据} .\;}
\tag{6.4}
$$

---

## §7 不许主张

1. **不主张** $t\beta$ "可以被定出"——**已否证**（§4.3）；
2. **不主张** $\beta=0$ "是唯一解"——本文只说**自然条件只给 $\beta=0$**（§3.4）；
3. **不主张**同齐性"意味着 $t$ 无意义"——它意味着 $t$ **不可单独观测**（§4.2）；
4. **不主张**任何组合命中是证据——**Q275 纪律**（§6.4）；
5. **不主张** $\beta$ "必须登记"已经完成——本文只是**建议**（§5.3）；
6. **不主张**本文改了 Z001 的参数账——它的表**未改**；本文只**指出**需补一项；
7. **不主张**本文改了任何正式篇——**暂存稿**，未落盘。

---

## §8 结论

$$
\boxed{\;
\begin{aligned}
&\textbf{同齐性障碍}:\ \text{相位}=t\beta n \Longrightarrow \textbf{只依赖组合} \Longrightarrow \textbf{标度不变} .\\[1mm]
&\textbf{故}:\ \text{任何相位侧的自洽条件}\ \textbf{都无法} \text{定}\ t\ \text{或}\ \beta .\\[1mm]
&\textbf{而分布侧}:\ \text{依赖}\ \beta\ \text{单独，但所有自然条件}\ \textbf{只给}\ \beta=0 .\\[1mm]
&\textbf{故 (甲)}\ \textbf{否定}:\ t\beta\ \text{不可由结构定};\quad t\ \text{是}\ \textbf{具名输入} .\\[1mm]
&\textbf{净}:\ \text{Zero 须多登记}\ \textbf{一条无量纲具名输入}\ \Theta=t\beta .
\end{aligned}
}
$$

**一句话：**

$$
\boxed{\;\textbf{$t$ 与 $\beta$ 只以乘积出现，所以结构只能定这个乘积的"地位"，定不出它的"值"。}\;}
$$

---

## 附录 A　核验

`E44_tbeta_check.py`：纯标准库（＋`numpy`）；不 import 项目脚本、不写文件。

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | 同齐性：$(1,1)$ 与 $(2,0.5)$ 同相位、异权重 | 10 |
| B | 分布侧条件：归一化恒成立；$\langle r\rangle=3\iff\beta=0$；最大熵 | 10 |
| C | 组合搜索（7 个命中）＋ Q275 纪律 | 8 |
| D | 能定／不能定的清单 | 6 |
| E | 语料引文（Z001 §6.1／§6.2；Z005 定理 8） | 8 |
| | **合计** | **42** |

**核验纪律**：数值只算旁证；**本文不主张任何组合命中**（§6.4）。

---

*检验稿 · 暂存 · 未改动任何公理 · 完*
