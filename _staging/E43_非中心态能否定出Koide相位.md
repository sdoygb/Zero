# E43：非中心态能否定出 Koide 的 $\theta$ —— **不能**，但它给出 $\theta$ 的精确形式（附一处自我更正）

**篇性**：**检验 ＋ 否证 ＋ 更正篇（暂存讨论稿）**。用户指示："(甲)"（检验 $p_\beta\propto e^{-\beta r}$ 的新相位能否定出 Koide 的 $\theta$）。
**结论**：**(一)** **能给出** $\theta$ 的**形式**：$\theta=n\,t\,\beta$（$n$ 为小整数）——即 $\theta$ **确实是一个相对相位**，与 Koide 的结构**相符**；**(二)** **但不能定出**其值：$t$（钟）与 $\beta$（倾角）**都不由结构定** $\Longrightarrow$ $\theta$ **仍是自由的**；**(三)** 故 (甲) 的检验结果是**否定**的，**但原因比预期的精确**；**(四)** **本文更正自己两处**（见 §6）。
**纪律**：**未改动任何公理**；**暂存稿**，未落盘到根目录。根目录 Z001–Z079 未改动。
**依据**：正式篇 Z001（公理③）；Z004（定理 5）；Z005（定理 8：钟是具名输入）；Z014（定理 27–32；Koide 走廊）；本系列 E38（反推上游）、E39（最小代价）、E42（最小破缺）。
**核验**：[`E43_theta_check.py`](E43_theta_check.py)。

---

## 摘要

$$
\boxed{\;\textbf{（一）}\ \text{非中心态}\ \textbf{能} \text{给出}\ \theta\ \text{的形式}:\quad \theta=n\,t\,\beta .\;}
\tag{0.1}
$$

$$
\boxed{\;\textbf{（二）}\ \textbf{但不能定出} \text{其值}:\ t\ \text{与}\ \beta\ \text{都不由结构定} .\;}
\tag{0.2}
$$

$$
\boxed{\;\textbf{（三）}\ \text{故 (甲) 的检验}\ \textbf{否定} \text{，但原因精确} .\;}
\tag{0.3}
$$

$$
\boxed{\;\textbf{（四）}\ \text{须把态扩展到}\ \textbf{词的算符代数} \text{（非对角相干）} \text{—— 那是}\ \textbf{新结构} .\;}
\tag{0.4}
$$

---

## §1 检验的设计

$$
\textbf{E42 的候选态}:\quad p_\beta(w)\propto e^{-\beta\,r(w)};\qquad \rho_\beta=\operatorname{diag}(p_\beta) .
\tag{1.1}
$$

$$
\textbf{问题}:\quad \text{它的新独立相位能否}\ \textbf{定出} \ \theta=2.3166247339\ \text{rad}\ ?
\tag{1.2}
$$

$$
\textbf{Koide 的复形式}:\quad m_i=A^2\Bigl(1+b\cos\Bigl(\theta+\frac{2\pi i}3\Bigr)\Bigr)^2
\Longrightarrow z_i=A\Bigl(1+b\,e^{i(\theta+2\pi i/3)}\Bigr) .
\tag{1.3}
$$

$$
\boxed{\;\Longrightarrow\quad \theta=\text{两个复项之间的}\ \textbf{相对相位} .\;}
\tag{1.4}
$$

$$
\textbf{模流给的相位（**核验 §A**）}:\quad \sigma_t\bigl(\lvert i\rangle\langle j\rvert\bigr)=e^{\,i\,t\,(K_j-K_i)}\lvert i\rangle\langle j\rvert .
\tag{1.5}
$$

$$
\Longrightarrow\quad \text{相位}=t\,(K_j-K_i) \Longrightarrow \textbf{形式相符} .
\tag{1.6}
$$

---

## §2 计算（**核验 §B**）

$$
L=4:\quad p_\beta\propto e^{-\beta r(w)};\qquad r\in\{2,3,4\}\ (\text{各 2 个}) .
\tag{2.1}
$$

$$
\textbf{核验 §B}:\quad \beta=0.3\Rightarrow K\in\{1.5215,\ 1.8215,\ 2.1215\},\ K\ \text{差}\in\{0.3,\ 0.6\} .
\tag{2.2}
$$

$$
\textbf{核验 §B}:\quad \beta=1.0\Rightarrow K\in\{1.1008,\ 2.1008,\ 3.1008\},\ K\ \text{差}\in\{1.0,\ 2.0\} .
\tag{2.3}
$$

$$
\Longrightarrow\quad K\ \text{差}=\beta\cdot n,\ n\in\{1,2\} \Longrightarrow
\boxed{\;\theta=t\,\beta\,n\;}\qquad(t=\text{模流时间}) .
\tag{2.4}
$$

---

## §3 为什么**定不出**值

$$
\theta=n\,t\,\beta\ \text{含}\ \textbf{两个} \text{自由量}:\quad
\begin{cases}
t & \textbf{钟} \text{—— Z005 定理 8：}\textbf{具名输入} \text{，不由公理推出}\\
\beta & \textbf{倾角} \text{—— E42：新连续参数，值不由结构定}
\end{cases}
\tag{3.1}
$$

$$
\Longrightarrow\quad \text{对任意}\ \theta,\ \text{都可取}\ t\beta=\theta/n \Longrightarrow \textbf{结构无法排除任何}\ \theta .
\tag{3.2}
$$

$$
\boxed{\;\textbf{故 (甲) 的检验}\ \textbf{否定}:\ \text{非中心态给}\ \theta\ \text{的}\ \textbf{形式} \text{，但不给}\ \textbf{值} .\;}
\tag{3.3}
$$

$$
\textbf{与 E38 一致}:\quad \text{缺的是}\ \textbf{一个无量纲数};\ \text{而此处它表现为}\ t\beta .
\tag{3.4}
$$

---

## §4 与 Zero 已有数的对照（**核验 §C**）

| 候选值 | 值 | $\theta/\text{候选}$ |
|:--|--:|--:|
| $1/2$ | $0.5$ | $4.633249$ |
| $5/9$ | $0.555556$ | $4.169925$ |
| $2/3$ | $0.666667$ | $3.474937$ |
| $\log2$ | $0.693147$ | $3.342183$ |
| $\log3$ | $1.098612$ | $2.108683$ |
| $\theta/2$ | $1.158312$ | $2.000000$ |
| $\pi/2$ | $1.570796$ | $1.474809$ |
| $\log6$ | $1.791759$ | $1.292933$ |

$$
\textbf{核验 §C}:\quad \textbf{没有任何} \text{候选给出整数或简单有理比} \Longrightarrow \theta\ \text{与 Zero 现有数}\ \textbf{无关} \text{（在现有结构下）} .
\tag{4.1}
$$

$$
\textbf{一处须写的警告}:\quad \text{低分母数的"接近"}\ \textbf{没有判别力} \text{（Q275 的纪律）} \Longrightarrow \text{本文}\ \textbf{不主张} \text{任何接近} .
\tag{4.2}
$$

---

## §5 要 $\theta$ 需要什么（精确要求）

$$
\textbf{缺口}:\quad \text{对角态（经典）}\ \textbf{给不出} \text{独立连续相位} .
\tag{5.1}
$$

$$
\textbf{理由（核验 §D）}:\quad \rho=\operatorname{diag}(p)\ \text{与词基}\ \textbf{同时对角} \Longrightarrow
\text{模流的"方向"就是词本身} \Longrightarrow \text{词是}\ \textbf{经典标签} .
\tag{5.2}
$$

$$
\textbf{须做的扩展}:\quad \text{把态放在}\ \textbf{词的算符代数} \text{上} \Longrightarrow \text{允许}\ \textbf{非对角相干} \text{（coherence）} .
\tag{5.3}
$$

$$
\boxed{\;\Longrightarrow\quad \text{"权重依赖}\ r\text{"}\ \textbf{不足} \text{；须}\ \textbf{相干} \text{—— 那是新结构，}\textbf{不是} \text{E42 的}\ \beta .\;}
\tag{5.4}
$$

$$
\textbf{与 E39 的关系}:\quad \text{E39 说"独立连续参数须块内非标量"};\ \text{本文说}\ \textbf{更进一步} \text{：须}\ \textbf{非对角} .
\tag{5.5}
$$

---

## §6 **本文更正自己两处**

$$
\textbf{更正一}:\quad \text{本文曾写"模流相位无观测意义，因词是经典标签"}\quad\textbf{—— 错} .
\tag{6.1}
$$

$$
\textbf{理由}:\quad \lvert i\rangle\langle j\rvert\ (i\ne j)\ \textbf{是} \ M_n\ \text{的元素} \Longrightarrow \text{其相位}\ \textbf{有} \text{观测意义} .
\tag{6.2}
$$

$$
\textbf{更正二}:\quad \text{本文曾写"$L=4$ 上第一层内}\ r{=}3\ \text{与}\ r{=}4\ \text{各 2 个，故层内比}=1"\quad\textbf{—— 表述错} .
\tag{6.3}
$$

$$
\textbf{正解}:\quad \text{层内比}\ n_{r=4}/n_{r=3}=2/2=1\ \text{确实为 1，}\textbf{但这不影响}\ K\ \text{差};\quad K\ \text{差}=\beta\cdot n\ \text{仍成立} .
\tag{6.4}
$$

$$
\Longrightarrow\quad \textbf{主结论（§3）不受影响} .
\tag{6.5}
$$

---

## §7 不许主张

1. **不主张** $\theta$ "已被解释"——**没有**（§3.3）；
2. **不主张** $\theta=n t\beta$ "就是答案"——它只是**形式相符**（§1.6）；
3. **不主张**任何"接近"是证据——**Q275 纪律**（§4.2）；
4. **不主张**"相干"能解决问题——**未验**（§5.3 只说它是**要求**）；
5. **不主张** E42 "错了"——它关于**中心性**的结论正确；本文只指出**它不足**（§5.4）；
6. **不主张**钟 $t$ "可以被定出"——Z005 定理 8 说它是**具名输入**；
7. **不主张**本文改了任何正式篇——**暂存稿**，未落盘。

---

## §8 结论

$$
\boxed{\;
\begin{aligned}
&\textbf{检验结果}:\ \textbf{否定} \text{—— 非中心态定不出}\ \theta\ \text{的值} .\\[1mm]
&\textbf{但它给出形式}:\quad \theta=n\,t\,\beta \Longrightarrow \theta\ \textbf{确实是相对相位} \text{（与 Koide 结构相符）} .\\[1mm]
&\textbf{两个自由量}:\ t\ (\text{钟，具名输入})\ \text{与}\ \beta\ (\text{倾角，新参数}) \Longrightarrow \textbf{结构无法排除任何}\ \theta .\\[1mm]
&\textbf{故须}:\ \text{把态扩展到}\ \textbf{词的算符代数} \text{上的}\ \textbf{相干} \text{—— 那是新结构} .\\[1mm]
&\textbf{与 E38 一致}:\ \text{缺的始终是}\ \textbf{一个无量纲数} .
\end{aligned}
}
$$

**一句话：**

$$
\boxed{\;\textbf{非中心态给出了 $\theta$ 的正确"形状"（一个相对相位），却给不出它的值 —— 因为值里含着钟，而钟是具名输入。}\;}
$$

---

## 附录 A　核验

`E43_theta_check.py`：纯标准库（＋`numpy`）；不 import 项目脚本、不写文件。

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | 模流的相位公式；$\theta$ 是相对相位 | 6 |
| B | $L=4$ 的 $K$ 与 $K$ 差（$\beta=0.3,1.0$）；$K$ 差 $=\beta n$ | 10 |
| C | 与 Zero 已有数的对照（无整数比） | 8 |
| D | 对角态与词基同时对角；须非对角相干 | 8 |
| E | 语料引文（Z005 定理 8；Z014 定理 27–32） | 8 |
| F | 两处自我更正的核验 | 6 |
| | **合计** | **48** |

**核验纪律**：数值只算旁证；**本文更正自己两处**（§6）；**不主张任何接近**（§4.2）。

---

*检验稿 · 暂存 · 未改动任何公理 · 完*
