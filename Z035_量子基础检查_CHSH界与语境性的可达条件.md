# 量子基础检查：CHSH 界、中心态的经典结构与语境性的可达条件

**篇性**：推导篇 $\cdot$ **两条完整证明 ＋ 一处否证 ＋ 一个判据**。本文检查缺口审计的**量子基础**一栏（`Bell`／`EPR`／`语境性`），逐条给出**完整证明**而非引用。

**Zero 项目 · 基础推导系列 · 第 Z035 号**

**编号**：本文为 $Z_{35}$。前篇 $Z_{34}$（量子信息与开放系统）。
**前件**：[量子运动学](Z004_量子运动学_代数与Born规则.md)（$\mathcal A_L=\bigoplus_hM_2(\mathbb C)$、Born 规则）；[纠缠关联与复合](Z008_纠缠关联与复合.md)（纠缠闭式律、CHSH 的硬结果）；[概率与测量](Z007_概率与测量_效应与退相干.md)（账本态）。
**层指标**：$\mathrm{L1}$（量子态）、$\mathrm{L1}'$（账本态）、$\mathcal R$（联合分布）。
**核验**：[`Z035_qfound_check.py`](_staging/Z035_qfound_check.py)（**26 条断言 / 0 不符**）。

---

## 摘要

审计指出："纠缠与 CHSH 有硬结果"，但**没有 Bell／EPR／语境性的系统检查**。本文补上，并把每一步**证完**。

$$
\boxed{\;\textbf{定理 23}:\ \text{CHSH 的经典界是 }2,\ \text{量子界是 }2\sqrt2;\ \text{本框架的代数}\textbf{允许}\text{后者}.\;}
\tag{0.1}
$$

| # | 内容 | 档位 |
|--:|:--|:--|
| **一** | **CHSH 经典界 $=2$**（完整证明：凸性 ＋ 确定性模型） | **已证** |
| **二** | **Tsirelson 界 $=2\sqrt2$**（上半：完整证明；达界：构造 ＋ 数值） | **已证** |
| **三** | **否证"中心态 $\Rightarrow$ 因子化"**：$\langle AB\rangle\ne\langle A\rangle\langle B\rangle$（中心态） | **排除** |
| **四** | **语境性判据**：块标量（中心）读数 $\Rightarrow$ 交换 $\Rightarrow$ **无语境性**；块内读数 $\Rightarrow$ 非交换 $\Rightarrow$ **可能** | **条件证成** |

$$
\boxed{\;\text{净效果}:\ \text{"量子基础缺席"不成立};\ \text{但}\textbf{语境性的可达性由读出的分辨率决定}.\;}
\tag{0.2}
$$

---

## §1 记号与前件

### 1.1 对象

| 符号 | 名称 | 定义 |
|:--|:--|:--|
| $\mathcal A_L=\bigoplus_{h=1}^{H}M_2(\mathbb C)$ | 可观察代数 | $H=\lfloor L/2\rfloor$ |
| $\mathcal Z_L=\bigoplus_h\mathbb C I_2$ | **中心** | 块标量算符 |
| $\rho$ | 量子态 | $\rho\ge0$，$\operatorname{tr}\rho=1$ |
| $\rho_{\rm led}=\bigoplus_h\frac{w_h}{2}I_2$ | **账本态** | $w_h\ge0$，$\sum_hw_h=1$ |
| $E(a,b)$ | 关联函数 | $\operatorname{tr}(\rho\,\sigma_a\otimes\sigma_b)$ |
| $S$ | **CHSH 量** | $E(a,b)-E(a,b')+E(a',b)+E(a',b')$ |

### 1.2 两条前件

> **前件 A（Born 规则）.** 对 $A\in\mathcal A_L$，$\langle A\rangle=\operatorname{tr}(\rho A)$。

> **前件 B（账本态是中心元）.** $\rho_{\rm led}\in\mathcal Z_L$，即 $[\rho_{\rm led},A]=0$ 对一切 $A\in\mathcal A_L$。

**本文的核验确认前件 B**（核验 §D2c）。

---

## §2 主结果一：CHSH 的经典界（完整证明）

### 2.1 陈述与证明

> **定理 23（经典界）.** 设 $\lambda$ 取值于任意测度空间，$A_1,A_2,B_1,B_2:\lambda\mapsto\{\pm1\}$ 为任意函数。则
> $$
> \bigl\lvert A_1B_1-A_1B_2+A_2B_1+A_2B_2\bigr\rvert\ \le\ 2 .
> $$

*证明（两行）.* 记 $s=A_1B_1-A_1B_2+A_2B_1+A_2B_2$。由 $B_j^2=1$：
$$
s=A_1(B_1-B_2)+A_2(B_1+B_2).
$$
对每个 $\lambda$，$B_1,B_2\in\{\pm1\}$ 恰使 $(B_1-B_2,B_1+B_2)$ 取 $(0,\pm2)$ 或 $(\pm2,0)$。故该 $\lambda$ 处**恰有一项非零**，且其绝对值 $=2$。因此 $\lvert s\rvert\le2$ 逐点成立；对 $\lambda$ 积分仍 $\le2$。$\blacksquare$

$$
\boxed{\;\text{关键}:\ \text{逐点}\ \lvert s\rvert\le2;\ \text{故对任意（含隐变量的）局域模型都成立}.\;}
\tag{2.1}
$$

**推论（凸性）**：局域模型的关联函数集合是**凸集**，其极点由确定性模型给出。由定理 23，凸集含于 $\{\lvert S\rvert\le2\}$ 的条带内。

$$
\Longrightarrow\quad \boxed{\;\text{经典界}\ \lvert S\rvert\le2\ \text{对一切局域实在论模型成立}.\;}
\tag{2.2}
$$

**核验 §A**：穷举 $2^4=16$ 个确定性模型，$\max\lvert S\rvert=2$（**精确达到，非估计**）。

### 2.2 乘积态不违反

> **命题 18.** 若 $\rho=\rho_1\otimes\rho_2$，则 $E(a,b)=\langle\sigma_a\rangle_1\langle\sigma_b\rangle_2$，且 $\lvert S\rvert\le2$。

*证明.* 期望值在张量积下因子化：
$$
E(a,b)=\operatorname{tr}\bigl[(\rho_1\otimes\rho_2)(\sigma_a\otimes\sigma_b)\bigr]
=\operatorname{tr}(\rho_1\sigma_a)\,\operatorname{tr}(\rho_2\sigma_b)=\alpha_a\beta_b ,
$$
其中 $\alpha_a:=\langle\sigma_a\rangle_1$，$\beta_b:=\langle\sigma_b\rangle_2$，且 $\lvert\alpha_a\rvert,\lvert\beta_b\rvert\le1$（因 $\sigma_a^2=I$ 给 $\lvert\langle\sigma_a\rangle\rvert\le1$）。

把 $S$ 按 $\beta$ 归并：
$$
S=\alpha_a\beta_b-\alpha_a\beta_{b'}+\alpha_{a'}\beta_b+\alpha_{a'}\beta_{b'}
=(\alpha_a+\alpha_{a'})\beta_b+(\alpha_{a'}-\alpha_a)\beta_{b'} .
$$
用 $\lvert\beta_b\rvert,\lvert\beta_{b'}\rvert\le1$ 与三角不等式：
$$
\lvert S\rvert\le\lvert\alpha_a+\alpha_{a'}\rvert+\lvert\alpha_{a'}-\alpha_a\rvert
=2\max\bigl(\lvert\alpha_a\rvert,\lvert\alpha_{a'}\rvert\bigr)\le2 ,
$$
最后一步用恒等式 $\lvert x+y\rvert+\lvert y-x\rvert=2\max(\lvert x\rvert,\lvert y\rvert)$（取 $x=\alpha_a$，$y=\alpha_{a'}$）。$\blacksquare$

$$
\boxed{\;\text{关键恒等式}:\ \lvert x+y\rvert+\lvert y-x\rvert=2\max(\lvert x\rvert,\lvert y\rvert);\ \text{这是}\ \mathbb R^2\ \text{中}\ \ell^1\ \text{范数的旋转不变性}.\;}
\tag{2.3}
$$

**核验 §A1**：乘积态 $\lvert00\rangle$ 给 $S=0$。

---

## §3 主结果二：Tsirelson 界（上半完整证明，达界构造）

### 3.1 上半：$\lvert S\rvert\le2\sqrt2$

> **定理 24（Tsirelson 界）.** 对任意态 $\rho$ 与任意方向 $a,a',b,b'$，$\lvert S\rvert\le2\sqrt2$。

*证明（用算子范数，四步）.*

**第一步：把 $S$ 写成算子的期望。** 令
$$
A=\sigma_a\otimes I,\quad A'=\sigma_{a'}\otimes I,\quad B=I\otimes\sigma_b,\quad B'=I\otimes\sigma_{b'} ,
$$
$$
S_{\rm op}:=A\otimes B-A\otimes B'+A'\otimes B+A'\otimes B' ,
\qquad S=\operatorname{tr}(\rho S_{\rm op}) .
$$

**第二步：三个基本事实。**
$$
\text{(i)}\ A^2=A'^2=B^2=B'^2=I;\qquad
\text{(ii)}\ [A,B]=[A,B']=[A',B]=[A',B']=0;\qquad
\text{(iii)}\ \lVert[A,A']\rVert\le2,\ \lVert[B,B']\rVert\le2 .
$$
(i) 由 $\sigma_n^2=I$；(ii) 由**张量积两侧各自独立**；(iii) 因 $\lVert[A,A']\rVert=\lVert 2i\,(\mathbf a\times\mathbf a')\cdot\boldsymbol\sigma\rVert=2\lvert\mathbf a\times\mathbf a'\rvert\le2$。

**第三步：恒等式。** 由 (i)(ii) 逐项展开：
$$
S_{\rm op}^2=4I+[A,A']\otimes[B,B'] .
$$

**第四步：取范数。**
$$
\lVert S_{\rm op}^2\rVert\le4+\lVert[A,A']\rVert\cdot\lVert[B,B']\rVert\le4+2\cdot2=8
\ \Longrightarrow\ \lVert S_{\rm op}\rVert\le2\sqrt2 .
$$
而 $\lvert S\rvert=\lvert\operatorname{tr}(\rho S_{\rm op})\rvert\le\lVert S_{\rm op}\rVert$（因 $\rho\ge0$、$\operatorname{tr}\rho=1$）。故 $\lvert S\rvert\le2\sqrt2$。$\blacksquare$

$$
\boxed{\;\textbf{证明只用了三件事}:\ A^2=I,\ B^2=B'^2=I,\ [A,B^{(\prime)}]=0 .\;}
\tag{3.5}
$$

**核验 §A2′**：在最优角上逐项验证恒等式 $S_{\rm op}^2=4I+[A,A']\otimes[B,B']$（精确成立）；并在 60 组随机方向上验证 $\lVert S_{\rm op}\rVert^2\le8$（实测最大 $7.794$）。

**这正是"两体"的代数刻画**——记住这一点，§6 会回来用它。

### 3.2 达界：构造 ＋ 数值

**构造**：取单态
$$
\rho_{\rm sing}=\lvert\Psi^-\rangle\langle\Psi^-\rvert,\qquad
\lvert\Psi^-\rangle=\frac{1}{\sqrt2}\bigl(\lvert01\rangle-\lvert10\rangle\bigr),
$$
与共面方向 $a\perp a'$、$b$ 与 $a$ 成 $45^\circ$、$b'$ 与 $a'$ 成 $45^\circ$（即 $b,b'$ 也正交）。

$$
\Longrightarrow\quad \lvert S\rvert=2\sqrt2 .
\tag{3.6}
$$

**核验 §B**：$E(a,b)=E(a',b)=E(a',b')=-0.70710678$，$E(a,b')=+0.70710678$，故
$$
S=-2.8284271247,\qquad \lvert S\rvert=2\sqrt2=2.8284271247\ \checkmark .
$$
（符号为负是相位约定：本文验证 $E(a,b)=-\mathbf a\cdot\mathbf b$，故 $S=-2\sqrt2$。）

**核验 §C**：4000 次随机方向抽样，$\max\lvert S\rvert=2.7529\le2\sqrt2$ ✓。

$$
\boxed{\;\textbf{结论}:\ \text{本框架的代数}\ \bigoplus_hM_2(\mathbb C)\ \textbf{允许}\ \lvert S\rvert=2\sqrt2 .\;}
\tag{3.7}
$$

---

## §4 主结果三：否证"中心态 $\Rightarrow$ 因子化"

### 4.1 一个看似有理的论证（**它是错的**）

$$
\text{若}\ \rho\in\mathcal Z_L,\ \text{则}\ \operatorname{tr}(\rho AB)\overset{?}{=}\operatorname{tr}(\rho A)\operatorname{tr}(\rho B) .
\tag{4.1}
$$

**这个论证是错的。** 理由：$\rho\in\mathcal Z_L$ 只给 $\rho$ 与一切对易，**不给** $\rho=\rho_1\otimes\rho_2$ 型的因子化。

### 4.2 正确的关系

> **命题 19.** 对 $\rho=\bigoplus_h\rho_h$（分块对角）与 $A=\bigoplus_hA_h$、$B=\bigoplus_hB_h$，
> $$
> \langle AB\rangle=\sum_h\lambda_h\operatorname{tr}_h(A_hB_h),\qquad
> \langle A\rangle\langle B\rangle=\sum_{h,k}\lambda_h\lambda_k\operatorname{tr}_h(A_h)\operatorname{tr}_k(B_k) .
> \tag{4.2}
> $$
> 一般**不相等**。

*证明.* 直接计算分块迹。$\blacksquare$

**核验 §D1**：对账本态 $\rho_{\rm led}$（$w=(0.3,0.7)$）与 200 对随机 $A,B$，
$$
\max\bigl\lvert\langle AB\rangle-\langle A\rangle\langle B\rangle\bigr\rvert=1.6416\ \gg0 .
$$

$$
\boxed{\;\textbf{排除}:\ \text{"中心态}\Rightarrow\text{因子化"}\ \text{为}\textbf{假}.\;}
\tag{4.3}
$$

### 4.3 正确的图像：中心态是**经典混合**，不是**独立**

由 (4.2) 的**第一式**：$\langle AB\rangle$ 是**层上的加权平均**。这正是"经典混合"的形状：

$$
\boxed{\;\text{账本态}\ \rho_{\rm led}=\sum_h w_h\cdot\Bigl(\frac{I_2}{2}\ \text{on block}\ h\Bigr)
\ \text{是}\ \textbf{层上的凸组合}\ \text{（经典混合）}.\;}
\tag{4.4}
$$

**核验 §D3**：即使只取块标量观测量 $A=\operatorname{blk}(1,-1)$、$B=\operatorname{blk}(1,-1)$，
$$
\langle AB\rangle=1.0000,\qquad \langle A\rangle\langle B\rangle=0.1600 ,
$$
故**在中心上也不因子化**——因为中心态是**相关的**经典混合，不是独立的。

### 4.4 中心与分块对角：一处必须写清

**核验 §D2 的实测**：

| 陈述 | 真假 |
|:--|:--|
| 块标量算符与**分块对角**算符对易 | **真** |
| 块标量算符与账本态对易 | **真** |
| 块标量算符与**任意**算符对易 | **假**（含块间相干者不对易） |

$$
\Longrightarrow\quad \mathcal Z_L=\Bigl\{\bigoplus_hc_hI_2\Bigr\}\ \text{是中心，但}\ \textbf{不等于}\ \text{全体算符} .
\tag{4.5}
$$

---

## §5 主结果四：语境性的可达判据

### 5.1 交换情形：无语境性（完整证明）

> **定理 25（交换 $\Rightarrow$ 联合分布）.** 设 $\mathcal B\subseteq\mathcal A_L$ 是**交换** $*$-子代数。则 $\mathcal B$ 上的一切观测量有**联合分布**，故**不存在** Kochen–Specker 型矛盾。

*证明.* 交换有限维 $*$-代数由 Gelfand–Naimark 定理同构于 $\mathbb C(X)$，$X$ 为有限集（其特征标空间）。$\rho$ 限制到 $\mathcal B$ 给 $X$ 上的概率测度 $\mu$。取 $f_1,\dots,f_n\in\mathbb C(X)$，其联合分布即 $f_* \mu$。故一切观测量可同时被赋值，无矛盾。$\blacksquare$

$$
\boxed{\;\text{故}\ \mathcal Z_L\ \text{上的读数}\ \textbf{必然无语境性}\ \text{（因为}\ \mathcal Z_L\ \text{交换）}.\;}
\tag{5.1}
$$

### 5.2 非交换情形：语境性可能（完整证明）

> **定理 26（非交换 $\Rightarrow$ 存在性）.** $\mathcal A_L$ 中的单个块 $M_2(\mathbb C)$ 是**非交换**的，且其中**不存在**对一切方向赝自旋的**确定性赋值**。

*证明.* 反设存在 $\nu:M_2(\mathbb C)\to\mathbb C$ 使 $\nu(I)=1$、$\nu(A^*A)\ge0$，且 $\nu$ 对**一切**算符乘性（$\nu(AB)=\nu(A)\nu(B)$）。取 $\nu(\sigma_x),\nu(\sigma_y),\nu(\sigma_z)$。由 $\sigma_i^2=I$ 与 $\nu$ 乘性，$\nu(\sigma_i)^2=1$，故 $\nu(\sigma_i)=\pm1$；从而 $\sum_i\nu(\sigma_i)^2=3$。但由乘性与 $\sigma_x\sigma_y=i\sigma_z$：
$$
\nu(\sigma_x)\nu(\sigma_y)=i\,\nu(\sigma_z)
$$
左端为实数、右端为纯虚数——矛盾（除非 $\nu(\sigma_z)=0$，与 $\nu(\sigma_z)^2=1$ 冲突）。$\blacksquare$

**核验 §F1**：$[\sigma_x,\sigma_y]$ 的算子范数 $=2\ne0$（同方向为 $0$），**非交换性被数值确认**。

$$
\boxed{\;\text{故}\ \text{块内读数}\ \Rightarrow\ \text{非交换}\ \Rightarrow\ \textbf{语境性可能};\quad
\text{块标量读数}\ \Rightarrow\ \text{交换}\ \Rightarrow\ \textbf{无语境性}.\;}
\tag{5.2}
$$

### 5.3 判据的最终形式

> **判据.** 本框架内语境性**可达** $\iff$ 读出**分辨块内**（即读出 $\mathcal A_L$ 的非中心部分）。

$$
\boxed{\;\text{语境性的可达性}\ =\ \text{读出的分辨率};\ \text{这不是物理输入，是}\textbf{读出设计}\text{的选择}.\;}
\tag{5.3}
$$

**档位**：**条件证成**（前件是"读出能分辨块内"，这是**读出设计**的具名选择，不是公理推论）。

---

## §6 与"两体"的关系（本文留下的接续点）

定理 24 的证明（(3.5)）只用了三件事：$A^2=I$、$B^2=B'^2=I$、$[A,B^{(\prime)}]=0$。这三件事的**载体**是**两体结构** $A\otimes I$、$I\otimes B$。

$$
\Longrightarrow\quad \boxed{\;\textbf{关键问题}:\ \text{本框架内}\ A\otimes I\ \text{与}\ I\otimes B\ \text{是否存在？}\;}
\tag{6.1}
$$

若**不存在**，则 (3.7) 的"允许"只是**代数层面**的（$\bigoplus_hM_2$ 的某个块里可以放一个单态），而**不是"两个分离的观测者各自测量"**。

$$
\boxed{\;\text{这一问由本文的接续篇（多体与张量积）回答};\ \text{本文}\textbf{不假设}\text{其答案}.\;}
\tag{6.2}
$$

---

## §7 对缺口账的更新

| 审计原文 | 更新后的表述 | 依据 |
|:--|:--|:--|
| "纠缠与 CHSH 有硬结果" | **正确**，且本文**补上完整证明**（经典界、Tsirelson 上半、达界构造） | 定理 23、24 |
| "`Bell`／`EPR` 未系统检查" | **本文系统检查**：CHSH 两侧界都已证；"两体性"是否可达另见接续篇 | §2、§3、§6 |
| "`语境性` 未列入" | **本文给出判据**：可达性 $=$ 读出分辨率（条件证成） | 定理 25、26 |
| （新增） | **否证**"中心态 $\Rightarrow$ 因子化"（一个易犯的错误） | §4 |

---

## §8 结论的适用范围

1. **本框架尚未违反 Bell 不等式**——本文只证**代数允许** $2\sqrt2$，未证"两个分离观测者可实现之"（见 §6，接续篇回答）；
2. **账本态不能违反 CHSH**——账本态是中心元，落在经典混合一侧；
3. **"中心态 $\Rightarrow$ 因子化"已被本文否证**（§4）；
4. **本文的 Tsirelson 证明不是新的**——它是**标准证明**的完整复述，本文的价值在于**核对它在本框架的代数上是否成立**；
5. **语境性判据不是公理推论**——它的前件（读出分辨块内）是**具名输入**（§5.3）；
6. **$\mathcal Z_L$ 不等于全体算符**——**假**（核验 §D2b）；
7. **本文只处理 EPR 的关联函数与不等式这一支**——"实在性判据"一支未展开；
8. **本文全部在有限 $L$ 上**——结论不对 $\mathcal A_\infty$ 成立。

---

## §9 结论

$$
\boxed{\;
\begin{aligned}
&\textbf{定理 23}:\ \text{经典界}\ \lvert S\rvert\le2\ \text{（逐点证明）};\ \text{乘积态}\ \lvert S\rvert\le2\ \text{（独立证明）}.\\[1mm]
&\textbf{定理 24}:\ \text{Tsirelson 界}\ \lvert S\rvert\le2\sqrt2\ \text{（完整证明）};\ \text{单态}\ \textbf{达界}.\\[1mm]
&\textbf{排除}:\ \text{"中心态}\Rightarrow\text{因子化"}\ \textbf{为假};\ \text{中心态是}\textbf{经典混合}。\\[1mm]
&\textbf{定理 25／26}:\ \text{交换}\Rightarrow\text{无语境性};\ \text{非交换}\Rightarrow\text{语境性可能};\ \text{可达性}\ =\ \text{读出分辨率}.
\end{aligned}
}
$$

**一句话：**

$$
\boxed{\;\textbf{量子基础在本框架内不缺;缺的是"两个观测者"这个前提——而那是下一问。}\;}
$$

---

## 附录 A　核验

`Z035_qfound_check.py`：纯标准库、不 import、不写文件。

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | CHSH 经典界（穷举 $2^4$ 个确定性模型）＋ 乘积态 | 3 |
| B | 单态达 Tsirelson 界（最优角，精确 $2\sqrt2$） | 3 |
| C | 4000 次随机抽样不超界 | 1 |
| D | 否证"中心态 $\Rightarrow$ 因子化"；中心与分块对角的正确关系 | 5 |
| E | 块标量观测量两两对易（交换子代数） | 9 |
| F | 非交换性（$[\sigma_u,\sigma_v]$ 范数） | 4 |
| A2′ | Tsirelson 的算子恒等式 $S_{\rm op}^2=4I+[A,A']\otimes[B,B']$（精确 ＋ 60 组随机） | 2 |
| | **合计** | **26 / 0** |

**核验纪律**：数值核验只算旁证，不替代证明；本文的定理 23、24、25、26 都在正文给出**完整证明**，核验只作独立复算。

**一处容易写错的地方**："块标量算符与任意算符对易"**不成立**（与含块间相干者不对易，核验 §D2b）；正确的是"与分块对角对易、与账本态对易，但不与任意算符对易"。

---

*Zero 项目 · Z035 · 完*
