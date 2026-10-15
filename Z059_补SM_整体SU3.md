# 补 SM（二）：整体 $SU(3)$ —— **可补，但代价是一条结构公设**

**Zero 项目 · 基础推导系列 · 第 Z059 号**

**篇性**：**推导 ＋ 定位篇**。用户指示："补第 5 项"（整体 $SU(3)$）。
**结论**：**可补** —— $\mathfrak{su}(3)$ 的 8 个生成元与李括号**一次得到**；但代价是承认**色空间 $\mathbb C^3$** 为新的基本对象（一条结构公设）。
**纪律**：**未改动任何公理**；本文落盘时未改动 Z001–Z036 任何一篇。
**依据**：正式篇 Z004（$\mathbb C$-代数）、Z021（Heawood 定理 8）；本系列 Z043（$M_2$ 被逼出）、Z057（夸克）、Z058（三色可区分）。
**核验**：[`Z059_su3_check.py`](_staging/Z059_su3_check.py)。

---

## 摘要

$$
\boxed{\;\textbf{（一）关键困难}:\ \mathbb Z_3\ \text{与}\ \mathrm{Aut}(\text{Fano})\ \text{都}\ \textbf{不给} \ 8\ \text{维李代数} \ (\text{阿贝尔／有限}).\;}
\tag{0.1}
$$

$$
\boxed{\;\textbf{（二）出路}:\ \text{一旦有}\ \mathbb C^3,\ \text{则}\ \operatorname{End}=\mathbb M_3(\mathbb C)\ \Longrightarrow\ \mathfrak{su}(3)\ \textbf{自动} \ (\text{无需额外结构}).\;}
\tag{0.2}
$$

$$
\boxed{\;\textbf{（三）代价}:\ \text{承认"色空间}\ \mathbb C^3\ "\ \text{为基本对象} = \text{一条}\ \textbf{结构公设} .\;}
\tag{0.3}
$$

$$
\boxed{\;\textbf{（四）剩余差别}:\ \text{此}\ \mathfrak{su}(3)\ \text{是}\ \textbf{整体} \text{的};\ \text{QCD 需}\ \textbf{局域} \text{（第 6 项）} .\;}
\tag{0.4}
$$

---

## §1 为什么难：两个候选都不给李代数

### 1.1 $\mathbb Z_3$ 是阿贝尔的

$$
\mathfrak{Lie}(\mathbb Z_3)=0\ \Longrightarrow\ \text{无升降算符，无李括号} .
\tag{1.1}
$$

$$
\mathbb Z_3\ \text{只给}\ 3\ \text{个特征标（Z058）} \Longrightarrow \textbf{3 个可区分的态}，\text{不给}\ 8\ \text{个生成元} .
\tag{1.2}
$$

### 1.2 $\mathrm{Aut}(\text{Fano})=PSL(2,7)$ 是**有限**群

$$
\lvert PSL(2,7)\rvert=168=2^3\cdot3\cdot7\ \Longrightarrow\ \text{李代数}\ \textbf{0 维} .
\tag{1.3}
$$

**核验 §A**：阶分布 $1^{(1)},2^{(21)},3^{(56)},4^{(42)},7^{(48)}$。

$$
9\nmid168\ \Longrightarrow\ \text{（拉格朗日）}\ \textbf{不含}\ \mathbb Z_3\times\mathbb Z_3\ \Longrightarrow\ \text{rank-2 Cartan}\ \textbf{缺} .
\tag{1.4}
$$

### 1.3 但 $PSL(2,7)$ **有** 3 维不可约表示

$$
\text{不可约维数}:\ 1,\ \mathbf 3,\ \mathbf 3,\ 6,\ 7,\ 8 .
\tag{1.5}
$$

$$
\Longrightarrow\quad \textbf{色三重态可以承载};\ \text{但它们}\ \textbf{不生成} \ \mathfrak{su}(3)\ \text{（有限群的像是有限的）} .
\tag{1.6}
$$

---

## §2 出路：$\mathfrak{su}(3)$ 是 $\mathbb C^3$ 的**自动**推论

### 2.1 一步到位

$$
\text{设}\ V=\mathbb C^3\ (\text{色空间}) \Longrightarrow
\operatorname{End}(V)=\mathbb M_3(\mathbb C),\ \dim_{\mathbb C}=9 .
\tag{2.1}
$$

$$
\mathfrak{su}(3)=\{X\in\mathbb M_3(\mathbb C):X^\dagger=-X,\ \operatorname{tr}X=0\}
\ \Longrightarrow\
\dim_{\mathbb R}\mathfrak{su}(3)=18-1-\underbrace{9}_{X^\dagger=-X\ \text{的实维数}}=8 .
\tag{2.2}
$$

**核验 §B**：$\dim_{\mathbb R}\{X:X^\dagger=-X\}=9$；减 $\operatorname{tr}=0$ 给 $8$ ✓。

$$
\boxed{\;\text{故}\ 8\ \text{个生成元 ＋ 李括号}\ \textbf{一次得到};\ \text{不需要}\ \mathbb Z_3\times\mathbb Z_3 .\;}
\tag{2.3}
$$

### 2.2 为什么这**解决了** §1.2 的困难

$$
\text{§1.2 的困难是"有限群不给李代数"};\ \text{但它}\ \textbf{不需要} \text{给} .
\tag{2.4}
$$

$$
\text{李代数来自}\ \operatorname{End}(V)\ \text{（自同态代数）};\ \text{群只给}\ \textbf{表示} .
\tag{2.5}
$$

$$
\boxed{\;\text{这是一个}\textbf{标准的数学事实}:\ \mathfrak{su}(n)=\mathfrak{sl}(n)\cap\mathfrak u(n);\ \text{只要有}\ \mathbb C^n\ \text{就有它} .\;}
\tag{2.6}
$$

### 2.3 显式基（Gell-Mann）

$$
\mathfrak{su}(3)\ \text{的基}:\ \lambda_1,\dots,\lambda_8;\qquad
[\lambda_a,\lambda_b]=2if_{abc}\lambda_c .
\tag{2.7}
$$

$$
\text{Cartan}\ \mathfrak h=\langle\lambda_3,\lambda_8\rangle\ (\text{rank}\ 2);\qquad
\underbrace{3\ \text{对}}_{\lambda_{1,2},\lambda_{4,5},\lambda_{6,7}}\ \text{升降算符}=6 .
\tag{2.8}
$$

$$
\Longrightarrow\quad 2+6=8\ \checkmark .
\tag{2.9}
$$

---

## §3 色三重态从哪里来：Fano 线

### 3.1 Fano 平面的关联性质（**本文核验**）

$$
7\ \text{点},\ 7\ \text{线},\ \text{每线 3 点},\ \text{每点 3 线} .
\tag{3.1}
$$

$$
\textbf{核验 §C}:\quad \text{任两线}\ \textbf{恰交一点}\ (0\ \text{反例});\qquad
\text{任两点}\ \textbf{恰在一条线上}\ (0\ \text{反例}) .
\tag{3.2}
$$

$$
\Longrightarrow\quad \text{取一条线}\ L_0\ \Longrightarrow\ L_0\cong\{\text{3 点}\}\ \Longrightarrow\ \mathbb C^3 .
\tag{3.3}
$$

### 3.2 与 $L0$ 的接口：$M_2\hookrightarrow M_3$

$$
\text{Zero 现有的矩阵因子是}\ \mathbb M_2(\mathbb C)\ (\text{每层},\ \text{Z043}) .
\tag{3.4}
$$

$$
\text{固定一个点}\ p\in L_0;\ \text{过}\ p\ \text{有 3 条线}\ L_1,L_2,L_3 .
\tag{3.5}
$$

$$
\Longrightarrow\quad \mathbb M_3(\mathbb C)\ \text{作用在}\ \{L_1,L_2,L_3\}\ \text{上};\quad
\mathbb M_2(\mathbb C)\hookrightarrow\mathbb M_3(\mathbb C)\ (\text{左上块}) .
\tag{3.6}
$$

$$
\boxed{\;\text{故}\ L0\ \text{的}\ M_2\ \textbf{可以嵌入} \ \text{色代数}\ M_3;\ \text{且}\ 7\ \text{个点给}\ 7\ \text{个}\ \textbf{共轭} \text{的嵌入} .\;}
\tag{3.7}
$$

---

## §4 代价与收益（**必须写清**）

### 4.1 收益

| # | 得到 | 来源 |
|--:|:--|:--|
| $1$ | **8 个生成元** | $\operatorname{End}(\mathbb C^3)$（§2.1） |
| $2$ | **李括号** | 矩阵交换子（§2.7） |
| $3$ | **rank-2 Cartan** | $\lambda_3,\lambda_8$（§2.8） |
| $4$ | **色三重态** $\mathbf 3$ | Fano 线 3 点（§3.1） |
| $5$ | **色单态** | $\epsilon_{ijk}$（Z057） |
| $6$ | **三色可区分** | $\mathbb Z_3$ 特征标（Z058） |

### 4.2 代价

$$
\boxed{\;\textbf{代价}:\ \text{承认}\ \mathbb C^3\ (\text{色空间})\ \text{是}\ \textbf{新的基本对象} .\;}
\tag{4.1}
$$

$$
\text{它}\ \textbf{不是} \text{现有公理①②③的推论};\ \text{是一条}\ \textbf{结构公设} .
\tag{4.2}
$$

$$
\Longrightarrow\quad \text{记作}\ \mathrm{C1}'\ (\text{与 Z039 的}\ \mathrm{C1}\ \text{不同}):\quad
\boxed{\;\text{每条 Fano 线}\ L\ \text{张成一个}\ \mathbb C^3\ \text{色空间} .\;}
\tag{4.3}
$$

### 4.3 为什么**不是**"内部矛盾"

$$
\text{Zero 已有的矩阵结构是}\ \mathbb M_2\ (\text{每层}) ;\ \text{本构造加}\ \mathbb M_3 .
\tag{4.4}
$$

$$
\mathbb M_2\hookrightarrow\mathbb M_3\ (\text{§3.2}) \Longrightarrow \textbf{相容} \text{（不是替换）} .
\tag{4.5}
$$

$$
\boxed{\;\text{故}\ \mathrm{C1}'\ \text{是}\ \textbf{扩充} \text{，不是修改};\ \text{与 Z039 的"改代数型"是}\textbf{两件事} .\;}
\tag{4.6}
$$

---

## §5 剩下的差别：整体 vs 局域

$$
\text{本构造给}\ \mathfrak{su}(3)\ \text{作用在}\ \textbf{一个} \ \mathbb C^3\ \text{上} \Longrightarrow \textbf{整体} \text{对称性} .
\tag{5.1}
$$

$$
\text{QCD 需}\ A_\mu^a(x)\ \text{带}\ x \Longrightarrow \textbf{局域} .
\tag{5.2}
$$

$$
\boxed{\;\text{局域化需要}\ \textbf{位置};\ \text{而位置}\ \textbf{不可导出} \ (\text{D14}=\text{Z028}=\text{Z012}) .\;}
\tag{5.3}
$$

$$
\Longrightarrow\quad \text{第 5 项}\ \textbf{可补};\ \text{第 6 项}\ \text{是}\ \textbf{同一个老缺口} .
\tag{5.4}
$$

---

## §6 对"算质子"的意义（**必须诚实**）

$$
\text{即使}\ \mathfrak{su}(3)\ \text{是整体的，也可能足够}:\ \text{禁闭}\ \text{发生在}\ \textbf{动力学} \text{，不在规范性} .
\tag{6.1}
$$

$$
\Longrightarrow\quad \text{要算质子，真正需要的仍是}\ \text{Z056 的}\ \Lambda_{\rm QCD}\ (\text{非幂律跑动}) .
\tag{6.2}
$$

$$
\boxed{\;\text{故本文补了}\ \textbf{代数} \text{，没补}\ \textbf{动力学};\ \text{Z056 的判定}\ \textbf{不变} .\;}
\tag{6.3}
$$

$$
\textbf{但一处改善}:\ \text{现在有}\ 8\ \text{个生成元可用来算}\ \textbf{一圈系数} \ b_0 .
\tag{6.4}
$$

$$
b_0=\frac{11N_c-2N_f}{12\pi};\qquad N_c=3\ \Longrightarrow\
b_0=\frac{33-2N_f}{12\pi} .
\tag{6.5}
$$

$$
\boxed{\;\text{故}\ \mathrm{C1}'\ \text{使}\ b_0\ \textbf{可写};\ \text{但}\ \text{Z056 §4.2 的"RGT 非自治"仍挡住普适}\ \beta .\;}
\tag{6.6}
$$

---

## §7 结论

$$
\boxed{\;
\begin{aligned}
&\textbf{难在哪}:\ \mathbb Z_3\ \text{阿贝尔};\ PSL(2,7)\ \text{有限};\ \text{两者都不给 8 维李代数}.\\[1mm]
&\textbf{出路}:\ \text{取}\ \mathbb C^3\ \Longrightarrow\ \mathfrak{su}(3)\ \textbf{自动} \ (\operatorname{End}\ \text{的迹零反厄米部分}).\\[1mm]
&\textbf{来源}:\ \text{Fano 线 3 点}\ \Longrightarrow\ \mathbb C^3;\ \text{且}\ M_2\hookrightarrow M_3\ \text{相容}.\\[1mm]
&\textbf{代价}:\ \mathrm{C1}'\ \text{``每条 Fano 线张成一个}\ \mathbb C^3\text{''}\ \text{——一条结构公设}.\\[1mm]
&\textbf{剩余}:\ \text{整体 vs 局域} = \text{"位置"缺口（第 6 项）}.
\end{aligned}
}
$$

**一句话：**

$$
\boxed{\;\textbf{8 个生成元不用找——只要有}\ \mathbb C^3\ \text{它就在那儿；要找的是承认}\ \mathbb C^3\ \text{的理由。}\;}
$$

---

## §8 不许主张

1. **不主张** $\mathrm{C1}'$ 可从公理①②③**推出**——**不能**（§4.2）；
2. **不主张** Fano 线"**就是**"色空间——它是**本文的构造**；
3. **不主张** $M_2\hookrightarrow M_3$ 的嵌入是**唯一**的——$7$ 个点给 $7$ 个共轭嵌入（§3.2）；
4. **不主张**整体 $SU(3)$ "够用"——它是否足以禁闭**未查**（§6.1）；
5. **不主张** $b_0$ 因此可算——Z056 §4.2 的非自治**仍**挡住（§6.6）；
6. **不主张**本文改了 Z056／Z057／Z058 的其余结论；
7. **不主张** $\mathbb M_3$ 是 Zero **已有**的——Zero 现有的是 $\mathbb M_2$（每层）；
8. **不主张**第 6 项（局域化）因此变容易——它是**同一个**老缺口。

---

## 附录 A　核验

`Z059_su3_check.py`：纯标准库（＋`numpy`）；不 import 项目脚本、不写文件。

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | $\mathbb Z_3$ 阿贝尔；$PSL(2,7)$ 有限、阶分布、$9\nmid168$ | 8 |
| B | $\dim\mathfrak{su}(3)=9-1=8$；Gell-Mann 计数 $2+6$ | 6 |
| C | Fano 关联：任两线恰交一点、任两点恰一线 | 6 |
| D | $M_2\hookrightarrow M_3$ 的维数核算 | 4 |
| E | $b_0$ 的写法（$N_c=3$，$N_f=3,4,5,6$） | 6 |
| | **合计** | **30** |

**核验纪律**：数值只算旁证；**本文的核心是一步标准数学**（$\mathfrak{su}(n)$ 是 $\operatorname{End}$ 的迹零反厄米部分）。

---

*推导稿 · 已落盘 · 未改动任何公理 · 完*
