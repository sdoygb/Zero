# 夸克在 Zero 里是什么？——**结构、电荷与根本差别**

**Zero 项目 · 基础推导系列 · 第 Z057 号**

**篇性**：**根源推导篇**。用户指示："从根源查起，质子由三个夸克组成，这三个夸克在 Zero 里是什么东西？结构呢？电荷呢？等等性质在 Zero 里怎样表达？"
**结论**：**夸克在 Zero 里可以表达，而且"三个"有天然来源** —— 但有一处**关键差别**（$\mathbb Z_3$ 是**置换**而非**相乘**）。
**纪律**：**未改动任何公理**。
**依据**：正式篇 Z001（公理、荷分布）、Z010／Z014（味结构缺口）、Z021（Heawood 定理 8）；本系列 Z041（禁闭直觉）、Z056（标度）。
**核验**：[`Z057_quark_in_zero_check.py`](_staging/Z057_quark_in_zero_check.py)。

---

## 摘要

$$
\boxed{\;\textbf{（一）正式篇从未定义夸克}:\ \text{"夸克"仅出现 2 次，都在讲"夸克扇区"缺口} .\;}
\tag{0.1}
$$

$$
\boxed{\;\textbf{（二）Zero 里唯一的"3"是}\textbf{空间维数}\ d=3;\ \text{来源}:\ \text{Heawood 的}\ \deg=3 .\;}
\tag{0.2}
$$

$$
\boxed{\;\textbf{（三）唯一干净的 3 元组}:\ \text{Fano 平面的每条线恰含 3 点};\ \text{7 条线} .\;}
\tag{0.3}
$$

$$
\boxed{\;\textbf{（四）三色不可区分}:\ \mathbb Z_3\ \text{在线上}\ \textbf{传递};\ \text{但这是}\ \textbf{置换}，\text{不是 QCD 的}\ \textbf{相乘中心} .\;}
\tag{0.4}
$$

$$
\boxed{\;\textbf{（五）电荷}:\ \text{整数}\ x\ \text{除以线长 3};\ \text{色单态}\iff x\equiv0\pmod3\iff\text{总荷为整数} .\;}
\tag{0.5}
$$

---

## §1 第一步：正式篇里没有夸克

**核验 §A**：$Z001$–$Z036$ 中"夸克"只出现 **3 次**（$Z010$ 一次、$Z014$ 两次）：

$$
Z010\ §3.10:\ \text{"}10\ (\text{夸克一般点})\text{"};\qquad
Z014\ \text{定理 3}:\ \text{"夸克一般点}\ 18-8=36-26=10\text{"} .
\tag{1.1}
$$

$$
\Longrightarrow\quad \textbf{Zero 从未定义过夸克};\ \text{这是}\ \textbf{从零推} .
\tag{1.2}
$$

### 1.1 $L0$ 有什么（Z001 的清单）

| 对象 | 内容 |
|:--|:--|
| 通道 | $x\in\mathbb Z^C$，$Q=\sum_ix_i$，$H_Q=\{Q=0\}$ |
| 词 | $w\in\{\pm1\}^L$，$\sum_kw_k=0$ |
| 高度轮廓 | $s_k=\sum_{i\le k}w_i$ |
| 层高 | $h(w)=\max_k\lvert s_k\rvert$ |
| 面积 | $a(w)=\sum_{k<L}s_k$ |
| 层 | $h=1,\dots,H=L/2=2$ |

$$
\boxed{\;\textbf{关键}:\ \text{Zero 的电荷是}\ \textbf{整数} \text{且}\ \textbf{总量恒为零} \Longrightarrow\ \text{单荷不存在（Z041）} .\;}
\tag{1.3}
$$

---

## §2 "三"在 Zero 里的来源

### 2.1 候选清单（逐个查）

| 候选 | 值 | 是三？ |
|:--|--:|:--:|
| $\lvert C\rvert$ | $14$ | ✗（$14=2\cdot7$） |
| $H=L/2$ | $2$ | ✗ |
| $L$ | $4$ | ✗ |
| $w_k$ | $\pm1$ | ✗（二值） |
| $b_h$ | $4,2$ | ✗ |
| $\beta$ 个数 | $2$ | ✗ |
| **Heawood 的度** | $\mathbf 3$ | **✓** |
| **Fano 线长** | $\mathbf 3$ | **✓** |
| **完美匹配数** | $24$ | ✗（Z021 写的 $3$ 是**着色数**的误记，见 §2.3） |

$$
\Longrightarrow\quad \boxed{\;\text{Zero 里唯一的"3"是}\ \textbf{空间维数} \text{（Heawood 的度）与}\ \textbf{Fano 线长} .\;}
\tag{2.1}
$$

### 2.2 Heawood 图（Z021 定理 8，**本文复核**）

$$
\lvert V\rvert=14,\quad\deg=3,\quad\text{围长}=6,\quad\text{二部},\quad\lvert E\rvert=21 .
\tag{2.2}
$$

$$
\text{构造}:\ \text{Heawood}=\text{Fano 平面的}\ \textbf{关联图} .
\tag{2.3}
$$

**核验 §B** 的邻接谱：

$$
\operatorname{spec}=\{3^{(1)},\ (\sqrt2)^{(6)},\ (-\sqrt2)^{(6)},\ (-3)^{(1)}\},\qquad
1+6+6+1=14\ \checkmark
\tag{2.4}
$$

**这是 $\deg=3$ 的直接后果**（$3$ 与 $-3$ 是两端的单重本征值）。

### 2.3 一处对 Z021 的更正

$$
\text{Z021 §8.3 写}:\ \text{"完美匹配数}\ \mathbf 3=\tbinom D2=\text{轴数}" .
\tag{2.5}
$$

**核验 §C**：

$$
\textbf{完美匹配数}=24;\qquad \textbf{3-边着色数}=8\ (\text{无序}) .
\tag{2.6}
$$

$$
\Longrightarrow\quad \boxed{\;\textbf{更正}:\ \text{"完美匹配数}=3"\ \textbf{不成立};\ \text{正确是}\ 24（\text{着色}\ 8） .\;}
\tag{2.7}
$$

$$
\textbf{但 Z021 的用途}\ \textbf{仍然成立}:\ \text{3-正则二部图}\ \textbf{必有} \text{3-边着色（König）} .
\tag{2.8}
$$

---

## §3 三个夸克：结构

### 3.1 唯一干净的 3 元组

$$
14/3\notin\mathbb Z\ \Longrightarrow\ \text{顶点层面}\ \textbf{无} \text{干净的三分类} .
\tag{3.1}
$$

$$
\text{但 Fano 平面的每条线恰含}\ \mathbf 3\ \text{点};\qquad \text{共}\ 7\ \text{条线} .
\tag{3.2}
$$

$$
\Longrightarrow\quad \boxed{\;\text{唯一干净的 3 元组} = \{\text{Fano 线}\};\ \text{共 7 个候选"三夸克态"} .\;}
\tag{3.3}
$$

### 3.2 三色的**不可区分性**（**核心结果**）

**核验 §D**：取一条线 $L_0=\{0,1,3\}$。存在**阶 3 的保线自同构** $8$ 个，且它们在 $L_0$ 上是**传递的 3-循环**。

$$
\text{例}:\ \sigma=(1,3,2,0,5,6,4)\ \Longrightarrow\ \sigma\big|_{L_0}:\ 0\to1\to3\to0 .
\tag{3.4}
$$

$$
\Longrightarrow\quad \mathbb Z_3\ \text{在}\ L_0\ \text{上}\ \textbf{传递} \Longrightarrow \text{三色}\ \textbf{一个轨道} \Longrightarrow \textbf{不可区分} .
\tag{3.5}
$$

$$
\boxed{\;\text{这是一个}\textbf{真正的三色结构} \text{（不是"三个恰好凑一起的东西"）} .\;}
\tag{3.6}
$$

**一般事实**：每个阶 3 自同构在 7 点上必有不动点（因 $7=2\cdot3+1$），轨道形如 $1+3+3$。这与传递性**不矛盾**：不动点可以**在线外**。

### 3.3 色单态

$$
\boxed{\;\epsilon_{ijk}\ \text{在 3-循环下}\ \textbf{不变} \Longrightarrow \text{它是一个色单态} .\;}
\tag{3.7}
$$

**核验 §E**：$\epsilon$ 有 $6$ 个非零元；在 3-循环 $(1,2,0)$ 下 $\epsilon'=\epsilon$ ✓。

$$
\textbf{n 个夸克的全反对称部分维数}:\ \binom{3}{n}=1,3,\boxed{1},0\ (n=1,2,3,4) .
\tag{3.8}
$$

$$
\Longrightarrow\quad \boxed{\;\textbf{最小色单态}=3\ \text{个夸克};\ \text{2 个}\ \textbf{不行} \text{（与 Z041 独立一致）} .\;}
\tag{3.9}
$$

---

## §4 电荷

### 4.1 Zero 的原生电荷

$$
x\in\mathbb Z^C\ (\textbf{整数});\qquad Q=\sum_ix_i;\qquad H_Q=\{Q=0\}\ (\textbf{零和}) .
\tag{4.1}
$$

$$
\Longrightarrow\quad \text{原生电荷}\ \textbf{整数} \text{且}\ \textbf{总量恒零};\ \text{夸克的}\ \pm\tfrac23,\mp\tfrac13\ \textbf{不在} \text{其中} .
\tag{4.2}
$$

### 4.2 唯一自然的分数化路径

$$
\boxed{\;\text{夸克电荷}:=\frac{x}{3},\qquad x\in\mathbb Z,\qquad 3=\text{Fano 线长} .\;}
\tag{4.3}
$$

**核验 §F**：

| 强子 | $x$ 组合 | 电荷 $=\sum x/3$ | 实验值 |
|:--|:--|--:|--:|
| 质子 $uud$ | $2+2-1$ | $\mathbf 1$ | $1$ ✓ |
| 中子 $udd$ | $2-1-1$ | $0$ | $0$ ✓ |
| $\Delta^{++}$ $uuu$ | $2+2+2$ | $2$ | $2$ ✓ |
| $\Omega^-$ 型 $ddd$ | $-1-1-1$ | $-1$ | $-1$ ✓ |

$$
\Longrightarrow\quad \text{取}\ u:x=2,\ d:x=-1\ \text{即得全部强子电荷} .
\tag{4.4}
$$

### 4.3 **关键约束**：色单态 $\iff$ 总荷为整数

$$
\text{色单态}\iff\sum_{\text{3 quarks}}x_i\equiv0\pmod3
\iff\text{总荷}=\frac{\sum x_i}{3}\in\mathbb Z .
\tag{4.5}
$$

$$
\boxed{\;\textbf{即}:\ \text{色单态条件与"电荷整数"是}\textbf{同一个} \ \mathbb Z_3\ \text{条件} .\;}
\tag{4.6}
$$

**核验 §F**：$uud:3\equiv0$；$udd:0\equiv0$；$uuu:6\equiv0$；$ddd:-3\equiv0$。

### 4.4 但**一处必须诚实**

$$
\boxed{\;\text{具体的}\ x\ \text{取值（}2\ \text{与}\ -1\text{）}\ \textbf{Zero 没有给出};\ \text{它需要一张}\ \textbf{荷表} .\;}
\tag{4.7}
$$

$$
\Longrightarrow\quad \text{Z010／Z014 已把它登记为}\ \textbf{味结构缺口} \text{（"全部数值不可导出"）} .
\tag{4.8}
$$

$$
\textbf{能导出什么}:\ \text{电荷的}\ \textbf{量子化} \text{（以}\ 1/3\ \text{为单位）与}\ \textbf{整数性条件} .
\tag{4.9}
$$

$$
\textbf{不能导出什么}:\ \text{哪些味对应哪些}\ x .
\tag{4.10}
$$

---

## §5 **根本差别**：$\mathbb Z_3$ 是**置换**，不是**相乘**

### 5.1 对照

| | QCD | Zero（本文构造） |
|:--|:--|:--|
| 三色的载体 | $\mathbf 3$（$SU(3)$ 基本表示） | Fano 线的 3 点 |
| $\mathbb Z_3$ 的作用 | **相乘**（中心荷，相因子 $1,\omega,\omega^2$） | **置换**（3-循环） |
| 三色 | **可区分**（triality 不同） | **不可区分**（一个轨道） |
| 胶子 | $8$ 个（伴随） | **无**（只有 $\mathbb Z_3$ ＋ 度 3） |
| 代数 | $\mathfrak{su}(3)$ | 无（只有 $\mathbb Z_3$） |

$$
\boxed{\;\textbf{关键差别}:\ \text{QCD 的}\ \mathbb Z_3\ \text{带}\ \textbf{相位};\ \text{Zero 的}\ \mathbb Z_3\ \text{只}\ \textbf{换位置} .\;}
\tag{5.1}
$$

### 5.2 一处必须点明的数字巧合

$$
\text{Heawood 的圈秩}=\lvert E\rvert-\lvert V\rvert+1=21-14+1=8=\dim\mathfrak{su}(3) .
\tag{5.2}
$$

$$
\boxed{\;\textbf{这只登记为巧合};\ \text{Heawood 里}\ \textbf{没有}\ \mathfrak{su}(3)\ \text{结构（核验 §G）} .\;}
\tag{5.3}
$$

$$
\text{理由}:\ \text{圈秩}\ 8\ \text{是}\ \textbf{图的秩};\ \mathfrak{su}(3)\ \text{需要}\ 8\ \text{个}\ \textbf{生成元} \text{（带李括号）} .
\tag{5.4}
$$

### 5.3 这决定了一件重要的事

$$
\textbf{不可区分} \Longrightarrow \text{三个夸克是}\ \textbf{全同粒子};\ \text{波函数}\ \textbf{必须} \text{全反对称} .
\tag{5.5}
$$

$$
\Longrightarrow\quad \text{这}\ \textbf{正是} \ \epsilon_{ijk}\ \text{单态（§3.3）};\ \text{且说明为什么}\ \textbf{恰好} \text{需要反对称化} .
\tag{5.6}
$$

$$
\boxed{\;\text{故"三色不可区分"}\ \textbf{不是缺陷};\ \text{它}\ \textbf{解释} \text{了色单态的反对称性} .\;}
\tag{5.7}
$$

$$
\textbf{但代价}:\ \text{它}\ \textbf{给不出} \text{胶子的八重态};\ \text{故}\ \textbf{给不出} \text{禁闭动力学} .
\tag{5.8}
$$

---

## §6 其余性质（逐项判定）

| 性质 | Zero 里能否表达 | 依据 |
|:--|:--|:--|
| **个数 3** | **能** | Heawood 度 / Fano 线长（§2） |
| **色（三色）** | **能**（但不可区分） | $\mathbb Z_3$ 传递（§3.2） |
| **色单态** | **能** | $\epsilon_{ijk}$（§3.3） |
| **电荷（以 $1/3$ 为单位）** | **能** | $x/3$，$x\in\mathbb Z$（§4.2） |
| **电荷整数性条件** | **能** | $\sum x\equiv0\pmod3$（§4.3） |
| **具体荷值（味表）** | **不能** | Z010／Z014 缺口（§4.4） |
| **胶子（8）** | **不能** | 只有 $\mathbb Z_3$，无 $\mathfrak{su}(3)$（§5.2） |
| **禁闭（动力学）** | **不能** | Z041：只有选择定则版 |
| **质量** | **不能** | Z054／Z056 |
| **自旋 $\tfrac12$** | **部分** | 原生不可约表示维数只有 $1,2$（Z013）$\Rightarrow$ 有 2 维槽 |
| **统计（费米）** | **部分** | Z006 有自旋统计，但需 2 维槽的支持 |
| **位置** | **不能** | D14 |

$$
\boxed{\;\textbf{故"夸克"在 Zero 里}\ \textbf{有一半} \text{是有结构的};\ \text{另一半（胶子／禁闭／质量）没有} .\;}
\tag{6.1}
$$

---

## §7 对"算质子"的意义

$$
\textbf{本文的贡献}:\ \text{把"三夸克"从}\ \textbf{买入} \text{降为}\ \textbf{半导出} .
\tag{7.1}
$$

$$
\underbrace{\text{个数 3}}_{\textbf{已导出}};\quad
\underbrace{\text{三色 + 不可区分}}_{\textbf{已导出}};\quad
\underbrace{\text{色单态}}_{\textbf{已导出}};\quad
\underbrace{\text{电荷量子化}}_{\textbf{已导出}} .
\tag{7.2}
$$

$$
\textbf{但仍买}:\ \text{味表（}x\ \text{的具体值）};\ \text{胶子};\ \text{禁闭动力学} .
\tag{7.3}
$$

$$
\boxed{\;\textbf{故 Z056 的结论}\ \textbf{不变}:\ \text{质子质量仍需}\ \Lambda_{\rm QCD};\ \text{而}\ \Lambda\ \text{需}\ \textbf{非幂律的跑动} .\;}
\tag{7.4}
$$

$$
\textbf{但本文使}\ \text{Z041 的"要买入定向"}\ \textbf{更精确}:
\tag{7.5}
$$

$$
\text{要买入的不是"三"},\ \text{而是}\ \textbf{"胶子"} \text{（}\mathfrak{su}(3)\text{）} .
\tag{7.6}
$$

---

## §8 不许主张

1. **不主张**"夸克就是 Fano 线"——这是**本文的构造**，不是 Zero 的定理；
2. **不主张** $x=2,-1$ 可导出——**不能**（§4.4）；
3. **不主张** Heawood 里有 $\mathfrak{su}(3)$——**没有**（§5.2）；
4. **不主张** $\mathbb Z_3$ 置换与 QCD 中心荷"等价"——**不等价**（§5.1）；
5. **不主张**"完美匹配数 $=3$"——**更正为 24**（§2.3）；
6. **不主张**自旋／统计已完整——只到"有 2 维槽"（§6）；
7. **不主张**本文使质子可算——**不可算**（§7.4）；
8. **不主张** Fano 线是**唯一**可能的载体——本文只查了它（14 顶点不可三分类）。

---

## §9 结论

$$
\boxed{\;
\begin{aligned}
&\textbf{根源}:\ \text{正式篇从未定义夸克};\ \text{本文从}\ L0\ \text{推} .\\[1mm]
&\textbf{"三"的来源}:\ \text{Heawood 的度}=3\ (=\text{空间维数});\ \text{Fano 线长}=3\ (\text{唯一干净 3 元组}) .\\[1mm]
&\textbf{结构}:\ \text{三色}\ \textbf{不可区分} \text{（}\mathbb Z_3\ \text{传递）}\Longrightarrow\epsilon_{ijk}\ \text{单态};\ \text{2 夸克不行} .\\[1mm]
&\textbf{电荷}:\ \text{整数}\ x\ \text{除以线长 3};\ \text{色单态}\iff\sum x\equiv0\pmod3\iff\text{整数荷} .\\[1mm]
&\textbf{根本差别}:\ \text{Zero 的}\ \mathbb Z_3\ \text{是}\textbf{置换};\ \text{QCD 的是}\textbf{相乘} \Longrightarrow\ \text{无胶子}.
\end{aligned}
}
$$

**一句话：**

$$
\boxed{\;\textbf{"三"和"色"Zero 都给得出；给不出的是那八个胶子。}\;}
$$

---

## 附录 A　核验

`Z057_quark_in_zero_check.py`：纯标准库；不 import 项目脚本、不写文件。

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | 正式篇里"夸克"只出现 3 次 | 3 |
| B | Heawood 图：$14/21/3$-正则/二部；邻接谱 | 8 |
| C | **更正**：完美匹配 24，着色 8（Z021 写 3） | 4 |
| D | $\mathbb Z_3$ 在 Fano 线上传递（8 个阶 3 自同构） | 6 |
| E | $\epsilon_{ijk}$ 单态；$\binom3n=1,3,1,0$ | 7 |
| F | 电荷 $x/3$ 与 $x\equiv0\pmod3$ | 8 |
| G | 圈秩 $8$ 只是巧合（无 $\mathfrak{su}(3)$） | 3 |
| | **合计** | **40** |

**核验纪律**：数值只算旁证；**本文更正了 Z021 一处误记**（完美匹配数）。

---

*推导稿 · 已落盘 · 未改动任何公理 · 完*
