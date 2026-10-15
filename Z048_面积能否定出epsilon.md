# 面积能否定出 $\epsilon_h$？——**不能**，但面积确实**再分**了层

**Zero 项目 · 基础推导系列 · 第 Z048 号**

**篇性**：**推导 ＋ 否证篇**。按用户指示"查" Z047 §6.3 留下的线索：$\epsilon_h$ 是否**等于面积符号**？
**结论**：**线索不成立** —— 但查证过程给出一个**新的精确结构**（面积在层内的谱），以及一条**对 $L$ 的苛刻约束**。
**纪律**：**未改动任何公理**。
**依据**：正式篇 Z001（平衡词、层高）；本系列 Z042（$b_h$）、Z044（面积）、Z047（$\mathbb Z_2$、$\epsilon_h$）。
**核验**：[`Z048_area_sign_check.py`](_staging/Z048_area_sign_check.py)。

---

## 摘要

$$
\boxed{\;\textbf{线索判定}:\ \epsilon_h=\operatorname{sgn}(a)\ \textbf{不成立}\ \text{——因框架自定的}\ L=4\ \text{上有【零面积】词}.\;}
\tag{0.1}
$$

| # | 结果 | 内容 |
|--:|:--|:--|
| **一** | $\mathbb Z_2$ 配对保面积**大小** | $a(-w)=-a(w)\Rightarrow$ 层内面积谱关于 $0$ 对称 |
| **二** | **零面积词** | $L\equiv0\pmod4$ 时存在；$L\equiv2\pmod4$ 时**不存在** |
| **三** | **面积秩** $A_h$ | $L\equiv2\pmod4\Rightarrow A_h$ **恒偶**；$L\equiv0\pmod4\Rightarrow A_h$ **可奇** |
| **四** | **面积跨层** | 面积**不是**层内标签（同一面积跨多层） |
| **五** | 层内结构 | 面积把每层再分成 $A_h$ 个**子层**（对角算符） |
| **六** | **否定** | 框架自定 $L=4$ 给 $A_1=3$（奇）$\Rightarrow$ 面积符号**不干净** |

$$
\boxed{\;\textbf{净}:\ \text{面积确实再分层，但它给的是}\ A_h\ \text{个子层（很多）};\ \text{不是}\ 2\ \text{个}.\;}
\tag{0.2}
$$

---

## §1 面积的基本性质

### 1.1 $\mathbb Z_2$-奇

$$
a(w)=\sum_{k<L}s_k;\qquad a(-w)=\sum_{k<L}(-s_k)=-a(w) .
\tag{1.1}
$$

**核验 §A**：$L=8,10$ 上逐词成立。

$$
\Longrightarrow\quad \boxed{\;\text{层内的面积谱关于}\ 0\ \textbf{对称};\ \text{非零值成对出现}.\;}
\tag{1.2}
$$

### 1.2 零面积词

$$
a(w)=0\ \text{的词}:\quad \text{在}\ \mathbb Z_2\ \text{下是}\ \textbf{不动点} .
\tag{1.3}
$$

**核验 §B**：

| $L$ | $4$ | $6$ | $8$ | $10$ | $12$ | $14$ | $16$ | $18$ | $20$ |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| $L\bmod4$ | $0$ | $2$ | $0$ | $2$ | $0$ | $2$ | $0$ | $2$ | $0$ |
| 零面积词数 | $2$ | $\mathbf 0$ | $8$ | $\mathbf 0$ | $58$ | $\mathbf 0$ | $526$ | $\mathbf 0$ | $5448$ |

$$
\boxed{\;L\equiv2\pmod4\ \Longrightarrow\ \textbf{无零面积词};\qquad
L\equiv0\pmod4\ \Longrightarrow\ \textbf{有零面积词}.\;}
\tag{1.4}
$$

$$
\Longrightarrow\quad \boxed{\;\text{故面积符号在}\ L\equiv2\pmod4\ \text{上是}\textbf{干净的二分类};\ \text{在}\ L\equiv0\pmod4\ \text{上}\textbf{不是}.\;}
\tag{1.5}
$$

---

## §2 面积秩 $A_h$（层内的面积谱）

### 2.1 定义与结果

$$
A_h:=\#\{a(w):h(w)=h\}\qquad(\text{层}\ h\ \text{内的面积取值数}) .
\tag{2.1}
$$

**核验 §C**：

| $L$ | $h=1$ | $2$ | $3$ | $4$ | $5$ | $6$ | $7$ | $8$ |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| $4$ | $\mathbf 3$ | $2$ | | | | | | |
| $6$ | $4$ | $6$ | $2$ | | | | | |
| $8$ | $\mathbf 5$ | $\mathbf{11}$ | $8$ | $2$ | | | | |
| $10$ | $6$ | $14$ | $16$ | $10$ | $2$ | | | |
| $12$ | $\mathbf 7$ | $\mathbf{17}$ | $\mathbf{25}$ | $20$ | $12$ | $2$ | | |
| $14$ | $8$ | $20$ | $30$ | $32$ | $24$ | $14$ | $2$ | |
| $16$ | $\mathbf 9$ | $\mathbf{23}$ | $\mathbf{35}$ | $\mathbf{45}$ | $38$ | $28$ | $16$ | $2$ |

$$
\boxed{\;L\equiv2\pmod4\Rightarrow A_h\ \textbf{恒偶};\qquad
L\equiv0\pmod4\Rightarrow A_h\ \textbf{可奇}\ (\text{低层}).\;}
\tag{2.2}
$$

$$
\text{奇偶性的来源}:\ A_h\ \text{的奇偶}\ \longleftrightarrow\ \text{层}\ h\ \text{内是否有零面积词} .
\tag{2.3}
$$

### 2.2 框架自定 $L=4$ 的完整结构

| 层 | $b_h$ | 面积谱（值$\times$重数） | $A_h$ |
|--:|--:|:--|--:|
| $1$ | $4$ | $-2\times1,\ \mathbf 0\times2,\ +2\times1$ | $\mathbf 3$ |
| $2$ | $2$ | $-4\times1,\ +4\times1$ | $2$ |

$$
\boxed{\;\text{故在}\ L=4\ \text{上，}\ \epsilon_h=\operatorname{sgn}(a)\ \textbf{不能定义}\ \text{——层 1 有两条零面积词}.\;}
\tag{2.4}
$$

---

## §3 线索**不成立**（本文的主要否证）

### 3.1 线索的原样

$$
\text{Z047 §6.3 的线索}:\quad \epsilon_h\ \overset{?}{=}\ \operatorname{sgn}\bigl(a(w)\bigr) .
\tag{3.1}
$$

### 3.2 为什么它失败

$$
\boxed{\;\text{三条独立的失败理由}:}
\tag{3.2}
$$

**(i) 零面积词。** $L\equiv0\pmod4$ 时有 $a(w)=0$ 的词。
$$
\Longrightarrow\quad \operatorname{sgn}(0)\ \text{无定义};\quad \text{层 1 有 2 条这样的词（}L=4\text{）} .
\tag{3.3}
$$

**(ii) 面积跨层。** 面积**不是**层内标签。

**核验 §D**：$L=8$ 有 $17$ 个面积取值，其中 $9$ 个**跨多层**；$L=10$：$26$ 个中 $20$ 个跨层。

$$
\Longrightarrow\quad \text{面积是}\ \textbf{整个}\ W_L\ \text{上的算符};\ \text{要当层内结构，必须先}\ \textbf{条件化} \text{到层} .
\tag{3.4}
$$

**(iii) 面积给的是 $A_h$ 个类，不是 2 个。**

$$
A_h\ \text{一般远大于}\ 2\ (L=16,\ h=4:\ A_h=45) .
\tag{3.5}
$$

$$
\Longrightarrow\quad \text{面积}\ \textbf{再分} \text{了层，但分成}\ A_h\ \text{块，不是 2 块} .
\tag{3.6}
$$

### 3.3 判定

$$
\boxed{\;\textbf{否证}:\ \epsilon_h=\operatorname{sgn}(a)\ \text{在框架自定的}\ L=4\ \text{上}\textbf{不成立};\ \text{在}\ L\equiv2\pmod4\ \text{上才勉强成立}.\;}
\tag{3.7}
$$

$$
\Longrightarrow\quad \boxed{\;\text{层内输入}\ \{\epsilon_h\}\ \textbf{未能归零}.\;}
\tag{3.8}
$$

---

## §4 但面积确实给了**更有价值**的东西

### 4.1 面积把层再分成子层

$$
\text{层}\ h\ \longrightarrow\ \{h\}\times\{\text{面积值}\};\qquad \text{共}\ A_h\ \text{个子层} .
\tag{4.1}
$$

$$
\Longrightarrow\quad \text{层的总数}:\ H\ \longrightarrow\ \sum_{h=1}^{H}A_h .
\tag{4.2}
$$

**核验 §E**：$L=4$ 给 $A_1+A_2=3+2=5$；$L=8$ 给 $5+11+8+2=26$。

$$
\boxed{\;\text{这是一个}\textbf{二阶粗粒化}:\ \text{词}\to\text{层}\to\text{面积};\ \text{且每一级都是纯组合}.\;}
\tag{4.3}
$$

### 4.2 二阶粗粒化能给出什么

$$
\text{若用（层, 面积）当宏观类}:\quad
\omega_{(h,a)}=\frac{\#\{w:h(w)=h,\ a(w)=a\}}{\lvert W_L\rvert} .
\tag{4.4}
$$

$$
\Longrightarrow\quad K_{(h,a)}=-\log\omega_{(h,a)}=-\log\lvert W_L\rvert+\log\#\{w:h(w)=h,a(w)=a\} .
\tag{4.5}
$$

$$
\boxed{\;\text{故}\ K\ \text{的谱}:\ \text{层间（由}\ b_h\text{）}\ +\ \text{层内（由面积谱）};\ \textbf{两者都是纯组合}.\;}
\tag{4.6}
$$

### 4.3 这比 Z047 的 $\mathbb Z_2$ 结构更**细**

$$
\underbrace{\mathbb Z_2\ (w\mapsto-w)}_{\text{给 }2\text{ 类，}\epsilon_h}\qquad\text{vs}\qquad
\underbrace{\text{面积}}_{\text{给 }A_h\text{ 类}} .
\tag{4.7}
$$

$$
\Longrightarrow\quad \boxed{\;\text{面积的分类}\textbf{比}\ \mathbb Z_2\ \text{更细};\ \text{但"更细"}\ne\text{"更小"}——\text{它带来的输入更多}.\;}
\tag{4.8}
$$

**这是一处必须点明的反直觉**：更细的粗粒化**不一定**减少输入，因为每一类都要一个标签。

---

## §5 对 $L$ 的一条苛刻约束（**本文的附加收获**）

$$
\text{干净的二分类}\ \Longleftrightarrow\ L\equiv2\pmod4 .
\tag{5.1}
$$

$$
\text{框架自定}\ L=4\ \equiv\ 0\pmod4\ \Longrightarrow\ \textbf{不干净} .
\tag{5.2}
$$

$$
\boxed{\;\text{故"用面积符号定}\ \epsilon_h"\ \text{这条路}\textbf{要求}\ L\equiv2\pmod4;\ \text{与框架自定的}\ L=4\ \textbf{冲突}.\;}
\tag{5.3}
$$

**注意**：$L$ 不改。$L=4$ 是框架的结论（四维窗 ＋ 容量），**不因本文而改**。

$$
\Longrightarrow\quad \text{这条约束的价值在于}:\ \text{它}\textbf{再次}\text{把"缺什么"定位到}\ L\ \text{之外的输入} .
\tag{5.4}
$$

---

## §6 那么 $\{\epsilon_h\}$ 到底还需不需要输入

### 6.1 两条路的对照

| 路 | 层内的类数 | 输入 | 状态 |
|:--|:--|:--|:--|
| $\mathbb Z_2$（Z047） | $2$ | $\{\epsilon_h\}\in\{\pm1\}^H$ | **仍需输入**（本文确认） |
| 面积（本文） | $A_h$ | 每个子层一个标签 | **更贵** |

$$
\boxed{\;\text{两条路都}\textbf{没有}\text{把层内输入归零};\ \text{而}\ \mathbb Z_2\ \text{那条}\textbf{更便宜}.\;}
\tag{6.1}
$$

### 6.2 建议的取向

$$
\boxed{\;\text{保留 Z047 的}\ \mathbb Z_2\ \text{结构（}\epsilon_h\text{），}\textbf{放弃}\text{"}\epsilon_h=\text{面积符号"}\ \text{这条线索}.\;}
\tag{6.2}
$$

$$
\Longrightarrow\quad \text{当前状态}:\ \text{层内输入} = \{\epsilon_h\}\ (\textbf{二元符号},\ H\ \text{个}) .
\tag{6.3}
$$

$$
\text{相比 Z046 的}\ 2H\ \text{个连续参数}:\ \text{已是}\textbf{实质简化};\ \text{但}\textbf{未归零} .
\tag{6.4}
$$

---

## §7 一处必须报告的复核结论

Z047 §6.3 的线索是我自己写下的：

> 若把 $\epsilon_h$ 认同为"面积符号"，则 $\epsilon_h$ 与面积一一对应 $\Longrightarrow$ **又回到纯组合**？

**本文否证了它。** 两处错：

$$
\textbf{(a)}\ \text{我假设了面积符号是干净的二分类};\ \text{但}\ L\equiv0\pmod4\ \text{时有零面积词} .
\tag{7.1}
$$

$$
\textbf{(b)}\ \text{我假设了面积是层内标签};\ \text{但核验 §D 显示面积}\textbf{跨层} .
\tag{7.2}
$$

$$
\boxed{\;\text{§6.3 的线索}\ \textbf{不成立} \text{（被否证的提案）}.\;}
\tag{7.3}
$$

---

## §8 结论的适用范围

1. **面积符号不是"$\epsilon_h$ 的来源"**——**否证**（§3.3）；
2. **零面积词的个数没有闭式**——本文只数值列出（§1.2 表）；
3. **$A_h$ 没有闭式**——本文只数值列出（§2.1 表）；
4. **"$A_h$ 奇偶 $\leftrightarrow$ 零面积词存在"是观察**，不是定理（核验 §C 支持）；
5. **不应改 $L$**——$L=4$ 是框架的结论（§5 注）；
6. **面积分类并不"更好"**——它更细，但输入更多（§4.3）；
7. **本文只查了面积与 $\mathbb Z_2$**——未穷尽层内结构；
8. **$\{\epsilon_h\}$ 未必永远不能归零**——本文只否证了**一条**路径。

---

## §9 结论

$$
\boxed{\;
\begin{aligned}
&\textbf{线索判定}:\ \epsilon_h=\operatorname{sgn}(a)\ \textbf{不成立};\ \text{三条独立理由（零面积词／跨层／类数}\gg2\text{）}.\\[1mm]
&\textbf{零面积的规律}:\ L\equiv2\pmod4\ \text{无};\ L\equiv0\pmod4\ \text{有};\ \text{框架自定}\ L=4\ \text{属后者}.\\[1mm]
&\textbf{面积的真正作用}:\ \text{把层再分成}\ A_h\ \text{个子层（}\textbf{二阶粗粒化}）;\ \text{全部纯组合}.\\[1mm]
&\textbf{当前状态}:\ \text{层内输入仍是}\ \{\epsilon_h\}\ (\text{二元},\ H\ \text{个});\ \text{比 Z046 便宜，但}\textbf{未归零}.\\[1mm]
&\textbf{判定}:\ \text{§6.3 的线索}\textbf{不成立}.
\end{aligned}
}
$$

**一句话：**

$$
\boxed{\;\textbf{面积确实再分了层，但它分成很多块而不是两块——}\epsilon_h\ \text{还是要买。}\;}
$$

---

## 附录 A　核验

`Z048_area_sign_check.py`：纯标准库；不 import 项目脚本、不写文件。

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | 面积是 $\mathbb Z_2$-奇（$L=8,10$ 逐词） | 4 |
| B | 零面积词的规律（$L=4,\dots,20$；$L\bmod4$ 判定） | 14 |
| C | 面积秩 $A_h$ 的奇偶与零面积词的对应 | 33 |
| D | 面积**跨层**（不是层内标签） | 6 |
| E | 二阶粗粒化的子层总数 $\sum_hA_h$ | 5 |
| F | $L=4$ 的完整结构（$A_1=3$ 含零） | 5 |
| | **合计** | **67** |

**核验纪律**：数值只算旁证；本文的否证只用**反例**（$L=4$ 的零面积词）。

---

*推导稿 · 已落盘 · 未改动任何公理 · 完*
