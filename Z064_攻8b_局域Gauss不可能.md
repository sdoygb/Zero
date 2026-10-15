# 攻 $8b$：加"局域 Gauss 定律" —— **不可能成功**（一条定理）

**Zero 项目 · 基础推导系列 · 第 Z064 号**

**篇性**：**否证篇**。用户指示："先攻 $8b$"（加局域 Gauss 定律以打开禁闭）。
**结论**：**表述本身不可能成功** —— 在连通图上"逐顶点零和"与"全局零和"**横截**，合起来**只剩平凡解**。这不是禁闭，是**真空**。
**纪律**：**未改动任何公理**。
**依据**：正式篇 Z024（引理 1：$x(\tau)-x(0)=Bf$）；本系列 Z041（线索）、Z060（局域群）、Z061（质量隙）、Z063（九条判据）。
**核验**：[`Z064_local_gauss_check.py`](_staging/Z064_local_gauss_check.py)。

---

## 摘要

$$
\boxed{\;\textbf{定理 28}:\ \text{在}\ \textbf{连通} \text{图上},\quad
\ker B^{\mathsf T}=\operatorname{span}\{(1,\dots,1)\};\quad
\ker B^{\mathsf T}\cap H_Q=\{0\} .\;}
\tag{0.1}
$$

$$
\boxed{\;\textbf{推论}:\ \text{"局域 Gauss} ＋ \text{全局零和"}\ \textbf{什么都不许};\ \text{不是禁闭} .\;}
\tag{0.2}
$$

$$
\boxed{\;\textbf{且}:\ \text{Zero 的零和}\ \textbf{就是} \text{"荷由流决定"（}x=Bf\text{）};\ \text{它}\ \textbf{自动} \text{给全局中性} .\;}
\tag{0.3}
$$

$$
\boxed{\;\textbf{故}\ 8b\ \textbf{关闭};\ \text{但}\ 8a\ \text{（质量隙）}\ \textbf{仍然成立} .\;}
\tag{0.4}
$$

---

## §1 把"局域 Gauss 定律"写**精确**

### 1.1 两个空间

$$
C^0=\mathbb R^{14}\ (\textbf{顶点}/\text{荷});\qquad C^1=\mathbb R^{21}\ (\textbf{边}/\text{流}) .
\tag{1.1}
$$

$$
B:C^1\to C^0\ (\textbf{散度});\qquad B^{\mathsf T}:C^0\to C^1\ (\textbf{梯度}) .
\tag{1.2}
$$

**核验 §A**：

$$
\operatorname{rank}B=\mathbf{13};\qquad \dim\ker B=21-13=\mathbf 8\ (=\text{圈秩}) .
\tag{1.3}
$$

$$
\operatorname{rank}B^{\mathsf T}=\mathbf{13};\qquad \dim\ker B^{\mathsf T}=14-13=\mathbf 1 .
\tag{1.4}
$$

### 1.2 两条"零和"的**精确形式**

$$
\textbf{（甲）全局零和}:\quad H_Q=\Bigl\{x:\sum_ix_i=0\Bigr\};\qquad \dim H_Q=13 .
\tag{1.5}
$$

$$
\textbf{（乙）局域 Gauss}:\quad G=\{x:B^{\mathsf T}x=0\};\qquad \dim G=1 .
\tag{1.6}
$$

$$
\textbf{核验 §A}:\ G=\operatorname{span}\{(1,1,\dots,1)\}\ (\text{常数向量}) .
\tag{1.7}
$$

**理由**：$B^{\mathsf T}x=0$ 即"每条边上两端点的荷差为零"；图**连通** $\Longrightarrow$ 全部 $x_i$ 相等。$\blacksquare$

---

## §2 定理 28：两者**横截**

$$
\boxed{\;\ker B^{\mathsf T}\cap H_Q=\{0\} .\;}
\tag{2.1}
$$

*证明.* 设 $x\in\ker B^{\mathsf T}$，则 $x=q(1,\dots,1)$（式 1.7）。再设 $x\in H_Q$，则
$$
0=\sum_ix_i=q\cdot14\ \Longrightarrow\ q=0\ \Longrightarrow\ x=0 .
$$
$\blacksquare$

### 2.1 所以不是"更强"，是"**正交方向**"

$$
(1,\dots,1)\notin H_Q\ (\text{其和为}\ 14\ne0)\ \Longrightarrow\ G\ \textbf{与}\ H_Q\ \textbf{横截} .
\tag{2.2}
$$

$$
\boxed{\;\textbf{故}\ G+H_Q=C^0\ (\text{整个空间})\ \text{而}\ G\cap H_Q=\{0\}\ \Longrightarrow\ \text{合起来}\ \textbf{什么都不约束} .\;}
\tag{2.3}
$$

### 2.2 逐条物理判定

| 组态 | 满足全局零和？ | 满足局域 Gauss？ | 加了 $G$ 后 |
|:--|:--:|:--:|:--|
| 单荷 $+q$（一点） | ✗ | ✗ | 禁 |
| 荷-反荷 $(+q,-q)$ | ✓ | ✗ | **禁** |
| 真空 $x=0$ | ✓ | ✓ | 允许 |

$$
\boxed{\;\text{加}\ G\ \text{会禁掉}\ \textbf{一切} \text{非零组态} \Longrightarrow\ \text{这不是禁闭，是}\ \textbf{真空} .\;}
\tag{2.4}
$$

---

## §3 正确的做法：让规范场携带通量 —— **而 Zero 已经这么做了**

### 3.1 真正的 Gauss 定律

$$
\textbf{（真）}:\quad (\operatorname{div}E)(v)=\rho(v)\quad\Longleftrightarrow\quad B\,a=x\quad(a=\text{边上的场}) .
\tag{3.1}
$$

$$
\Longrightarrow\quad \sum_ix_i=\sum_v(Ba)_v=0\ \textbf{自动} .
\tag{3.2}
$$

### 3.2 Zero 的零和**就是**这一条

$$
\textbf{Z024 引理 1}:\quad x(\tau)-x(0)=Bf .
\tag{3.3}
$$

$$
\boxed{\;\text{故 Zero 的零和} = \text{"荷由流决定"}\ \text{的陈述};\ \text{它}\ \textbf{自动} \text{给全局中性} .\;}
\tag{3.4}
$$

$$
\Longrightarrow\quad \text{Zero}\ \textbf{已经} \text{有正确的 Gauss 定律};\ \text{它}\ \textbf{不是} \text{"缺"这一条} .
\tag{3.5}
$$

$$
\textbf{而它给的是}:\ \text{全局中性（净荷恒零）};\quad \textbf{不给}:\ \text{弦张力／线性势} .
\tag{3.6}
$$

---

## §4 为什么"加 $G$"必然失败（**一般论证**）

$$
\text{规范不变性要求}:\ \text{物理态} = \text{规范不变子空间} .
\tag{4.1}
$$

$$
\textbf{（甲）}:\ \text{全局中性} = \text{"净荷 }=0\text{"} \Longrightarrow \textbf{规范不变} \checkmark .
\tag{4.2}
$$

$$
\textbf{（乙）}:\ \text{逐顶点零} = \text{"每处荷 }=0\text{"} \Longrightarrow \textbf{不是} \text{规范不变} .
\tag{4.3}
$$

$$
\boxed{\;\textbf{故（乙）连规范不变性都不满足};\ \text{它是一个}\ \textbf{过强} \text{的外源要求} .\;}
\tag{4.4}
$$

$$
\Longrightarrow\quad \text{规范理论里}\ \textbf{从不} \text{要求"每点荷为零"（那会禁掉一切源）} .
\tag{4.5}
$$

---

## §5 那么 $8b$ 的真正内容是什么

### 5.1 拆开 $8b$

$$
8b:\quad \text{"无带色渐近态"} = \text{"谱中只有规范单态"} .
\tag{5.1}
$$

$$
\textbf{（甲）运动学部分}:\ \text{规范不变性} \Longrightarrow \text{物理 Hilbert 空间}\ \mathcal H_{\rm phys} .
\tag{5.2}
$$

$$
\textbf{（乙）动力学部分}:\ \text{弦张力}\ \sigma>0 \Longrightarrow \text{线性势} .
\tag{5.3}
$$

### 5.2 逐部分判定

$$
\textbf{（甲）}\ \text{在 Zero 里}\ \textbf{可得}:\ \text{规范不变态即"净荷为零"的态（由零和保证）} .
\tag{5.4}
$$

$$
\textbf{（乙）}\ \text{需}\ \sigma\ \text{从小格点到大格点的}\ \textbf{标度行为} \Longrightarrow\ \textbf{需热力学极限} .
\tag{5.5}
$$

$$
\boxed{\;\textbf{故}\ 8b\ \text{的（甲）可得、（乙）不可得};\ \text{与}\ 8c\ \textbf{同一障碍} .\;}
\tag{5.6}
$$

### 5.3 一条具体的检查

$$
\text{在}\ \textbf{有限} \text{格点上，静态势}\ V(R)\ \textbf{不} \text{线性增长};\ \text{它}\ \textbf{饱和} \text{（因有质量隙）} .
\tag{5.7}
$$

$$
\text{故即使做}\ \mathbb Z_3\ \text{规范理论的格点计算，}\ V(R)\ \text{也}\ \textbf{给不出}\ \sigma .
\tag{5.8}
$$

$$
\Longrightarrow\quad \text{与 Z061 的}\ P^k\to1/7\ne0\ \textbf{是同一件事} .
\tag{5.9}
$$

---

## §6 对第 8 项的最终判定（**更新**）

| 子问题 | Z063 的状态 | **本文之后的最终状态** |
|:--|:--|:--|
| $8a$ 质量隙 $>0$ | **已得**（$m=0.7520386984$） | **已得** |
| $8b$ 无带色渐近态 | 需局域 Gauss | **（甲）可得；（乙）不可得** |
| $8c$ 面积律／弦张力 | 关闭 | **关闭** |

$$
\boxed{\;\text{第 8 项的净结果}:\ 8a\ \textbf{已得};\ 8b\ \textbf{拆开后（乙）与}\ 8c\ \textbf{同一障碍};\ 8c\ \textbf{关闭} .\;}
\tag{6.1}
$$

$$
\Longrightarrow\quad \text{三条门}\ \textbf{全部检查过};\ \text{没有一条能打开弦张力} .
\tag{6.2}
$$

---

## §7 剩下的路（**诚实**）

| # | 路 | 代价 | 现状 |
|--:|:--|:--|:--|
| $1$ | 加局域 Gauss | **不可能**（定理 28） | **关闭** |
| $2$ | 用 $\mathbb Z_3$ 中心对称 | 需 Polyakov 环；且它**也不需**热力学极限 | 未算 |
| $3$ | 松开 $L=4$ | 放弃"读出四维"（Z062） | 关闭 |
| $4$ | **在 $L=4$ 上另找质量来源** | 未知 | **未查** |
| $5$ | 承认弦张力不可得 | 无 | **当前结论** |

$$
\boxed{\;\text{路 1 与 3 关闭};\ \text{路 2 未算};\ \text{路 4 未查};\ \text{路 5 是结论} .\;}
\tag{7.1}
$$

### 7.1 一处**必须澄清**的事（关于质子质量）

$$
\textbf{质子质量}\ m_p\approx(3\text{–}4)\Lambda_{\rm QCD};\qquad
\Lambda_{\rm QCD}\ \text{来自}\ \textbf{量纲嬗变} .
\tag{7.2}
$$

$$
\Longrightarrow\quad \text{算质子真正缺的是}\ \textbf{标度生成} \text{（Z056），}\ \textbf{不是} \text{禁闭本身} .
\tag{7.3}
$$

$$
\boxed{\;\text{故}\ 8b\ \text{与}\ 8c\ \text{关闭}\ \textbf{不直接} \text{决定质子质量};\ \text{Z056 的障碍（线性}\ \beta\text{）才是} .\;}
\tag{7.4}
$$

$$
\textbf{这是一个重要的重新定位}:\ \text{质子质量的障碍}\ \textbf{不} \text{在禁闭侧} .
\tag{7.5}
$$

---

## §8 不许主张

1. **不主张**"加局域 Gauss 能打开禁闭"——**不能**（定理 28）；
2. **不主张** Zero "缺 Gauss 定律"——它**已经**有（$x=Bf$，§3.2）；
3. **不主张** $G$ 与 $H_Q$ 有包含关系——它们**横截**（§2.1）；
4. **不主张** $8b$ 的（甲）部分"不可得"——它**可得**（§5.2）；
5. **不主张**本文算了 Polyakov 环判据——**未算**（路 2）；
6. **不主张**"在 $L=4$ 上另找质量来源"不存在——**未查**（路 4）；
7. **不主张**质子质量的障碍"只"在量纲嬗变——本文只说它**不直接**由禁闭决定；
8. **不主张**定理 28 是**新**的——它是线性代数的初等事实；新的是它的**物理判定**。

---

## §9 结论

$$
\boxed{\;
\begin{aligned}
&\textbf{定理 28}:\ \ker B^{\mathsf T}=\operatorname{span}\{(1,\dots,1)\};\quad \ker B^{\mathsf T}\cap H_Q=\{0\}.\\[1mm]
&\textbf{故}:\ \text{"局域 Gauss}＋\text{全局零和"}\ \textbf{什么都不许};\ \text{不是禁闭} .\\[1mm]
&\textbf{且}:\ \text{Zero 的}\ x=Bf\ \textbf{就是} \text{正确的 Gauss 定律（荷由流决定）} .\\[1mm]
&\textbf{重定位}:\ \text{质子质量的障碍}\ \textbf{不} \text{在禁闭侧，而在}\ \textbf{标度生成} \text{（Z056）}.\\[1mm]
&\textbf{净}:\ 8a\ \text{已得};\ 8b\ \text{拆开后与}\ 8c\ \text{同障碍};\ \text{三路检查完，无一条打开弦张力}.
\end{aligned}
}
$$

**一句话：**

$$
\boxed{\;\textbf{"每点荷为零"连规范不变性都不满足——它不是更强的禁闭，是真空。}\;}
$$

---

## 附录 A　核验

`Z064_local_gauss_check.py`：纯标准库（＋`numpy`）；不 import 项目脚本、不写文件。

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | $\operatorname{rank}B=13$；$\dim\ker B=8$；$\dim\ker B^{\mathsf T}=1$；$\ker B^{\mathsf T}=\operatorname{span}\{(1,\dots,1)\}$ | 8 |
| B | 定理 28：$\ker B^{\mathsf T}\cap H_Q=\{0\}$；$G+H_Q=C^0$ | 6 |
| C | 三个组态的逐条判定 | 6 |
| D | $x=Bf$ 是 Zero 的零和；$\sum(Bf)=0$ 自动 | 5 |
| E | 静态势在有限格点饱和（与 Z061 一致） | 4 |
| | **合计** | **29** |

**核验纪律**：数值只算旁证；定理 28 由**连通性 + 线性代数**给出。

---

*推导稿 · 已落盘 · 未改动任何公理 · 完*
