# E31：从 Zero 构造"场" —— 一个**成立一半**的构造

**篇性**：**构造 ＋ 限定篇（暂存讨论稿）**。用户的构造提案（照录）：
> "我有个想法，场是一种势，带方向的，那我可不可这么从 Zero 构建出场：刚开始零是乱动的，但随着这个动分化出一些闭合的动，一些比较有规律的动，然后很多一些同一方向的动和有规律的动，就形成了场？"

**结论**：**(一)** 用户的四步里，**前三步 Zero 都有**；**(二)** 第四步**部分成立** —— Zero 的层结构**自发地**诱导出一个**局部、有方向、非平凡**的场；**(三)** 但它止步于**阿贝尔、离散、无量纲**，到不了非阿贝尔规范场；**(四)** 本文更正了我自己此前一处**过严**的判定。
**纪律**：**未改动任何公理**。
**依据**：正式篇 Z001（公理③）、Z024（引理 1）、Z026（无热力学极限）、Z027／Z031（定向被禁）；本系列 E20、E24、E28；旧库 `U15`（往零走）。
**核验**：[`E31_field_from_zero_check.py`](E31_field_from_zero_check.py)。

---

## 摘要

$$
\boxed{\;\textbf{（一）}\ \text{用户的四步里}\ 1,2,3\ \textbf{成立};\ 4\ \text{部分成立} .\;}
\tag{0.1}
$$

$$
\boxed{\;\textbf{（二）}\ \text{Zero 的层结构}\ \textbf{自发地} \text{给出一个局部方向场}\ p_\uparrow(h,k)\ne\tfrac12 .\;}
\tag{0.2}
$$

$$
\boxed{\;\textbf{（三）}\ \text{它}\ \textbf{不需要} \text{定向输入（用无向计数）}\Longrightarrow\ \text{墙"定向被禁"在此}\textbf{不适用} .\;}
\tag{0.3}
$$

$$
\boxed{\;\textbf{（四）}\ \text{但它}\ \textbf{阿贝尔、离散、无量纲} \Longrightarrow\ \text{到不了}\ SU(3)\ \text{规范场} .\;}
\tag{0.4}
$$

---

## §1 用户的构造：四步对照

| 步 | 用户的说法 | Zero 里对应 | 判定 |
|--:|:--|:--|:--|
| $1$ | 零开始乱动 | 词 $w\in\{\pm1\}^L$（公理③：一切移动被实例化） | **✓** |
| $2$ | 分化出闭合的动 | 零和 $\sum_kw_k=0$（游走回到起点） | **✓** |
| $3$ | 有规律的动 | 层高 $h(w)$、面积 $a(w)$、旋转类（可粗粒化） | **✓** |
| $4$ | 许多同方向的动形成场 | ← 本文要查的 | **部分** |

---

## §2 一处**过严**的判定（本文的自我更正）

### 2.1 我曾判定

$$
\textbf{（我此前的判定）}:\quad \text{"方向"需要}\ \textbf{定向} \text{输入，而定向}\ \textbf{被禁} \text{（Z027／Z031）} \Longrightarrow 4\ \text{不成立} .
\tag{2.1}
$$

### 2.2 为什么过严

$$
\textbf{有向流}:\ f_e\in\mathbb Z\ (\text{带正负})\ \Longrightarrow\ \text{需先给每条边选一个}\ \textbf{定向} .
\tag{2.2}
$$

$$
\textbf{无向计数}:\ a_e\ge0\ (\text{走过几次})\ \Longrightarrow\ \textbf{不需要} \text{定向} .
\tag{2.3}
$$

$$
\boxed{\;\Longrightarrow\quad \text{用}\ a_e\ \text{就绕过了"定向被禁"这面墙} .\;}
\tag{2.4}
$$

### 2.3 而"局部三重结构"是**自带的**

$$
\text{平衡游走在每个位置有}\ 3\ \text{种走向}:\ \uparrow,\ \to\ (\text{同层}),\ \downarrow .
\tag{2.5}
$$

$$
\Longrightarrow\quad \text{局部方向信息天然是}\ \textbf{3 分量的} .
\tag{2.6}
$$

---

## §3 载体的修正：不是图顶点，是 $(h,k)$ 格点

$$
\textbf{关键}:\ \text{"高度"是}\ \textbf{词} \text{的属性，}\textbf{不是} \text{图的属性} .
\tag{3.1}
$$

$$
\text{Heawood 的顶点没有原生"高度"} \Longrightarrow F\ \text{不能活在图上} .
\tag{3.2}
$$

$$
\boxed{\;\text{载体}=\Bigl\{(h,k):1\le h\le H,\ 0\le k\le h\Bigr\};\qquad \lvert\text{载体}\rvert=\frac{H(H+3)}2 .\;}
\tag{3.3}
$$

**核验 §A**：

| $L$ | $4$ | $8$ | $16$ |
|:--|--:|--:|--:|
| $H=L/2$ | $2$ | $4$ | $8$ |
| 载体点数 | $\mathbf 5$ | $14$ | $44$ |

---

## §4 场的定义与 $S_3$ 分解

### 4.1 定义

$$
F(h,k):=\bigl(p_\uparrow,\ p_\to,\ p_\downarrow\bigr)(h,k)\quad(\text{3 条出边的权重}) .
\tag{4.1}
$$

$$
p_\uparrow(h,k):=\Pr\bigl[w_j=+1\ \big|\ h(w)=h,\ s_j=k\bigr] .
\tag{4.2}
$$

### 4.2 $S_3$ 分解（**核验 §B**）

$$
F\ \text{是}\ S_3\ \text{的}\ \textbf{3 维置换表示};\qquad
\chi_F(\sigma)=\#\operatorname{Fix}(\sigma) .
\tag{4.3}
$$

$$
\langle\chi_F,\mathbf 1\rangle=1;\qquad
\langle\chi_F,\operatorname{sgn}\rangle=0;\qquad
\langle\chi_F,\mathbf 2\rangle=1 .
\tag{4.4}
$$

$$
\boxed{\;F\ =\ \mathbf 1\ \oplus\ \mathbf 2\qquad(3=1+2) .\;}
\tag{4.5}
$$

$$
\Longrightarrow\quad \text{方向信息的}\ \textbf{不可约内容}\ \text{是}\ \mathbf 1\ (\text{强度})\ \text{与}\ \mathbf 2\ (\text{方向}) .
\tag{4.6}
$$

$$
\boxed{\;\text{而}\ \dim\mathbf 2=2=\operatorname{rank}\mathfrak{su}(3)\ \text{（与 E23 的 Cartan 同维数）} .\;}
\tag{4.7}
$$

---

## §5 决定性检验：那个场是**非平凡**的吗

### 5.1 方法

$$
\text{把每一步的方向按}\ (h,k)\ \text{条件化，看}\ p_\uparrow .
\tag{5.1}
$$

### 5.2 结果（**核验 §C**，$L=16$）

| 层 $h$ | 高度 $k$ | 上 | 下 | $p_\uparrow$ | $\lvert p_\uparrow-\tfrac12\rvert$ |
|--:|--:|--:|--:|--:|--:|
| $2$ | $1$ | $5103$ | $11369$ | $0.3098$ | $\mathbf{0.1902}$ |
| $3$ | $1$ | $6321$ | $9335$ | $0.4037$ | $0.0963$ |
| $3$ | $2$ | $3776$ | $6321$ | $0.3740$ | $0.1260$ |
| $4$ | $1$ | $2601$ | $3497$ | $0.4265$ | $0.0735$ |
| $4$ | $2$ | $2274$ | $2601$ | $0.4665$ | $0.0335$ |
| $4$ | $3$ | $1700$ | $2274$ | $0.4278$ | $0.0722$ |
| $5$ | $2$ | $681$ | $712$ | $0.4889$ | $0.0111$ |
| $5$ | $4$ | $544$ | $663$ | $0.4507$ | $0.0493$ |
| $6$ | $1$ | $138$ | $168$ | $0.4510$ | $0.0490$ |
| $6$ | $3$ | $136$ | $136$ | $0.5000$ | $\mathbf{0.0000}$ |
| $6$ | $5$ | $119$ | $135$ | $0.4685$ | $0.0315$ |

$$
\textbf{核验 §C 的统计}:\quad \text{内部点}\ 15\ \text{个};\quad
\lvert p_\uparrow-\tfrac12\rvert>0.001\ \text{者}\ \mathbf{14/15};\quad
>0.01\ \text{者}\ 11/15;\quad \text{最大偏离}\ \mathbf{0.1902} .
\tag{5.2a}
$$

$$
\boxed{\;p_\uparrow\ne\tfrac12\ \text{对绝大多数}\ (h,k) \Longrightarrow \text{方向分布}\ \textbf{非平凡} .\;}
\tag{5.2}
$$

### 5.3 两个**结构性**特征

$$
\textbf{(a) 壁效应}:\quad k=h\ \Longrightarrow\ p_\uparrow=\mathbf 0\ (\text{只能向下}) .
\tag{5.3}
$$

$$
\textbf{(b) 内部偏置}:\quad p_\uparrow<\tfrac12\ \text{（有向下的偏置）} .
\tag{5.4}
$$

$$
\boxed{\;\text{零和下的}\ \textbf{最大层约束} \text{（}h(w)=\max_k\lvert s_k\rvert\text{）}\ \textbf{自发地} \text{产生了方向偏置} .\;}
\tag{5.5}
$$

### 5.4 与公理③不冲突

$$
\text{偏置来自}\ \textbf{度量} \text{（低层可完成数多），}\textbf{不是} \text{动力学偏好} .
\tag{5.6}
$$

$$
\Longrightarrow\quad \text{公理③（无偏好）}\ \textbf{未被违反} .
\tag{5.7}
$$

$$
\textbf{偏离的依赖（新）}:\quad
\begin{cases}
\text{层越小（约束越紧）}\Longrightarrow\text{偏离越大}: & \lvert\text{dev}(2,1)\rvert=0.1902>\lvert\text{dev}(6,1)\rvert=0.0490;\\
\text{唯一完全对称点}: & (h,k)=(6,3),\ p_\uparrow=\tfrac12\ \text{（恰在中线）} .
\end{cases}
\tag{5.7a}
$$

$$
\textbf{与旧库对照}:\ \text{`U15` 已算过同一效应（向零／离零}\approx2.047\text{）};\ \text{但那时裁定"不是力、时间对称"} .
\tag{5.8}
$$

$$
\Longrightarrow\quad \textbf{本文的新点}:\ \text{同一个效应}\ \textbf{可以} \text{读成一个局部方向场的值} .
\tag{5.9}
$$

---

## §6 一处**有意思**的副产物：壁处的自发破缺

$$
k=h\ \text{时}\ p_\uparrow=0 \Longrightarrow \text{3 条出边里有一条被}\ \textbf{完全禁掉} .
\tag{6.1}
$$

$$
\Longrightarrow\quad \text{局部}\ S_3\ \text{对称性破缺到}\ S_2\ (\text{两条向下边的交换}) .
\tag{6.2}
$$

$$
\boxed{\;\text{一个}\ \textbf{自发破缺} \text{的局部对称性}\ \Longrightarrow\ \text{有 Goldstone 型结构的可能} .\;}
\tag{6.3}
$$

$$
\textbf{但}:\ \text{这是}\ \textbf{线索};\ \text{本文}\ \textbf{未查证} \text{（Z013 的 Higgs 代数核可能接上）} .
\tag{6.4}
$$

---

## §7 这个场**到不了**哪里（**必须写清**）

| 想要的（QCD） | 本文构造给出的 |
|:--|:--|
| 向量场 $A_\mu^a(x)$ | 标量型漂移分布 $p_\uparrow(h,k)$ |
| 群 $SU(3)$（8 生成元，连续） | $S_3$（6 元素，**离散**） |
| 非阿贝尔 | **阿贝尔**（漂移是实数，可交换） |
| 自相互作用 $f^{abc}$ | **无** |
| 渐近自由 | **无** |
| 线性势／弦张力 | **无**（面积有界，E25／E29） |
| 热力学极限 | **无**（Z026） |

$$
\boxed{\;\text{得到的是}\ \textbf{阿贝尔、离散、标量型} \text{的场};\ \text{不是非阿贝尔规范场} .\;}
\tag{7.1}
$$

$$
\textbf{（载体的规模）}:\ L=4\ \Longrightarrow\ \lvert\text{载体}\rvert=5 .
\tag{7.2}
$$

$$
\Longrightarrow\quad \text{"许多动形成场"里的}\ \textbf{"许多"}\ \text{在 Zero 里是}\ O(1) .
\tag{7.3}
$$

---

## §8 对本轮讨论的净贡献

### 8.1 回到"约束能不能产力"

$$
\textbf{本文的答案}:\quad \text{约束}\ \textbf{确实} \text{产生一个}\ \textbf{方向场};\quad \text{但它}\ \textbf{不产生力} .
\tag{8.1}
$$

$$
\textbf{差别}:\quad
\text{力}\ \text{需要}\ \underbrace{\text{势的梯度}}_{\text{Zero 有（}E=f^{\mathsf T}Lf\text{）}}\ +\ \underbrace{\text{时间演化律}}_{\text{Zero}\ \textbf{无}（\text{Z005／U12}）} .
\tag{8.2}
$$

$$
\boxed{\;\text{约束给了"方向"，没有给"力"};\ \text{缺的那一半仍是}\ \textbf{演化律} .\;}
\tag{8.3}
$$

### 8.2 对"定向墙"的修正

$$
\text{Z027／Z031}:\ \text{四道排除归约为}\ \textbf{定向};\ \text{定向是被禁的输入} .
\tag{8.4}
$$

$$
\textbf{本文}:\ \text{在}\ \textbf{这个用途上} \text{不需要定向（用无向计数）} .
\tag{8.5}
$$

$$
\Longrightarrow\quad \boxed{\;\text{定向墙仍然存在，但}\ \textbf{不挡} \text{"从计数构造方向场"这条路} .\;}
\tag{8.6}
$$

---

## §9 不许主张

1. **不主张**本文"构造出了规范场"——**没有**（§7）；
2. **不主张** $p_\uparrow(h,k)$ 是"力"——它是**静态漂移分布**（§8.1）；
3. **不主张**壁处的 $S_3\to S_2$ 破缺已与 Higgs 接上——**未查证**（§6.4）；
4. **不主张**本文的载体选择是**唯一**的——本文取 $(h,k)$ 是因为"高度"只在词上有定义（§3）；
5. **不主张** $S_3$ 的 $\mathbf 2$ 与 $\mathfrak{su}(3)$ 的 Cartan "是同一个东西"——只是**同维数**（§4.7）；
6. **不主张**本文更正了 Z027／Z031——它们**仍然成立**（§8.2）；
7. **不主张**"许多动"能形成宏观场——**不能**（Z026，§7.3）；
8. **不主张**本文穷尽了"从 Zero 构造场"的路——本文只做了用户的这一条。

---

## §10 结论

$$
\boxed{\;
\begin{aligned}
&\textbf{用户的四步}:\ 1,2,3\ \textbf{成立};\ 4\ \textbf{部分成立} .\\[1mm]
&\textbf{自我更正}:\ \text{"方向需定向输入"}\ \textbf{过严};\ \text{无向计数绕过它} .\\[1mm]
&\textbf{载体}:\ (h,k)\ \text{格点，}\lvert\cdot\rvert=H(H+3)/2\ (\text{而非图顶点}) .\\[1mm]
&\textbf{场的存在}:\ F=\mathbf 1\oplus\mathbf 2;\quad p_\uparrow\ne\tfrac12\ \textbf{非平凡};\ \text{壁处}\ p_\uparrow=0 .\\[1mm]
&\textbf{到不了}:\ \text{非阿贝尔、自相互作用、渐近自由、线性势、热力学极限}.
\end{aligned}
}
$$

**一句话：**

$$
\boxed{\;\textbf{你的直觉对了一大半：Zero 的层结构}\textbf{自发地}\text{给出一个方向场，且不需要定向输入；但它止步于阿贝尔、离散、无量纲。}\;}
$$

---

## 附录 A　核验

`E31_field_from_zero_check.py`：纯标准库；不 import 项目脚本、不写文件。

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | 载体点数 $H(H+3)/2$（$5/14/44$） | 5 |
| B | $S_3$ 的置换表示 $=1\oplus2$（内积核验 $1,0,1$；$\langle\mathbf 2,\mathbf 2\rangle=1$） | 8 |
| C | $p_\uparrow(h,k)$ 的精确表；壁处 $=0$；内部 $14/15$ 偏离 $>0.001$，最大 $0.1902$ | 16 |
| D | 公理③未被违反（偏置来自度量） | 4 |
| E | 到不了哪里（离散／阿贝尔／规模 $O(1)$） | 6 |
| | **合计** | **45** |

**核验纪律**：数值只算旁证；**本文更正了自己一处过严的判定**（§2）。

---

*构造稿 · 未落盘 · 未改动任何公理 · 完*
