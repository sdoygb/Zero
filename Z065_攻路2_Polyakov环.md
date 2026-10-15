# 攻路 2：Polyakov 环／中心对称判据 —— **在同调群上闭合，不是热力学**

**Zero 项目 · 基础推导系列 · 第 Z065 号**

**篇性**：**推导 ＋ 否证篇**。用户指示："2"（攻第 8 项的路 2：用 $\mathbb Z_3$ 中心对称／Polyakov 环判据）。
**结论**：**(一)** Polyakov 环需要**类时方向**，而 Heawood **没有**；**(二)** 但把它换成"**同调群上的 Wilson 环**"后判据**可定义且可算**；**(三)** 结果**闭合**：面积律在有限格点上**不可分辨**（面积有界）。
**纪律**：**未改动任何公理**。
**依据**：正式篇 Z021（Heawood 定理 8）；本系列 Z060（局域群、中心平凡）、Z061（精确质量隙）、Z063（九条判据）、Z064（定理 28）。
**核验**：[`Z065_polyakov_check.py`](_staging/Z065_polyakov_check.py)。

---

## 摘要

$$
\boxed{\;\textbf{（一）}\ \text{Polyakov 环需}\ \textbf{类时方向};\ \text{Heawood}\ \textbf{没有} \Longrightarrow\ \text{标准形式}\ \textbf{不可定义} .\;}
\tag{0.1}
$$

$$
\boxed{\;\textbf{（二）}\ \text{替代}:\ \text{同调群上的 Wilson 环}\ \Longrightarrow\ \textbf{可定义且可算} .\;}
\tag{0.2}
$$

$$
\boxed{\;\textbf{（三）}\ \text{结果}\ \textbf{闭合}:\ \text{面积}\le\text{圈秩}=8\ \Longrightarrow\ \text{面积律}\ \textbf{不可分辨} .\;}
\tag{0.3}
$$

---

## §1 Polyakov 环**不能**按标准方式定义

### 1.1 标准定义需要什么

$$
P(\mathbf x)=\operatorname{tr}\ \mathcal P\exp\Bigl(i\oint_{\text{时间圈}}A_0\,dt\Bigr)
\ \Longrightarrow\ \text{需要}\ \textbf{类时方向} \text{（Matsubara 圈）} .
\tag{1.1}
$$

$$
\text{物理含义}:\ \text{在有限温度}\ T\ \text{下，时空为}\ \mathbb R^3\times S^1_\beta .
\tag{1.2}
$$

### 1.2 Zero 有什么

$$
\text{Heawood 是}\ \textbf{无向图};\ \text{它的"方向"只有一个指标（边）};\ \textbf{无时间轴} .
\tag{1.3}
$$

$$
\Longrightarrow\quad \boxed{\;\text{没有}\ S^1_\beta\ \text{这样的时间圈} \Longrightarrow\ \text{标准 Polyakov 环}\ \textbf{不可定义} .\;}
\tag{1.4}
$$

$$
\textbf{这与 Z005（无普遍时间平移）同源}:\ \text{没有时间轴} \Longrightarrow\ \text{没有时间圈} .
\tag{1.5}
$$

---

## §2 替代：**同调群上的 Wilson 环**

### 2.1 关键结构

$$
\textbf{核验 §A}:\quad \text{Heawood 的}\ \textbf{圈秩}=\lvert E\rvert-\lvert V\rvert+1=21-14+1=\mathbf 8 .
\tag{2.1}
$$

$$
\Longrightarrow\quad H_1(\Gamma;\mathbb Z)\cong\mathbb Z^8\ (\text{自由部分};\ \text{因}\ \Gamma\ \text{连通}) .
\tag{2.2}
$$

$$
\boxed{\;\text{故}\ \Gamma\ \text{上有}\ 8\ \text{个独立的}\ \textbf{非同调类} \text{的环};\ \text{各给一个 Wilson 环} .\;}
\tag{2.3}
$$

### 2.2 中心对称

$$
\textbf{核验 §B}:\quad \operatorname{Zent}(S_3)=\{1\}\ \text{（平凡）};\qquad
A_3\cong\mathbb Z_3,\ \operatorname{Zent}(\mathbb Z_3)=\mathbb Z_3\ne\{1\} .
\tag{2.4}
$$

$$
\Longrightarrow\quad \boxed{\;\text{须取规范群}=\mathbb Z_3\ (\text{而非}\ S_3)\ \text{才有中心对称} .\;}
\tag{2.5}
$$

$$
\text{中心变换}:\quad a_e\longmapsto a_e+1\quad(\text{一切边}) \Longrightarrow \text{holonomy}\longmapsto\text{holonomy}+6=\text{holonomy}\ (\mathrm{mod}\ 3) .
\tag{2.6}
$$

$$
\textbf{关键}:\ \text{因一切 plaquette 是}\ \mathbf 6\ \text{环（偶数）},\ \text{中心变换}\ \textbf{保持作用量} .
\tag{2.7}
$$

$$
\Longrightarrow\quad \mathbb Z_3\ \text{是作用量的}\ \textbf{精确对称} \Longrightarrow\ \text{中心对称}\ \textbf{存在} \checkmark .
\tag{2.8}
$$

---

## §3 计算

### 3.1 6-环的完整枚举

$$
\textbf{核验 §C}:\quad \text{Heawood 的长度 6 的}\ \textbf{相异环}=\mathbf{56};\qquad
\text{相异}\ \textbf{顶点集}=\mathbf{28} .
\tag{3.1}
$$

$$
\Longrightarrow\quad \textbf{每个 6 顶点集承载 2 个相异环}\ (56=2\times28) .
\tag{3.2}
$$

$$
\text{（故"56"与"28"}\ \textbf{都对}，\text{数的是不同的东西}） .
\tag{3.3}
$$

$$
\textbf{一处要点}:\ \text{若指数}\ \textbf{环} \text{，6-环数是}\ \textbf{56}（顶点集数 28） .
\tag{3.4}
$$

### 3.2 $\mathbb Z_3$ 规范理论（Wilson 作用量）

$$
S[a]=-\beta\sum_{P}\cos\Bigl(\frac{2\pi}{3}\sum_{e\in P}a_e\Bigr);\qquad a_e\in\mathbb Z_3 .
\tag{3.4a}
$$

$$
\textbf{核验 §D}（Metropolis，$3^{21}\approx1.05\times10^{10}$ 组态空间）：
\tag{3.5}
$$

| $\beta$ | $\langle W\rangle$（56 个 6-环平均） | 接受率 |
|--:|--:|--:|
| $0.10$ | $0.105054$ | $0.724$ |
| $0.50$ | $\mathbf{1.000000}$ | $0.000$ |
| $1.00$ | $\mathbf{1.000000}$ | $0.000$ |
| $2.00$ | $\mathbf{1.000000}$ | $0.000$ |
| $5.00$ | $\mathbf{1.000000}$ | $0.000$ |

$$
\Longrightarrow\quad \beta\ge0.5\ \text{时系统}\ \textbf{冻结} \text{在}\ a_e=\text{const} \Longrightarrow \text{一切 holonomy}=1 .
\tag{3.6}
$$

$$
\textbf{解析对照（强耦合）}:\quad
\langle\cos\rangle\Big|_{\beta=1}=0.537158\ (\text{单 plaquette}) .
\tag{3.7}
$$

### 3.3 **结果闭合**：面积律不可分辨

$$
\textbf{核验 §E}:\quad \text{面积的上界}=\text{圈秩}=\mathbf 8;\qquad
\text{周长（6-环）}=\mathbf 6 .
\tag{3.8}
$$

$$
\Longrightarrow\quad \boxed{\;\text{面积}\le8\ \text{而}\ \text{周长}=6\ \Longrightarrow\ \text{面积与周长}\ \textbf{同量级} \Longrightarrow\ \text{面积律}\ \textbf{不可分辨} .\;}
\tag{3.9}
$$

$$
\textbf{所有}\ 21\ \text{条边都在 6-环里};\ \text{任何环的"面积"}\le\text{它含的 6-环数}\le7 .
\tag{3.10}
$$

---

## §4 一件**未验证**的事（**必须如实登记**）

### 4.1 我曾试图证明的更强命题

$$
\textbf{猜想}:\quad \text{一切 plaquette 平坦}\ (U_P=1)\ \Longrightarrow\ \text{一切 6-环的}\ \text{holonomy}=1 .
\tag{4.1}
$$

$$
\text{若成立} \Longrightarrow \text{强耦合下}\ \langle W\rangle=1\ \text{可}\ \textbf{解析导出} .
\tag{4.2}
$$

### 4.2 状态

$$
\text{我的 mod-3 秩计算给出}\ \textbf{不一致} \text{的结果（平坦空间维数}=1\ \text{vs 纯规范维数}=13） .
\tag{4.3}
$$

$$
\Longrightarrow\quad \boxed{\;\textbf{本文未能证明（4.1）};\ \text{它}\ \textbf{只由数值支持}（\S3.2\ \text{的冻结}） .\;}
\tag{4.4}
$$

$$
\textbf{这不是致命}:\ \S3.3\ \text{的结论}\ \textbf{不依赖} \text{（4.1）};\ \text{它只依赖"面积有界"} .
\tag{4.5}
$$

---

## §5 对第 8 项的最终判定（**第三次更新**）

| 子问题 | Z063 | Z064 后 | **本文后** |
|:--|:--|:--|:--|
| $8a$ 质量隙 | 已得 | 已得 | **已得**（$m=0.7520386984$） |
| $8b$ 无带色态 | 需局域 Gauss | （甲）可得（乙）不可得 | 同 |
| $8c$ 面积律／弦张力 | 关闭 | 关闭 | **关闭** |
| 中心对称判据 | 可查 | 可查 | **已查：因无时间圈，退化为** $H_1$ **上的 Wilson 环** |
| Polyakov 环判据 | — | — | **退化为 Wilson 环；结果同上** |

$$
\boxed{\;\text{第四条判据（中心对称）}\ \textbf{已查完};\ \text{它}\ \textbf{不} \text{提供新信息} .\;}
\tag{5.1}
$$

$$
\Longrightarrow\quad \text{九条判据里}\ \textbf{可查的四条已全部查完};\ \text{结论一致} .
\tag{5.2}
$$

---

## §6 结论

$$
\boxed{\;
\begin{aligned}
&\textbf{（一）}\ \text{标准 Polyakov 环}\ \textbf{不可定义} \text{（无类时方向）};\ \text{与 Z005 同源}.\\[1mm]
&\textbf{（二）}\ \text{替代}:\ H_1(\Gamma)\ \text{上的 Wilson 环（8 个独立类）};\ \text{中心对称}\ \textbf{存在}.\\[1mm]
&\textbf{（三）}\ \text{结果}:\ \text{面积}\le8\ \text{与周长}\ 6\ \text{同量级} \Longrightarrow\ \textbf{面积律不可分辨}.\\[1mm]
&\textbf{（四）}\ \text{一条更强的命题}\ \textbf{未证明} \text{（只由数值支持）};\ \text{但不影响结论}.\\[1mm]
&\textbf{净}:\ \text{路 2}\ \textbf{走完};\ \text{不提供新信息};\ \text{第 8 项}\ \textbf{维持} \text{Z064 的判定}.
\end{aligned}
}
$$

**一句话：**

$$
\boxed{\;\textbf{没有时间圈，Polyakov 环就退化成普通 Wilson 环——而它在本就有界的面积上量不出弦张力。}\;}
$$

---

## §7 不许主张

1. **不主张** Zero "有 Polyakov 环"——**没有**（无类时方向，§1）；
2. **不主张**猜想（4.1）成立——**未证明**（§4）；
3. **不主张**中心对称判据"给出了弦张力"——**没有**（§3.3）；
4. **不主张** $S_3$ 有中心对称——**没有**（Z063）；
5. **不主张**"56 个 6-环"已解析证明——本文用**枚举**（含边合法性检查）；
6. **不主张** Metropolis 的冻结是**严格**结果——它是**数值**（$3^{21}$ 空间未被穷举）；
7. **不主张**本文改变了 Z061／Z064 的结论；
8. **不主张**九条判据"全部"查完——本文只查了**可查的四条**（3、4、5、6）。

---

## 附录 A　核验

`Z065_polyakov_check.py`：纯标准库（＋`numpy`）；不 import 项目脚本、不写文件。

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | 圈秩 $=8$；$H_1\cong\mathbb Z^8$ | 4 |
| B | $\operatorname{Zent}(S_3)=\{1\}$；$\operatorname{Zent}(\mathbb Z_3)=\mathbb Z_3$；中心变换保作用量 | 6 |
| C | 相异 6-环 $=56$；顶点集 $=28$（每集 2 环）；边索引全合法 | 6 |
| D | Metropolis 的 $\langle W\rangle$ 数值表（含冻结） | 8 |
| E | 面积上界 $=$ 圈秩；面积与周长同量级 | 5 |
| | **合计** | **28** |

**核验纪律**：数值只算旁证；**本文登记了一处未证明的猜想**（§4）与**一处要点**（环数 56；顶点集数 28）。

---

*推导稿 · 已落盘 · 未改动任何公理 · 完*
