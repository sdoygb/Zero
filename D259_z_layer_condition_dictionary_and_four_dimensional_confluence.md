# D259 · Z 层条件字典与四维时空的汇流构造

**日期**：本轮对话 · **依赖**：`AXIOMS.md` §0.3、§1、§4、`D143`、`D144`、`D155`、`D193`、`D194`、`D195`、`D196`、`D211`、`D248`、`D251`、`D252`、`D253`、`D256`、`D257`、`D258`
**测试模型**：Z 层结果类与局部历史分离、五通道零和秩、完整交换各向同性、时间线提升、体积尺度、区域字典汇流、局部标号置换和局部时间线非唯一性。它们不是 `U1-U4` 的推论。
**预先结构**：局部闭合单元、五通道零和证书、完整交换型、局部时间线、时间定向、体积尺度、区域覆盖、转移映射和细化规则。它们不是 `U1-U4` 的推论。
**核验**：[`verify/d259_z_layer_condition_dictionary_and_four_dimensional_confluence.py`](verify/d259_z_layer_condition_dictionary_and_four_dimensional_confluence.py) —— **通过 / 不符**见运行输出
**v0.5 定位**：`D211` 让 Z 保存闭合结果类，让局部历史保存在 P 层；`D194-D196` 给出零和秩、时间线和条件四维洛伦兹度规；`D257-D258` 给出交换权重、局部各向同性和外边反作用缺口。本文把这些片段整理成一张 Z 层条件字典，并明确跨区域历史如何汇流到同一个四维时空类。

> **【路线本地状态｜[`STATUS.md`](STATUS.md)】** 本文的“未导出／未解输入”是 **D259 借入路线本地**的状态，不是主 `Z/G` 路线的当前总账。主路线已经在 E1–E4 下条件恢复 GR；两条路线尚未证明等价，进度不能相加。全项目当前状态见 [`STATUS.md`](STATUS.md) §3、§7。

$$
\boxed{
\text{Z 保存“结果”，P/E 保存“怎样得到结果”。}
}
$$

$$
\boxed{
\begin{aligned}
&\text{五通道零和局域证书}\\
&+\text{完整交换各向同性}\\
&+\text{无向时间线与时间定向}\\
&+\text{体积／尺度}\\
&+\text{跨区域字典汇流}\\
&\Longrightarrow_{\rm cond}
\text{四维洛伦兹图册}.
\end{aligned}
}
$$

$$
\boxed{
\text{Z 的无时间性}
\not\Longrightarrow
\text{五通道证书、时间线、尺度或跨区域汇流自动存在}.
}
$$

本文登记恢复层结构 `R-Z-INVARIANT-HISTORY-SPLIT`、`R-Z-FIVE-CHANNEL-RANK-CERTIFICATE`、`R-Z-CLOSURE-TIME-LINE-CERTIFICATE`、`R-Z-CROSS-REGION-CONFLUENCE`、`R-Z-CONDITIONAL-4D-ATLAS` 与缺口 `R-Z-DICTIONARY-SOURCE-GAP`。

本文不修改 `U1-U4+C1`，不新增 `U5`。

---

## §0 推导过程

**第 1 步｜问题。**

现有层结构可写成

$$
E_i
\longrightarrow_{\rm closure}
P_i+\mathcal Z_\ast
\longrightarrow_{\rm seed}
E_i'.
$$

其中 $\mathcal Z_\ast$ 保存闭合结果类，局部 $P_i$ 保存形成该结果的历史。一个闭合类可以有多个不同实现：

$$
\pi_i(w_i)=\alpha,
\qquad
\pi_j(w_j)=\alpha,
\qquad
w_i\ne w_j .
$$

这符合“定理住在 Z 层，而证明住在历史层”的解释。要推进四维时空，必须回答：

> 哪些局部闭合不变量足以条件组装四维洛伦兹几何，并让不同区域的实现汇流到同一个时空类？

本文不把这个问题的答案写成 Z 的无时间性结论。本文给出一张显式条件字典。

**第 2 步｜结果与历史分离。**

把一个 Z 闭合类记为 $\alpha$。局部历史写成 $w_i$，并取两个读出：

$$
\pi_i(w_i)=\alpha\in[\mathcal C],
\qquad
\Pi_i(w_i)=\mathfrak d_i .
$$

这里 $\mathfrak d_i$ 是局部几何字典候选，而不是历史本身。

若

$$
w_i\ne w_j,
\qquad
\Pi_i(w_i)=\Pi_j(w_j),
$$

则两份历史在低能几何读出上不可区分。若它们在重叠区域还要给同一时空，则必须有汇流条件：

$$
\boxed{
\pi_i(w_i)=\pi_j(w_j)=\alpha,
\qquad
\Pi_i(w_i)|_{U_i\cap U_j}
\sim
\Pi_j(w_j)|_{U_i\cap U_j}.
}
$$

这里 $\sim$ 表示允许由局部转移映射联系的同一几何数据。

$$
\boxed{
\text{同一闭合类允许多种历史；}
\text{同一时空必需重叠一致性。}
}
$$

这登记为 `R-Z-INVARIANT-HISTORY-SPLIT`。

**第 3 步｜五通道零和秩证书。**

`D194` 已证明零和约束的局域方向数为

$$
\dim H_Q
=
\dim
\left\{
x\in\mathbb R^m:
\sum_{i=1}^m x_i=0
\right\}
=
m-1 .
$$

因此若局部闭合单元给出五个基本通道和一个零和约束，则局域方向数为

$$
m=5
\quad\Longrightarrow\quad
\dim H_Q=4 .
$$

相反，四通道只给

$$
m=4
\quad\Longrightarrow\quad
\dim H_Q=3 .
$$

所以四维局部方向的最小零和证书是：

$$
\boxed{
\text{五通道零和局部闭合单元}
\Longrightarrow_{\rm cond}
\text{四个独立局域方向}.
}
$$

这不是“四维时空已经导出”。它只把维数输入从裸写的 $d=4$ 改写为一张更接近闭合事件语言的证书。

$$
\boxed{
\text{五通道证书本身仍是恢复层输入，}
\text{只是它把维数输入压成秩输入。}
}
$$

这登记为 `R-Z-FIVE-CHANNEL-RANK-CERTIFICATE`。

**第 4 步｜完整交换给各向同性正定型。**

在 $H_Q$ 上取完整交换二次型

$$
\mathcal Q(p)
=
\sum_{1\le i<j\le5}
(p_j-p_i)^2 .
$$

因为 $\sum_i p_i=0$，有恒等式

$$
\mathcal Q(p)
=
5\sum_{i=1}^5p_i^2 .
$$

因此完整交换型在四个局域方向上各向同性。

`D258` 给出更一般的交换投影：

$$
A_E^0
=
\sum_e K_e d_e^0\otimes d_e^0
=
\operatorname{proj}_{\operatorname{Ran}L}.
$$

当局部闭合并只包含完整五通道交换图时，范围维数为四，故

$$
A_E^0=I_4
$$

在四个局域方向上成立。

$$
\boxed{
\text{完整局域交换}
\Longrightarrow_{\rm cond}
\text{各向同性四维正定型 }g_R.
}
$$

对四顶点空间单元，`D258` 的 $K_4$ 正例给出三维各向同性空间截面。对外部边进入局部交换图的情况，`R-Z-LOCAL-EXCHANGE-BACKREACTION-GAP` 仍然有效。

**第 5 步｜时间线证书。**

四维正定型 $g_R$ 本身没有洛伦兹号差。按 `D195`，还需要一条无向时间线

$$
[\ell],
\qquad
\dim \ell_x=1 .
$$

取单位截面 $u$，定义

$$
g_L
=
2u^\flat\otimes u^\flat-g_R .
$$

在 $u$ 方向与 $u^\perp$ 上分别有

$$
g_L(u,u)=1,
\qquad
g_L|_{u^\perp}=-g_R|_{u^\perp}.
$$

因此

$$
\boxed{
\operatorname{signature}(g_L)=(1,3).
}
$$

但 `D195` 已证明零和置换对称不能自然选出这条线。本文也保留这一点：

$$
\boxed{
\text{无向时间线是独立闭合证书，不是无时间 Z 层的自动结果。}
}
$$

这登记为 `R-Z-CLOSURE-TIME-LINE-CERTIFICATE`。

**第 6 步｜时间定向与体积尺度。**

为了把无向线分成未来和过去，还需要时间定向

$$
o:\ell\to\{\text{未来},\text{过去}\}.
$$

为了把共形类 $ [g_R]$ 固定为具体代表，还需要体积尺度 $\mu$。按 `D143`／`D196`，给定

$$
[g_R],\qquad \mu,\qquad [\ell],\qquad o,
$$

可以条件得到唯一四维洛伦兹度规代表 $g_L$，并保持

$$
\operatorname{vol}_{g_L}=\mu .
$$

于是局部字典中的第一条桥为

$$
\boxed{
(\text{五通道秩},\text{各向同性 }g_R,[\ell],o,\mu)
\Longrightarrow_{\rm cond}
(U_\sigma,g_{L,\sigma}).
}
$$

它只给局部图表，不给跨区域汇流。

**第 7 步｜跨区域字典汇流。**

设区域 $U_i,U_j$ 有非空重叠 $U_i\cap U_j$。若两个局部读出的历史不同，但仍要描述同一时空，则需存在转移映射

$$
\Phi_{ji}:U_i\cap U_j\longrightarrow U_j\cap U_i
$$

使得

$$
\Phi_{ji}^*g_{L,j}=g_{L,i},
\qquad
\Phi_{ji}^*[\ell_j]=[\ell_i],
\qquad
\Phi_{ji}^*o_j=o_i,
\qquad
\Phi_{ji}^*\mu_j=\mu_i .
$$

若这些条件在全部重叠上成立并满足 cocycle 相容性

$$
\Phi_{ki}
=
\Phi_{kj}\circ\Phi_{ji},
$$

则局部图表可粘成一张条件四维洛伦兹图册：

$$
\boxed{
\{(U_i,g_{L,i},[\ell_i],o_i,\mu_i,\Phi_{ji})\}
\Longrightarrow_{\rm cond}
(M,g_L).
}
$$

这登记为 `R-Z-CROSS-REGION-CONFLUENCE` 与 `R-Z-CONDITIONAL-4D-ATLAS`。

如果两个区域给出同一个 Z 闭合类，但在重叠上不能由 $\Phi_{ji}$ 相容，那么不能汇流为单一四维时空。此时 Z 只保存共同结果类，区域历史仍给不同局部几何候选。

$$
\boxed{
\text{同一 Z 结果}
\not\Longrightarrow
\text{自动同一局部几何}.
}
$$

**第 8 步｜结果先于历史还是历史先于结果。**

在活动层和历史层中，闭合过程有局部前缀顺序：

$$
w_i\prec_{\rm loc}\alpha .
$$

这里“$\prec_{\rm loc}$”只表示局部历史的闭合写入，不表示 Z 内部的时间。

在 Z 内部，$\alpha$ 没有内部先后：

$$
\alpha\not\prec_Z\alpha',
\qquad
\text{Z 内部无时间推进}.
$$

因此本文支持的表述是：

$$
\boxed{
\text{局部历史中，过程先闭合成结果；}
\quad
\text{Z 内部，结果没有先后属性。}
}
$$

本文不声称 $\alpha$ 先在 Z 中出现并因果选择历史。若需要这种解释，必须额外加入律在先或读出选择器。

**第 9 步｜无时间 Z 层不能单独给四维。**

仅有闭合类多重集

$$
\mathcal Z_\ast=\mathbb N^{([\mathcal C])}
$$

时，缺少以下数据：

1. 局部五通道秩证书；
2. 各向同性正定型；
3. 时间线与时间定向；
4. 体积尺度；
5. 区域覆盖与转移映射；
6. 细化与局部化规则。

其中任一项缺失，四维洛伦兹图册都不能由 $\mathcal Z_\ast$ 单独恢复。

$$
\boxed{
\mathcal Z_\ast
\not\Longrightarrow
(M,g_L).
}
$$

这不是否定 Z 层，而是说明 Z 层保存的是四维时空的“结果类”，不是构造四维时空的完整字典。

**第 10 步｜条件字典的净结果。**

本步把目标结构压成：

$$
\boxed{
\mathfrak D
=
\left(
m=5,\;
g_R,\;
[\ell],\;
o,\;
\mu,\;
\Phi_{ij}
\right)
}
$$

并给出条件链：

$$
\boxed{
\mathfrak D
\Longrightarrow_{\rm cond}
(M,g_L).
}
$$

它不导出 $\mathfrak D$。但它的价值是：只要后续能从 Z 层闭合事件、历史层与支持映射中逐项生成 $\mathfrak D$，四维时空就不必作为裸输入塞进来。

当前最小下一步应是逐项生成：

1. 五通道零和局部证书；
2. 无向时间线；
3. 体积尺度；
4. 跨区域汇流映射；
5. 外部边局部隔离。

---

## §1 核验内容

核验脚本验证：

1. D259 与恢复结构已登记；
2. Z 结果类与历史实现分离；
3. 五通道零和秩为四，四通道秩为三；
4. 完整交换型的各向同性；
5. 时间线提升给 $(1,3)$ 号差；
6. 体积与共形尺度可条件归一；
7. 正确转移映射可汇流，错误映射会破坏同一度规；
8. 零和置换对称不能选出唯一时间线；
9. 文档登记条件字典输入与缺口；
10. 文档不把 Z 层无时间性写成四维时空导出；
11. 上游边界保持；
12. 文档不引用项目外体系名或外部路径。

---

## §2 恢复层结构

| 标识 | 含义 | 状态 |
|:--|:--|:--|
| `R-Z-INVARIANT-HISTORY-SPLIT` | Z 保存闭合结果类，P/E 保存不同历史实现；同一结果可有多个历史 | 条件结构 |
| `R-Z-FIVE-CHANNEL-RANK-CERTIFICATE` | 五通道加零和约束给四个局域方向，是四维的最小秩证书 | 恢复层输入 |
| `R-Z-CLOSURE-TIME-LINE-CERTIFICATE` | 局部四条方向中选出无向时间线；零和对称不能自动选择 | 恢复层输入 |
| `R-Z-CROSS-REGION-CONFLUENCE` | 不同区域历史若给同一 Z 类，局部字典还须在重叠上由转移映射相容 | 条件约束 |
| `R-Z-CONDITIONAL-4D-ATLAS` | 秩、各向同性、时间线、定向、体积与汇流条件给四维洛伦兹图册 | 条件构造 |
| `R-Z-DICTIONARY-SOURCE-GAP` | 五通道证书、时间线、体积尺度和汇流映射的来源仍未导出 | 未解输入 |

这些结构不修改 `U1-U4+C1`，也不新增 `U5`。

---

## §3 输入预算

| 输入 | 作用 | 状态 |
|:--|:--|:--|
| 五通道零和局域证书 | 给四个局域方向 | 未导出 |
| 完整局域交换图 | 给各向同性四维正定型 | `D258` 条件正例；局部化仍缺 |
| 无向时间线 | 把正定型提升为洛伦兹型 | 未导出 |
| 时间定向 | 区分未来与过去 | 未导出 |
| 体积／尺度 | 从共形类选出具体度规 | 未导出 |
| 区域覆盖与转移映射 | 把局部图粘成全局图册 | 未导出 |
| 细化规则 | 给共同连续极限 | 未导出 |
| 外部边局部隔离 | 防止局部各向同性被全局交换破坏 | `R-Z-LOCAL-EXCHANGE-BACKREACTION-GAP` |

---

## §4 结论

Z 层可以作为无时间的闭合结果层，保存类似

$$
1+1=2
$$

这样的结果类；其局部证明、演化和区域差异放在 P/E 层。四维时空若要在这一框架中出现，不应被写成 Z 层的裸预设，而应写成跨区域闭合历史共同汇流到的闭合时空类。

$$
\boxed{
\text{Z 保存结果类；}
\quad
\text{历史保存局部实现；}
\quad
\text{四维时空由条件字典与区域汇流恢复。}
}
$$

本文不声称四维时空已经无条件导出，也不声称 Z 的无时间性自动给出四维。本文给出的是一张可逐项攻击的条件字典：下一步应优先从闭合事件和局部支持生成五通道证书、时间线与跨区域汇流。
