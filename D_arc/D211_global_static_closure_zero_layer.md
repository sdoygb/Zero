# D211 · 全局静态闭合零层：多重零记录、局部历史与全局读回边界

**日期**：本轮对话 · **依赖**：`AXIOMS.md` §0.3、§3、§6、§7、`D175`、`D178`、`D179`、`D180`、`D187`、`D191`、`D210`
**测试模型**：有限零和闭合历史、循环等价类、多重集记录、局部历史纤维、可逆周期与吸收终端补全。
**预先结构**：闭合谓词、循环等价、全局多重集、吸收写入规则、局部历史划分与局部重建规则。
**核验**：[`verify/d211_global_static_closure_zero_layer.py`](../verify/d211_global_static_closure_zero_layer.py) —— **通过 / 不符**见运行输出
**数值模拟**：[`simulations/zero_sum_closure_exit.py`](../simulations/zero_sum_closure_exit.py)、[`simulations/zero_sum_closure_exit.md`](../zero_sum_closure_exit.md)、[`simulations/zero_sum_open_reservoir.py`](../simulations/zero_sum_open_reservoir.py)、[`simulations/zero_sum_open_reservoir.md`](../zero_sum_open_reservoir.md)、[`simulations/zero_sum_periodic_destruction.py`](../simulations/zero_sum_periodic_destruction.py)、[`simulations/zero_sum_periodic_destruction.md`](../zero_sum_periodic_destruction.md)
**v0.5 定位**：`D191` 的零层是满足零约束的动态状态层；本文新增的 $\mathcal Z\_\ast$ 是已经闭合的振动模式停驻其上的全局静态记录层。前者承载运动，后者承载完成的闭合结果。

$$

\text{动态零层 } \mathcal Z_{\rm dyn}\text{ 承载未完成的零约束运动；}
\quad
\text{闭合零层 } \mathcal Z_\ast\text{ 承载已经完成的零闭合记录。}

$$

---

## §0 推导过程

**第 1 步｜要检验的命题。**

设局部区域为 $i$。候选层级为

$$
\mathcal Z_\ast
\quad\longrightarrow\quad
P_i
\quad\longrightarrow\quad
E_i
\quad\longrightarrow\quad
D_i .
$$

其中：

| 层 | 作用 | 范围 |
|:--|:--|:--|
| $\mathcal Z\_\ast$ | 已经闭合的振动模式及其多重零记录 | 全局 |
| $P\_i$ | 第 $i$ 区的精确闭合历史 | 局部 |
| $E\_i$ | 第 $i$ 区的有限寿命活动过程 | 局部 |
| $D\_i$ | 第 $i$ 区销毁的开放路径 | 局部 |

要检验的不是“所有东西都是零”，而是：

1. 全局零层能否容纳多个不同的闭合模式？
2. “停在零”是记录静止，还是动力学吸收？
3. 全局零层会不会抹掉局部历史和区域差异？

**第 2 步｜闭合模式必须先有区分。**

取闭合历史集合

$$
\mathcal C
=
\{w:\ \mathsf{Closed}(w)\}.
$$

同一个闭合模式允许存在不同相位起点。定义循环等价

$$
w\sim w'
\quad\Longleftrightarrow\quad
w'\ \text{是}\ w\ \text{的循环移位}.
$$

于是闭合模式类为

$$
[\mathcal C]
=
\mathcal C/{\sim}.
$$

如果 $\mathcal Z\_\ast$ 只是字面数值 $0$，则所有类都必须映射到同一个点：

$$
\pi([w])=0_*
\qquad
\text{对所有 }[w]\in[\mathcal C].
$$

此时

$$
\ker\pi=[\mathcal C],
$$

不同闭合模式无法再被区分。

$$

\text{“很多零”不能是同一个标量 }0\text{；}
\text{每个零必须带着自己的闭合类标签。}

$$

**第 3 步｜全局零层是多重的。**

令

$$
\mathcal Z_\ast
=
\mathbb N^{([\mathcal C])}.
$$

一个元素可以写成有限形式和

$$
Z
=
\sum_{\alpha\in[\mathcal C]}
n_\alpha\,0_\alpha,
\qquad
n_\alpha\in\mathbb N,
$$

其中 $0\_\alpha$ 表示“闭合类 $\alpha$ 的静态零记录”。

所以 $0\_\alpha\ne0\_\beta$ 当且仅当两个闭合模式类不同。若同一模式闭合两次，则 $n\_\alpha$ 增加。

$$

\mathcal Z_\ast\text{ 是闭合零记录的多重集，不是单点零集合。}

$$

如果物理上不需要记录重复次数，可以退化为集合；但那样会丢失谱系数量和闭合多重数。

**第 4 步｜“静态”只描述层内动力学。**

层内没有从一个闭合模式到另一个闭合模式的转变：

$$
T_{\mathcal Z_\ast}(0_\alpha)=0_\alpha.
$$

当局部过程完成一次闭合时，只允许向层内写入新记录：

$$
Z_{t+1}
=
Z_t
+
\sum_i\Delta_i(t),
\qquad
\Delta_i(t)\in\mathbb N^{([\mathcal C])}.
$$

若不允许删除闭合记录，则每个坐标满足

$$
n_\alpha(t+1)\ge n_\alpha(t).
$$

$$

\text{层内静态}\ne\text{总层不变；}
\text{静态指没有内部模式跃迁，写入仍可增加记录。}

$$

**第 5 步｜“停在零”有两种解释。**

第一种是时间遗忘商：

$$
w\longmapsto [w],
$$

物理模式完成一个周期后仍可继续运动。进入 $\mathcal Z\_\ast$ 的不是带时间的轨道，
而是轨道在“忘掉相位与时间参数”后的闭合类。这个商把整条周期轨道压成一个静态零记录：

$$
\mathcal Q_T:
\{\text{闭合周期轨道}\}
\longrightarrow
\mathcal Z_\ast,
\qquad
\gamma\longmapsto0_{[\gamma]}.
$$

$\mathcal Q\_T$ 不是时间中的状态跃迁，而是把 $S^1$ 型周期轨道识别为一个点的层间读出。
因此它不需要阻尼，也不要求物理分支在 $E\_i$ 中消失。

零层内部没有时间坐标，所以同一条闭合轨道上的不同相位在层内不可区分：

$$
\gamma(t_1)\sim_{\mathcal Z}\gamma(t_2)
\qquad
\text{对所有 }t_1,t_2.
$$

这个不可区分关系把整条闭合轨道压成一个静态零记录。若引入一个仍能读出相位的
外部时间观察者，则看到的仍是运动中的周期轨道，而不是静止零。

第二种是吸收式停止：

$$
T(w)=0_{[w]},
\qquad
T(0_{[w]})=0_{[w]}.
$$

这是非可逆终端映射。零和约束本身不会在两种补全之间自动选择；`D175-D179` 已证明吸收方向与终端记录需要额外规则。

$$

\text{闭环可以只产生静态记录；}
\text{若要求物理振动真正停在零，则还必须加入非可逆吸收规则。}

$$

本文采用时间遗忘商作为最小版本，把第二种登记为额外输入。

**第 5A 步｜零作用量为什么不能单独推出物理停止。**

取谐振子

$$
L
=
\frac12m\dot x^2
-
\frac12m\omega^2x^2,
\qquad
x(t)=A\sin(\omega t).
$$

在一个完整周期 $T=2\pi/\omega$ 内，位置返回零：

$$
x(0)=x(T)=0.
$$

总作用量为

$$
S
=
\int_0^T
\left(
\frac12m\dot x^2
-
\frac12m\omega^2x^2
\right)
dt
=
0.
$$

但在返回点速度并不为零：

$$
\dot x(0)=\dot x(T)=A\omega\ne0.
$$

所以

$$

S=0
\quad\text{和}\quad
x(0)=x(T)=0
\quad\not\Longrightarrow\quad
(x,\dot x)=(0,0).

$$

零作用量表示一整圈的净作用量抵消，不表示轨迹在终点失去速度。若要求物理运动真正停止，还必须额外要求

$$
\dot x(T)=0
$$

并给出使终点成为稳定或吸收终端的机制。

$$

\text{闭合给出周期等价类；}
\text{时间遗忘商把它变成静态零；}
\text{零作用量只保证净作用量为零。}

$$

**第 5B 步｜闭合分支退出活动层。**

你现在补的机制可以写成一条分支终止规则。设活动层 $E\_i$ 只保存尚未闭合的
路径。若某个历史 $w$ 已经闭合，则它在这一层不再有后继：

$$
w\in\mathcal C_i
\quad\Longrightarrow\quad
\text{Succ}_{E_i}(w)=\varnothing.
$$

闭合事件不继续延长 $w$，而是把它送到记录层：

$$
w
\longrightarrow
(w,[w])
\in
P_i\times\mathcal Z_\ast,
\qquad
w\notin E_i'.
$$

于是活动层中的每条路径只有三种去向：

$$
E\to E
\quad\text{(继续未闭合演化)},
$$

$$
E\to P+Z
\quad\text{(闭合，演化分支结束)},
$$

$$
E\to D
\quad\text{(层寿命耗尽，开放路径销毁)}.
$$

这不需要物理轨迹速度为零，也不需要阻尼。它只要求：**闭合是活动层中的终点。**

$$

\text{闭合分支退出 }E_i\text{；}
\text{精确历史进入 }P_i\text{，闭合类进入 }\mathcal Z_\ast.

$$

这条规则是分支局部的。一个分支闭合并退出，不会让其他开放分支停止；
全局演化仍在其他 $E\_j$ 中继续。

若还要求一个在外部时间中本可继续运动的物理状态真正停止，则那仍然是
`R-Z-ZERO-STOPPING-GAP`。若只要求演化分支结束并留下记录，则本步已经足够。

**第 5C 步｜活动层的定期毁灭与再生。**

设毁灭周期为 $T$。在第 $n$ 个周期末，

$$
t=nT.
$$

活动层中尚未闭合的路径全部进入 $D\_i$：

$$
E_i^{(n)}\longrightarrow D_i.
$$

此时

$$
E_i=\varnothing,
$$

但历史层和全局闭合零层保留：

$$
P_i\ \text{保留},
\qquad
\mathcal Z_\ast\ \text{保留}.
$$

下一活动层由闭合历史重新播种：

$$
E_i^{(n+1)}
=
\text{Seed}(P_i,\mathcal Z_\ast).
$$

在最小程序中，每条闭合历史 $w$ 生成两个开放种子：

$$
w\longmapsto\{w+,w-\}.
$$

这条规则形成：

$$

E
\overset{\ \text{周期毁灭}\ }{\longrightarrow}
D,
\qquad
P+\mathcal Z_\ast
\overset{\ \text{重新播种}\ }{\longrightarrow}
E'.

$$

毁灭周期 $T$ 和播种映射 $\text{Seed}$ 都仍是恢复层输入。
程序只验证它们能形成自洽的毁灭和再生循环，不从零和约束唯一导出它们。

**第 5D 步｜$D\_i$ 是终端汇，不是返回通道。**

周期毁灭时，进入 $D\_i$ 的只能是尚未闭合的活动路径：

$$
w\in E_i^{(n)}
\quad\text{且}\quad
\mathsf{Closed}(w)=\mathsf{false}
\quad\Longrightarrow\quad
w\in D_i.
$$

已经闭合的 $w$ 在更早一步已经退出活动层，写入
$(w,[w])\in P\_i\times\mathcal Z\_\ast$，所以它不会进入 $D\_i$。
因此 $D\_i$ 在最小模型中只有入边：

$$
E_i\longrightarrow D_i,
\qquad
\text{Succ}_{D_i}=\varnothing.
$$

更具体地说，本文要求：

1. $D\_i$ 不向 $E\_i$ 回写；
2. $D\_i$ 不向 $P\_i$ 回写；
3. $D\_i$ 不向 $\mathcal Z\_\ast$ 回写；
4. $D\_i$ 不参与 $\text{Seed}(P\_i,\mathcal Z\_\ast)$；
5. $D\_i$ 只保存已经发生的毁灭记录，供审计，不充当新活动层的源。

因此毁灭不是“物质从活动层转化为另一种可继续产生后代的物质”，而是
**一条开放分支失去继续资格并进入终端账本**。它也不自动提供时间箭头：
终端方向来自周期毁灭规则本身，不是由 $D\_i$ 自动导出。

零和方面还有一条必须保留的边界。令

$$
Q_{D_i}
=
\sum_{w\in D_i}
\sum_t w_t.
$$

最小程序验证的是

$$
\sum_i Q_{D_i}=0,
$$

而不是

$$
Q_{D_i}=0
\qquad\text{对每个 }i.
$$

事实上模拟得到

$$
Q_D
=
[8,-8,16,-16,24,-24].
$$

单个局部毁灭层可以带走非零平衡，补偿发生在不同区域之间。因此若物理上
要求每个区域独立零和，必须新增成对毁灭规则

$$
\ell+(-\ell)\to D_i,
$$

或者给每个 $D\_i$ 增加局部补偿记录。仅凭“总作用量为零”与当前毁灭
规则，不能推出局部 $D\_i$ 自守恒。

$$

\begin{aligned}
&D_i\text{ 是只进不出的局部终端汇；}\\
&D_i\text{ 不参与重播，也不是新的活动层；}\\
&\sum_iQ_{D_i}=0\text{ 已核验，但 }Q_{D_i}=0\text{ 未被推出。}
\end{aligned}

$$

**第 6 步｜全局范围和全局读回不是一回事。**

$\mathcal Z\_\ast$ 全局表示所有局部区域共同写入同一个闭合记录层。
每个区域先产生一个局部贡献

$$
\Delta Z_i
\in
\mathbb N^{([\mathcal C])},
$$

全局层是这些贡献的和：

$$
\mathcal Z_\ast
=
\sum_i\Delta Z_i.
$$

但精确历史仍留在局部：

$$
\pi_i:P_i\to[\mathcal C],
\qquad
\pi_i(w)=[w].
$$

同一个闭合类可以来自不同本地路径。若只保留 $Z$，精确路径已经丢失：

$$
\pi_i(a)=\pi_i(b)=\alpha
\quad
\not\Longrightarrow
\quad
a=b.
$$

因此全局可见的 $Z$ 不能自动恢复局部 $P\_i$。

$$

\text{全局零层保存闭合结果类；局部历史层保存它怎样被组装。}

$$

**第 7 步｜为什么局部重建必须留给 $P\_i$。**

若允许所有区域从全局零层直接重建活动历史：

$$
\widehat P_i
=
\rho(\mathcal Z_\ast),
$$

则所有区域的 $\widehat P\_i$ 相同。于是收敛为

$$
\text{overlap}(\widehat P_i,\widehat P_j)
\to 1.
$$

前一轮有限模拟给出：

| 重建来源 | 最终站点历史重叠 | 最终主导结果数 |
|:--|--:|--:|
| 局部 $P\_i$ | 0.011 | 6 |
| 全局历史 | 0.9995 | 1 |
| 全局代表历史 | 0.9990 | 1 |

$$

\text{全局 }Z\text{ 可以是共享记录，但不能代替局部 }P_i\text{ 充当本地历史的唯一源。}

$$

**第 8 步｜零层不是动态零态层。**

`D191` 中的动态零层为

$$
\mathcal Z_{\rm dyn}
=
\{x:Q(x)=0\}.
$$

它是仍可继续运动的状态约束层。

本文的闭合零层为

$$
\mathcal Z_\ast
=
\mathbb N^{([\mathcal C])}.
$$

它是已经闭合并被记为静态零的的结果层。

二者不能混同：

$$
\mathcal Z_{\rm dyn}\ne\mathcal Z_\ast.
$$

$$

\text{能运动的是动态零约束层；}
\text{已经停在零的记录进入全局闭合零层。}

$$

**第 9 步｜层级更新。**

原来的 $R/P/E/D$ 需要相应改写。把 $R$ 提升为全局闭合零层 $\mathcal Z\_\ast$，得到

$$
E_i
\longrightarrow
P_i+\mathcal Z_\ast
\longrightarrow
E_i',
$$

以及失败路径

$$
E_i\longrightarrow D_i.
$$

更准确地说：

1. 未闭合路径继续留在 $E\_i$；
2. 闭合路径 $w$ 不再延长，而是退出 $E\_i$；
3. $w$ 的精确历史写入 $P\_i$；
4. 它的循环闭合类 $[w]$ 写入全局 $\mathcal Z\_\ast$；
5. 其他开放路径继续演化；若另起新层，则新层是重新播种的分支，不是 $w$ 的继续；
6. 不允许把 $\mathcal Z\_\ast$ 单独当作所有 $E\_i$ 的全局初始化器。

**第 10 步｜这一改动的收益与代价。**

收益：

1. 零层能够保存多个闭合模式，而不是把所有零压缩成一个点。
2. “静态零”有了精确含义：没有层内跃迁，只有闭合写入。
3. 全局性和局部历史分工清楚，不自动消灭区域差异。
4. 闭合路径退出活动层，不需要人为让物理速度归零。

代价：

1. 零层需要闭合类标签和重数，不能再把它看作纯标量零。
2. 闭合分支是活动层的终止节点；若要求外部时间中的物理状态也停止，仍需非可逆吸收。
3. 全局零层若被当成本地重建源，仍会造成全局同化。
4. 局部区域 $i$ 的存在与边界仍是输入，不由零层自动产生。

---

## §1 恢复层结构

| 结构 | 内容 | 当前地位 |
|:--|:--|:--|
| `R-Z-STATIC-CLOSURE-LAYER` | 全局闭合零记录的多重集、静态层内规则与写入语义 | 条件构造 |
| `R-Z-CLOSURE-TIME-QUOTIENT` | 用时间遗忘商把闭合周期轨道压成静态零记录 | 条件构造 |
| `R-Z-CLOSURE-EXIT-RULE` | 闭合路径退出活动层，写入 $P\_i+\mathcal Z\_\ast$ | 条件构造 |
| `R-Z-PERIODIC-WASHOUT-RESEED` | 活动层按周期毁灭，并由 $P\_i+\mathcal Z\_\ast$ 重新播种 | 条件构造 |
| `R-Z-DEATH-SINK-ZERO-SUM` | $D\_i$ 是只接收未闭合路径的局部终端汇；全局零和已核验，局部零和未推出 | 条件构造与否定性边界 |
| `R-Z-ZERO-STOPPING-GAP` | 让外部时间中的物理轨迹也真正停止所需的非可逆吸收规则 | 未解选择器 |

新增输入：

1. 闭合谓词 $\mathsf{Closed}$；
2. 循环等价 $\sim$；
3. 全局闭合零层 $\mathcal Z\_\ast$ 与其重数；
4. 局部历史层 $P\_i$；
5. 闭合事件写入 $\mathcal Z\_\ast$ 的规则；
6. 本地重建只使用 $P\_i$ 的规则；
7. 闭合是活动层终止节点的规则；
8. 若要求物理停止而不是分支退出，还需非可逆吸收规则。
9. 毁灭周期 $T$ 与播种映射 $\text{Seed}$ 的定义。
10. $D\_i$ 只增加记录、不回写活动层且不参与重播的终端语义。

未被本文导出：

1. 为什么闭合模式必须有静态记录；
2. 为什么全局层是多重集而不是单点零；
3. 局部区域划从哪里来；
4. 为什么闭合必须等价于活动层终止；
5. 外部时间中真正停止所需的非可逆方向；
6. 为什么毁灭周期必须是 $T$，以及为什么只能按 $P\_i+\mathcal Z\_\ast$ 播种；
7. 为什么毁灭路径必须进入 $D\_i$ 而不是重新进入活动层；
8. 为什么局部 $D\_i$ 必须自守恒，或者是否可以允许跨区域补偿；
9. 全局零层怎样接到连续时空。

---

## §2 结论

若“零层”指已经闭合的振动模式集合，那么正确的全局对象不是单个标量 $0$，而是

$$

\mathcal Z_\ast
=
\mathbb N^{(\mathcal C/{\sim})}.

$$

它的层内动力学是静态的，但允许新闭合记录写入。这里有两步：

1. 闭合路径先退出活动层 $E\_i$；
2. 时间遗忘商再把闭合轨道压成没有时间参数的静态零记录。

精确历史必须留在局部 $P\_i$。这样，之前“$R$ 全局、$P$ 局部”的猜想获得更清楚的物理含义：

$$

\text{全局保存一切闭合结果，}
\quad
\text{局部保存每个结果怎样形成。}

$$

**后续状态｜`D222` 把单层历史改为有限记忆。**
`D222` 允许 $E\_i,P\_i,\mathcal Z\_\ast$ 都有亚层，并规定：活动亚层在局部寿命到达时全清；历史层只保留最高两层；零层压缩不能合并不同闭合类标签。局部寿命 $\tau\_i$ 可以随区域不同，因此本文的单一全局周期只是最小版本，不是必需结构。
