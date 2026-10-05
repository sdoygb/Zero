# Zero 系列 · 繁殖转移定理：五条规则、谱半径与幂零灭绝

**日期**：本轮（由程序写成的文章） · **来源程序**：[`zero_sum_reproduction_audit.py`](zero_sum_reproduction_audit.py)
**性质**：把程序里的**线性代数定理**写成文章；并给出 [`G90`](G90_zero_series_reference_triage.md) 判为【打问号】的那个对象（"零和保持的繁殖规则"）的**部分答案**。
**等级标签**：【定义】/【定理】/【数值核验】/【no-go】/【开放】。
**核验**：[`zero_sum_reproduction_transition_theorems_check.py`](zero_sum_reproduction_transition_theorems_check.py) —— **独立实断言 30 / 结论行 0 / 不符 0**，退出码 `0`
（全部数字取自该程序的结果 JSON（`matrix_audit` / `sequences`），本版逐项复读并制表）

---

## §0 本文写什么

程序 `zero_sum_reproduction_audit.py` 的自我定位是**审计**：先前沙盒用了 $M=1+r$、$M=2^r$ 这类乘数，
把**"可能的描述数／分支程序数"**当成了**"新物理圈数"**。本文把这次审计的**定理内容**写出来：

$$
\ \text{闭合}\ \Rightarrow\ \text{持续};\qquad \text{持续}\ \ne\ \text{繁殖};\qquad \text{只有显式复制给指数增长。}\ 
$$

---

## §1 定义：五条繁殖规则

模式（modes）$=$ 周期 $\le 12$ 的全部循环平衡词（旋转类）。对父模式 $p$，程序给出五条**子代生成规则**：

| 规则 | 含义 |
|:--|:--|
| `persist` | 只保留父（不生成子） |
| `copy` | 显式**复制**父 |
| `split_parent` | 在内部归零点**切分**，父保留 |
| `split_only` | 切分，但**移除父** |
| `split_all` | 切分（全部切点） |

---

## §2 定理 T1（繁殖转移矩阵）

**定义 D1**：在模式集上定义**繁殖转移矩阵**

$$
T_{c,p}\;:=\;\#\{\ \text{父 }p\text{ 的某条分支中所含的 }c\ \},
$$

即"父 $p$ 每代产生多少个 $c$"（程序 `transition_matrix`）。种群向量按

$$
v_{k+1}=T\,v_k,\qquad N_k=\mathbf 1^{\!\top}T^{k}v_0
$$

演化（程序 `sequence`，18 代）。

---

## §3 定理 T2（谱半径与渐近行为，五规则全表）

**核验（逐项取自结果 JSON）**：

**两个不同的指数不可混为一谈**（本版更正：首版把它们并成一列）：

- $\rho(T)=1$ 时 $T$ 是**幂单**（$T=I+N$，$N$ 幂零）⇒ 报的是 $\;I-T$ 的幂零指数；
- $\rho(T)=0$ 时 $T$ **本身幂零** ⇒ 报的是 $\;T$ 的幂零指数；
- $\rho(T)=2$（`copy`）两者皆非 ⇒ **不报指数**。

| 规则 | 谱半径 $\rho(T)$ | 零特征值数 | $I-T$ 的幂零指数（幂单时） | $T$ 的幂零指数 | 观测到的序列行为 |
|:--|--:|--:|--:|--:|:--|
| `persist` | $1.0$ | 0 | 1 | —（非幂零） | **constant**（常数） |
| `copy` | $\mathbf{2.0}$ | 0 | — | — | **exponential**（指数） |
| `split_parent` | $1.0$ | 0 | 6 | —（非幂零） | constant ／ polynomial |
| `split_only` | $\mathbf{0.0}$ | **123** | — | **6** | **finite extinction**（有限灭绝） |
| `split_all` | $1.0$ | 0 | 6 | —（非幂零） | constant ／ polynomial |

$$
\ \text{五条规则中，}\textbf{只有显式复制（copy）给指数增长};\ \text{其余至多常数／多项式。}\ 
$$

---

## §4 定理 T3（幂零 $\iff$ 有限灭绝）

`split_only` 的 $T$ 满足

$$
T^{6}=0,\qquad \rho(T)=0,\qquad \#\{\text{零特征值}\}=123 ,
$$

故 $N\_k=\mathbf 1^\top T^k v\_0$ 在**有限代内归零**——程序把这一指数直接算作 `nilpotent_index`：

$$
\ \rho(T)=0\ \Longleftrightarrow\ T\ \text{幂零}\ \Longleftrightarrow\ \text{谱系在有限代内灭绝};\quad\text{灭绝时间}=\text{幂零指数}.
$$

**这条是**"闭合不等于永续"**的定量形式**：切分而**不复制**，父一撤，谱系就死。

---

## §5 定理 T4（$M=1+r$ 与 $M=2^r$ **不是**物理后代数）【no-go】

程序把先前被混为一谈的四个对象**分开**（其自述）：**closure（闭合）／persistence（持续）／reproduction（繁殖）／branch-program count（分支程序计数）**。审计结论逐条：

| 结论（程序自述，逐字） | 汉译 |
|:--|:--|
| *"among the five tested rules, only the explicit copy rule gives exponential growth; multi-type autocatalytic cycles remain a separate possibility"* | 五条规则里只有显式复制给指数增长；**多类型自催化环仍是另一种可能**（未测） |
| *"closure alone gives persistence, not reproduction"* | **闭合本身只给持续，不给繁殖** |
| *"splitting into non-reproducing subloops gives at most polynomial accumulation; if the parent is removed, the process can die in finite generations"* | 切成不复制的子圈至多给多项式积累；父一撤可有限代灭绝 |
| *"M=1+r and M=2\*\*r are branch-program counts. They are not automatically physical offspring counts and do not by themselves imply exponential lineage growth"* | $M=1+r$、$M=2^r$ 是**分支程序计数**，不自动等于物理后代数，也不自动蕴含谱系指数增长 |

$$
\Longrightarrow\ \textbf{先前的 } M=1+r\ /\ 2^{r}\ \textbf{含一个隐藏的"复制假设"（本版判为 no-go）。}
$$

---

## §6 与体系的关系

| 本文对象 | 在体系里的位置 |
|:--|:--|
| `persist` $\Rightarrow$ 常数 | 对应 $Z0$ 的"持续"款；[`zero_sum_persistence_theorems.md`](zero_sum_persistence_theorems.md) |
| `copy` $\Rightarrow$ 指数、$\rho=2$ | 这正是 [`G90`](G90_zero_series_reference_triage.md) 判【打问号】的那个**待造对象**：**零和保持的繁殖／自催化规则**——本文给出它的**五规则分类**，但仍**没有**给出物理形态 |
| `split_only` $\Rightarrow$ 幂零灭绝 | 与 [`zero_sum_closure_exit.md`](zero_sum_closure_exit.md) 的"活动层熄灭"同型：**不补充就死** |
| $M=1+r$／$2^r$ 的 no-go | 修正 [`zero_sum_rotation_class_algebra.md`](zero_sum_rotation_class_algebra.md) §4 的重数表：那三套**只是代数**，不是繁殖律 |

---

## §7 诚实边界与开放项

**边界**

1. 五条规则是**被测规则的枚举**，不是"全部可能的繁殖律"；程序自己写明**多类型自催化环未测**。
2. 谱半径／幂零指数都在**有限模式集（周期 $\le12$）**上算，**未证**对一般周期成立；$123$ 个零特征值也是该规模下的数。
2b. **一处程序侧隐患（本版查出）**：幂零判定用 `int64` 累乘。若误把它用于 $\rho=2$ 的 `copy` 规则，$2^k$ 会**溢出归零**，从而**假报**"幂零指数 $=64$"。程序本体只在 $\rho\in\{0,1\}$ 时调用它，故未触发；但任何复用该函数的地方必须带此警告。
3. 18 代的序列分类用**启发式阈值**（比值 $>1.5$ 判指数），不是渐近证明。

**开放项**

1. **多类型自催化环**：程序点名但未测——这可能是"繁殖"的真正形态。
2. **一般周期的幂零指数**是否有闭式（现只知 $T\le12$ 时 `split_only` 给 $6$）。
3. **$\rho(T)$ 与周期的关系**：五个 $\rho$ 里只有 `copy` 给 $2$，其余为 $1$ 或 $0$——是否存在一般判据？
