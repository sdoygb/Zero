# Z142 · 连续化的三步：翘曲乘积的纤维 Ricci、Gauss 缩并与面积展开

> 篇号：Z142
> 状态：正式篇（推导篇）
> 前置：Z002（连续极限的 Γ-收敛条件定理与嵌入读法）、Z033（离散泊松方程）；本篇为三步的独立复算
> 关联脚本：`ZeroGPU/src/gr_chain_continuum_check.py`、`ZeroGPU/src/a2_bridge_check.py`
> 数据：`ZeroGPU/results/gr_chain_continuum_check.json`
> 登记日期：2026-10-16

## 摘要

把"离散图 → 连续几何"的读法拆成三步，本篇对三步各自给出**独立复算**（不引用任何既有数字），7 条判据全 PASS：

$$
\begin{aligned}
&\textbf{第一步（曲率）}:\ \text{翘曲乘积的纤维 Ricci 含}\ \textbf{交叉项}\ -\tfrac{n_bR_a'R_b'}{R_aR_b};\\
&\textbf{第二步（缩并）}:\ R_{ij}{}^{ij}={}^{(3)}R\ \textbf{精确};\quad 2G_{00}={}^{(3)}R\ \textbf{一阶};\\
&\textbf{第三步（面积）}:\ \text{测地球面积展开的}\ \ell^2\ \text{系数}=\mathrm{Ric}(\dot\gamma,\dot\gamma)/6 .
\end{aligned}
\tag{0.1}
$$

**关键点（第一步的机制）**：**一个"界"的演化率出现在另一个"界"的场方程里**——即曲率非零来源于**界与界的耦合**，而不是单界的自身项。

**本篇给出作用量桥的常数核验**（定理 142.4：二次型 → Einstein–Hilbert 由热核系数 $a_2$ 给出，其中 $R$ 的系数 $1/6$ 由独立脚本核算），**并给出它的入口条件**——舞台／嵌入那一笔具名输入（§5.3）。

---

## §1 记号：翘曲乘积

取 $D=n_1+n_2+1$ 维翘曲乘积

$$
ds^2=d\rho^2+R_1(\rho)^2\,d\Omega_{n_1}^2+R_2(\rho)^2\,d\Omega_{n_2}^2 ,
\qquad n_a=\dim\Omega_{n_a}.
\tag{1.1}
$$

本文的数值例取 $n_1=n_2=2$（$D=5$）。

> **定义 142.0（"界"）.** 翘曲乘积 (1.1) 里的**每一个因子块**称为一个**界**。第 $a$ 个界带两样东西：**维数** $n_a$（$n_a$ 维单位球因子的维数）与**演化率** $R_a(\rho)$（尺度因子，其对数导数 $h_a:=R_a'/R_a$ 是"膨胀率"）。**"界与界的耦合"**即 (2.1) 的交叉项 $-(R_a'/R_a)\sum_{b\ne a}n_bR_b'/R_b$——**一个界的演化率出现在另一个界的场方程里**。（"界"不是层，也不是亚层：层是 [分层架构与层纪律](Z016_分层架构与层纪律.md) 的 L0／L1／L2，界是度规的因子块。）

---

## §2 第一步：纤维 Ricci 与交叉项

> **定理 142.1（纤维 Ricci）.** 对 (1.1)，

$$
\frac{R^{(a)}_{ii}}{g^{(a)}_{ii}}
=\frac{(n_a-1)-R_aR_a''-(n_a-1)R_a'^2}{R_a^2}
\;\underbrace{-\;\frac{R_a'}{R_a}\sum_{b\ne a}n_b\frac{R_b'}{R_b}}_{\textbf{交叉项}} ,
\tag{2.1}
$$

> 且

$$
R_{\rho\rho}=-\sum_a n_a\frac{R_a''}{R_a} .
\tag{2.2}
$$

**机器核验**（A 组，sympy 精确符号）：自算 Christoffel → Riemann → Ricci（先自检缩并约定），

| 判据 | 结果 |
|:--|:--|
| A1a (2.1) 与独立推导**逐项相同** | 差 $=0$ |
| A1b $R_2'=0$ 时交叉项消失 | 差 $=0$ |
| A1c (2.2) | 差 $=0$ |

$n_a=n_b=2$ 时交叉项 $=-\dfrac{2R_1'R_2'}{R_1R_2}$。

> **读法（机制）.** 交叉项**只由两侧的 $R'$ 组成**。故 $R_2'\equiv0$（单因子）时交叉项为零，而**自身项仍可非零**——即

$$
\boxed{\;\text{耦合}\ \textbf{既非曲率的必要条件};\quad \text{交叉项是"一个界的演化率进入另一个界的场方程"} .\;}
\tag{2.3}
$$

---

## §3 第二步：Gauss 缩并与 Hamiltonian 约束

静态弱场：$g_{00}=-(1+2\Phi)$，$g_{ij}=(1-2\Phi)\delta_{ij}$。

> **定理 142.2.** 在 $K_{ij}=0$（无外曲率）切片上，

$$
\sum_{\text{空间}}g^{ik}g^{jl}R_{ijkl}={}^{(3)}R\quad\textbf{（精确）};\qquad
2G_{00}={}^{(3)}R+O(\Phi^2)\quad\textbf{（一阶）} .
\tag{3.1}
$$

**它在 ADM 里的位置（任意叶化）.** 设 lapse $N$、shift $\beta^i$，外曲率

$$
K_{ij}=\frac{1}{2N}\bigl(\partial_th_{ij}-\nabla_i\beta_j-\nabla_j\beta_i\bigr),
$$

则法帧分量满足

$$
\boxed{\;2\,n^an^bG_{ab}={}^{(3)}R+K^2-K_{ij}K^{ij};\qquad
n^aG_{ai}=-\,D_j\bigl(K^j{}_i-\delta^j{}_iK\bigr) .\;}
$$

故 (3.1) 的 $2G_{00}={}^{(3)}R$ 是**该恒等式在 $K_{ij}=0$ 处的特例**：$K=0$ 等价于「$h$ 不随 $t$ 变且 $\beta=0$」（此时那 3 个混合分量 $n^aG_{ai}\equiv0$，是恒等式）。**一般情形那 3 个分量由 $N$ 与 $\beta^i$ 供给**——它们是 $4$ 个叶化数据（见 [Z002](Z002_从Zero条件退化出广义相对论.md) §4.4）。

**机器核验**（B 组）：

| 判据 | 结果 |
|:--|:--|
| B1a Gauss 缩并 $={}^{(3)}R$ | 差 $=0$（精确） |
| B1b $D:=2G_{00}-{}^{(3)}R$ 的 $\varepsilon$-标度 | $D/\varepsilon^2\to-39.265$（常数），$D/\varepsilon\to0$ ⟹ $D\propto\varepsilon^2$ |

**数值**（$\Phi=\varepsilon\times$解析高斯，取样点 $(0.3,0.2,0.1)$）：

| $\varepsilon$ | $D=2G_{00}-{}^{(3)}R$ | $D/\varepsilon$ | $D/\varepsilon^2$ |
|--:|--:|--:|--:|
| $10^{-2}$ | $-3.996316\times10^{-3}$ | $-3.996\times10^{-1}$ | $-39.9632$ |
| $10^{-3}$ | $-3.932773\times10^{-5}$ | $-3.933\times10^{-2}$ | $-39.3277$ |
| $10^{-4}$ | $-3.926495\times10^{-7}$ | $-3.926\times10^{-3}$ | $-39.2650$ |

> **一处必须写明的区分.** 精确的恒等式是 **$R_{ij}{}^{ij}$（空间指标的 Gauss 缩并）$={}^{(3)}R$**；而 **${}^{(3)}R=2G_{00}$ 是一阶**（Hamiltonian 约束在 $N=1$、$\beta=0$、$K=0$ 下的形式）。两者**不同阶**，不可写成同一个恒等式。

---

## §4 第三步：测地球面积展开的 $\ell^2$ 系数

> **定理 142.3.** 沿测地线的测地球面积展开

$$
\frac{\mathrm{Area}(\partial B_r)}{\mathrm{Area}(S^{n-1})}=r^{\,n-1}\Bigl(1-\frac{\mathrm{Ric}(\dot\gamma,\dot\gamma)}{6}\,r^2+O(r^4)\Bigr) ,
\tag{4.1}
$$

> 球面（Einstein）情形退化为系数 $\mathcal R/(6n)$。

**机器核验**（C 组）：

| 判据 | 结果 |
|:--|:--|
| C1a $S^n$（$n=2,3,4$）：$\ell^2$ 系数 $=\mathcal R/(6n)$ | 最大相对差 $6.77\times10^{-11}$ |
| C2a 一般情形：系数 $=\mathrm{Ric}(\dot\gamma,\dot\gamma)/6$（用 (2.2) 的 $\rho\rho$ 分量） | 例：$R(0)=1,R'(0)=0.3,R''(0)=0.1$ ⟹ $\mathrm{Ric}_{\rho\rho}=-0.2$，系数 $=-0.033333$ |

> **读法.** (4.1) 给出一条**从面积读出 Ricci**的通道——面积展开的 $\ell^2$ 系数就是 Ricci 沿测地线的分量；球面上它是 $\mathcal R/(6n)$。

---

## §5 三处接口

### 5.1 三步各自的位置

| 步 | 给什么 | 用于 |
|:--|:--|:--|
| 一（(2.1)） | 交叉项 $-\frac{R_a'}{R_a}\sum_b n_b\frac{R_b'}{R_b}$ | 曲率非零的**来源**（界间耦合） |
| 二（(3.1)） | $R_{ij}{}^{ij}={}^{(3)}R$（精确）；$2G_{00}={}^{(3)}R$（一阶） | 面积型量与 Einstein 张量的对接 |
| 三（(4.1)） | 面积展开首两项；$\ell^2$ 系数 $=\mathrm{Ric}/6$ | 从面积读曲率 |

### 5.2 与 $\Gamma$-收敛条件定理的关系

Z002 的连续极限是**条件定理**：在正则假设（周期格、$k$ 固定、标度场 $c\in C^{1,1}$、$c\ge c_{\min}>0$）下，紧性／$\Gamma$-liminf／恢复列／唯一极限**已证**，字典误差界显式且速率尖锐（$C^{1,1}$ 给一致 $O(a^2)$，间断场处该速率不成立）。其**嵌入读法**（"类 → 物理点"）是一笔**具名输入**（Z002 §7.2／§11）。**这笔输入的精确形式**是：

$$
\text{嵌入输入}=\textbf{一个局部标度场}\ c(x)\ +\ \textbf{单元形状};
\qquad
w_e=P_{k,d}(c_e)\bigl(1+O(a^2)\bigr),\qquad
P_{k,d}(c)=\sum_{m\le k}m\,W_d(m-1)\,c^{\,m},
\tag{5.3}
$$

> 其中 $W_d$ 是 $d$ 维格的游走计数。**即**：给定 $c(x)$ 与单元形状，边上权重 $w_e$ 由 (5.3) 显式给出（误差 $O(a^2)$）。

本篇的三步是在**坐标层面**把该条件定理用到的几何量显式化（纤维 Ricci、Gauss 缩并、面积展开）。

### 5.3 作用量桥：由**热核系数的 $a_2$** 给出

$\Gamma$-收敛的极限泛函是**二阶椭圆型**（(2.1) 的自身项 ＋ 交叉项即其系数）。把它接到**曲率型作用量**由**热核展开的第二系数 $a_2$** 给出：

> **定理 142.4（二次型 → Einstein–Hilbert 作用量；外部标准结果）.** 对 Laplace 型算子，Seeley–DeWitt／Gilkey 展开的第二系数为

$$
a_2=\int\sqrt{g}\ \sum_s n_s\Bigl(\frac16R-E_s\Bigr),
\tag{5.1}
$$

> 即**二次型的 $a_2$ 系数含 Einstein–Hilbert 项**（$R$ 的系数为 $\frac16$）。**这一条的来源是外部标准结果**（与 [Z002](Z002_从Zero条件退化出广义相对论.md) §3.10 引理 9 的 Lovelock 同类：嫁接的经典数学）；把格上的二次型送到该算子的二次型那一步，是 Z002 §11 第 1 项的外部条件定理。
>
> **常数 $1/6$ 的独立核验**（`ZeroGPU/src/a2_bridge_check.py`）：直接从**谱**算 heat trace、拟合 $a_1/a_0$ 与 $R/6$ 比较，**不引用**上面的公式文本。球面 $a=1,2,3$：相对差 $\le1.3\times10^{-9}$；平坦环面对照（$R=0$）：$|a_1/a_0|\le1.2\times10^{-12}$。5 条判据全 PASS。

**配套的两环**：

| 环 | 内容 | 状态 |
|:--|:--|:--|
| **误差界** | $\max_{e}\bigl\lvert\widehat W_{a,e}/P_{k,d}(c(m_e))-1\bigr\rvert\le Ca^2$（$a=1/N$）；一阶项被"关于边中点的**反射**＋**游程反转**"配对面**精确消去**，故 $O(a)\to O(a^2)$ | 显式、速率**尖锐** |
| **紧性·唯一极限** | 紧性／$\Gamma$-liminf／恢复列／唯一极限 | 在四条正则假设下成立 |

**速率表**（正则性 ↓）：$c\in C^{1,1}\Rightarrow O(a^2)$；$c\in C^{1,\alpha}\Rightarrow O(a^{1+\alpha})$；$c\in C^{0,\alpha}\Rightarrow O(a^{\alpha})$；仅连续仍 $\to0$；**$c$ 间断 $\Rightarrow O(1)$**。

$$
\boxed{\;\text{故}:\quad
\underbrace{\text{二次型}}_{\text{已显式（本三步）}}\ +\ \underbrace{a_2\ \text{（定理 142.4，外部标准结果；常数 }1/6\text{ 已独立核验）}}_{\text{曲率项的来源}}\ +\ \underbrace{\text{误差界／紧性}}_{\text{显式}} .\;}
\tag{5.2}
$$

**入口条件**：三环的前提是"**周期立方环面 $\mathbb T^d$ ＋ 格距 $a=1/N$**"，即**舞台／嵌入这一笔具名输入**（Z002 §11、[Z145](Z145_内容图与格点：有限维几何可长性的判据与两族的对照.md)）。故桥本身已建，**它的入口是那一笔输入**。

---

## §6 结论的适用范围

1. **第一步的精确性对一切 $R_a(\rho)>0$ 成立**；$n_a=n_b=2$ 的数值例是 $D=5$；
2. **第二步的"精确"部分限 $K_{ij}=0$ 切片**；${}^{(3)}R=2G_{00}$ 只在 $\Phi$ 的一阶成立（(3.1) 分列两阶）；一般叶化见该节的 ADM 形式（Hamiltonian 含 $K^2-K_{ij}K^{ij}$，混合分量含 $K$ 的散度）；
3. **第三步的 $\mathcal R/(6n)$ 仅对 Einstein（球面）情形**；一般情形用 $\mathrm{Ric}(\dot\gamma,\dot\gamma)/6$；
4. **C1a 的相对差 $6.77\times10^{-11}$ 是数值核验**（$n=2,3,4$），不是解析界；
5. **§5.3 的作用量桥由定理 142.4 给出**，其来源是外部标准结果、$R$ 系数 $1/6$ 由 `a2_bridge_check.py` 独立核验；入口条件是舞台／嵌入这一笔具名输入。

---

## §7 结论

$$
\begin{aligned}
&\textbf{曲率的来源}:\ \text{交叉项}\ -\tfrac{R_a'}{R_a}\sum_{b\ne a}n_b\tfrac{R_b'}{R_b};\quad \text{单侧坍缩时它为零};\quad\textbf{已证（符号）}.\\
&\textbf{缩并}:\ R_{ij}{}^{ij}={}^{(3)}R\ \textbf{精确};\ 2G_{00}={}^{(3)}R\ \textbf{一阶};\quad\textbf{已证（符号＋标度）}.\\
&\textbf{面积}:\ \ell^2\ \text{系数}=\mathrm{Ric}(\dot\gamma,\dot\gamma)/6;\quad \text{球面}\ \mathcal R/(6n);\quad\textbf{已证（含数值核验）}.\\
&\textbf{作用量桥}:\ \text{二次型}\to\text{Einstein--Hilbert 由 }a_2\text{ 给出（定理 142.4，外部标准结果；}R\text{ 系数 }1/6\text{ 经独立核验）}.
\end{aligned}
$$

**一句话：**

$$
\boxed{\;\text{曲率来自界与界的耦合；面积展开的}\ \ell^2\ \text{系数即 Ricci；二次型到 Einstein–Hilbert 的作用量桥由 }a_2\ \text{给出}.\;}
$$

---

## 附录 A　核验

`gr_chain_continuum_check.py`（纯标准库 ＋ `numpy`／`sympy`；不 import 项目脚本、不写文件）：

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | 翘曲乘积纤维 Ricci、交叉项的消失条件、$R_{\rho\rho}$ | 3 |
| B | Gauss 缩并精确性、$2G_{00}-{}^{(3)}R$ 的 $\varepsilon$-标度 | 2 |
| C | $S^n$ 的 $\ell^2$ 系数、一般情形的 $\mathrm{Ric}/6$ | 2 |
| D | $a_2$ 的 $R$ 系数 $=1/6$（`a2_bridge_check.py`：球面 $a=1,2,3$ ＋平坦环面对照） | 5 |
| E | ADM Hamiltonian 约束（含 $K\ne0$）＋动量约束（$K=0$ 恒满足／$K\ne0$ 时 $-2\epsilon\partial_t\partial_x\Phi$） | 2 |
| | **合计** | **14** |

**核验纪律**：数值只算旁证；A 组的比较是 sympy **精确符号**（差 $=0$），B1a 亦为精确；B1b／C1a／D 组是**数值**（标度、相对差、拟合系数）。D 组从**谱**算 heat trace，与 §5.3 引用的公式文本相互独立。

复跑：

```bash
cd ZeroGPU/src && /usr/bin/python3 gr_chain_continuum_check.py   # A/B/C 组 7 条
cd ZeroGPU/src && /usr/bin/python3 a2_bridge_check.py            # D 组 5 条
cd _staging/GR缺口重评 && /usr/bin/python3 a2_adm_completion_check.py   # E 组（Hamiltonian）
cd _staging/GR缺口重评 && /usr/bin/python3 a2_momentum_check.py         # E 组（动量）
```

---

*Zero 项目 · Z142 · 完*
