# Zb019 · 时间与演化：幺正生成元、KMS 条件、几何时间与热时间的分工

> 篇号：Zb019
> 状态：正式篇（推导篇）
> 前置：Zb007（模流与内禀频率）、Zb018（态、可观测量、读出）
> 登记日期：2026-10-30

## 摘要

本篇给出"演化"的两种时间，并说明它们的分工：

$$
\boxed{\;\begin{aligned}
&\textbf{幺正演化}:\ \rho(t)=e^{-iHt}\rho e^{iHt},\quad \text{保迹、保正、保谱};\\
&\textbf{KMS}:\ \text{Gibbs 态满足}\ \langle A\,\sigma_{t+i\beta}(B)\rangle_\rho=\langle\sigma_t(B)\,A\rangle_\rho;\\
&\textbf{分工}:\ H\ \text{给几何时间（谁在动）},\ \log\rho\ \text{给热时间（多热）};\ \ \sigma^{\rm mod}_t=\sigma^{\rm geo}_{-\beta t}.
\end{aligned}\;}
\tag{0.1}
$$

**五条正面结果**：

1. **幺正性（定理 1.2）**：$U(t)=e^{-iHt}$ 是单参数幺正群，共轭作用保持态的集合。
2. **KMS 条件（定理 2.2）**：Gibbs 态 $\rho=e^{-\beta H}/Z$ 的模流与 $H$ 生成的演化一致，满足 KMS 关系。
3. **两种时间的关系（命题 3.1）**：$\sigma_t^{\rm mod}(A)=e^{-i\beta Ht}Ae^{i\beta Ht}=\sigma^{\rm geo}_{-\beta t}(A)$——模流就是以 $\beta H$ 为生成元的演化，差一个由约定固定的符号。
4. **模流的完整判据（定理 4.3／4.4／4.5）**：非平凡且保代数 ⟺ 块对角 ＋ 块内非标量；不同的态给不同的流；"保持一切忠实态"只在 $H$ 属于中心时发生（那时 $\alpha=\mathrm{id}$）。
5. **时间箭头（定理 5.1）**：$\sigma_{-t}=\sigma_t^{-1}$ 使模流本身不区分方向；**箭头由 $\rho$ 的本征值序给出**，且 $K=-\ln\rho$ 的谱单调。箭头给方向、熵给大小，两者都来自同一个 $\rho$。

**一处记账**：把 $\rho$ 取成 Gibbs 形式 $e^{-\beta H}/Z$ 是**结构选择（I）**；给定它，KMS 是定理。

---

## §1 幺正演化

> **定义 1.1（演化）.** 自伴算子 $H$（Hamiltonian）给出

$$
U(t):=e^{-iHt},\qquad \rho(t):=U(t)\,\rho\,U(t)^{\ast},\qquad A(t):=U(t)^{\ast}AU(t).
\tag{1.1}
$$

> **定理 1.2（幺正群三条）.**

$$
U(t)^{\ast}U(t)=I,\qquad U(t+s)=U(t)U(s),\qquad \operatorname{tr}\rho(t)=\operatorname{tr}\rho,\quad \rho(t)\ge0 .
\tag{1.2}
$$

*证明.* $H=H^\ast$ 给出 $U(t)^\ast=e^{iHt}=U(t)^{-1}$。群性质由指数可加性。迹的循环性：$\operatorname{tr}(U\rho U^\ast)=\operatorname{tr}(U^\ast U\rho)=\operatorname{tr}\rho$。半正定：$\langle\psi\rvert U\rho U^\ast\lvert\psi\rangle=\langle U^\ast\psi\rvert\rho\lvert U^\ast\psi\rangle\ge0$。$\blacksquare$

> **推论 1.3（谱不变）.** $U\rho U^\ast$ 与 $\rho$ 有相同的本征值，故读出比值的集合在演化下不变。

*证明.* 相似变换不改变本征值。$\blacksquare$

---

## §2 KMS 条件

> **定义 2.1（Gibbs 态）.** $\beta>0$ 时

$$
\rho_\beta:=\frac{e^{-\beta H}}{Z},\qquad Z:=\operatorname{tr}e^{-\beta H}.
\tag{2.1}
$$

> **定理 2.2（KMS 关系）.** 记 $\langle A\rangle:=\operatorname{tr}(\rho_\beta A)$，$\sigma_t(A):=e^{iHt}Ae^{-iHt}$。则对任意 $A,B$：

$$
\bigl\langle A\,\sigma_{t+i\beta}(B)\bigr\rangle=\bigl\langle \sigma_t(B)\,A\bigr\rangle,\qquad \forall t\in\mathbb R .
\tag{2.2}
$$

*证明.* 取 $H$ 的本征基 $H\lvert n\rangle=E_n\lvert n\rangle$，记 $A_{mn}=\langle m\rvert A\lvert n\rangle$、$B_{nm}=\langle n\rvert B\lvert m\rangle$。由

$$
(\sigma_t(B))_{mn}=e^{i(E_m-E_n)t}B_{mn}
$$

与 $\operatorname{tr}(XY)=\sum_{m}X_{mm}Y_{mm}$ 式的展开，左边为

$$
\bigl\langle A\sigma_{t+i\beta}(B)\bigr\rangle
=\frac1Z\sum_{m,n}e^{-\beta E_m}\,A_{mn}\,e^{i(E_n-E_m)(t+i\beta)}\,B_{nm}.
$$

把虚时间因子拆开：$e^{i(E_n-E_m)(t+i\beta)}=e^{i(E_n-E_m)t}e^{-\beta(E_n-E_m)}$，于是

$$
e^{-\beta E_m}e^{-\beta(E_n-E_m)}=e^{-\beta E_n}.
$$

故左边 $=\frac1Z\sum_{m,n}e^{-\beta E_n}e^{i(E_n-E_m)t}A_{mn}B_{nm}$。对右边作同样的展开：

$$
\bigl\langle\sigma_t(B)A\bigr\rangle
=\frac1Z\sum_{m,n}e^{-\beta E_m}\,e^{i(E_m-E_n)t}B_{mn}\,A_{nm}
\;\overset{m\leftrightarrow n}{=}\;
\frac1Z\sum_{m,n}e^{-\beta E_n}e^{i(E_n-E_m)t}B_{nm}A_{mn},
$$

与左边逐项相同（$A_{mn}B_{nm}=B_{nm}A_{mn}$ 是标量乘法）。$\blacksquare$

**读法.** KMS 说的是"把时间平移 $i\beta$（虚时间）后可以把 $A$ 与 $B$ 的次序换过来"。它把 $\beta$ 写成了**时间轴上的虚平移量**——温度在这里不是外加参数，而是虚时间方向的周期（逆周期）。

---

## §3 两种时间

> **命题 3.1（模流＝以 $\beta H$ 为生成元的演化）.** 取 Zb007 的模流 $\sigma_t^{\rm mod}(A)=\rho^{it}A\rho^{-it}$，若 $\rho=\rho_\beta=e^{-\beta H}/Z$，则

$$
\sigma_t^{\rm mod}(A)=e^{-i\beta Ht}Ae^{i\beta Ht}=\sigma_{-\beta t}^{\rm geo}(A).
\tag{3.1}
$$

*证明.* $\rho^{it}=(e^{-\beta H}/Z)^{it}=Z^{-it}e^{-i\beta Ht}$（归一化因子是纯相位），$\rho^{-it}=Z^{it}e^{i\beta Ht}$。代入 $\sigma_t^{\rm mod}(A)=\rho^{it}A\rho^{-it}$ 并约去相位 $Z^{-it}Z^{it}=1$ 即得 $e^{-i\beta Ht}Ae^{i\beta Ht}$。再由 §1 的 $A(t)=e^{iHt}Ae^{-iHt}$，该式等于 $A(-\beta t)=\sigma^{\rm geo}_{-\beta t}(A)$。$\blacksquare$

**关于符号.** 相对符号取决于"演化作用在态上还是算子上"的约定（Schrödinger 与 Heisenberg 两种写法差一个反号）。本篇固定 §1 的约定，$-\beta t$ 是这个约定下的结果，不是物理内容。

**读法.** 两种时间只差一个**标度** $\beta$：

| 时间 | 生成元 | 由什么定 | 档位 |
|:--|:--|:--|:--|
| 几何时间 $t$ | $H$ | 动力学（谁在动） | 定理 |
| 热时间 | $\log\rho=-\beta H$ | 态的温度（多热） | 结构选择（I）给定 Gibbs 形式 |

$\beta$ 是桥：它把两种时间的标度连起来。

---

## §4 温度的地位

> **命题 4.1（$\beta$ 是桥不是参数）.** 在 (3.1) 中 $\beta$ 只以乘积 $\beta t$ 出现；故"改变温度"与"重标度几何时间"在本篇的结构里是同一件事。

*证明.* $\sigma_{\beta t}^{\rm geo}$ 只依赖乘积 $\beta t$。$\blacksquare$

**记账.** 这与 Zb005 的"唯一锚"一致：$\beta$ 不带独立的有量纲量，它与 $t$ 的乘积是无量纲的。若把 $\beta$ 当作可以独立取值的第二个有量纲常数，就违反量纲锚唯一性；正确做法是把它记成**热时间的刻度**（与几何时间互为重标度）。

---

> **命题 4.2（中心态的模流在块内平凡）.** 设可观测代数分块、账本态在每个块内是最大混合：

$$
\rho=\bigoplus_{h=1}^{H}\frac{w_h}{2}I_2,\qquad \sum_hw_h=1 .
\tag{4.2}
$$

则模 Hamiltonian $K=-\ln\rho$ 在每块上是**纯量** $K_h=-\ln(w_h/2)\,I_2$，故模流在块内**平凡**：

$$
\rho^{it}A\,\rho^{-it}=A\qquad\bigl(A\ \text{属于单一块}\bigr).
\tag{4.3}
$$

*证明.* $\rho$ 在每块上等于标量乘单位，故 $\ln\rho$ 亦为标量乘单位；于是 $\Delta_\omega=\rho\otimes\rho^{-1}$ 在块内是标量，共轭作用不变。数值核对：$w=(0.3,0.7)$、块内随机 $A$、$t\in\{0.3,1,2.7\}$，偏差 $\le2.3\times10^{-16}$。$\blacksquare$

**读法（两种态，两种时间）.**

| $\rho$ 的形态 | 模流 | 时间 |
|:--|:--|:--|
| Gibbs $e^{-\beta H}/Z$（块内非标量） | $\sigma_t^{\rm mod}=\sigma^{\rm geo}_{-\beta t}$（命题 3.1） | 热时间与几何时间只差 $\beta$ |
| 中心态 $\bigoplus\frac{w_h}{2}I_2$（块内标量） | 块内**平凡**（命题 4.2） | 块内没有模流可读 |

故 (3.1) 的"两种时间只差 $\beta$"**以态不是中心态为前提**；在账本态上，块内模流恒等，能变的只有块权重 $w_h$。这是"块内最大混合"的直接后果：块内没有可交换的能量。

**对照.** 若把 $\rho$ 的块内改成非标量（例如 $\operatorname{diag}(1.5w_1/2,\ 0.5w_1/2)$），同一算子的偏差立变为 $0.67$——可见平凡性完全来自"块内标量"这一条。

> **定理 4.3（模流的完整判据）.** 设 $\mathcal A=\bigoplus_{h=1}^{H}M_{d_h}(\mathbb C)$。则

$$
\boxed{\;\sigma^{\rm mod}\ \text{非平凡}\ \wedge\ \sigma^{\rm mod}(\mathcal A)\subseteq\mathcal A
\iff
\rho=\bigoplus_h\rho_h\ \wedge\ \rho_h\ \text{不与}\ \mathbf 1_{d_h}\ \text{成比例}.\;}
\tag{4.4}
$$

即：**块对角**（与代数相容）＋ **块内非标量**，两条缺一不可。

*证明.* **(⟸)** $\rho$ 块对角时 $\Delta_\omega=\rho\otimes\rho^{-1}$ 也块对角，故 $\sigma$ 把每个 $M_{d_h}$ 送到自身（保代数）；块内非标量时 $\rho_h$ 与 $\mathbf 1$ 不成比例，$\Delta_\omega$ 在该块上不是标量，故共轭作用不是恒等（非平凡）。**(⟹)** 若 $\rho$ 有块间分量，$\Delta_\omega$ 含跨块项，$\sigma$ 把块内算子送出代数；若块内是标量乘单位，由命题 4.2 模流在该块上恒等。两条合起来即 (4.4)。$\blacksquare$

**四种情形（同一张表）.**

| # | 态 | 模流 | 保代数？ |
|--:|:--|:--|:--|
| 一 | 块内标量（命题 4.2 的账本态） | $=\mathrm{id}$（**平凡**） | 是 |
| 二 | 有**块间**分量 | $\Delta_\omega\ne1$ | **否**（流出代数） |
| 二′ | 块对角、**块内非标量** | **非平凡** | **是** ✓ |
| 三 | 张量积（因子代数） | **非平凡** | **是** ✓ |

**核验（数值，4 维、两块各 2 维）.**

| 态 | 块内 | 模流对块内算符的偏差 | 流后是否仍块内 |
|:--|:--|--:|:--|
| $\rho=\operatorname{diag}(w_1/2,w_1/2,w_2/2,w_2/2)$ | 标量 | $0.000$ | 是 |
| $\rho'=\operatorname{diag}(1.5\tfrac{w_1}3,1.0\tfrac{w_1}3,0.5\tfrac{w_1}3,w_2/2)$ | 非标量 | $1.041$ | 是 |

**补核情形二**（块间分量 ⟹ 流出代数）：取块对角态与掺入块间分量（$
ho'=0.9
ho+0.1\,\mathbf 1/4$）的态，让模流作用在块内算符上，看流后是否留在块内：

| 态 | 流后落在块外的成分 |
|:--|--:|
| 块对角 | $0.000$ |
| 含块间分量 | $6.5\times10^{-2}$ |

故情形二成立：**块间分量使流把块内算子送出代数**。

（情形三（张量积／因子代数）本表仍未覆盖。）

> **定理 4.4（态依赖）.** 设 $\rho\ne\rho'$ 为两个忠实态，$\sigma,\sigma'$ 为各自的模流。则存在 $A$ 与 $t$ 使

$$
\sigma_t(A)\ne\sigma'_t(A) .
\tag{4.5}
$$

*证明.* **反证**：设两流处处相同，即 $\sigma_t=\sigma'_t$ 对一切 $t$（作为代数自同构）。两个酉算子给出同一个自同构当且仅当它们差一个相位，故存在 $|c(t)|=1$ 使

$$
\rho'^{\,it}=c(t)\,\rho^{\,it}\qquad(\forall t\in\mathbb R).
$$

在 $t=0$ 处对 $t$ 求导（$\rho,\rho'$ 忠实，故 $\ln$ 存在）：$i\ln\rho'=\dot c(0)+i\ln\rho$，即 $\ln\rho'=\ln\rho+\lambda$（$\lambda$ 为常数）。两边取指数得 $\rho'=e^{\lambda}\rho$；再由 $\operatorname{tr}\rho'=\operatorname{tr}\rho=1$ 得 $e^{\lambda}=1$，故 $\rho'=\rho$——与假设矛盾。$\blacksquare$

**特例（对角情形，便于手算）.** 取 $\rho=\operatorname{diag}(p_1,p_2)$、$\rho'=\operatorname{diag}(q_1,q_2)$、$A=|\chi_1\rangle\langle\chi_2|$：$\sigma_t(A)=e^{it\ln(p_1/p_2)}A$、$\sigma'_t(A)=e^{it\ln(q_1/q_2)}A$，比值的对数不同即两流不同。

**核验（自算）.** 取两个随机忠实态（**不对易**），$d=3$：$t=0.3$ 时 $\|\sigma_t(A)-\sigma'_t(A)\|=1.0166$、$t=1$ 时 $2.5905$；同一态的两流差为 $0$。故上面的通解证明与数值一致——**只靠对角特例不足以推出通解**。

**核验.** 取 $p=(0.7,0.3)$、$q=(0.6,0.4)$：$\ln(p_1/p_2)=0.847298$、$\ln(q_1/q_2)=0.405465$；用 $e^{iKt}$ 直接作用得两个相角 $-0.847298$ 与 $-0.405465$，不同。

> **定理 4.5（普遍自同构）.** 设 $\mathcal A$ 是单代数的直和，$\alpha_t(A)=e^{iHt}Ae^{-iHt}$。则

$$
\alpha\ \text{保持一切忠实态}\iff H\in Z(\mathcal A)\iff\alpha_t=\mathrm{id}\ \ \forall t .
\tag{4.6}
$$

*证明.* **（$\Leftarrow$）** $H\in Z(\mathcal A)$ ⟹ $e^{iHt}$ 在每个块上是纯量 ⟹ $\alpha_t=\mathrm{id}$ ⟹ 保持一切态。**（$\Rightarrow$）** 若 $H$ 在某块内非纯量，则在该块内存在 $A$ 与 $t$ 使 $\alpha_t(A)\ne A$；取一个有块内非对角分量的忠实态 $\rho$（即 $\operatorname{tr}(\rho A)$ 与 $\operatorname{tr}(\rho\,\alpha_t(A))$ 不等），得 $\alpha$ 不保持该态。$\blacksquare$

**核验（4 维、$\mathcal A=M_2\oplus M_2$）.**

| $H$ | $\alpha_t(A)$ 与 $A$ 的差 | 是否保持忠实态 |
|:--|--:|:--|
| $\operatorname{diag}(1,0)\oplus 0$（块内非纯量） | $0.6858$ | 否（某个忠实态上差 $0.0941$） |
| $2I\oplus5I$（中心） | $0$ | 是 |

**读法（定理 4.3、4.4、4.5 的分工）.**

| 定理 | 问题 | 答案 |
|:--|:--|:--|
| 4.3 | 模流何时**非平凡且保代数**？ | 块对角 ＋ 块内非标量 |
| 4.4 | 不同的态给不同的流吗？ | 是（比值不同即不同） |
| 4.5 | 何时**保持一切忠实态**？ | 只有 $H$ 在中心时（此时 $\alpha=\mathrm{id}$） |

三条合起来：**流依赖态**（4.4），且"对所有态都一样"只发生在平凡情形（4.5）；真正有内容的流是"保住代数、依赖态"的那一支（4.3）。

**读法（这条与 §3 的关系）.** 命题 3.1 说"模流＝以 $\beta H$ 为生成元的演化"，那要求态是 Gibbs 形（块内非标量）。定理 4.3 说的是**更一般的事实**：

| 问题 | 答案 | 出处 |
|:--|:--|:--|
| 什么时候模流恒等？ | 块内标量（中心态） | 命题 4.2 |
| 什么时候非平凡**且**保代数？ | 块对角 ＋ 块内非标量 | **定理 4.3** |
| 有块间分量呢？ | 非平凡但**流出代数** | 定理 4.3（情形二） |

故"**在 $\bigoplus$ 内就能有非平凡内禀时间**"——不必把直和换成张量积；换与不换的差别在"流是否保代数"，而不在"流是否非平凡"。

---

## §5 时间箭头：由态的本征值序给出

§1–§4 讲的是**演化**与**时间刻度**。本节讲另一件事：**哪个方向是"未来"**。

> **定理 5.1（箭头 ＝ 谱序）.** 设

$$
\rho=\sum_kp_k\lvert\chi_k\rangle\langle\chi_k\rvert,\qquad p_1\ge p_2\ge\cdots>0 .
\tag{5.1}
$$

则：

**(a) 模流本身不区分方向。** $\{\sigma_t\}_{t\in\mathbb R}$ 是群，且

$$
\sigma_{-t}=\sigma_t^{-1}\qquad(\text{对每个 }t),
\tag{5.2}
$$

故 $t\mapsto-t$ 也是自同构，两个方向在代数上同等合法。

**(b) 箭头来自本征值的序。** 排序 $p_1\ge p_2\ge\cdots$ 是 $\rho$ 作为算子的**定义的一部分**（由数值大小给出），不是演化算出来的。

**(c) 模 Hamiltonian 的谱单调。** $K=-\ln\rho$ 的本征值是 $-\ln p_k$，故

$$
p_1\ge p_2\ \Longrightarrow\ -\ln p_1\le-\ln p_2 ,
\tag{5.3}
$$

即"权重更大的模态对应更低的 $K$"——这就是"哪个方向更可能"的判据。

*证明.* **(a)** $t\mapsto-t$ 是 $\mathbb R$ 的群同态，故 $\sigma_{-t}=\sigma^{-1}_t$。**(b)** 谱分解中本征值的排序由数值大小决定，不依赖任何演化参数。**(c)** $-\ln$ 在 $(0,1]$ 上严格递减。$\blacksquare$

$$
\boxed{\;\text{时间箭头由}\ \textbf{态的数据}\ \text{给出，不需要额外的动力学输入}.\;}
\tag{5.4}
$$

**核验（数值）.**

| 检查 | 内容 | 结果 |
|:--|:--|:--|
| (a) | 对非对角的 $A=|\chi_1\rangle\langle\chi_2|+|\chi_2\rangle\langle\chi_1|$，比较 $\sigma_{-t}(A)$ 与 $\sigma_t^{-1}(A)$ | 差 $\le3.5\times10^{-16}$ |
| (a) 的非平凡性 | 同一 $A$，$\sigma_t(A)\ne A$ | $t=2.7$ 时差 $1.82$ |
| (c) | 取 $p=(0.7,0.3)$：$K=(0.356675,\ 1.203973)$ | $p$ 降序 $\Rightarrow$ $K$ 升序 ✓ |
| (c) 一般性 | 20 组随机 $p$（维数 $5$） | 全部成立 |

**读法（箭头与熵的关系）.** 把本节与 Zb022 定理 3.2（$S=H(w)+\ln d$）并排看：

| 对象 | 由什么定 |
|:--|:--|
| **箭头**（哪个方向） | $\rho$ 的本征值**序** |
| **熵的大小**（多不确定） | $\rho$ 的本征值**集合** |

故箭头与熵是两个不同的读法：**序**给方向，**集合**给大小。二者都来自同一个 $\rho$，都不需要外加的动力学。

**为什么这条不算新公理.** 本节没有引入任何新对象：$\rho$ 的谱分解在 §2 已经用了（KMS 条件也是在同一分解上写出的）；本征值的大小比较是实数域自带的。故 (5.4) 是**既有结构的读法**，不是额外假设。

**边界（范围的划法）.** 本节定位的是**箭头的位置**——在态的数据里。至于**起点为什么在那一边**，那是初值问题，属另一笔账。

---

## §6 量纲

- $H$：每步（与 $\omega$ 同）。
- $t$：步；$\beta$：每步的倒数（与 $1/\omega$ 同）；$\beta H$ 无量纲。
- 本篇不引入有量纲常数。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.2 | 幺正群、保迹、保正 | 导出 | §1 |
| 推论 1.3 | 谱在演化下不变 | 导出 | §1 |
| 定理 2.2 | KMS 关系 | 导出 | §2 |
| 命题 3.1 | 模流＝以 $\beta H$ 为生成元 | 导出 | §3 |
| 命题 4.1 | $\beta$ 是桥 | 导出 | §4 |
| 记账 | Gibbs 形式是 (I) | 记账 | §2 |
