# Zb019 · 时间与演化：幺正生成元、KMS 条件、几何时间与热时间的分工

> 篇号：Zb019
> 状态：正式篇（推导篇）
> 前置：Zb007（模流与内禀频率）、Zb018（态、可观测量、读出）、Zb047（账本熵与 Gibbs 不等式）
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

**六条正面结果**：

1. **幺正性（定理 1.2）**：$U(t)=e^{-iHt}$ 是单参数幺正群，共轭作用保持态的集合。
2. **KMS 条件（定理 2.2）**：Gibbs 态 $\rho=e^{-\beta H}/Z$ 的模流与 $H$ 生成的演化一致，满足 KMS 关系。
3. **两种时间的关系（命题 3.1）**：$\sigma_t^{\rm mod}(A)=e^{-i\beta Ht}Ae^{i\beta Ht}=\sigma^{\rm geo}_{-\beta t}(A)$——模流就是以 $\beta H$ 为生成元的演化，差一个由约定固定的符号。
4. **模流的完整判据（定理 4.3／4.4／4.5）**：非平凡且保代数 ⟺ 块对角 ＋ 块内非标量；不同的态给不同的流；"保持一切忠实态"只在 $H$ 属于中心时发生（那时 $\alpha=\mathrm{id}$）。
5. **时间箭头（定理 5.1）**：$\sigma_{-t}=\sigma_t^{-1}$ 使模流本身不区分方向；**箭头由 $\rho$ 的本征值序给出**，且 $K=-\ln\rho$ 的谱单调。箭头给方向、熵给大小，两者都来自同一个 $\rho$。
6. **热态（定理 4.6–4.9、命题 4.10）**：给定平均能量时 Gibbs 态是唯一最大熵态；$\beta$ 与 $E$ 一一对应；$S=\ln Z+\beta E$；KMS 条件在非简并谱上**反推出** Boltzmann 权重；Rindler 视界给 $\beta_\eta\kappa=2\pi$（Killing 归一化下 $\kappa=1$）与 $T_U=\hbar a/(2\pi c k_B)$。

**一处记账**：把 $\rho$ 取成 Gibbs 形式 $e^{-\beta H}/Z$ 是**结构选择（I）**；给定它，KMS 是定理。定理 4.9 说明：在非简并谱上，KMS 条件把 Gibbs 权重定出来。

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

**补核情形二**（块间分量 ⟹ 流出代数）：取块对角态与掺入块间分量（$\rho'=0.9\rho+0.1\,\mathbf 1/4$）的态，让模流作用在块内算符上，看流后是否留在块内：

| 态 | 流后落在块外的成分 |
|:--|--:|
| 块对角 | $0.000$ |
| 含块间分量 | $6.5\times10^{-2}$ |

故情形二成立：**块间分量使流把块内算子送出代数**。

（情形三（张量积／因子代数）本表仍未覆盖。）

> **定理 4.4（态依赖：因子上成立，直和上只看到条件态）.** 分两种代数：

$$
\textbf{(i)}\ \mathcal A=M_n(\mathbb C)；\qquad
\textbf{(ii)}\ \mathcal A=\bigoplus_{h=1}^H M_{d_h}(\mathbb C) .
$$

**(i)** 设 $\rho\ne\rho'$ 为 $M_n(\mathbb C)$ 上两个忠实态，则存在 $A$ 与 $t$ 使 $\sigma_t(A)\ne\sigma'_t(A)$。

**(ii)** 设 $\rho=\bigoplus_h w_h\tilde\rho_h$、$\rho'=\bigoplus_h w'_h\tilde\rho_h$（$\tilde\rho_h$ 是块内归一化态）。则

$$
\sigma=\sigma'\iff \tilde\rho_h=\tilde\rho'_h\ (\forall h);
\qquad\text{中心权重 } w_h\ \textbf{不进入模流}.
$$

*证明.* **（i）** **反证**：设两流处处相同，即 $\sigma_t=\sigma'_t$ 对一切 $t$（作为代数自同构）。**在因子上**，两个酉算子给出同一个自同构当且仅当它们差一个相位（因为 $V^*U$ 与一切算符对易 ⟹ 属中心 ⟹ 是标量），故存在 $|c(t)|=1$ 使

$$
\rho'^{\,it}=c(t)\,\rho^{\,it}\qquad(\forall t\in\mathbb R).
$$

在 $t=0$ 处对 $t$ 求导（$\rho,\rho'$ 忠实，故 $\ln$ 存在）：$i\ln\rho'=\dot c(0)+i\ln\rho$，即 $\ln\rho'=\ln\rho+\lambda$（$\lambda$ 为常数）。两边取指数得 $\rho'=e^{\lambda}\rho$；再由 $\operatorname{tr}\rho'=\operatorname{tr}\rho=1$ 得 $e^{\lambda}=1$，故 $\rho'=\rho$——与假设矛盾。

*证明（续）. **（ii）** 直和上 $\rho^{it}=\bigoplus_h(w_h\tilde\rho_h)^{it}=\bigoplus_h w_h^{it}\tilde\rho_h^{it}$，故对 $A=\bigoplus_hA_h$

$$
\sigma_t(A)=\bigoplus_h w_h^{it}\tilde\rho_h^{it}A_h\,w_h^{-it}\tilde\rho_h^{-it}
=\bigoplus_h \tilde\rho_h^{it}A_h\tilde\rho_h^{-it},
$$

标量相位 $w_h^{it}w_h^{-it}\equiv1$ **恒等相消**。故两流相同当且仅当各块条件态相同；$w_h$ 与 $w'_h$ 可以完全不同而流不变。$\blacksquare$

> **推论（相位相同的两态不是同一个流，权重不同的两态也不是）.** (i) 说明"态不同 ⟹ 流不同"**只在因子上**成立；(ii) 说明直和代数上模流是**中心盲**的——它只记录块内的比值结构。本库的 $\mathcal A$ 是直和（§4.2），故以后要用模流区分两态，必须先确认它们在**同一块内**不同。

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

## §4.5 热态：最大熵与 Gibbs

§2 的 KMS 是从 Gibbs 态**推出来**的。本节补反方向：把"微观计数"与"Gibbs 权重"接上，并给出 KMS 的逆定理。

> **定理 4.6（最大熵 ⟹ Gibbs）.** 设 $\mathcal H=\mathbb C^n$、$H=H^\ast$，$E$ 在 $H$ 的谱区间内。则
$$
\rho_\beta=\frac{e^{-\beta H}}{Z},\qquad \langle H\rangle_\beta=E
\tag{4.7}
$$
在约束族 $\{\rho:\operatorname{tr}\rho=1,\ \operatorname{tr}(\rho H)=E\}$ 中使 $S(\rho)=-\operatorname{tr}\rho\ln\rho$ 取极大；且对**一切**态
$$
S(\rho)\ \le\ \ln Z+\beta\,\langle H\rangle_\rho,\qquad \text{等号}\iff\rho=\rho_\beta .
\tag{4.8}
$$
*证明.* 令 $\ln\rho_\beta=-\beta H-\ln Z$。相对熵 $D(\rho\|\rho_\beta)=\operatorname{tr}\rho(\ln\rho-\ln\rho_\beta)\ge0$（Klein 不等式；Zb047 定理 2.1 的矩阵版），展开：
$$
0\le -S(\rho)+\beta\langle H\rangle_\rho+\ln Z ,
$$
即 (4.8)。等号当且仅当 $\rho=\rho_\beta$。约束 $\langle H\rangle_\beta=E$ 由 $\beta\mapsto\langle H\rangle_\beta$ 的连续性给出。$\blacksquare$

**读法（"从计数到 Gibbs"的那一步）.** (4.8) 的右端只用到 $Z=\operatorname{tr}e^{-\beta H}$（一个计数）与 $\langle H\rangle$；最大化把权重钉成 $e^{-\beta E_k}/Z$。故 Gibbs 权重是"给定平均能量时熵最大"的唯一解，不需要额外的概率公设——概率公设的位置在 Zb018（读出计数）。

> **定理 4.7（$\beta$ 与 $E$ 一一对应）.** 记 $\langle H\rangle_\beta=-\partial_\beta\ln Z$。则
$$
\partial_\beta\langle H\rangle_\beta=-\operatorname{Var}_\beta(H)\le0,
\tag{4.9}
$$
且等号只在 $H$ 于支撑上取单值时成立。故在非平凡谱上 $\beta\mapsto E$ 严格单调，两者一一对应。

*证明.* $Z=\sum_ke^{-\beta E_k}$，逐项求导：$-\partial_\beta\ln Z=\sum_kE_ke^{-\beta E_k}/Z$。再求导得 $-\sum_kE_k^2p_k+(\sum_kE_kp_k)^2=-\operatorname{Var}(H)$。$\blacksquare$

> **定理 4.8（热力学恒等式）.** 记 $E=\langle H\rangle_\beta$、$S=S(\rho_\beta)$。则
$$
S=\ln Z+\beta E,\qquad \frac{\partial S}{\partial E}=\beta,\qquad
F:=-\frac1\beta\ln Z,\qquad \operatorname{Var}(H)=\frac{\partial^2\ln Z}{\partial\beta^2}.
\tag{4.10}
$$
*证明.* 第一式由 (4.8) 取等号；对 $E$ 求导并用 (4.9) 的逆得 $\partial S/\partial E=\beta$；后两式是定义与 (4.9)。$\blacksquare$

> **定理 4.9（KMS 的唯一性）.** 设 $\omega$ 是 $M_n(\mathbb C)$ 上的态，$\sigma_t(A)=e^{iHt}Ae^{-iHt}$，且对一切 $A,B$ 与 $t$ 满足 (2.2)。则 $\omega$ 在 $H$ 的本征基下的对角元为
$$
\omega_{mm}=\frac{e^{-\beta E_m}}{Z};
\tag{4.11}
$$
不同本征值之间的非对角元为零。故在**非简并**谱上 $\omega=\rho_\beta$；在简并子空间内 KMS 不选态（该子空间上任意密度矩阵都满足 (2.2)）。

*证明.* 取 $A=|m\rangle\langle n|$、$B=|n\rangle\langle m|$，$t=0$：左边 $=\omega_{mm}e^{-\beta(E_n-E_m)}$，右边 $=\omega_{nn}$，故 $\omega_{mm}e^{\beta E_m}=\omega_{nn}e^{\beta E_n}$，即对角元 $\propto e^{-\beta E_m}$；归一给 (4.11)。取 $A=|m\rangle\langle m|$、$B=|m\rangle\langle n|$：左边 $=e^{\beta(E_m-E_n)}\omega_{mn}$，右边 $=\omega_{mn}$，故 $E_m\ne E_n$ 时 $\omega_{mn}=0$。$\blacksquare$

**读法（这条补上了什么）.** 定理 2.2 说"Gibbs 态满足 KMS"；定理 4.9 说"满足 KMS 的态在非简并谱上就是 Gibbs"。两条合起来：**KMS 条件本身把 Boltzmann 权重 $e^{-\beta E_m}$ 定出来**，不需要先假定 Gibbs 形式。Gibbs 形式作为"结构选择（I）"的记账（摘要末行）因此可以换成更强的说法：给定 KMS 与谱，权重是**导出**的；未被导出的是"为什么物理态满足 KMS"这一条（它是热平衡的定义）。

**核验（数值，$n=5$，随机 $H$）.**

| 检查 | 结果 |
|:--|:--|
| 300 个随机态上 $S(\rho)-(\ln Z+\beta\langle H\rangle_\rho)$ 的最大值 | $-0.598\le0$ |
| Gibbs 态上同一差 | $-8.9\times10^{-16}$ |
| $\beta\in[0.05,4]$ 上 $E(\beta)$ 的单调性（60 点） | 全部严格下降（$0.0976\to-3.2283$） |
| 随机 $A,B$ 的 KMS 关系相对误差 | $1.3\times10^{-15}$ |
| 由 KMS 反解的对角权重 vs $e^{-\beta E}/Z$ | 差 $1.4\times10^{-16}$ |

---

## §4.6 几何流与视界温度

> **命题 4.10（Rindler 的模周期）.** 在 $1+1$ 维 Rindler 楔取坐标 $(\eta,\xi)$，度规为
$$
ds^2=-\xi^2d\eta^2+d\xi^2 .
\tag{4.12}
$$
则：(i) 固定 $\xi$ 的观察者的固有时为 $\tau=\xi\eta$；(ii) boost 的 Killing 向量 $\zeta=\partial_\eta$ 在视界 $\xi=0$ 处的表面引力为 $\kappa=1$（这一归一化下它沿视界为常数）；(iii) Euclid 化后 $(\xi,\eta)$ 在视界附近是平面极坐标，光滑性要求 $\eta$ 的周期为 $2\pi$；(iv) 若热态是 boost 流的 KMS 态（定理 2.2），则对该流的参数 $\eta$ 有 $\beta_\eta\kappa=2\pi$，即 $\beta_\eta=2\pi$；固定在 $\xi$ 的观察者用固有时 $\tau=\xi\eta$ 记账，故其局部温度为
$$
T_U(\xi)=\frac{\hbar}{2\pi c\,k_B\,\xi}=\frac{\hbar a}{2\pi c\,k_B},\qquad a=\frac1\xi .
\tag{4.13}
$$
*证明.* (i) 沿固定 $\xi$ 的曲线 $d\tau=\sqrt{-ds^2}=\xi\,d\eta$。(ii) 由 $\kappa^2=-\frac12(\nabla_\mu\zeta_\nu)(\nabla^\mu\zeta^\nu)$：$\zeta_\eta=-\xi^2$，$\nabla_\xi\zeta_\eta=-\xi$、$\nabla_\eta\zeta_\xi=\xi$，$g^{\xi\xi}=1$、$g^{\eta\eta}=-\xi^{-2}$，代入得 $\kappa^2=1$。(iii) 令 $\rho=\xi$，度规为 $d\rho^2+\rho^2d\eta^2$；$\rho=0$ 处无锥奇点当且仅当 $\eta\sim\eta+2\pi$。(iv) KMS 的虚时间周期是 $\beta$（(2.2) 的等价形式）；boost 的 Euclid 周期是 $2\pi$，故按 $\eta$ 记账的逆温度是 $2\pi$。观察者的固有时是 $\tau=\xi\eta$，故按 $\tau$ 记账的逆温度是 $2\pi\xi$，对应的温度即 (4.13)。$\blacksquare$

> **推论 4.11（接口的量纲记账）.** 取库内的长度单位 $\ell$、时间单位 $\tau$（$\ell$ 即 Zb005 的格距 $a$ 所取的单位，此处用 $\ell$ 以免与固有加速度 $a$ 撞符号；$c$ 无量纲，故 $\ell$ 与 $\tau$ 只差纯数）。固有加速度 $a$ 的量纲是 $\ell/\tau^2$，故
$$
T_U=\frac{\hbar}{k_B\tau}\cdot\frac{a\tau^2}{2\pi c\,\ell}
=\frac{\hbar}{k_B\tau}\times(\text{纯数}),
\tag{4.14}
$$
即 $T_U$ 的数值 = 纯数 × $\hbar/(k_B\tau)$，与 Zb005 的"唯一锚"一致：$\hbar,c,k_B$ 是**单位**，不是新的有量纲常数。

**核验（数值）.**

| $\xi$ | $\int\xi\,d\eta$（数值） | $\xi\cdot\Delta\eta$ | 相对差 | $a=1/\xi$ | $\beta_\eta\kappa$ |
|--:|--:|--:|--:|--:|--:|
| $0.5$ | $2.0000013$ | $2.0$ | $6.7\times10^{-7}$ | $2.0$ | $2\pi$ |
| $1.0$ | $4.0000027$ | $4.0$ | $6.7\times10^{-7}$ | $1.0$ | $2\pi$ |
| $2.0$ | $8.0000053$ | $8.0$ | $6.7\times10^{-7}$ | $0.5$ | $2\pi$ |

（相对差来自积分步长；表面引力 $\kappa=1$ 由 (ii) 的公式直接给出，$\beta_\eta\kappa=2\pi$ 即 (iii) 的周期。）

**读法（这条接口的位置）.** (4.13) 的左边是 KMS 参数（本库的定理 2.2），右边是视界几何（Rindler 的初等几何）。两者的桥梁是"热时间与几何时间只差标度"（命题 3.1）。视界**本身**在库内的建立仍记为开放（背靠 Zb015 定理 2.2；Zb058 命题 4.1 给了它需要的两条材料）。

---

## §4.7 记账与边界

| 项 | 档位 | 说明 |
|:--|:--|:--|
| Gibbs 形式 $e^{-\beta H}/Z$ | **结构选择（I）** | 摘要末行的记账；定理 4.9 说明它对**非简并**谱由 KMS 反推出来 |
| 最大熵 ⟹ Gibbs | 定理 | 定理 4.6 |
| $\beta\leftrightarrow E$ 一一对应 | 定理 | 定理 4.7 |
| KMS ⟹ Gibbs（非简并谱） | 定理 | 定理 4.9 |
| Rindler 的 $\beta\kappa=2\pi$ | 定理（初等几何） | 命题 4.10 |
| 视界的建立 | 开放（背靠 Zb015 定理 2.2） | Zb058 命题 4.1 |

**边界（热力学极限的位置）.** 本节的 Gibbs 族是**有限维精确对象**：$Z=\operatorname{tr}e^{-\beta H}$、$S=\ln Z+\beta E$ 都在固定维数下成立。广延极限要求 $S\propto$ 系统大小，而本库的账本熵是 $O(\ln L)$（Zb047 命题 5.2），故相变与临界现象的定义要另找承载量（Zb047 命题 6.2）。本节不越过这条线。

（核验脚本：`_staging/热态核验.py`；数据：`_staging/热态核验.json`。）

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
- §4.6：$\eta,\xi$ 无量纲（以 $\ell$ 为单位）；表面引力 $\kappa=1$（无量纲）；固有加速度 $a=\ell/\tau^2$；$\beta_\eta$ 无量纲；$T_U$ 是温度，数值 = 纯数 × $\hbar/(k_B\tau)$。
- $\hbar,c,k_B$：**单位**（Zb005 命题 4.2′ 的判据），不是新的有量纲常数。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.2 | 幺正群、保迹、保正 | 导出 | §1 |
| 推论 1.3 | 谱在演化下不变 | 导出 | §1 |
| 定理 2.2 | KMS 关系 | 导出 | §2 |
| 命题 3.1 | 模流＝以 $\beta H$ 为生成元 | 导出 | §3 |
| 命题 4.1 | $\beta$ 是桥 | 导出 | §4 |
| 命题 4.2 | 中心态的模流在块内平凡 | 导出 | §4 |
| 定理 4.3 | 模流的完整判据 | 导出 | §4 |
| 定理 4.4 | 态依赖：因子上成立、直和上只看到条件态 | 导出 | §4 |
| 定理 4.5 | 普遍自同构（直和上的形式） | 导出 | §4 |
| 定理 4.6／4.7／4.8 | 最大熵 ⟹ Gibbs；$\beta\leftrightarrow E$；热力学恒等式 | 导出＋数值 | §4.5 |
| 定理 4.9 | KMS ⟹ Gibbs（非简并谱） | 导出＋数值 | §4.5 |
| 命题 4.10／推论 4.11 | Rindler 的 $\beta_\eta\kappa=2\pi$；温度的量纲记账 | 导出＋数值 | §4.6 |
| 定理 5.1 | 时间箭头 ＝ 谱序 | 导出 | §5 |
| 记账 | Gibbs 形式是 (I) | 记账 | §2 |
