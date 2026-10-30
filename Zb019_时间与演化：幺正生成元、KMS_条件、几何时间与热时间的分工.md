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

**三条正面结果**：

1. **幺正性（定理 1.2）**：$U(t)=e^{-iHt}$ 是单参数幺正群，共轭作用保持态的集合。
2. **KMS 条件（定理 2.2）**：Gibbs 态 $\rho=e^{-\beta H}/Z$ 的模流与 $H$ 生成的演化一致，满足 KMS 关系。
3. **两种时间的关系（命题 3.1）**：$\sigma_t^{\rm mod}(A)=e^{-i\beta Ht}Ae^{i\beta Ht}=\sigma^{\rm geo}_{-\beta t}(A)$——模流就是以 $\beta H$ 为生成元的演化，差一个由约定固定的符号。

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

## §5 量纲

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
