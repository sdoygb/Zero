# Zb018 · 运动学：态、可观测量与读出计数——Born 比值不是新公理

> 篇号：Zb018
> 状态：正式篇（推导篇）
> 前置：Zb001（两条公理：零不停留、不设概率）、Zb003（层与词层）、Zb007（模流与内禀频率）
> 登记日期：2026-10-30

## 摘要

本篇建立量子运动学，并把"概率"的位置说清：

$$
\boxed{\;\begin{aligned}
&\textbf{可观测量}:\ \text{自伴算子}\ A=A^\ast\ \text{（有限维、实谱）};\\
&\textbf{读出计数}:\ N_k=\operatorname{tr}(\rho P_k);\\
&\textbf{Born 比值}:\ p_k=\frac{\operatorname{tr}(\rho P_k)}{\operatorname{tr}\rho}=\operatorname{tr}(\rho P_k)\quad(\operatorname{tr}\rho=1);\\
&\textbf{记账}:\ p_k\ \text{是}\textbf{计数比}，\text{由公理②（不设概率）留给} \text{读出}.
\end{aligned}\;}
\tag{0.1}
$$

**三条正面结果**：

1. **谱定理（定理 1.2）**：有限维自伴算子有实谱与正交本征基，可观测量的"读数"因此是实数。
2. **Born 比值定理（定理 2.2）**：把读数定义为 $\operatorname{tr}(\rho P_k)$ 后，比值自动非负、归一、对正交分解可加——**不需要额外公设**。
3. **纯态的两次方（推论 2.4）**：$\rho=\lvert\psi\rangle\langle\psi\rvert$ 时 $p_k=\lvert\langle k\rvert\psi\rangle\rvert^2$，即"振幅的模方"由迹的定义直接给出。

**一处记账**：把"读出"实现为投影 $P_k$ 上的迹，是**结构选择（I）**；它的后果（非负、归一、可加）是定理。

---

## §1 态与可观测量

> **定义 1.1（有限维量子系统）.** 系统由一个复 Hilbert 空间 $\mathcal H=\mathbb C^n$ 描述。

- **态**：$\rho\in M_n(\mathbb C)$，满足 $\rho=\rho^\ast$、$\rho\ge0$（半正定）、$\operatorname{tr}\rho=1$。
- **可观测量**：自伴算子 $A=A^\ast$。
- **纯态**：$\rho^2=\rho$ 的态；此时 $\rho=\lvert\psi\rangle\langle\psi\rvert$ 对某个单位向量 $\psi$。

> **定理 1.2（谱定理，有限维）.** 自伴 $A$ 有分解

$$
A=\sum_{k}\lambda_k P_k,\qquad \lambda_k\in\mathbb R,\qquad P_kP_{k'}=\delta_{kk'}P_k,\qquad \sum_kP_k=I .
\tag{1.1}
$$

*证明.* 有限维复空间上的自伴算子有本征值：特征多项式在 $\mathbb C$ 上有根，设 $A\psi=\lambda\psi$，则 $\lambda\|\psi\|^2=\langle\psi,A\psi\rangle=\overline{\langle A\psi,\psi\rangle}=\bar\lambda\|\psi\|^2$，故 $\lambda\in\mathbb R$。不同本征值的本征子空间互相正交：$(\lambda-\lambda')\langle u,v\rangle=\langle Au,v\rangle-\langle u,Av\rangle=0$。取各子空间的正交基即得 (1.1)。$\blacksquare$

**读法.** 可观测量的"读数"是实数，这一点由自伴性给出（不是由公设加上的）。

---

## §2 读出与 Born 比值

> **结构选择（I）（读出）.** 对可观测量的谱分解 (1.1)，**读出**定义为：给定量子态 $\rho$，第 $k$ 个结果的出现**计数**是

$$
N_k(\rho):=\operatorname{tr}(\rho P_k).
\tag{2.1}
$$

> **定理 2.2（Born 比值）.** 记 $N:=\sum_kN_k(\rho)$。则

$$
p_k:=\frac{N_k}{N}=\operatorname{tr}(\rho P_k),\qquad p_k\ge0,\qquad \sum_kp_k=1,\qquad p_k+p_{k'}=p_{k\cup k'},
\tag{2.2}
$$

其中最后一条指 $P_k+P_{k'}$ 仍是投影时的合并读数。

*证明.* **归一**：$\sum_kN_k=\sum_k\operatorname{tr}(\rho P_k)=\operatorname{tr}\bigl(\rho\sum_kP_k\bigr)=\operatorname{tr}(\rho I)=\operatorname{tr}\rho=1$，故 $N=1$、$p_k=N_k$。
**非负**：$\rho\ge0$ 与 $P_k\ge0$ 给出 $\rho^{1/2}P_k\rho^{1/2}\ge0$，故 $\operatorname{tr}(\rho P_k)=\operatorname{tr}(\rho^{1/2}P_k\rho^{1/2})\ge0$。
**可加**：若 $P_kP_{k'}=0$，则 $P_k+P_{k'}$ 是到直和的投影，$\operatorname{tr}(\rho(P_k+P_{k'}))=\operatorname{tr}(\rho P_k)+\operatorname{tr}(\rho P_{k'})$。$\blacksquare$

**读法.** 三条性质都由 $\operatorname{tr}$ 的线性与 $\rho$ 的半正定性给出。**没有第三条公理**说"概率等于迹的模方"——那正是 (2.1) 的定义加上 (2.2) 的推论。

---

## §3 与公理②的关系

公理②说"不设概率"。本篇的处理是：

| 层 | 内容 | 档位 |
|:--|:--|:--|
| 动力学 | 幺正演化（确定性的） | 定理（见下篇） |
| 读出 | 投影 $P_k$ 上的迹计数 (2.1) | 结构选择（I） |
| 比值 | $p_k=\operatorname{tr}(\rho P_k)$ | **定理** 2.2 |
| 频率解释 | 大量独立读出的相对频率趋近 $p_k$ | 定理（大数定律，§4） |

> **命题 3.1（比值与公理的分工）.** 公理②禁止的是"把概率当作基本的、不可再解释的东西"；它允许"计数比"存在。(2.2) 的 $p_k$ 是计数比，其地位与掷骰子的相对频率相同。

*证明.* 定义 2.1 与定理 2.2 只用到迹与半正定性，未引入任何概率测度或随机跳变。$\blacksquare$

---

## §4 纯态：两次方的来源

> **定理 4.1（纯态）.** 若 $\rho=\lvert\psi\rangle\langle\psi\rvert$（$\|\psi\|=1$），则对任意投影 $P_k$：

$$
p_k=\operatorname{tr}(\rho P_k)=\langle\psi\rvert P_k\lvert\psi\rangle .
\tag{4.1}
$$

若进一步 $P_k=\lvert k\rangle\langle k\rvert$ 是某个正交基的投影，则

$$
p_k=\lvert\langle k\rvert\psi\rangle\rvert^2 .
\tag{4.2}
$$

*证明.* 用迹的循环性：$\operatorname{tr}(\lvert\psi\rangle\langle\psi\rvert P_k)=\operatorname{tr}(\langle\psi\rvert P_k\lvert\psi\rangle)=\langle\psi\rvert P_k\lvert\psi\rangle$（内积是标量，迹等于它自身）。取 $P_k=\lvert k\rangle\langle k\rvert$ 得 $\langle\psi\rvert k\rangle\langle k\rvert\psi\rangle=\lvert\langle k\rvert\psi\rangle\rvert^2$。$\blacksquare$

**读法.** "振幅的模方"是**算子迹**的直接结果，不是独立公设。这也解释了为什么 $p_k$ 必须是非负的：它是 $\langle\psi\rvert P_k\lvert\psi\rangle$，一个半正定二次型。

---

## §5 频率趋近（引用）

大量**独立**读出的相对频率趋近 $p_k$，这是概率论的大数定律（外部定理）。它的前提是"读出彼此独立"，这在本库是**结构选择（I）**：把每次读出当作同分布且互不影响。

**读法.** 大数定律把"计数比"与"观测到的频率"接起来；本库不为它付新的物理公设。

---

## §6 量纲

- $\rho$、$P_k$、$A$ 无量纲。
- 本征值 $\lambda_k$ 的量纲由具体可观测量决定。
- 本篇不引入有量纲常数。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.2 | 有限维谱定理 | 导出 | §1 |
| 定理 2.2 | Born 比值：非负、归一、可加 | 导出 | §2 |
| 命题 3.1 | 计数比与公理②的分工 | 记账 | §3 |
| 定理 4.1 | 纯态 ⟹ 模方 | 导出 | §4 |
| §5 | 频率趋近引大数定律 | 引用＋(I) | §5 |
