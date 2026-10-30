# Zb025 · 开放系统：完全正、Kraus 表示与 Lindblad 生成元的保迹性

> 篇号：Zb025
> 状态：正式篇（推导篇）
> 前置：Zb019（演化与 KMS）、Zb021（退相干与效应）、Zb023（门与耦合）
> 登记日期：2026-10-30

## 摘要

本篇把开放系统的三条代数事实写清：

$$
\boxed{\;\begin{aligned}
&\textbf{Kraus}:\ \Phi(\rho)=\sum_kM_k\rho M_k^\ast\ \Longrightarrow\ \Phi\ \text{完全正};\\
&\textbf{保迹}:\ \operatorname{tr}\Phi(\rho)=\operatorname{tr}\rho\iff\textstyle\sum_kM_k^\ast M_k=I;\\
&\textbf{Lindblad}:\ \dot\rho=-i[H,\rho]+\sum_k\Bigl(A_k\rho A_k^\ast-\tfrac12\{A_k^\ast A_k,\rho\}\Bigr)\ \text{保迹}.
\end{aligned}\;}
\tag{0.1}
$$

**三条正面结果**：

1. **完全正（定理 1.2）**：任何 Kraus 形式都给出完全正映射，且 Kraus 形式对"系统-环境耦合 + 部分迹"完全等价（Zb023 定理 2.2 的逆命题）。
2. **保迹判据（定理 2.1）**：$\sum_kM_k^\ast M_k=I$ ⟺ 保迹。
3. **Lindblad 保迹（定理 3.2）**：GKSL 生成元逐项保迹，故 $\operatorname{tr}\rho(t)$ 是常数。

**一处记账**：取 GKSL 形式（而不取任意保迹正映射）是**结构选择（I）**——它限制在"马尔可夫、时间齐次"的一支；由它推出的保迹性是定理。

---

## §1 Kraus 表示与完全正

> **定义 1.1（完全正）.** 线性映射 $\Phi:\mathcal B(\mathcal H)\to\mathcal B(\mathcal H)$ 称为**完全正**，若对任意辅助维数 $d$，$\Phi\otimes I_d$ 都是正映射（把正算子映到正算子）。

> **定理 1.2（Kraus ⟹ 完全正）.** 若 $\Phi(\rho)=\sum_kM_k\rho M_k^\ast$，则 $\Phi$ 完全正。

*证明.* 对任意 $d$ 与任意 $\rho'\ge0$（在 $\mathcal H\otimes\mathbb C^d$ 上），

$$
(\Phi\otimes I_d)(\rho')=\sum_k(M_k\otimes I_d)\,\rho'\,(M_k\otimes I_d)^\ast .
$$

每一项都是 $BXX^\ast B^\ast$ 型（半正定）：对任意 $\lvert v\rangle$，$\langle v\rvert (M\otimes I)\rho'(M\otimes I)^\ast\lvert v\rangle=\langle (M\otimes I)^\ast v\rvert\rho'\lvert(M\otimes I)^\ast v\rangle\ge0$。有限个半正算子之和半正。$\blacksquare$

> **命题 1.3（等价性）.** Kraus 形式与"耦合 + 部分迹"是同一类映射：每个 Kraus 形式都可由 Zb023 定理 2.2 的酉耦合实现，反之亦然。

*证明.* "Kraus ⟹ 耦合"见 Zb023 定理 2.2；"耦合 ⟹ Kraus"：对 $U(\rho\otimes\lvert0\rangle\langle0\rvert)U^\ast$ 取 $\mathcal H_E$ 的基 $\{\lvert e_k\rangle\}$，得

$$
\operatorname{tr}_E(\cdots)=\sum_k(I\otimes\langle e_k\rvert)U(\rho\otimes\lvert0\rangle\langle0\rvert)U^\ast(I\otimes\lvert e_k\rangle)=\sum_kM_k\rho M_k^\ast,
$$

其中 $M_k:=(I\otimes\langle e_k\rvert)U(I\otimes\lvert0\rangle)$。$\blacksquare$

---

## §2 保迹判据

> **定理 2.1（保迹 ⟺ Kraus 完备）.**

$$
\operatorname{tr}\Phi(\rho)=\operatorname{tr}\rho\ \ \forall\rho\iff\sum_kM_k^\ast M_k=I .
\tag{2.1}
$$

*证明.* **（⟸）**：$\operatorname{tr}\bigl(\sum_kM_k\rho M_k^\ast\bigr)=\sum_k\operatorname{tr}(M_k\rho M_k^\ast)=\sum_k\operatorname{tr}(M_k^\ast M_k\rho)=\operatorname{tr}\bigl((\sum_kM_k^\ast M_k)\rho\bigr)=\operatorname{tr}\rho$。
**（⟹）**：取 $\rho=\lvert v\rangle\langle v\rvert$，得 $\langle v\rvert\sum_kM_k^\ast M_k\lvert v\rangle=\langle v\rvert v\rangle=1$ 对一切 $v$。故 $S:=\sum_kM_k^\ast M_k$ 满足 $\langle v\rvert S\lvert v\rangle=\langle v\rvert v\rangle$ 对一切 $v$，即 $S=I$（极化恒等式）。$\blacksquare$

> **命题 2.2（正但不完全正的例子）.** 转置映射 $T(\rho)=\rho^{\mathsf T}$ 正、保迹，但**不**完全正。

*证明.* 正性：$\rho\ge0\Rightarrow\rho^{\mathsf T}\ge0$；保迹：$\operatorname{tr}\rho^{\mathsf T}=\operatorname{tr}\rho$。不完全正：取 $\lvert\Omega\rangle=\frac{1}{\sqrt2}(\lvert00\rangle+\lvert11\rangle)$，则 $(T\otimes I)(\lvert\Omega\rangle\langle\Omega\rvert)$ 有一个负本征值 $-1/2$（Bell 态的转置）。$\blacksquare$

**读法.** 这条是"为什么要求完全正而不是正"的判据：正性说的是"对系统本身"，完全正说的是"对系统加任意辅助系统"。开放系统必须用后者——因为环境随时可能纠缠进来。

---

## §3 Lindblad 生成元

> **定义 3.1（GKSL 形式）.**

$$
\dot\rho=\mathcal L(\rho):=-i[H,\rho]+\sum_k\Bigl(A_k\rho A_k^\ast-\tfrac12\{A_k^\ast A_k,\rho\}\Bigr),
\tag{3.1}
$$

其中 $H=H^\ast$，$A_k$ 任意，$\{X,\rho\}:=X\rho+\rho X$。

> **定理 3.2（保迹）.** $\operatorname{tr}\mathcal L(\rho)=0$ 对一切 $\rho$。

*证明.* 逐项：

$$
\operatorname{tr}(-i[H,\rho])=-i\bigl(\operatorname{tr}(H\rho)-\operatorname{tr}(\rho H)\bigr)=0 ;
$$

$$
\operatorname{tr}\bigl(A\rho A^\ast\bigr)=\operatorname{tr}(A^\ast A\rho),\qquad
\operatorname{tr}\bigl(\tfrac12\{A^\ast A,\rho\}\bigr)=\tfrac12\bigl(\operatorname{tr}(A^\ast A\rho)+\operatorname{tr}(\rho A^\ast A)\bigr)=\operatorname{tr}(A^\ast A\rho).
$$

两项相消。对 $k$ 求和仍为零。$\blacksquare$

> **推论 3.3.** $\operatorname{tr}\rho(t)=\operatorname{tr}\rho(0)$ 对一切 $t$。

*证明.* 由定理 3.2，$\frac{d}{dt}\operatorname{tr}\rho=\operatorname{tr}\mathcal L(\rho)=0$。$\blacksquare$

> **命题 3.4（振幅阻尼的显式解）.** 取 $H=0$、单个 $A=\sqrt\gamma\,\sigma_-$（$\sigma_-=\lvert0\rangle\langle1\rvert$）。则

$$
\rho(t)=\sum_{k=0,1}M_k(t)\,\rho(0)\,M_k(t)^\ast,\qquad
M_0=\begin{pmatrix}1&0\\0&\sqrt{1-p}\end{pmatrix},\quad
M_1=\begin{pmatrix}0&\sqrt p\\0&0\end{pmatrix},
\tag{3.2}
$$

其中 $p=1-e^{-\gamma t}$，且 $\sum_kM_k^\ast M_k=I$。

*证明.* 直接计算 $M_0^\ast M_0+M_1^\ast M_1=\operatorname{diag}(1,1-p)+\operatorname{diag}(0,p)=I$。把它对 $t$ 求导并与 (3.1)（$A=\sqrt\gamma\sigma_-$）比较，两项逐项相同：$\dot p=\gamma(1-p)$。$\blacksquare$

**读法.** 振幅阻尼是 GKSL 形式的一个**可显式求解**的例子；它的 Kraus 完备性（$\sum M^\ast M=I$）同时验证了定理 2.1 与定理 3.2。

---

## §4 记账

| 项 | 档位 | 说明 |
|:--|:--|:--|
| 取 GKSL 形式 | 结构选择（I） | 限制在马尔可夫、时间齐次的一支 |
| $H$、$A_k$ | 由耦合给定 | 来自 Zb023 的耦合模型 |
| 保迹、完全正（Kraus） | 定理 | 定理 1.2、2.1、3.2 |

---

## §5 量纲

- $\rho$ 无量纲；$H$：每步；$A_k$：每步$^{1/2}$（使 $A^\ast A$ 与 $H$ 同量纲）。
- 本篇不引入有量纲常数。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.2 | Kraus ⟹ 完全正 | 导出 | §1 |
| 命题 1.3 | Kraus 与"耦合＋部分迹"等价 | 导出 | §1 |
| 定理 2.1 | 保迹 ⟺ Kraus 完备 | 导出 | §2 |
| 命题 2.2 | 转置：正但不完全正 | 导出 | §2 |
| 定理 3.2／推论 3.3 | GKSL 保迹 | 导出 | §3 |
| 命题 3.4 | 振幅阻尼显式解 | 导出 | §3 |
