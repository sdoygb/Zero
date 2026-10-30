# Zb030 · 谱、破缺与尺度：真空流形、Goldstone 计数与自由度守恒

> 篇号：Zb030
> 状态：正式篇（推导篇）
> 前置：Zb020（自旋）、Zb028（载体）、Zb029（规范原理与荷）
> 登记日期：2026-10-30

## 摘要

本篇给出破缺的三件事：

$$
\boxed{\;\begin{aligned}
&\textbf{真空流形}:\ V(\phi)=-\mu^2|\phi|^2+\lambda|\phi|^4\ \Longrightarrow\ |\phi_0|=\frac{v}{\sqrt2},\ v=\frac{\mu}{\sqrt\lambda};\\
&\textbf{Goldstone 计数}:\ \dim(G/H)=\dim G-\dim H;\\
&\textbf{自由度守恒}:\ \underbrace{2}_{\text{无质量规范}}+\underbrace{2}_{\text{复标量}}=\underbrace{3}_{\text{有质量规范}}+\underbrace{1}_{\text{Higgs}} .
\end{aligned}\;}
\tag{0.1}
$$

**三条正面结果**：

1. **真空（定理 1.2）**：$V$ 的极小点构成一圈（$U(1)$ 情形），半径 $v/\sqrt2$。
2. **Goldstone 计数（定理 2.2）**：破缺的生成元个数等于真空流形的维数。
3. **自由度守恒（定理 3.2）**：规范场"吃掉"Goldstone 后，总自由度在破缺前后**相等**——这是 Higgs 机制的代数内容。

**一处记账**：把 $\mu^2>0$ 作为输入（"负质量平方"）是**结构选择（I）**；由它推出的真空位置与计数是定理。

---

## §1 真空流形

> **定义 1.1（Abelian Higgs 势）.** 复标量 $\phi$ 的势

$$
V(\phi)=-\mu^2|\phi|^2+\lambda|\phi|^4,\qquad \mu^2>0,\ \lambda>0 .
\tag{1.1}
$$

> **定理 1.2（极小点）.** $V$ 的极小值为 $-\dfrac{\mu^4}{4\lambda}$，取到于

$$
|\phi_0|^2=\frac{\mu^2}{2\lambda}\equiv\frac{v^2}{2},\qquad v:=\frac{\mu}{\sqrt\lambda}.
\tag{1.2}
$$

即真空流形是半径 $v/\sqrt2$ 的圆周（把 $\phi$ 看作 $\mathbb R^2$）。

*证明.* 令 $r=|\phi|^2\ge0$，则 $V=-\mu^2r+\lambda r^2$，这是 $r$ 的二次函数，顶点在 $r=\dfrac{\mu^2}{2\lambda}$，值 $V=-\dfrac{\mu^4}{4\lambda}$。二阶导 $2\lambda>0$ 故为极小。$|\phi|^2=v^2/2$ 的解集是圆周。$\blacksquare$

> **定义 1.3（破缺）.** 势在 $\phi\mapsto e^{i\alpha}\phi$ 下不变，但极小点**不是**单个点——它在变换下被移到别的极小点。称对称性被**自发破缺**。

**读法.** "破缺"的含义是：**方程**有 $U(1)$ 对称，**解**没有。这不是把对称性去掉，而是解集比单点大。

---

## §2 Goldstone 计数

> **定义 2.1（真空流形）.** 设对称群 $G$、未破缺子群 $H$（保持某个选定 $\phi_0$ 的变换）。真空流形是

$$
\mathcal M:=\{g\phi_0:g\in G\}\cong G/H .
\tag{2.1}
$$

> **定理 2.2（Goldstone 计数）.** 若 $G$ 是连通李群、$H$ 是其闭子群，则

$$
\dim\mathcal M=\dim G-\dim H,
\tag{2.2}
$$

即有 $\dim G-\dim H$ 个无能隙模（Goldstone 玻色子）。

*证明.* 商流形的维数公式 $\dim(G/H)=\dim G-\dim H$（$H$ 是闭子群，$\mathfrak g=\mathfrak h\oplus\mathfrak m$ 作为向量空间分解）。每个 $\mathfrak m$ 方向对应一个沿真空流形的平坦方向，故能量在该方向上不涨——即无质量模。$\blacksquare$

> **推论 2.3（$U(1)$ 情形）.** $G=U(1)$、$H=\{1\}$：$\dim=1-0=1$，恰一个 Goldstone 模。

*证明.* 代入 (2.2)。$\blacksquare$

**读法.** 计数只用到群的维数，与势的具体形状无关——这是它可用的原因。

---

## §3 自由度守恒

> **定理 3.1（无质量矢量的偏振数）.** 四维中无质量矢量有 $2$ 个横偏振；有质量矢量有 $3$ 个。

*证明.* 无质量：$k^2=0$，规范等价类把 $4$ 个分量减到 $2$（Zb011 命题 2.4 的 $D=4$ 情形）。有质量：$k^2=m^2$，横向 $2$ 个加上纵向 $1$ 个。$\blacksquare$

> **定理 3.2（自由度守恒）.** 破缺前后总自由度相等：

$$
\underbrace{2}_{\text{无质量规范场}}+\underbrace{2}_{\text{复标量}}=\underbrace{3}_{\text{有质量规范场}}+\underbrace{1}_{\text{实标量}} .
\tag{3.3}
$$

两边都是 $4$。

*证明.* 左边：$2+2=4$。右边：$3+1=4$。复标量有 $2$ 个实分量（$\operatorname{Re}\phi,\operatorname{Im}\phi$）；破缺后一个变成 Higgs（沿径向），另一个被规范场吃掉变成纵向偏振。$\blacksquare$

**读法.** Higgs 机制**不产生**自由度也不**消灭**自由度：它把 Goldstone 模**重排**成纵向偏振。这解释了为什么 $W$、$Z$ 有质量——它们获得了第三个偏振态。

> **命题 3.3（质量公式）.** 规范场获得的质量为

$$
m_A=q\,v .
\tag{3.4}
$$

*证明.* 把 $\phi=\frac{1}{\sqrt2}(v+h)e^{i\theta/v}$ 代入 $|D_\mu\phi|^2$（$D_\mu=\partial_\mu+iqA_\mu$，Zb029 定理 1.2），取 $h=0$：

$$
|D_\mu\phi|^2\supset\frac{v^2}{2}q^2A_\mu A^\mu ,
$$

这是质量为 $m_A=qv$ 的 Proca 项。$\blacksquare$

---

## §4 尺度的地位

| 量 | 档位 | 说明 |
|:--|:--|:--|
| $\mu^2$ | 初值／结构选择（I） | 取正值是选择；其数值属初值 |
| $\lambda$ | 初值 | 无量纲耦合 |
| $v=\mu/\sqrt\lambda$ | 派生 | 由上式定 |
| $m_A=qv$ | 定理 | 命题 3.3 |

**读法.** 破缺把"一个负的质量平方"变成一个**正的质量**——这是 SM 里质量的来源之一。$v$ 是本篇唯一的新尺度，与 Zb005 的锚分开记账（它是耦合乘幂，不是单位）。

---

## §5 量纲

- $[\phi]=$ 质量；$\mu$：质量；$\lambda$ 无量纲；$v$：质量。
- $m_A=qv$：$q$ 无量纲，故 $m_A$ 与 $v$ 同量纲。
- 本篇不引入有量纲常数。

---

## 自检

| 项 | 内容 | 档位 | 依据 |
|:--|:--|:--|:--|
| 定理 1.2 | 极小点 $|\phi_0|=v/\sqrt2$ | 导出 | §1 |
| 定理 2.2／推论 2.3 | Goldstone 计数 $\dim G-\dim H$ | 导出 | §2 |
| 定理 3.1／3.2 | 偏振数；自由度守恒 $4=4$ | 导出 | §3 |
| 命题 3.3 | $m_A=qv$ | 导出 | §3 |
