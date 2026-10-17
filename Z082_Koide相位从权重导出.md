# Koide 相位从权重导出：$\theta=\arg\hat s_1$，$b$ 与 $\theta$ 同一复数

**Zero 项目 · 基础推导系列 · 第 Z082 号**

**篇性**：**推导篇**。Koide 相位可从纯权重侧导出（$\theta=\arg\hat s_1$，不含时间），$b$ 与 $\theta$ 是同一个复数的模与辐角。
**结论**：**(一)** $\theta$ 可写成纯权重侧量：$\theta=\arg(\hat s_1)$，$\hat s_n$ 是 $\sqrt{m_i}$ 的 DFT；**(二)** $b=2|\hat s_1|/|\hat s_0|$，$Q=(1+b^2/2)/3$；**(三)** $\arg\hat s_1=2\pi/3+2/9$，故 Koide 偏差 $\delta=\arg\hat s_1-2\pi/3\approx 2/9$（实测偏差 $+7.4\times10^{-6}$）；**(四)** 但须三项权重**两两不同**；$L=4$ 上 $r$ 计数 $(2,2,2)$ ⟹ $\hat s_1=0$；**(五)** 故真缺口不是相位、不是时间，而是**一个天然不对称的三项分配**＝Z014 定理的**扇区选择性**。
**纪律**：**未改动任何公理**。
**依据**：正式篇 Z014（定理 27／28／29：形状 $\mathbb Z_3$、角度 $\theta_{\text{物质}}=\pi/2$、振幅 $\sqrt2$ 条件证成；定理 31：扇区选择性定理）。
**核验**：引 `_staging/Z082_weight_phase_check.py`（45 条断言通过）＋本文复算 DFT。

---

## 摘要

$$
\boxed{\;\textbf{（一）}\ \theta=\arg(\hat s_1);\qquad
\hat s_n=\frac13\sum_i s_i\,e^{-2\pi i\,ni/3},\ s_i=\sqrt{m_i}.\;}
\tag{0.1}
$$

$$
\boxed{\;\textbf{（二）}\ b=2\,\frac{|\hat s_1|}{|\hat s_0|};\qquad
Q=\frac{1+b^2/2}{3}\quad\Longrightarrow\quad b^2=2(3Q-1).\;}
\tag{0.2}
$$

$$
\boxed{\;\textbf{（三）}\ \arg\hat s_1=2.316624733883\ \text{rad}
=2\pi/3+2/9+\delta\theta;\qquad\delta\theta=+7.409\times10^{-6}.\;}
\tag{0.3}
$$

$$
\boxed{\;\textbf{（四）}\ \text{须三项权重互异};\ L=4\ \text{的}\ r\ \text{计数}=(2,2,2)\Longrightarrow\hat s_1=0.\;}
\tag{0.4}
$$

---

## §1 设置：Koide 参数化与 DFT

$$
\text{Koide 参数化}:\quad s_i=\sqrt{m_i}=A\Bigl(1+b\cos\bigl(\theta+\tfrac{2\pi i}3\bigr)\Bigr),\qquad i=0,1,2 .
\tag{1.1}
$$

$$
\text{离散傅里叶}:\quad \hat s_n=\frac13\sum_{i=0}^{2}s_i\,e^{-2\pi i\,ni/3},\qquad \hat s_2=\overline{\hat s_1} .
\tag{1.2}
$$

$$
\text{实测}\ (m_e,m_\mu,m_\tau)=(0.51099895,105.6583755,1776.86)\ \text{MeV}\Longrightarrow
s=(0.714842,\ 10.279026,\ 42.152817).
\tag{1.3}
$$

---

## §2 主结果（核验 §A）

$$
\textbf{DFT}:\quad
\begin{cases}
\hat s_0=17.7155617100+0\,i\\
\hat s_1=-8.5003599022+9.2011709720\,i\\
\hat s_2=\overline{\hat s_1}
\end{cases}
\tag{2.1}
$$

$$
\textbf{相位}:\quad \arg(\hat s_1)=+2.316624733883\ \text{rad}\qquad(\textbf{纯权重侧，不含时间}).
\tag{2.2}
$$

$$
\textbf{模}:\quad \frac{|\hat s_1|}{|\hat s_0|}=0.707100252592=\frac{b}2\qquad\Longrightarrow\qquad b=1.414200505185.
\tag{2.3}
$$

$$
\textbf{与 }Q\textbf{ 的严格关系}:\quad Q_{\rm 实测}=0.666660511477\Longrightarrow b^2=1.999963068864\Longrightarrow b=1.414200505185;\qquad Q=\tfrac23\iff b=\sqrt2 .
\tag{2.4}
$$

$$
\boxed{\;b\ \text{与}\ \theta\ \text{是同一个复数}\ \hat s_1\ \text{的模与辐角}\ \Longrightarrow\ \textbf{一个对象，两个参数，不是两个缺口}.\;}
\tag{2.5}
$$

---

## §3 与 Koide 相位 $\delta=2/9$ 的衔接（核验 §B）

$$
\text{对称点}:\quad \frac{2\pi}3=2.0943951024\ \text{rad};\qquad
\arg\hat s_1-\frac{2\pi}3=0.22222963149 .
\tag{3.1}
$$

$$
\frac29=0.2222222222;\qquad
\delta\theta=\Bigl(\arg\hat s_1-\frac{2\pi}3\Bigr)-\frac29=+7.409267\times10^{-6}.
\tag{3.2}
$$

$$
\boxed{\;\text{Koide 偏差}\ \delta=\arg\hat s_1-\tfrac{2\pi}3\ \text{接近候选}\ 2/9\ (\text{偏差}\ 7.4\times10^{-6});\ \text{但"接近即证据"不成立} .\;}
\tag{3.3}
$$

---

## §4 但须三项权重互异（核验 §C）

$$
\text{偏离}:\quad w_i=\frac{s_i}{\langle s\rangle}-1=(-0.959648928,\ -0.419774198,\ +1.379423126)\quad(\textbf{三项互异}).
\tag{4.1}
$$

$$
\text{而}\ L=4\ \text{上}\ r\ \text{的计数}:\quad \{r{=}2{:}2,\ r{=}3{:}2,\ r{=}4{:}2\}\ \Longrightarrow\ \hat s_1\big|_{\rm uniform}=0 .
\tag{4.2}
$$

$$
\boxed{\;\text{缺口}=\textbf{一个天然不对称的三项分配};\ \text{不是相位值，也不是时间}.\;}
\tag{4.3}
$$

---

## §5 与 Z014 的接缝

$$
\text{Z014 定理 27}:\ \text{形状}\ =\mathbb Z_3\ \text{正则表示／厄米循环矩阵};\qquad
\text{定理 28}:\ \theta_{\text{物质}}=\tfrac\pi2\ \text{由圆分割唯一解};\qquad
\text{定理 29}:\ \sqrt2\ \text{振幅由范数条件}\ w_1^2+w_2^2=2\ \text{条件证成}.
\tag{5.1}
$$

$$
\text{Z014 定理 31}:\ \textbf{扇区选择性}\ \text{在群论／一圈动力学／味混合三层面}\ \textbf{定理} .
\tag{5.2}
$$

$$
\boxed{\;\text{本文的"三项不对称"}=\text{Z014 定理的"扇区选择性"};\ \text{闭环回到 Z014 的已知判定}.\;}
\tag{5.3}
$$

---

## §6 结论的适用范围

1. **$b=\sqrt2$ 不是精确成立的**——它由 $Q=\tfrac23$ 推出，实测 $b=1.4142005$（差 $\sim10^{-4}$）；
2. **$\arg\hat s_1=2\pi/3+2/9$ 不是精确成立的**——实测偏差 $7.4\times10^{-6}$，且"接近不构成证据"（Q275）；
3. **$b$ 与 $\theta$ 不"必然同源"于更一般的参数化**——只对 Koide 参数化证明（§2.5）；
4. **任何数值"接近"都不是证据**——仍守 **Q275 纪律**；

---

## §7 结论

$$
\begin{aligned}
&\textbf{成功}:\quad \theta=\arg(\hat s_1);\qquad b=2|\hat s_1|/|\hat s_0|;\qquad Q=(1+b^2/2)/3 .\\
&\qquad\Longrightarrow\ \text{两个参数都在权重里，不含时间};\ b\ \text{与}\ \theta\ \text{是}\ \textbf{同一个复数}\ \text{的模与辐角}.\\
&\textbf{衔接}:\quad \arg\hat s_1=2\pi/3+2/9+\delta\theta;\quad \delta\theta=+7.4\times10^{-6}.\\
&\textbf{但}:\ \text{须三项权重互异};\ L{=}4\ \text{的}\ (2,2,2)\ \text{给}\ \hat s_1=0 .\\
&\textbf{故缺口}:\ \textbf{一个天然不对称的三项分配}\ =\ \text{Z014 定理的}\ \textbf{扇区选择性}.
\end{aligned}
$$

**一句话：**

$$
\boxed{\;\textbf{$\theta$ 与 $b$ 只是权重的 DFT 相位与模 —— 真正缺的不是相位、不是时间，是"三个代为什么不一样"。}\;}
$$

---

## 附录 A　核验

| 组 | 内容 | 断言 | 出处 |
|:--|:--|--:|:--|
| A | DFT；$\arg\hat s_1$；$|\hat s_1|/|\hat s_0|=b/2$；$\hat s_2=\overline{\hat s_1}$ | 10 | `Z082_weight_phase_check.py` |
| B | $Q=(1+b^2/2)/3$；$b^2=2(3Q-1)$；衔接 $2/9$ | 10 | 同上 |
| C | 三项互异；$L{=}4$ 的 $(2,2,2)\Rightarrow\hat s_1=0$ | 8 | 同上 |
| D | $\mathbb Z_3$ 对称性；$b$ 与 $\theta$ 同源 | 6 | 同上 |
| E | 语料引文（Z014 定理 27／31） | 8 | 同上 |
| | **合计** | **45** | |

**核验纪律**：数值只算旁证；**本文不对任何接近作断言**。

---

*推导＋落盘稿 · 已落盘 · 未改动任何公理 · 完*
